#!/home/obsuser/miniconda3/bin/python

# Copyright 2021 Daniel Estevez <daniel@destevez.net>
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

import argparse
import json
import pathlib

from astropy.coordinates import Angle, SkyCoord

import numpy as np
from astropy.time import Time, TimeDelta

import phasing
import utils
import uvfits
import astropy


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Convert x-engine raw with JSON metadata '
                     'correlator output to UVFITS'))
    parser.add_argument(
        'coordinates', metavar='COORDINATES', help='antenna ECEF coordinates')
    parser.add_argument(
        'delays', metavar='DELAYS', help='antenna delays JSON file (in ns)')
    parser.add_argument(
        'json', metavar='JSON', help='observation JSON metadata')
    parser.add_argument(
        't_avg', metavar='T_AVG', type=int, help='time average')
    parser.add_argument(
        'f_avg', metavar='F_AVG', type=int, help='frequency average')
    parser.add_argument(
        'output', metavar='OUTPUT', help='output directory')
    parser.add_argument(
        '-i', '--infiles', nargs="+", help='input files', required=True)
    return parser.parse_args()


def read_sync_time(metadata):
    # Both snap_sync_time and sync_timestamp are
    # allowed keywords for the SNAP sync timestamp
    try:
        sync = metadata['snap_sync_time']
    except KeyError:
        sync = metadata['sync_timestamp']
    if type(sync) is str:
        return Time(sync)
    elif type(sync) is int:
        return Time(sync, format='unix')
    else:
        raise ValueError('Invalid snap_sync_time format')


def main():
    args = parse_args()

    ant_coordinates = utils.read_antenna_coordinates(args.coordinates)
    with open(args.json) as f:
        obs_metadata = json.load(f)
    with open(args.delays) as f:
        delays = json.load(f)

    scans = args.infiles

    output_path = pathlib.Path(args.output)

    
    source_name = obs_metadata['object_name']
    ants = obs_metadata['antenna_names'].split(",")
    n_pols = obs_metadata['polarizations']
    n_chans = obs_metadata['channels']
    chan_fs = obs_metadata['channel_width']
    t_sync = obs_metadata['sync_timestamp']
    n_ants_xgpu = obs_metadata['n_ants_xgpu']

    source = SkyCoord.from_name(source_name)

    n_ants = len(ants)
    n_baselines_valid = n_ants * (n_ants + 1) // 2
    n_baselines_total = n_ants_xgpu * (n_ants_xgpu + 1) // 2
    n = n_baselines_total * n_pols * n_chans
    f_first = obs_metadata['first_channel_center_freq']
    T = obs_metadata['ntime'] / chan_fs

    t0  = Time(t_sync, format='unix')
    t0 += TimeDelta(obs_metadata['first_seq_num'] / chan_fs,
                            format='sec')

    for scan in scans:
        print(f'Processing {scan}')

        x = np.fromfile(scan, dtype='complex64')

        x = x[:x.size//n*n].reshape((-1, n_chans, n_baselines_total, n_pols**2))
        x = x[:, :, :n_baselines_valid, :]

        f = n_chans/2*chan_fs + f_first

        baselines = [f'{a}-{b}' for j, a in enumerate(ants)
                     for b in ants[:j+1]]


        uvw = phasing.apply_phasing(
            phasing.midpoint_timestamps(t0, x.shape[0], T),
            x, source, f, chan_fs,
            utils.filter_by_antennas(ants, ant_coordinates),
            utils.filter_by_antennas(ants, delays) * 1e-9)

        uvfits.export_uvfits(
            t0, T, x, uvw, source_name, source, f, chan_fs, baselines,
            ant_coordinates, args.t_avg, args.f_avg, output_path)

        t0 += TimeDelta(T*x.shape[0], format='sec')


if __name__ == '__main__':
    main()
