#!/usr/bin/env python3

# Copyright 2021 Daniel Estevez <daniel@destevez.net>
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

import argparse
import json
import pathlib

import numpy as np
from astropy.time import Time, TimeDelta

import phasing
import utils
import uvfits


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
    parser.add_argument('input', metavar='INPUT', help='input directory')
    parser.add_argument(
        't_avg', metavar='T_AVG', type=int, help='time average')
    parser.add_argument(
        'f_avg', metavar='F_AVG', type=int, help='frequency average')
    parser.add_argument(
        'output', metavar='OUTPUT', help='output directory')
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
    sources = utils.read_sources(obs_metadata)

    input_path = pathlib.Path(args.input)
    output_path = pathlib.Path(args.output)
    scans = [s.name.replace('.json', '')
             for s in input_path.glob('*.json')]

    for scan in scans:
        print(f'Processing {scan}')

        x = np.fromfile(input_path / scan, dtype='complex64')

        with open(input_path / f'{scan}.json') as f:
            metadata = json.load(f)

        ants = metadata['antenna_names']
        n_ants = len(ants)
        n_baselines = n_ants * (n_ants + 1) // 2
        n_pols = metadata['polarizations']
        n_chans = metadata['channels']
        chan_fs = metadata['channel_width']
        n = n_baselines * n_pols * n_chans
        x = x[:x.size//n*n].reshape((-1, n_chans, n_baselines, n_pols**2))

        t_sync = read_sync_time(metadata)
        t0 = t_sync + TimeDelta(metadata['first_seq_num'] / chan_fs,
                                format='sec')
        f_first = metadata['first_channel_center_freq']
        f = n_chans/2*chan_fs + f_first
        T = metadata['ntime'] / chan_fs

        source_name = metadata['object_name']
        source = sources[source_name]
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


if __name__ == '__main__':
    main()
