#!/usr/bin/env python3

# Copyright 2021 Daniel Estevez <daniel@destevez.net>
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

import argparse
import json
import pathlib

from astropy.time import Time, TimeDelta
from gnuradio.blocks import parse_file_metadata
import pmt
import numpy as np

import phasing
import utils
import uvfits


def read_meta(path):
    with open(str(path) + '.hdr', 'rb') as f:
        header = f.read(parse_file_metadata.HEADER_LENGTH)
        header = pmt.deserialize_str(header)
        header = parse_file_metadata.parse_header(header, False)
        extra = f.read(header['extra_len'])
    extra = pmt.deserialize_str(extra)
    parse_file_metadata.parse_extra_dict(extra, header, False)
    return header


def total_timestamps(r, path, nfft):
    size = (path / (r + f'power_0')).stat().st_size
    return size // (4 * nfft)  # 4 is sizeof(float32)


def load_correlations(r, path, nfft, first_timestamp=0, n_timestamps=None):
    N = (n_timestamps if n_timestamps is not None
         else total_timestamps(r, path, nfft))
    offset_c64 = first_timestamp * nfft * 8
    offset_f32 = first_timestamp * nfft * 4
    # Order is (time, frequency, baseline, polarization)
    x = np.zeros((N, nfft, 3, 4), dtype='complex64')
    min_size = N
    for j in range(4):
        baseline = 2 if j >= 2 else 0
        polarization = 0 if j % 2 == 0 else 3
        a = np.fromfile(path / (r + f'power_{j}'), dtype='float32',
                        offset=offset_f32, count=N*nfft)
        a = a[:a.size//nfft*nfft]
        a = a.reshape((-1, nfft))
        x[:a.shape[0], :, baseline, polarization] = a
        min_size = min(a.shape[0], min_size)
    for j in range(1, 5):
        # Swap polarizations XY and YX because we are swapping the baseline
        # order below also.
        polarization = [0, 2, 1, 3][j - 1]
        a = np.fromfile(path / (r + f'cross_{j}'), dtype='complex64',
                        offset=offset_c64, count=N*nfft)
        a = a[:a.size//nfft*nfft]
        a = a.reshape((-1, nfft))
        # Complex conjugate here is because the baseline in the correlation
        # file is 0-1 rather than 1-0 (which is required for lower triangular
        # format)
        x[:a.shape[0], :, 1, polarization] = np.conjugate(a)
        min_size = min(a.shape[0], min_size)
    for j in (0, 5):
        baseline = 0 if j == 0 else 2
        a = np.fromfile(path / (r + f'cross_{j}'), dtype='complex64',
                        offset=offset_c64, count=N*nfft)
        a = a[:a.size//nfft*nfft]
        a = a.reshape((-1, nfft))
        x[:a.shape[0], :, baseline, 1] = a  # XY
        x[:a.shape[0], :, baseline, 2] = np.conjugate(a)  # YX
        min_size = min(a.shape[0], min_size)
    return x[:min_size]


def parse_args():
    parser = argparse.ArgumentParser(
        description='Convert GNU Radio raw metadata '
        'correlator output to UVFITS')
    parser.add_argument(
        'coordinates', metavar='COORDINATES', help='antenna ECEF coordinates')
    parser.add_argument(
        'delays', metavar='DELAYS', help='antenna delays JSON file (in ns)')
    parser.add_argument(
        'json', metavar='JSON', help='observation JSON metadata')
    parser.add_argument(
        'input', metavar='INPUT', help='input directory')
    parser.add_argument(
        't_avg', metavar='T_AVG', type=int, help='time average')
    parser.add_argument(
        'f_avg', metavar='F_AVG', type=int, help='frequency average')
    parser.add_argument('output', metavar='OUTPUT', help='output directory')
    parser.add_argument('--max_timestamps', type=int,
                        help='maximum number of timestamps in RAM')
    return parser.parse_args()


def pol_swap(x, baseline):
    swapped = {'1c', '2b', '3l'}
    ant0, ant1 = baseline.split('-')
    pol_array = np.arange(4)
    if ant0 in swapped:
        pol_swap_ant0(x)
    if ant1 in swapped:
        pol_swap_ant1(x)


def pol_swap_ant0(x):
    # autocorrelations
    x[:, :, 0, :] = x[:, :, 0, [3, 2, 1, 0]]
    # crosscorrelations on baseline 1-0
    x[:, :, 1, :] = x[:, :, 1, [1, 0, 3, 2]]


def pol_swap_ant1(x):
    # autocorrelations
    x[:, :, 2, :] = x[:, :, 0, [3, 2, 1, 0]]
    # crosscorrelations on baseline 1-0
    x[:, :, 1, :] = x[:, :, 1, [2, 3, 0, 1]]


def main():
    args = parse_args()

    ant_coordinates = utils.read_antenna_coordinates(args.coordinates)
    with open(args.json) as f:
        metadata = json.load(f)
    with open(args.delays) as f:
        delays = json.load(f)
    sources = utils.read_sources(metadata)

    nfft = metadata['channels']
    T = metadata['t_int']
    samp_rate = metadata['samp_rate']

    input_path = pathlib.Path(args.input)
    output_path = pathlib.Path(args.output)
    scans = [s.name.replace('cross_0', '')
             for s in input_path.glob('*cross_0')]
    for scan in scans:
        print(f'Processing {scan}*')
        t0 = Time(read_meta(input_path / f'{scan}cross_0')['rx_time'],
                  format='unix')
        if t0.value == 0:
            # Some files have invalid timestamps in the metadata.
            # We read the timestamp from the filename in those cases.
            t0 = Time(np.datetime64(scan.split('_')[1]))
        freq = pmt.to_double(
            read_meta(input_path / f'{scan}cross_0')['lo_freq']) * 1e6

        source_name = scan.split('_')[0]
        source = sources[source_name]
        baseline = str(read_meta(input_path / f'{scan}cross_0')['antennas'])
        ants = [baseline.split('-')[0], baseline.split('-')[1]]
        baselines = ['-'.join([ants[0], ants[0]]), baseline,
                     '-'.join([ants[1], ants[1]])]

        timestamps = total_timestamps(scan, input_path, nfft)
        offsets = ([0] if args.max_timestamps is None
                   else range(0, timestamps, args.max_timestamps))
        if len(offsets) > 1:
            print(f'File is too large, processing in {len(offsets)} chunks')
        for offset in offsets:
            x = load_correlations(
                scan, input_path, nfft, offset,
                timestamps if len(offsets) == 1 else args.max_timestamps)

            # Swap polarizations if needed
            pol_swap(x, baseline)

            if x.size == 0:
                print('Skipping empty file')
                break

            t = t0 + TimeDelta(T * offset, format='sec')
            uvw = phasing.apply_phasing(
                phasing.midpoint_timestamps(t, x.shape[0], T),
                x, source, freq, samp_rate / nfft,
                utils.filter_by_antennas(ants, ant_coordinates),
                utils.filter_by_antennas(ants, delays) * 1e-9)

            uvfits.export_uvfits(
                t, T, x, uvw, source_name, source, freq, samp_rate / nfft,
                baselines, ant_coordinates, args.t_avg, args.f_avg,
                output_path)


if __name__ == '__main__':
    main()
