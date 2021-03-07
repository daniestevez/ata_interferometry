#!/usr/bin/env python3

import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u
from gnuradio.blocks import parse_file_metadata
import pmt
import pyuvdata

import argparse
import json
import pathlib

import fringe_stop

def read_meta(path):
    with open(str(path) + '.hdr', 'rb') as f:
        header = f.read(parse_file_metadata.HEADER_LENGTH)
        header = pmt.deserialize_str(header)
        header = parse_file_metadata.parse_header(header, False)
        extra = f.read(header['extra_len'])
    extra = pmt.deserialize_str(extra)
    parse_file_metadata.parse_extra_dict(extra, header, False)
    return header

def observed_powers(r, path, nfft):
    p = [np.fromfile(path / (r + f'power_{j}'), dtype = 'float32').reshape((-1,nfft))
                    for j in range(4)]
    l = np.min([pp.shape[0] for pp in p])
    return np.array([pp[:l] for pp in p])

def observed_crossspectra(r, path, nfft):
    cr = [np.fromfile(path / (r + f'cross_{j}'), dtype = 'complex64').reshape((-1,nfft))
                    for j in range(1,5)]
    l = np.min([c.shape[0] for c in cr])
    return np.array([c[:l] for c in cr])

def parse_args():
    parser = argparse.ArgumentParser(description =
                'Convert GNU Radio raw metadata correlator output to UVFITS')
    parser.add_argument('coordinates', metavar = 'COORDINATES', help = 'antenna ECEF coordinates')
    parser.add_argument('delays', metavar = 'DELAYS', help = 'antenna delays JSON file (in ns)')
    parser.add_argument('json', metavar = 'JSON', help = 'observation JSON metadata')
    parser.add_argument('input', metavar = 'INPUT', help = 'input directory')
    parser.add_argument('t_avg', metavar = 'T_AVG', type = int, help = 'time average')
    parser.add_argument('f_avg', metavar = 'F_AVG', type = int, help = 'frequency average')
    parser.add_argument('output', metavar = 'OUTPUT', help = 'output directory')
    return parser.parse_args()

def skycoord_from_source(source):
    if type(source) is str:
        return SkyCoord.from_name(source)
    elif type(source) is list and len(source) == 2:
        return SkyCoord(Angle(source[0]), Angle(source[1]))
    else:
        raise ValueError(f'Invalid source description: {source}')

def read_sources(metadata):
    d = metadata['sources']
    return {k : skycoord_from_source(v) for k,v in d.items()}

def export_uvfits(t0, T, y, uvw, source_name, source, fc, samp_rate,
                  baseline, ant_coordinates,
                  t_avg, f_avg, out_path,
                  telescope_location = np.array([40.8174, -121.472, 1043]),
                  refant = '1h'):
    # time average
    # TODO: do not drop the last fractional average
    uv = y[:y.shape[0]//t_avg*t_avg].reshape((-1, t_avg, y.shape[1], y.shape[2]))
    uv = np.average(uv, axis = 1)
    # frequency average
    uv = uv.reshape((uv.shape[0], -1, f_avg, uv.shape[2]))
    uv = np.average(uv, axis = 2)
    t = t0 + T*t_avg*u.s*np.arange(uv.shape[0])
    uvw = uvw[:,::t_avg]

    UV = pyuvdata.UVData()
    UV.telescope_location_lat_lon_alt_degrees = telescope_location
    UV.instrument = 'gnuradio1'
    UV.telescope_name = 'ATA'
    UV.object_name = source_name
    UV.history = ''
    UV.vis_units = 'UNCALIB'
    UV._set_phased()
    UV.phase_center_ra = source.ra.rad
    UV.phase_center_dec = source.dec.rad
    UV.phase_center_epoch = 2000.0
    UV.Nants_data = 2
    UV.Nants_telescope = len(ant_coordinates)
    UV.Ntimes = uv.shape[0]
    ant0, ant1 = baseline.split('-')
    UV.antenna_names = list(sorted(ant_coordinates))
    ant_numbers = {k : j + 1 for j, k in enumerate(sorted(ant_coordinates))}
    UV.antenna_numbers = np.arange(1, UV.Nants_telescope + 1)
    UV.ant_1_array = np.empty(UV.Ntimes, dtype = 'int')
    UV.ant_1_array[:] = ant_numbers[ant0]
    UV.ant_2_array = np.empty(UV.Ntimes, dtype = 'int')
    UV.ant_2_array[:] = ant_numbers[ant1]
    UV.antenna_positions = np.array([(ant_coordinates[a] - ant_coordinates[refant]).value
                                         for a in UV.antenna_names],
                                        dtype = 'float')
    UV.baseline_array = UV.antnums_to_baseline(UV.ant_1_array, UV.ant_2_array)
    UV.Nbls = len(np.unique(UV.baseline_array))
    UV.time_array = t.jd
    UV.integration_time = np.empty(UV.time_array.size)
    UV.integration_time[:] = T * t_avg
    UV.set_lsts_from_time_array()
    UV.freq_array = np.array([fc + np.fft.fftshift(np.fft.fftfreq(uv.shape[1], 1/samp_rate))])
    UV.spw_array = np.array([0])
    UV.channel_width = samp_rate/uv.shape[1]
    UV.polarization_array = np.array([-5,-6,-7,-8], dtype = 'int')
    UV.Nfreqs = UV.freq_array.size
    UV.Npols = UV.polarization_array.size
    UV.Nblts = UV.Nbls * UV.Ntimes
    UV.Nspws = UV.spw_array.size
    UV.uvw_array = np.zeros((UV.Nblts, 3), dtype = 'float')
    UV.uvw_array[:] = uvw.T[:UV.Nblts]

    UV.data_array = np.zeros((UV.Nblts, 1, UV.Nfreqs, UV.Npols), dtype = 'complex')
    UV.flag_array = np.zeros((UV.Nblts, 1, UV.Nfreqs, UV.Npols), dtype = 'bool')
    UV.nsample_array = np.empty((UV.Nblts, 1, UV.Nfreqs, UV.Npols), dtype = 'float')
    UV.nsample_array[:] = samp_rate * T * t_avg

    UV.data_array[:] = uv[:,np.newaxis,:,[0,3,1,2]]
    
    # correct for weird rotECEF convention
    UV.antenna_positions = pyuvdata.utils.ECEF_from_rotECEF(UV.antenna_positions,
                                                        UV.telescope_location_lat_lon_alt[1])
    UV.write_uvfits(str(out_path / f'{t0.jd}_{source_name}_{baseline}.uvfits'),
                    force_phase=True, spoof_nonessential=True)

def pol_swap(baseline):
    swapped = {'1c', '2b'}
    ant1, ant2 = baseline.split('-')
    pol_array = np.arange(4)
    if ant1 in swapped:
        pol_array = pol_array[np.array([2,3,0,1])]
    if ant2 in swapped:
        pol_array = pol_array[np.array([1,0,3,2])]
    return pol_array
    
def main():
    args = parse_args()

    ant_coordinates = fringe_stop.read_coordinates(args.coordinates)
    with open(args.json) as f:
        metadata = json.load(f)
    with open(args.delays) as f:
        delays = json.load(f)
    sources = read_sources(metadata)

    nfft = metadata['channels']
    T = metadata['t_int']
    samp_rate = metadata['samp_rate']

    input_path = pathlib.Path(args.input)
    output_path = pathlib.Path(args.output)
    scans = [s.name.replace('cross_0','')
               for s in input_path.glob('*cross_0')]
    for scan in scans:
        print(f'Processing {scan}*')
        t0 = Time(read_meta(input_path / f'{scan}cross_0')['rx_time'],
                  format = 'unix')
        if t0.value == 0:
            # some files have invalid timestamps in the metadata
            # we read the timestamp from the filename in those cases
            t0 = Time(np.datetime64(scan.split('_')[1]))
        freq = pmt.to_double(read_meta(input_path / f'{scan}cross_0')['lo_freq']) * 1e6

        source_name = scan.split('_')[0]
        source = sources[source_name]
        baseline = str(read_meta(input_path / f'{scan}cross_0')['antennas'])
        x = observed_crossspectra(scan, input_path, nfft)

        # hack to swap polarizations on 2b and 1c
        x = x[pol_swap(baseline)]
        
        if x.size == 0:
            print('Skipping empty file')
            continue
        # reshape to x-engine axis order
        # TODO include this into observed_crosspectra()
        x = np.einsum('ijk->jki', x)
        y, uvw = fringe_stop.stop_baseline(t0, x, freq, source, baseline, T,
                                           samp_rate / nfft, ant_coordinates,
                                           delays)
        export_uvfits(t0, T, y, uvw, source_name, source, freq, samp_rate,
                      baseline, ant_coordinates,
                      args.t_avg, args.f_avg, output_path)

if __name__ == '__main__':
    main()
