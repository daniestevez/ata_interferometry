#!/usr/bin/env python3

import numpy as np
from astropy.time import Time, TimeDelta
import astropy.units as u
import pyuvdata

import argparse
import json
import pathlib

import fringe_stop
import grmeta_to_uvfits

# this is different from grmeta_to_uvfits because it handles
# a list of antennas instead of a single baseline
def export_uvfits(t0, T, y, uvw, source_name, source, fc, samp_rate,
                  baselines, ant_coordinates,
                  t_avg, f_avg, out_path,
                  telescope_location = np.array([40.8174, -121.472, 1043]),
                  refant = '1h'):
    # time average
    # TODO: do not drop the last fractional average
    uv = y[:y.shape[0]//t_avg*t_avg].reshape((-1, t_avg, y.shape[1], y.shape[2], y.shape[3]))
    if uv.size == 0:
        return
    uv = np.average(uv, axis = 1)
    # frequency average
    uv = uv.reshape((uv.shape[0], -1, f_avg, uv.shape[2], uv.shape[3]))
    uv = np.average(uv, axis = 2)
    t = t0 + T*t_avg*u.s*np.arange(uv.shape[0])
    uvw = uvw[::t_avg]

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
    used_antennas = {b.split('-')[0] for b in baselines} | {b.split('-')[1] for b in baselines}
    UV.Nants_data = len(used_antennas)
    UV.Nants_telescope = len(ant_coordinates)
    UV.Ntimes = uv.shape[0]
    UV.Nbls = uv.shape[2]
    UV.Nblts = UV.Ntimes * UV.Nbls
    UV.antenna_names = list(sorted(ant_coordinates))
    ant_numbers = {k : j for j, k in enumerate(sorted(ant_coordinates))}
    UV.antenna_numbers = np.arange(UV.Nants_telescope)
    UV.ant_1_array = np.empty(UV.Nblts, dtype = 'int')
    UV.ant_1_array.reshape((-1, UV.Nbls))[:] = \
      [ant_numbers[b.split('-')[0]] for b in baselines]
    UV.ant_2_array = np.empty(UV.Nblts, dtype = 'int')
    UV.ant_2_array.reshape((-1, UV.Nbls))[:] = \
      [ant_numbers[b.split('-')[1]] for b in baselines]
    UV.antenna_positions = np.array([(ant_coordinates[a] - ant_coordinates[refant]).value
                                         for a in UV.antenna_names],
                                        dtype = 'float')
    UV.baseline_array = UV.antnums_to_baseline(UV.ant_1_array, UV.ant_2_array)
    UV.time_array = np.empty(UV.Nblts, dtype = 'float')
    UV.time_array.reshape((-1, UV.Nbls))[:] = t.jd[:,np.newaxis]
    UV.integration_time = np.empty(UV.Nblts, dtype = 'float')
    UV.integration_time[:] = T * t_avg
    UV.set_lsts_from_time_array()
    UV.freq_array = np.array([fc + np.fft.fftshift(np.fft.fftfreq(uv.shape[1], 1/samp_rate))])
    UV.spw_array = np.array([0])
    UV.channel_width = samp_rate/uv.shape[1]
    UV.polarization_array = np.array([-5,-6,-7,-8], dtype = 'int')
    UV.Nfreqs = UV.freq_array.size
    UV.Npols = UV.polarization_array.size
    UV.Nspws = UV.spw_array.size
    UV.uvw_array = np.zeros((UV.Nblts, 3), dtype = 'float')
    UV.uvw_array[:] = uvw.reshape((-1,3))[:UV.Nblts]

    UV.data_array = np.zeros((UV.Nblts, 1, UV.Nfreqs, UV.Npols), dtype = 'complex')
    UV.flag_array = np.zeros((UV.Nblts, 1, UV.Nfreqs, UV.Npols), dtype = 'bool')
    UV.nsample_array = np.empty((UV.Nblts, 1, UV.Nfreqs, UV.Npols), dtype = 'float')
    UV.nsample_array[:] = samp_rate * T * t_avg

    uv = np.einsum('ijkl->ikjl', uv)
    uv = uv.reshape((-1, uv.shape[2], uv.shape[3]))
    UV.data_array[:] = uv[:,np.newaxis,:,[0,3,1,2]]
    
    # correct for weird rotECEF convention
    UV.antenna_positions = pyuvdata.utils.ECEF_from_rotECEF(UV.antenna_positions,
                                                        UV.telescope_location_lat_lon_alt[1])
    UV.write_uvfits(str(out_path / f'{t0.jd}_{source_name}.uvfits'),
                    force_phase=True, spoof_nonessential=True)

def parse_args():
    parser = argparse.ArgumentParser(description =
                'Convert x-engine raw with JSON metadata correlator output to UVFITS')
    parser.add_argument('coordinates', metavar = 'COORDINATES', help = 'antenna ECEF coordinates')
    parser.add_argument('delays', metavar = 'DELAYS', help = 'antenna delays JSON file (in ns)')
    parser.add_argument('json', metavar = 'JSON', help = 'observation JSON metadata')
    parser.add_argument('input', metavar = 'INPUT', help = 'input directory')
    parser.add_argument('t_avg', metavar = 'T_AVG', type = int, help = 'time average')
    parser.add_argument('f_avg', metavar = 'F_AVG', type = int, help = 'frequency average')
    parser.add_argument('output', metavar = 'OUTPUT', help = 'output directory')
    return parser.parse_args()

def main():
    args = parse_args()

    ant_coordinates = fringe_stop.read_coordinates(args.coordinates)
    with open(args.json) as f:
        obs_metadata = json.load(f)
    with open(args.delays) as f:
        delays = json.load(f)
    sources = grmeta_to_uvfits.read_sources(obs_metadata)

    input_path = pathlib.Path(args.input)
    output_path = pathlib.Path(args.output)
    scans = [s.name.replace('.json','')
               for s in input_path.glob('*.json')]
    for scan in scans:
        print(f'Processing {scan}')

        x = np.fromfile(input_path / scan, dtype = 'complex64')
        
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
        t_sync = Time(metadata['snap_sync_time'])
        t0 = t_sync + TimeDelta(metadata['first_seq_num'] / chan_fs, format = 'sec')
        f_first = metadata['first_channel_center_freq']
        f = n_chans/2*chan_fs + f_first
        T = metadata['ntime'] / chan_fs
        source_name = metadata['object_name']
        source = sources[source_name]

        y = np.empty_like(x)
        # axes for uvw are (time, baseline, uvw)
        uvw = np.zeros((x.shape[0], x.shape[2], 3), dtype = 'float64')

        baselines = [f'{a}-{b}' for j,a in enumerate(ants) for b in ants[:j+1]]

        for j, baseline in enumerate(baselines):
            if baseline.split('-')[0] == baseline.split('-')[1]:
                # autocorrelation. no need to do anything
                y[:,:,j] = x[:,:,j]
            else:
                stop = fringe_stop.stop_baseline(t0, x[:, :, j], f,
                                                 source, baseline,
                                                 T, chan_fs, ant_coordinates, delays)
                y[:,:,j] = stop[0]
                uvw[:,j,:] = stop[1].T
        
        export_uvfits(t0, T, y, uvw,
                      source_name, source, f, n_chans * chan_fs,
                      baselines, ant_coordinates,
                      args.t_avg, args.f_avg, output_path)

if __name__ == '__main__':
    main()
