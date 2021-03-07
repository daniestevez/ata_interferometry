#!/usr/bin/env python3

import numpy as np
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, AltAz, ICRS, ITRS
import astropy.units as u
import astropy.constants as const

import argparse
import json

def read_coordinates(path):
    with open(path) as f:
        l = f.readlines()[1:]
    
    ants = [ll.split(',')[0] for ll in l]
    ecef = np.empty((len(ants), 3))
    for j in range(3):
        ecef[:,j] = np.array([float(ll.split(',')[j+1]) for ll in l])
    return {a.lower() : e * u.m for a, e in zip(ants, ecef)}

def stop_baseline(t0, cross, freq, source, baseline, T, ch_fs,
                      ant_coordinates, ant_delays_ns):
    ant1, ant2 = baseline.split('-')
    delay_offset = (ant_delays_ns[ant1] - ant_delays_ns[ant2]) * 1e-9
    baseline_itrs = ant_coordinates[ant1] - ant_coordinates[ant2]
    north_radec = [source.ra.deg, source.dec.deg + 90]
    if north_radec[1] > 90:
        north_radec[1] = 180 - north_radec[1]
        north_radec[0] = 180 + north_radec[0]
    north = SkyCoord(ra = north_radec[0]*u.deg, dec = north_radec[1]*u.deg)
    
    f_obs = freq * u.Hz
    ts = t0 + TimeDelta(T, format = 'sec') * np.arange(cross.shape[0])
    source_itrs = source.transform_to(ITRS(obstime = Time(ts))).cartesian
    north_itrs = north.transform_to(ITRS(obstime = Time(ts))).cartesian
    east_itrs = north_itrs.cross(source_itrs)
    ww = source_itrs.xyz.T.dot(baseline_itrs)
    vv = north_itrs.xyz.T.dot(baseline_itrs)
    uu = east_itrs.xyz.T.dot(baseline_itrs)
    w_cycles = (ww/const.c*f_obs).to(1).value
    w_seconds = (ww/const.c).to(u.s).value
    phase_corr = np.exp(-1j*2*np.pi*w_cycles)[:,np.newaxis,np.newaxis]
    nfft = cross.shape[1]
    ch_idx = np.arange(-nfft//2,nfft//2)[:,np.newaxis]
    delay_corr = np.exp(1j*2*np.pi*(delay_offset - w_seconds[:,np.newaxis,np.newaxis])*ch_idx*ch_fs)
    return (cross * phase_corr * delay_corr, 
            np.array([uu.value,vv.value,ww.value]))

def parse_args():
    parser = argparse.ArgumentParser(description='Perform fringe stopping on x-engine output')
    parser.add_argument('coordinates', metavar = 'COORDINATES', help = 'antenna ECEF coordinates')
    parser.add_argument('delays', metavar = 'DELAYS', help = 'antenna delays JSON file (in ns)')
    parser.add_argument('xengine', metavar = 'XENGINE', help = 'x-engine output file')
    parser.add_argument('json', metavar = 'JSON', help = 'observation JSON metadata')
    parser.add_argument('output', metavar = 'OUTPUT', help = 'correlations output file')
    parser.add_argument('uvw', metavar = 'UVW', help = 'UVW coordinates output file')
    return parser.parse_args()

def main():
    args = parse_args()
    ant_coordinates = read_coordinates(args.coordinates)
    
    x = np.fromfile(args.xengine, dtype = 'complex64')

    with open(args.json) as f:
        metadata = json.load(f)

    with open(args.delays) as f:
        delays = json.load(f)

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
    source = SkyCoord.from_name(metadata['object_name'])

    x_stop = np.empty_like(x)
    # axes for uvw are (time, baseline, uvw)
    uvw = np.zeros((x.shape[0], x.shape[2], 3), dtype = 'float64')

    baselines = [f'{a}-{b}' for j,a in enumerate(ants) for b in ants[:j+1]]

    for j, baseline in enumerate(baselines):
        if baseline.split('-')[0] == baseline.split('-')[1]:
            # autocorrelation. no need to do anything
            x_stop[:,:,j] = x[:,:,j]
        else:
            stop = stop_baseline(t0, x[:, :, j], f,
                                 source, baseline,
                                 T, chan_fs, ant_coordinates, delays)
            x_stop[:,:,j] = stop[0]
            uvw[:,j,:] = stop[1].T
        
    x_stop.tofile(args.output)
    uvw.tofile(args.uvw)

if __name__ == '__main__':
    main()
