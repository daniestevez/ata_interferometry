#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import polarimetric_interferometry as pol

from ATATools import ata_control as ac
from astropy.coordinates import SkyCoord, Angle, EarthLocation, AltAz
import astropy.units as u
from astropy.time import Time
import numpy as np

from gnuradio import gr

import sys
import signal
import time
import pmt
import logging
import collections
import functools

Source = collections.namedtuple('source', ['name', 'coord'])

# do_observation = True uses the array with ATA control
# do_observation = False doesn't touch the array and only uses
# the USRPs (useful for debugging when other people are using
# the array)
do_observation = True

# full_array = True selects all the antennas that are working
# and connected to the USRPs
# full_array = False selects a smaller subset of those antennas
# (the ones that we typically have available on weekend tests)
full_array = True

# List of antennas available at each USRP

if full_array:
    antennas_usrp = [['1a', '1c', '2a', '4g', '4j'],
                     ['1f', '1h', '1k', '2b', '2h', '5c']]
    best_antennas_usrp = ['1c', '1h'] # the antennas with better performance on each USRP
else:
    antennas_usrp = [['4g'], ['2b', '3c']]
    best_antennas_usrp = ['4g', '2b'] # the antennas with better performance on each USRP

antennas_all_usrps = [a for b in antennas_usrp for a in b]
        
# these assume there are only two USRPs
baselines_all = [[a, b] for a in antennas_usrp[0] for b in antennas_usrp[1]]
baselines_calib = [[best_antennas_usrp[0], b] for b in antennas_usrp[1]] \
                 + [[a, best_antennas_usrp[1]] for a in antennas_usrp[0]
                        if a != best_antennas_usrp[0]]

lo = 'd'
if_switch_att = 20

science_target = Source(name = 'CasA',
                        coord = SkyCoord.from_name('Cassiopeia A'))
phase_calibrators = [Source(name = '3C84',
                            coord = SkyCoord.from_name('3C84'))]
compact_sources = [Source(name = '3C84',
                            coord = SkyCoord.from_name('3C84')),
                   Source(name = '3C286',
                            coord = SkyCoord.from_name('3C286')),
                   Source(name = '3C48',
                            coord = SkyCoord.from_name('3C48')),
                   Source(name = '3C147',
                            coord = SkyCoord.from_name('3C147')),
                  ]
    
obs_freq = 4900

# Scan lengths for science target and calibrator (seconds)
science_scan_len = 30
cal_scan_len = 120
compact_sources_scan_len = 60

# number of runs on the science target before re-visiting a calibrator
science_runs = 1

ATA = EarthLocation(lat = 40.816410*u.deg, lon = -121.471828*u.deg, height = 1000*u.m)

def ac_wrap(f):
    """Wraps an ata_control function

    This is used to retry the function call if it fails and
    not run the function call if we are doing a fake observation"""
    def wrap(f, *args):
        if do_observation:
            while True:
                try:
                    return f(*args)
                except Exception as e:
                    logging.info('ATA control call failed')
                    #print(e, file = sys.stderr)
                    time.sleep(1)
                    logging.info('Retrying ATA control call')
        else:
            time.sleep(1)
    return functools.partial(wrap, f)

set_az_el = ac_wrap(ac.set_az_el)
move_ant_group = ac_wrap(ac.move_ant_group)
make_and_track_ra_dec = ac_wrap(ac.make_and_track_ra_dec)
try_on_lnas = ac_wrap(ac.try_on_lnas)
set_freq = ac_wrap(ac.set_freq)
set_rf_switch = ac_wrap(ac.set_rf_switch)
set_atten_thread = ac_wrap(ac.set_atten_thread)
autotune = ac_wrap(ac.autotune)
get_dets = ac_wrap(ac.get_dets)
get_pams = ac_wrap(ac.get_pams)

def park_antennas():
    logging.info('Parking antennas')
    set_az_el(antennas_all_usrps, 180, 18)
    logging.info('Antennas parked')
    logging.info('Freeing antennas')
    move_ant_group(antennas_all_usrps, 'atagr', 'none')

def slew(target):
    logging.info(f'Slewing to {target.name}')
    make_and_track_ra_dec(target.coord.ra.hour,
                          target.coord.dec.deg,
                          antennas_all_usrps)
    logging.info(f'On source {target.name}')

def choose_calibrator():
    """Choose the calibrator with closest elevation to the science target"""
    alt_calibrators = np.array([cal.coord.transform_to(AltAz(obstime = Time.now(),
                                                             location = ATA)).alt.deg
                                for cal in phase_calibrators])
    alt_science = science_target.coord.transform_to(AltAz(obstime = Time.now(),
                                                          location = ATA)).alt.deg
    return phase_calibrators[np.argmin(np.abs(alt_calibrators - alt_science))]

def source_has_set(source, elev_mask_deg = 18):
    alt = source.coord.transform_to(AltAz(obstime = Time.now(),
                                          location = ATA)).alt.deg
    return alt < elev_mask_deg

def set_baseline(b):
    # the baseline b is a list of antennas (as many as USRPs)
    logging.info(f'Setting baseline {b} in IF switch')
    set_rf_switch(b)
    set_atten_thread([[f'{ant}x',f'{ant}y'] for ant in b],
                     [[if_switch_att,if_switch_att] for _ in b])
    logging.info(f'Baseline {b} set in IF switch')

class flowgraph(pol.flowgraph):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.first_record = True

    def start_record(self, source, baseline):
        params = {'src_name' : source.name,
                  'antennas' : '-'.join(baseline),
                  }
        logging.info(f'Starting recording of {params["src_name"]} with {params["antennas"]}')
        if self.first_record:
            self.first_record = False
            self.make_sinks(**params)
            self.start()
        else:
            self.reload_sinks(**params)
            self.unlock()

    def stop_record(self):
        logging.info(f'Stopping recording')
        self.lock()

    def do_run(self, baselines, source, scan_len):
        for baseline in baselines:
            set_baseline(baseline)
            self.start_record(source, baseline)
            time.sleep(scan_len)
            self.stop_record()
            if source_has_set(source):
                return 'source set'

def setup_array():
    logging.info('Reserving antennas')
    try:
        ac.move_ant_group(antennas_all_usrps, 'none', 'atagr')
    except Exception:
        logging.info('Could not reserve antennas. Continuing nevertheless...')
    
    logging.info('Setting up feeds and LO')
    try_on_lnas(antennas_all_usrps)
    set_freq(obs_freq, antennas_all_usrps, lo)
    logging.info('Feeds and LO set')

    # move to calibrator for autotune and first calibration
    cal = choose_calibrator()
    slew(cal)

    logging.info('Running autotune')
    autotune(antennas_all_usrps)
    logging.info(f'Detectors: {get_dets(antennas_all_usrps)}')
    logging.info(f'PAMs: {get_pams(antennas_all_usrps)}')
    logging.info('Autotune done')

    return cal
            
def main():
    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)
    options = pol.argument_parser().parse_args()
    if gr.enable_realtime_scheduling() != gr.RT_OK:
        logging.warning("Error: failed to enable real-time scheduling.")
    # hardcoded parameters regardless of options
    nfft = 2048
    samp_rate = 40.96e6
    gaindB = 45
    logging.info(f'Running with hardcoded parameters: nfft = {nfft}, samp_rate = {samp_rate}, gaindB = {gaindB}')
    tb = flowgraph(directory=options.directory, gaindB=gaindB,
                   samp_rate=samp_rate, src_name=science_target.name,
                   lo_freq = obs_freq, nfft = nfft)
    def sig_handler(sig=None, frame=None):
        tb.stop()
        tb.wait()
        park_antennas()
        logging.shutdown()
        sys.exit(0)

    signal.signal(signal.SIGINT, sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)

    cal = setup_array()

    # observe CasA until it sets
    source_visible = True
    while source_visible:
        # calibration run
        logging.info(f'Starting calibration run with {cal.name}')
        tb.do_run(baselines_calib, cal, cal_scan_len)
        logging.info('Calibration run finished')

        # slew to science target
        slew(science_target)
        for run in range(science_runs):
            logging.info(f'Starting science run {run + 1}/{science_runs}')
            if tb.do_run(baselines_all, science_target, science_scan_len) == 'source set':
                logging.info(f'Science target has set. Ending runs earlier.')
                source_visible = False
                break
            else:
                logging.info(f'Science run {run + 1}/{science_runs} finished')

        # choose new calibrator and slew
        cal = choose_calibrator()
        slew(cal)

    # observe compact calibrators
    logging.info('Now observing compact calibrators')
    while True:
        for source in compact_sources:
            if source_has_set(source):
                logging.info(f'{source.name} has set. Ignoring.')
                continue
            logging.info(f'Observing {source.name}')
            slew(source)
            logging.info(f'Starting run on {source.name}')
            tb.do_run(baselines_all, source, compact_sources_scan_len)
            logging.info(f'Run finished')
        
if __name__ == '__main__':
    main()
