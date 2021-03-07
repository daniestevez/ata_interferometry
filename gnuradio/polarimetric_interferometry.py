#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# SPDX-License-Identifier: GPL-3.0
#
# GNU Radio Python Flow Graph
# Title: Polarimetric Interferometry
# Author: Daniel Estévez <daniel@destevez.net>
# GNU Radio version: 3.8.2.0

from datetime import datetime
from gnuradio import blocks
from gnuradio import fft
from gnuradio.fft import window
from gnuradio import gr
from gnuradio.filter import firdes
import sys
import signal
from argparse import ArgumentParser
from gnuradio.eng_arg import eng_float, intx
from gnuradio import eng_notation
from gnuradio import gr, blocks
from gnuradio import uhd
import time
import pmt


class flowgraph(gr.top_block):

    def __init__(self, directory="/mnt/buf0", gaindB=35, samp_rate=30.72e6, src_name='unknown_source', lo_freq = -1, antennas = 'no_antennas',
                     delay = 0, nfft = 256, only_autocorr = False, parallel_ffts = 2):
        gr.top_block.__init__(self, "Polarimetric Interferometry")
        self.directory = directory
        self.src_name = src_name
        self.samp_rate = samp_rate
        self.lo_freq = lo_freq
        self.antennas = antennas
        self.extra_pmt = None
        print('parallel FFTs:', parallel_ffts)
        
        self.nfft = nfft
        nthreads = 1
        fc = 512e6
        output_rate = 0.1
        self.decim = decim = int(output_rate * samp_rate / nfft)
        interleave_block = decim
        self.nchans = nchans = 4
        self.npairs = npairs = nchans * (nchans-1) // 2
        self.delay = delay
        self.do_crossproducts = not only_autocorr

        self.uhd_usrp_source_0 = uhd.usrp_source(
            ",".join(("addr0=10.11.1.20, addr1=10.11.1.22", "")),
            uhd.stream_args(
                cpu_format="sc16",
                otw_format="sc16",
                args='',
                channels=list(range(nchans)),
            ),
        )
        for mb in range(2):
            self.uhd_usrp_source_0.set_time_source('external', mb)
            self.uhd_usrp_source_0.set_clock_source('external', mb)
        for chan in range(nchans):
            self.uhd_usrp_source_0.set_center_freq(fc, chan)
            self.uhd_usrp_source_0.set_gain(gaindB, chan)
            self.uhd_usrp_source_0.set_antenna('RX2', chan)
        self.uhd_usrp_source_0.set_samp_rate(samp_rate)
        self.uhd_usrp_source_0.set_time_unknown_pps(uhd.time_spec())
        self.uhd_usrp_source_0.set_time_next_pps(uhd.time_spec(int(time.time())+1))
        self.uhd_usrp_source_0.set_lo_export_enabled(True, "lo1", 0)
        self.uhd_usrp_source_0.set_rx_lo_dist(True, "LO_OUT_0")
        self.uhd_usrp_source_0.set_rx_lo_dist(True, "LO_OUT_1")
        for chan in range(nchans):
            self.uhd_usrp_source_0.set_lo_source("external", "lo1", chan)
        time.sleep(1)

        #self.delays = [blocks.delay(gr.sizeof_short*2, 0) for _ in range(nchans)]
        self.set_delays()
        if parallel_ffts > 1:
            self.deinterleaves = [blocks.deinterleave(gr.sizeof_short*2, nfft*interleave_block) for _ in range(nchans)]
        self.ishort_to_complexs = [[blocks.interleaved_short_to_complex(True, False) for _ in range(parallel_ffts)] for _ in range(nchans)]
        self.stream_to_vectors = [[blocks.stream_to_vector(gr.sizeof_gr_complex, nfft) for _ in range(parallel_ffts)] for _ in range(nchans)]
        self.ffts = [[fft.fft_vcc(nfft, True, window.rectangular(nfft), True, nthreads) for _ in range(parallel_ffts)] for _ in range(nchans)]

        self.powers = [[blocks.complex_to_mag_squared(nfft) for _ in range(parallel_ffts)] for _ in range(nchans)]
        self.power_avgs = [[blocks.integrate_ff(decim, nfft) for _ in range(parallel_ffts)] for _ in range(nchans)]
        self.interleave_powers = [blocks.interleave(gr.sizeof_float*nfft, interleave_block//decim) for _ in range(nchans)]

        if self.do_crossproducts:
            self.crossprods = [[blocks.multiply_conjugate_cc(nfft) for _ in range(parallel_ffts)] for _ in range(npairs)]
            self.crossprods_avgs = [[blocks.integrate_cc(decim, nfft) for _ in range(parallel_ffts)] for _ in range(npairs)]
            self.interleave_crossprods = [blocks.interleave(gr.sizeof_gr_complex*nfft, interleave_block//decim) for _ in range(npairs)]

        for chan in range(nchans):
            #self.connect((self.uhd_usrp_source_0, chan), self.delays[chan], self.deinterleaves[chan])
            if parallel_ffts > 1:
                self.connect((self.uhd_usrp_source_0, chan), self.deinterleaves[chan])
                for j in range(parallel_ffts):
                    self.connect((self.deinterleaves[chan], j), self.ishort_to_complexs[chan][j])
            else:
                 self.connect((self.uhd_usrp_source_0, chan), self.ishort_to_complexs[chan][0])
            for j in range(parallel_ffts):
                self.connect(self.ishort_to_complexs[chan][j], self.stream_to_vectors[chan][j],
                                 self.ffts[chan][j], self.powers[chan][j],
                                 self.power_avgs[chan][j], (self.interleave_powers[chan], j))

        if self.do_crossproducts:
            pairs = [(j,k) for j in range(nchans) for k in range(j+1,nchans)]
            for j, p in enumerate(pairs):
                for k in range(parallel_ffts):
                    for t in range(2):
                        self.connect(self.ffts[p[t]][k], (self.crossprods[j][k], t))
                    self.connect(self.crossprods[j][k], self.crossprods_avgs[j][k], (self.interleave_crossprods[j], k))

    def set_delays(self):
        return # function disabled
        # hardcoded for nchans = 4
        if self.nchans != 4:
            raise ValueError('Delays only implemented for 4 channels')
        if self.delay >= 0:
            a = self.delay
            b = 0
        else:
            a = 0
            b = -self.delay
        for j in range(2):
            self.delays[j].set_dly(a)
        for j in range(2,4):
            self.delays[j].set_dly(b)

    def create_sinks(self):
        timestamp = datetime.utcnow().isoformat()
        source_dict = pmt.dict_add(pmt.make_dict(), pmt.intern('source'), pmt.intern(self.src_name))
        extra_dict = pmt.dict_add(source_dict, pmt.intern('lo_freq'), pmt.from_double(self.lo_freq))
        extra_dict = pmt.dict_add(extra_dict, pmt.intern('antennas'), pmt.intern(self.antennas))
        extra_dict = pmt.dict_add(extra_dict, pmt.intern('delay'), pmt.from_long(self.delay))
        if self.extra_pmt is not None:
            extra_dict = pmt.dict_update(extra_dict, self.extra_pmt)
        self.power_filesinks = [blocks.file_meta_sink(gr.sizeof_float*self.nfft, f'{self.directory}/{self.src_name}_{timestamp}_{self.antennas}_power_{chan}',
                                                          self.samp_rate,
                                                          1.0/(self.decim*self.nfft),
                                                          blocks.GR_FILE_FLOAT,
                                                          False,
                                                          int(10*self.samp_rate/(self.decim*self.nfft)),
                                                          extra_dict,
                                                          True) for chan in range(self.nchans)]
        if self.do_crossproducts:
            self.crossprods_filesinks = [blocks.file_meta_sink(gr.sizeof_gr_complex*self.nfft, f'{self.directory}/{self.src_name}_{timestamp}_{self.antennas}_cross_{pair}',
                                                               self.samp_rate,
                                                               1.0/(self.decim*self.nfft),
                                                               blocks.GR_FILE_FLOAT,
                                                               True,
                                                               int(10*self.samp_rate/(self.decim*self.nfft)),
                                                               extra_dict,
                                                               True) for pair in range(self.npairs)]

    def delete_sinks(self):
        cross = self.crossrods_filesinks if self.do_crossproducts else []
        for f in self.power_filesinks + cross:
            f.close()
            del f
        del self.power_filesinks
        if self.do_crossproducts:
            del self.crossprods_filesinks
        
    def connect_sinks(self):
        for chan in range(self.nchans):
            self.connect(self.interleave_powers[chan], self.power_filesinks[chan])
        if self.do_crossproducts:
            for pair in range(self.npairs):
                self.connect(self.interleave_crossprods[pair], self.crossprods_filesinks[pair])

    def disconnect_sinks(self):
        for chan in range(self.nchans):
            self.disconnect(self.interleave_powers[chan], self.power_filesinks[chan])
        if self.do_crossproducts:
            for pair in range(self.npairs):
                self.disconnect(self.interleave_crossprods[pair], self.crossprods_filesinks[pair])

    def make_sinks(self, lo_freq = None, src_name = None, antennas = None, delay = None, extra_pmt = None):
        if lo_freq is not None:
            self.lo_freq = lo_freq
        if src_name is not None:
            self.src_name = src_name
        if antennas is not None:
            self.antennas = antennas
        if delay is not None:
            self.delay = delay
        if extra_pmt is not None:
            self.extra_pmt = extra_pmt 
        self.create_sinks()
        self.connect_sinks()
        self.set_delays()
            
    def reload_sinks(self, lo_freq = None, src_name = None, antennas = None, delay = None, extra_pmt = None):
        self.disconnect_sinks()
        self.delete_sinks()
        self.make_sinks(lo_freq, src_name, antennas, delay, extra_pmt)

def argument_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "--directory", dest="directory", type=str, default="/mnt/buf0",
        help="Set directory [default=%(default)r]")
    parser.add_argument(
        "-g", "--gaindB", dest="gaindB", type=eng_float, default="45.0",
        help="Set rx gain in decibels [default=%(default)r]")
    parser.add_argument(
        "-b", "--samp-rate", dest="samp_rate", type=eng_float, default="40.96M",
        help="Set sampling rate [Hz] [default=%(default)r]")
    parser.add_argument(
        "--src-name", dest="src_name", type=str, default='unknown_source',
        help="Set src_name [default=%(default)r]")
    parser.add_argument(
        "--lo-freq", dest="lo_freq", type=eng_float, default='-1.',
        help="Set lo_freq [default=%(default)r]")
    parser.add_argument(
        "--antennas", dest="antennas", type=str, default='no_antennas',
        help="Set antennas [default=%(default)r]")
    parser.add_argument(
        "--delay", dest="delay", type=int, default=0,
        help="Set delay [default=%(default)r]")
    parser.add_argument(
        "--nfft", dest="nfft", type=int, default=256,
        help="Set FFT size [default=%(default)r]")
    parser.add_argument(
        "--only-autocorr", dest="only_autocorr", type=bool, default=False,
        help="Compute only autocorrelations [default=%(default)r]")
    parser.add_argument(
        "--parallel-ffts", dest="parallel_ffts", type=int, default=2,
        help="Set parallel FFTs [default=%(default)r]")
    return parser

def main(top_block_cls=flowgraph, options=None):
    if options is None:
        options = argument_parser().parse_args()
    if gr.enable_realtime_scheduling() != gr.RT_OK:
        print("Error: failed to enable real-time scheduling.")
    tb = top_block_cls(directory=options.directory, gaindB=options.gaindB, samp_rate=options.samp_rate,
                           src_name=options.src_name, lo_freq=options.lo_freq, antennas=options.antennas,
                           delay=options.delay, nfft=options.nfft, only_autocorr=options.only_autocorr,
                           parallel_ffts=options.parallel_ffts)
    def sig_handler(sig=None, frame=None):
        tb.stop()
        tb.wait()

        sys.exit(0)

    signal.signal(signal.SIGINT, sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)

    tb.make_sinks()
    tb.start()

    tb.wait()


if __name__ == '__main__':
    main()
