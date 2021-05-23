#!/usr/bin/env python3

# Copyright 2021 Daniel Estevez <daniel@destevez.net>
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

import astropy.units as u
import numpy as np
import pyuvdata


ata_location = np.array([40.8174, -121.472, 1043])


def export_uvfits(t0, T, corrs, uvw,
                  source_name, source,
                  freq, ch_bw,
                  baselines,
                  ant_coordinates,
                  t_avg, f_avg,
                  out_path,
                  telescope_location=ata_location,
                  ref_coordinates=None,
                  instrument='gnuradio1',
                  telescope_name='ATA'):
    """Export visibility data as an UVFITS file using pyuvdata

    The output file name is chosen automatically according to the
    observation start JD and the name of the source.

    Args:
        t0: initial timestamp for the visibility data
        T: integration time
        corrs: visibility data, indexed as
            (time, channel, baseline, polarization)
        uvw: UVW baseline coordinates in metres, indexed as
            (time, baseline, uvw)
        source_name: source name, for the UVFITS metadata
        source: source as an SkyCoord
        freq: sky frequency corresponding to the centre frequency in Hz
        ch_bw: bandwidth or sample rate of each channel in Hz
        baselines: list of baselines names written as '1a'-'2b', etc.
        ant_coordinates: dictionary of antenna ECEF coordinates indexed by
            antenna name
        t_avg: blocks of time integrations to average over
        f_avg: blocks of frequency channels to average over
        out_path: output path, as a Path object
        telescope_location: telescope lat, long, height, for the UVFITS
            metadata
        ref_coordinates: reference ECEF coordinates for the antenna
            coordinates; by default, the average of the antenna positions is
            used
        instrument: instrument name, for the UVFITS header
        telescope_name: telescope name, for the UVFITS header
    """
    if ref_coordinates is None:
        ref_coordinates = np.average(np.array(list(ant_coordinates.values())),
                                     axis=0)

    # Time average
    # TODO: do not drop the last fractional average
    uv = (
        corrs[:corrs.shape[0]//t_avg*t_avg]
        .reshape((-1, t_avg, corrs.shape[1], corrs.shape[2], corrs.shape[3]))
        )
    if uv.size == 0:
        # Not enough data. Do nothing.
        return
    uv = np.average(uv, axis=1)
    # Frequency average
    uv = uv.reshape((uv.shape[0], -1, f_avg, uv.shape[2], uv.shape[3]))
    uv = np.average(uv, axis=2)

    t = t0 + T*t_avg*u.s*np.arange(uv.shape[0])
    uvw = uvw[::t_avg]

    # UVFITS general metadata
    UV = pyuvdata.UVData()
    UV.telescope_location_lat_lon_alt_degrees = telescope_location
    UV.instrument = instrument
    UV.telescope_name = telescope_name
    UV.object_name = source_name
    UV.history = ''
    UV.vis_units = 'UNCALIB'
    UV._set_phased()

    # UVFITS source data
    UV.phase_center_ra = source.ra.rad
    UV.phase_center_dec = source.dec.rad
    UV.phase_center_epoch = 2000.0

    # UVFITS time axis data
    UV.Ntimes = uv.shape[0]
    UV.Nbls = uv.shape[2]
    UV.Nblts = UV.Ntimes * UV.Nbls
    UV.time_array = np.empty(UV.Nblts, 'float')
    UV.time_array.reshape((-1, UV.Nbls))[:] = t.jd[:, np.newaxis]
    UV.integration_time = np.empty(UV.Nblts, 'float')
    UV.integration_time[:] = T * t_avg
    UV.set_lsts_from_time_array()

    # UVFITS antenna data
    used_antennas = ({b.split('-')[0] for b in baselines}
                     | {b.split('-')[1] for b in baselines})
    UV.Nants_data = len(used_antennas)
    UV.Nants_telescope = len(ant_coordinates)
    UV.antenna_names = list(sorted(ant_coordinates))
    ant_numbers = {k: j for j, k in enumerate(sorted(ant_coordinates))}
    UV.antenna_numbers = np.arange(UV.Nants_telescope)
    UV.ant_1_array = np.empty(UV.Nblts, 'int')
    UV.ant_1_array.reshape((-1, UV.Nbls))[:] = [ant_numbers[b.split('-')[0]]
                                                for b in baselines]
    UV.ant_2_array = np.empty(UV.Nblts, 'int')
    UV.ant_2_array.reshape((-1, UV.Nbls))[:] = [ant_numbers[b.split('-')[1]]
                                                for b in baselines]
    UV.antenna_positions = np.array(
        [ant_coordinates[a] - ref_coordinates for a in UV.antenna_names],
        'float')
    UV.baseline_array = UV.antnums_to_baseline(UV.ant_1_array, UV.ant_2_array)

    # UVFITS frequency data
    UV.channel_width = ch_bw * f_avg
    total_fs = ch_bw * corrs.shape[1]
    UV.freq_array = np.array(
        [freq + np.fft.fftshift(np.fft.fftfreq(uv.shape[1], 1/total_fs))])
    UV.Nfreqs = UV.freq_array.size
    UV.spw_array = np.array([0])  # only one spectral window supported
    UV.Nspws = UV.spw_array.size

    # UVFITS linear polarization
    UV.polarization_array = np.array([-5, -6, -7, -8], 'int')
    UV.Npols = UV.polarization_array.size

    # UVFITS UVW coordinates
    UV.uvw_array = np.zeros((UV.Nblts, 3), 'float')
    UV.uvw_array[:] = uvw.reshape((-1, 3))[:UV.Nblts]

    # UVFITS visibility data
    UV.data_array = np.zeros((UV.Nblts, UV.Nspws, UV.Nfreqs, UV.Npols),
                             'complex')
    UV.flag_array = np.zeros((UV.Nblts, UV.Nspws, UV.Nfreqs, UV.Npols),
                             'bool')
    UV.nsample_array = np.empty((UV.Nblts, UV.Nspws, UV.Nfreqs, UV.Npols),
                                'float')
    UV.nsample_array[:] = 1
    # Reorder axes to (time, baseline, freq, polarization)
    uv = np.einsum('ijkl->ikjl', uv)
    # Flatten along (time, baseline) axes
    uv = uv.reshape((-1, uv.shape[2], uv.shape[3]))
    # Assign, taking into account:
    #     - Only one spectral window
    #     - Polarization reordering
    UV.data_array[:] = uv[:, np.newaxis, :, [0, 3, 1, 2]]

    # Workaround CASA bug involving flipping baselines on concat():
    # we flip baselines if necessary so as to get ant_1 <= ant_2
    flip_blts = UV.ant_1_array > UV.ant_2_array
    UV.ant_1_array[flip_blts], UV.ant_2_array[flip_blts] = (
        UV.ant_2_array[flip_blts], UV.ant_1_array[flip_blts])
    UV.uvw_array[flip_blts] *= -1
    # Perform complex conjugate and exchange XY and YX polarizations
    UV.data_array[flip_blts] = np.conjugate(
        UV.data_array[flip_blts][..., [0, 1, 3, 2]])

    # Correct for weird rotECEF convention
    UV.antenna_positions = pyuvdata.utils.ECEF_from_rotECEF(
        UV.antenna_positions, UV.telescope_location_lat_lon_alt[1])

    # Write UVFITS file
    UV.write_uvfits(str(out_path / f'{t0.jd}_{source_name}.uvfits'),
                    force_phase=True, spoof_nonessential=True)
    return
