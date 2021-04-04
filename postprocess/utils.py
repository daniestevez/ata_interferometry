#!/usr/bin/env python3

# Copyright 2021 Daniel Estevez <daniel@destevez.net>
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

from astropy.coordinates import Angle, SkyCoord
import numpy as np


def read_antenna_coordinates(path):
    """Read the antenna coordinates file

    Args:
        path: path to the file
    Returns:
        A dictionary indexed by antenna name containing the ECEF coordinates
        of each antenna
    """
    with open(path) as f:
        lines = f.readlines()[1:]

    ants = [line.split(',')[0] for line in lines]
    ecef = np.empty((len(ants), 3))
    for j in range(3):
        ecef[:, j] = np.array([float(line.split(',')[j+1]) for line in lines])
    return {a.lower(): e for a, e in zip(ants, ecef)}


def skycoord_from_source(source):
    """Convert a source into a SkyCoord

    If the source is a string, we fetch the source from name.
    Otherwise, we assume it is a list with RADEC coordinates

    Args:
        source: the source to convert

    Returns:
        The source converted to a SkyCoord object
    """
    if type(source) is str:
        return SkyCoord.from_name(source)
    elif type(source) is list and len(source) == 2:
        return SkyCoord(Angle(source[0]), Angle(source[1]))
    else:
        raise ValueError(f'Invalid source description: {source}')


def read_sources(metadata):
    """Read the source from the JSON metadata

    Args:
        metadata: dictionary parsed from JSON metadata

    Returns:
        A dictionary of SkyCoords indexed by source name
    """
    d = metadata['sources']
    return {k: skycoord_from_source(v) for k, v in d.items()}


def filter_by_antennas(antenna_list, d):
    """Returns an array where only the data for the used antennas appears

    Given a dictionary d of numpy arrays indexed by antenna names and a list
    of antenna names, this returns a numpy array whose first axis corresponds
    to the antennas and only contains the information for the antennas that
    appear in antenna_list (and in that order).

    Args:
        antenna_list: a list of antenna names
        d: a dictionary indexed by antenna names containing numpy arrays of
           the same shape

    Returns:
        A numpy array with one dimension more than the arrays in d which
        contains only the information about the antennas appearing in
        antenna_list, and in that order
    """
    return np.stack([d[k] for k in antenna_list])
