import pytest

import numpy as np
from math import sin, cos
from copy import copy

from hyvr.geo.contact_surface import ContactSurface
from hyvr.geo.grid import Grid
from hyvr.geo.trough import Trough
from hyvr.geo.trough_ae import TroughAE
import hyvr.geo.trough_utils

class MockStratum:
    bg_facies = 100
    bg_azim = 0
    bg_dip = 0



def trough_test(type_params, trough_params, point, expected):
    """
    Run a test for given type params, trough params (a, b, c, alpha), a given
    point, and assert that the expected values are set
    
    Parameters
    ----------
    type_params : dict
        Dictionary of type params
    trough_params : list/tuple
        (a, b, c, alpha)
    point : list/tuple
        (x, y, z)
    expected : list/tuple
        (expected_facies, expected_azim, expected_dip)
    """
    grid = Grid(0., 0., 0., 0.1, 0.1, 0.1, 1.0, 1.0, 1.0, 0)
    bottom_surface = ContactSurface(grid, mode='flat', z=0)
    top_surface = ContactSurface(grid, mode='flat', z=1)

    a, b, c, alpha = trough_params
    type_params = copy(type_params)
    type_params['length'] = 2*a
    type_params['width'] = 2*b
    type_params['depth'] = c
    type_params['paleoflow'] = [alpha, alpha]

    trough_ae = TroughAE(bottom_surface,
                         top_surface,
                         'test_trough',
                         type_params,
                         MockStratum(),
                         grid)

    # setup containers
    geo_ids = np.zeros(3, dtype=np.int32)
    angles = np.zeros(2, dtype=np.float)

    # test
    x, y, z = point
    trough_ae.maybe_assign_points_to_object(0, geo_ids, angles, x, y, z, 0, 0, grid)
    assert geo_ids[0] == expected[0]
    if geo_ids[0] != -1:
        assert angles[0] == expected[1]
        assert angles[1] == expected[2]


def test_flat_trough():
    flat_type_params = {
        'te_xyz':[(0.5, 0.5, 0.5)],
        'lag_height': 0.0,
        'lag_facies': 20,
        'structure': 'flat',
        'facies': [0],
        'azimuth': [10, 10],
        'dip': [15, 15],
        'ae_id':0,
        'bg_facies':-1,
        'bg_azim': np.nan,
        'bg_dip': np.nan,
        'size_ztrend': None,
    }
    trough_params = (0.2, 0.1, 0.1, 0)

    # above/upper half of ellipsoid -> expected is -1 (not inside)
    trough_test(flat_type_params, trough_params, (0.5, 0.5, 0.51),
                (-1, 0, 0))
    # inside bounding box, but outside of ellipsoid -> expected is -1
    trough_test(flat_type_params, trough_params, (0.65, 0.59, 0.49),
                (-1, 0, 0))
    # inside -> expected is 0
    trough_test(flat_type_params, trough_params, (0.6, 0.55, 0.49),
                (0, 10+trough_params[-1], 15))

    # do the same thing with alpha = 90° and a and b switched (leading to the
    # same trough)
    trough_params = (0.1, 0.2, 0.1, 90)

    # above/upper half of ellipsoid -> expected is -1 (not inside)
    trough_test(flat_type_params, trough_params, (0.5, 0.5, 0.51),
                (-1, 0, 0))
    # inside bounding box, but outside of ellipsoid -> expected is -1
    trough_test(flat_type_params, trough_params, (0.65, 0.59, 0.49),
                (-1, 0, 0))
    # inside -> expected is 0
    trough_test(flat_type_params, trough_params, (0.6, 0.55, 0.49),
                (0, 10, 15))


# def test_dip_trough():

#     grid = Grid(0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 1.0, 1.0, 1.0, 0)

#     dip = 10
#     azim = 10
#     type_params = {
#             "structure":'dip',
#             "dipset_dist":0.01,
#             "facies":[5, 2],
#             "altfacies":[[2],[5]],
#             "dip":[dip,dip],
#             "azimuth":[azim,azim],
#             'lag_height':0.0,
#             'lag_facies':-1,
#     }

#     a = 1.0
#     b = 0.5
#     c = 0.1
#     alpha = 30.0
#     trough = Trough(type_params, x=0.0, y=0.0, z=0.0, a=a, b=b, c=c, alpha=alpha)
#     trough.num_ha = 0
#     facies = np.zeros(1, dtype=np.int32)
#     angles = np.zeros(2)
#     ids = np.zeros(3, dtype=np.int32)

#     trough.maybe_assign_facies_azim_dip(facies, angles, ids, -0.2, 0.1, -0.02, 0, 0, grid)
#     fac1 = facies[0]

#     sin_dip = sin(dip*np.pi/180)
#     cos_dip = cos(dip*np.pi/180)
#     sin_azim = sin(azim*np.pi/180)
#     cos_azim = cos(azim*np.pi/180)
#     normvec_x = -sin_dip*cos_azim
#     normvec_y = sin_dip*sin_azim
#     normvec_z = cos_dip
#     x = -0.2+0.01*normvec_x
#     y = -0.2+0.01*normvec_y
#     z = -0.2+0.01*normvec_z

#     trough.maybe_assign_facies_azim_dip(facies, angles, ids, x, y, z, 0, 0, grid)
#     fac2 = facies[0]

#     assert fac1 != fac2
#     assert angles[0] == alpha + 10
#     assert angles[1] == 10


# def test_bulb_trough():

#     grid = Grid(0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 1.0, 1.0, 1.0, 0)

#     dip = 90
#     azim = 0
#     type_params = {
#             "structure":'bulb',
#             "dipset_dist":0.01,
#             "facies":[5, 2],
#             "altfacies":[[2],[5]],
#             "dip":[dip,dip],
#             "azimuth":[azim,azim],
#             'lag_height':0.0,
#             'lag_facies':-1,
#     }

#     a = 1.0
#     b = 0.5
#     c = 0.1
#     alpha = 0.0
#     trough = Trough(type_params, x=0.0, y=0.0, z=0.0, a=a, b=b, c=c, alpha=alpha)
#     trough.num_ha = 0
#     facies = np.zeros(1, dtype=np.int32)
#     angles = np.zeros(2)
#     ids = np.zeros(3, dtype=np.int32)


#     x = 0.0
#     y = 0.0
#     z = 0.0
#     trough.maybe_assign_facies_azim_dip(facies, angles, ids, x, y, z, 0, 0, grid)
#     assert facies[0] == 2 or facies[0] == 5
#     assert angles[1] == 0

#     x = 0.0
#     y = 0.0
#     z = -c
#     trough.maybe_assign_facies_azim_dip(facies, angles, ids, x, y, z, 0, 0, grid)
#     assert facies[0] == 2 or facies[0] == 5
#     assert angles[1] == 0

#     x = a
#     y = 0.0
#     z = 0.0
#     trough.maybe_assign_facies_azim_dip(facies, angles, ids, x, y, z, 0, 0, grid)
#     assert facies[0] == 2 or facies[0] == 5
#     assert angles[1] == 90

#     x = a/2
#     y = 0.0
#     z = 0.0
#     trough.maybe_assign_facies_azim_dip(facies, angles, ids, x, y, z, 0, 0, grid)
#     assert facies[0] == 2 or facies[0] == 5
#     assert angles[0] == -90 or angles[0] == 90


# def test_bulb_sets_trough():

#     grid = Grid(0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 1.0, 1.0, 1.0, 0)

#     dip = 90
#     azim = 0
#     type_params = {
#             "structure":'bulb_sets',
#             "bulbset_dist":0.01,
#             "facies":[5, 2],
#             "altfacies":[[2],[5]],
#             "dip":[dip,dip],
#             "azimuth":[azim,azim],
#             'lag_height':0.0,
#             'lag_facies':-1,
#     }

#     a = 1.0
#     b = 0.5
#     c = 0.1
#     alpha = 0.0
#     trough = Trough(type_params, x=0.0, y=0.0, z=0.0, a=a, b=b, c=c, alpha=alpha)
#     trough.num_ha = 0
#     facies = np.zeros(1, dtype=np.int32)
#     angles = np.zeros(2)
#     ids = np.zeros(3, dtype=np.int32)

#     x = 0.0
#     y = 0.0
#     z = -c+0.01
#     trough.maybe_assign_facies_azim_dip(facies, angles, ids, x, y, z, 0, 0, grid)
#     fac1 = facies[0]

#     x = 0.0
#     y = 0.0
#     z = -c+2*0.01
#     trough.maybe_assign_facies_azim_dip(facies, angles, ids, x, y, z, 0, 0, grid)
#     fac2 = facies[0]

#     assert fac2 != fac1


# def test_lag_trough():

#     grid = Grid(0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 1.0, 1.0, 1.0, 0)

#     dip = 0
#     azim = 0
#     type_params = {
#             "structure":'flat',
#             "facies":[5],
#             "dip":[dip,dip],
#             "azimuth":[azim,azim],
#             'lag_height':0.05,
#             'lag_facies':2,
#     }

#     a = 1.0
#     b = 0.5
#     c = 0.1
#     alpha = 0.0
#     trough = Trough(type_params, x=0.0, y=0.0, z=0.0, a=a, b=b, c=c, alpha=alpha)
#     trough.num_ha = 0
#     facies = np.zeros(1, dtype=np.int32)
#     angles = np.zeros(2)
#     ids = np.zeros(3, dtype=np.int32)

#     trough.maybe_assign_facies_azim_dip(facies, angles, ids, 0.0, 0.0, -0.06, 0, 0, grid)
#     assert facies[0] == 2
#     trough.maybe_assign_facies_azim_dip(facies, angles, ids, 0.0, 0.0, -0.04, 0, 0, grid)
#     assert facies[0] == 5

