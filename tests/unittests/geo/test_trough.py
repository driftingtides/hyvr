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

grid = Grid(0., 0., 0., 0.1, 0.1, 0.1, 1.0, 1.0, 1.0, 0)

def setup_trough(type_params, trough_params):
    """
    Setup a trough AE.
    
    Parameters
    ----------
    type_params : dict
        Dictionary of type params
    trough_params : list/tuple
        (a, b, c, alpha)
    
    Returns
    -------
    trough : trough AE
    """
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
    return trough_ae

def trough_get_values_at_points(trough, point):
    """
    Run a test for given type params, trough params (a, b, c, alpha), a given
    point, and assert that the expected values are set
    
    Parameters
    ----------
    trough : trough AE object
    point : list/tuple
        (x, y, z)
    
    Returns
    -------
    facies, azim, dip : float
        Facies, azimuth and dip at given location
    """
    # setup containers
    geo_ids = np.zeros(3, dtype=np.int32)
    angles = np.zeros(2, dtype=np.float)

    # test
    x, y, z = point
    trough.maybe_assign_points_to_object(0, geo_ids, angles, x, y, z, 0, 0, grid)
    return geo_ids[0], angles[0], angles[1]

def trough_test_expected(trough, point, expected):
    """
    Run a test for given type params, trough params (a, b, c, alpha), a given
    point, and assert that the expected values are set
    
    Parameters
    ----------
    trough : trough AE object
    point : list/tuple
        (x, y, z)
    expected : list/tuple
        (expected_facies, expected_azim, expected_dip)
    """
    # setup containers
    geo_ids = np.zeros(3, dtype=np.int32)
    angles = np.zeros(2, dtype=np.float)

    # test
    x, y, z = point
    trough.maybe_assign_points_to_object(0, geo_ids, angles, x, y, z, 0, 0, grid)
    assert geo_ids[0] == expected[0]
    if geo_ids[0] != -1:
        assert angles[0] == expected[1]
        assert angles[1] == expected[2]



@pytest.mark.unit
def test_massive_trough():
    alpha = 10
    azim = 10
    dip = 15
    flat_type_params = {
        'te_xyz':[(0.5, 0.5, 0.5)],
        'lag_height': 0.0,
        'lag_facies': 20,
        'structure': 'massive',
        'facies': [0],
        'azimuth': [azim, azim],
        'dip': [dip, dip],
        'ae_id':0,
        'bg_facies':-1,
        'bg_azim': np.nan,
        'bg_dip': np.nan,
        'size_ztrend': None,
    }
    trough_params = (0.2, 0.1, 0.1, alpha)
    trough = setup_trough(flat_type_params, trough_params)

    # above/upper half of ellipsoid -> expected is -1 (not inside)
    trough_test_expected(trough, (0.5, 0.5, 0.51), (-1, 0, 0))
    # inside bounding box, but outside of ellipsoid -> expected is -1
    trough_test_expected(trough, (0.65, 0.59, 0.49), (-1, 0, 0))
    # inside -> expected is 0
    trough_test_expected(trough, (0.6, 0.55, 0.49), (0, azim+alpha, dip))

    # do the same thing with alpha = 90° and a and b switched (leading to the
    # same trough)
    alpha = 90
    trough_params = (0.1, 0.2, 0.1, alpha)
    trough = setup_trough(flat_type_params, trough_params)

    # above/upper half of ellipsoid -> expected is -1 (not inside)
    trough_test_expected(trough, (0.5, 0.5, 0.51), (-1, 0, 0))
    # inside bounding box, but outside of ellipsoid -> expected is -1
    trough_test_expected(trough, (0.65, 0.59, 0.49), (-1, 0, 0))
    # inside -> expected is 0
    trough_test_expected(trough, (0.6, 0.55, 0.49), (0, alpha+azim, dip))


@pytest.mark.unit
def test_dip_trough():

    dip = 10
    azim = 10
    alpha = 30
    dist = 0.01
    dip_type_params = {
        'te_xyz':[(0.5, 0.5, 0.5)],
        "structure":'dip',
        "dipset_dist":dist,
        "facies":[5, 2],
        "altfacies":[[2],[5]],
        "dip":[dip,dip],
        "azimuth":[azim,azim],
        'lag_height':0.0,
        'lag_facies':-1,
        'ae_id':0,
        'bg_facies':-1,
        'bg_azim': np.nan,
        'bg_dip': np.nan,
        'size_ztrend': None,
    }

    trough_params = (0.2, 0.2, 0.2, alpha)
    trough = setup_trough(dip_type_params, trough_params)

    # get first facies
    x0, y0, z0 = 0.5+0.1, 0.5+0.1, 0.5-0.1
    fac1, azim1, dip1 = trough_get_values_at_points(trough, (x0, y0, z0))
    assert dip1 == dip
    assert azim1 == azim + alpha

    # get second facies: first, find shift vector
    sin_dip = sin(dip*np.pi/180)
    cos_dip = cos(dip*np.pi/180)
    sin_azim = sin((azim+alpha)*np.pi/180)
    cos_azim = cos((azim+alpha)*np.pi/180)
    normvec_x = sin_dip*cos_azim
    normvec_y = -sin_dip*sin_azim
    normvec_z = cos_dip
    x = x0+dist*normvec_x
    y = y0+dist*normvec_y
    z = z0+dist*normvec_z
    fac2, azim2, dip2 = trough_get_values_at_points(trough, (x, y, z))
    assert dip2 == dip
    assert azim2 == azim + alpha

    assert fac1 != fac2
    assert fac1 in dip_type_params['facies']
    assert fac2 in dip_type_params['facies']

@pytest.mark.unit
def test_bulb_trough():

    dip = 90 # maximum possible dip
    azim = 0
    dist = 0.01
    type_params = {
        'te_xyz':[(0.5, 0.5, 0.5)],
        "structure":'bulb',
        "dipset_dist":dist,
        "facies":[5, 2],
        "altfacies":[[2],[5]],
        "dip":[dip,dip],
        "azimuth":[azim,azim],
        'lag_height':0.0,
        'lag_facies':-1,
        'ae_id':0,
        'bg_facies':-1,
        'bg_azim': np.nan,
        'bg_dip': np.nan,
        'size_ztrend': None,
    }

    trough_params = (0.2, 0.1, 0.05, 0.0)
    trough = setup_trough(type_params, trough_params)
    x0, y0, z0 = 0.5, 0.5, 0.5

    # in the center of the trough, the dip should be 0
    fac, azim, dip = trough_get_values_at_points(trough, (x0, y0, z0))
    assert fac in [5, 2]
    assert dip == 0

    # also at the bottom
    fac, azim, dip = trough_get_values_at_points(trough, (x0, y0, z0-0.05))
    assert fac in [5, 2]
    assert dip == 0

    # at the side in positive x-direction, the dip should be -90° and azim 0°
    fac, azim, dip = trough_get_values_at_points(trough, (x0+0.2, y0, z0))
    assert fac in [5, 2]
    assert azim == 0
    assert dip == -90

    # at the side in negative x-direction, the dip should be 90° and azim 0°
    fac, azim, dip = trough_get_values_at_points(trough, (x0-0.2, y0, z0))
    assert fac in [5, 2]
    assert azim == 0
    assert dip == 90

    # at the side in positive y-direction, the dip should be +-90° and azim -90°
    fac, azim, dip = trough_get_values_at_points(trough, (x0, y0+0.1, z0))
    assert fac in [5, 2]
    assert azim == -90
    assert abs(dip) == 90

    # at the side in negative y-direction, the dip should be +-90° and azim +90°
    fac, azim, dip = trough_get_values_at_points(trough, (x0, y0-0.1, z0))
    assert fac in [5, 2]
    assert azim == 90
    assert abs(dip) == 90


@pytest.mark.unit
def test_bulb_trough_limited_dip():

    dip = 45 # maximum possible dip
    azim = 0
    dist = 0.01
    type_params = {
        'te_xyz':[(0.5, 0.5, 0.5)],
        "structure":'bulb',
        "dipset_dist":dist,
        "facies":[5, 2],
        "altfacies":[[2],[5]],
        "dip":[dip,dip],
        "azimuth":[azim,azim],
        'lag_height':0.0,
        'lag_facies':-1,
        'ae_id':0,
        'bg_facies':-1,
        'bg_azim': np.nan,
        'bg_dip': np.nan,
        'size_ztrend': None,
    }

    trough_params = (0.2, 0.1, 0.05, 0.0)
    trough = setup_trough(type_params, trough_params)
    x0, y0, z0 = 0.5, 0.5, 0.5

    # at the side in x-direction, the dip should be 45° and azim 0°
    fac, azim, dip = trough_get_values_at_points(trough, (x0+0.2, y0, z0))
    assert fac in [5, 2]
    assert azim == 0
    assert abs(dip) == 45


@pytest.mark.unit
def test_bulb_sets_trough():

    dip = 90
    azim = 0
    type_params = {
        'te_xyz':[(0.5, 0.5, 0.5)],
        "structure":'bulb_sets',
        "bulbset_dist":0.01,
        "facies":[5, 2],
        "altfacies":[[2],[5]],
        "dip":[dip,dip],
        "azimuth":[azim,azim],
        'lag_height':0.0,
        'lag_facies':-1,
        'ae_id':0,
        'bg_facies':-1,
        'bg_azim': np.nan,
        'bg_dip': np.nan,
        'size_ztrend': None,
    }

    trough_params = (0.5, 0.2, 0.1, 0.0)
    trough = setup_trough(type_params, trough_params)

    fac1, azim, dip = trough_get_values_at_points(trough, (0.5, 0.5, 0.5-0.09))
    fac2, azim, dip = trough_get_values_at_points(trough, (0.5, 0.5, 0.5-0.08))

    assert fac2 != fac1


@pytest.mark.unit
def test_lag_trough():

    dip = 0
    azim = 0
    type_params = {
        'te_xyz':[(0.5, 0.5, 0.5)],
        "structure":'massive',
        "facies":[5],
        "dip":[dip,dip],
        "azimuth":[azim,azim],
        'lag_height':0.05,
        'lag_facies':2,
        'ae_id':0,
        'bg_facies':-1,
        'bg_azim': np.nan,
        'bg_dip': np.nan,
        'size_ztrend': None,
    }
    trough_params = (0.5, 0.2, 0.1, 0.0)
    trough = setup_trough(type_params, trough_params)

    fac, azim, dip = trough_get_values_at_points(trough, (0.5, 0.5, 0.44))
    assert fac == 2
    fac, azim, dip = trough_get_values_at_points(trough, (0.5, 0.5, 0.46))
    assert fac == 5
