from hyvr.classes.grid import Grid
from hyvr.classes.trough import Trough 
import numpy as np
from math import sin, cos


def test_flat_trough():

    grid = Grid(0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 1.0, 1.0, 1.0, 0)

    type_params = {
            "structure":'flat',
            "facies":[5],
            "dip":[10,10],
            "azimuth":[10,10],
            'lag_height':0.0,
            'lag_facies':-1,
    }

    a = 1.0
    b = 0.5
    c = 0.1
    alpha = 30.0
    trough = Trough(type_params, x=0.0, y=0.0, z=0.0, a=a, b=b, c=c, alpha=alpha)
    trough.num_ha = 0
    facies = np.zeros(1, dtype=np.int32)
    angles = np.zeros(2)
    ids = np.zeros(3, dtype=np.int32)

    # First test point: above
    trough.maybe_assign_facies_azim_dip(facies, angles, ids, 0.0, 0.0, 1.0, 0, 0, grid)
    assert facies[0] == -1

    # Second point: inside bounding box but not inside trough
    trough.maybe_assign_facies_azim_dip(facies, angles, ids, -0.5, 0.5, 0.0, 0, 0, grid)
    assert facies[0] == -1

    # Third point: inside trough
    trough.maybe_assign_facies_azim_dip(facies, angles, ids, -0.2, 0.1, -0.02, 0, 0, grid)
    assert facies[0] == 5
    assert angles[0] == 10
    assert angles[1] == 10


def test_dip_trough():

    grid = Grid(0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 1.0, 1.0, 1.0, 0)

    dip = 10
    azim = 10
    type_params = {
            "structure":'dip',
            "dipset_dist":0.01,
            "facies":[5, 2],
            "altfacies":[[2],[5]],
            "dip":[dip,dip],
            "azimuth":[azim,azim],
            'lag_height':0.0,
            'lag_facies':-1,
    }

    a = 1.0
    b = 0.5
    c = 0.1
    alpha = 30.0
    trough = Trough(type_params, x=0.0, y=0.0, z=0.0, a=a, b=b, c=c, alpha=alpha)
    trough.num_ha = 0
    facies = np.zeros(1, dtype=np.int32)
    angles = np.zeros(2)
    ids = np.zeros(3, dtype=np.int32)

    trough.maybe_assign_facies_azim_dip(facies, angles, ids, -0.2, 0.1, -0.02, 0, 0, grid)
    fac1 = facies[0]

    sin_dip = sin(dip*np.pi/180)
    cos_dip = cos(dip*np.pi/180)
    sin_azim = sin(azim*np.pi/180)
    cos_azim = cos(azim*np.pi/180)
    normvec_x = -sin_dip*cos_azim
    normvec_y = sin_dip*sin_azim
    normvec_z = cos_dip
    x = -0.2+0.01*normvec_x
    y = -0.2+0.01*normvec_y
    z = -0.2+0.01*normvec_z

    trough.maybe_assign_facies_azim_dip(facies, angles, ids, x, y, z, 0, 0, grid)
    fac2 = facies[0]

    assert fac1 != fac2
    assert angles[0] == alpha + 10
    assert angles[1] == 10


def test_bulb_trough():

    grid = Grid(0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 1.0, 1.0, 1.0, 0)

    dip = 90
    azim = 0
    type_params = {
            "structure":'bulb',
            "dipset_dist":0.01,
            "facies":[5, 2],
            "altfacies":[[2],[5]],
            "dip":[dip,dip],
            "azimuth":[azim,azim],
            'lag_height':0.0,
            'lag_facies':-1,
    }

    a = 1.0
    b = 0.5
    c = 0.1
    alpha = 0.0
    trough = Trough(type_params, x=0.0, y=0.0, z=0.0, a=a, b=b, c=c, alpha=alpha)
    trough.num_ha = 0
    facies = np.zeros(1, dtype=np.int32)
    angles = np.zeros(2)
    ids = np.zeros(3, dtype=np.int32)


    x = 0.0
    y = 0.0
    z = 0.0
    trough.maybe_assign_facies_azim_dip(facies, angles, ids, x, y, z, 0, 0, grid)
    assert facies[0] == 2 or facies[0] == 5
    assert angles[1] == 0

    x = 0.0
    y = 0.0
    z = -c
    trough.maybe_assign_facies_azim_dip(facies, angles, ids, x, y, z, 0, 0, grid)
    assert facies[0] == 2 or facies[0] == 5
    assert angles[1] == 0

    x = a
    y = 0.0
    z = 0.0
    trough.maybe_assign_facies_azim_dip(facies, angles, ids, x, y, z, 0, 0, grid)
    assert facies[0] == 2 or facies[0] == 5
    assert angles[1] == 90

    x = a/2
    y = 0.0
    z = 0.0
    trough.maybe_assign_facies_azim_dip(facies, angles, ids, x, y, z, 0, 0, grid)
    assert facies[0] == 2 or facies[0] == 5
    assert angles[0] == -90 or angles[0] == 90


def test_bulb_sets_trough():

    grid = Grid(0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 1.0, 1.0, 1.0, 0)

    dip = 90
    azim = 0
    type_params = {
            "structure":'bulb_sets',
            "bulbset_dist":0.01,
            "facies":[5, 2],
            "altfacies":[[2],[5]],
            "dip":[dip,dip],
            "azimuth":[azim,azim],
            'lag_height':0.0,
            'lag_facies':-1,
    }

    a = 1.0
    b = 0.5
    c = 0.1
    alpha = 0.0
    trough = Trough(type_params, x=0.0, y=0.0, z=0.0, a=a, b=b, c=c, alpha=alpha)
    trough.num_ha = 0
    facies = np.zeros(1, dtype=np.int32)
    angles = np.zeros(2)
    ids = np.zeros(3, dtype=np.int32)

    x = 0.0
    y = 0.0
    z = -c+0.01
    trough.maybe_assign_facies_azim_dip(facies, angles, ids, x, y, z, 0, 0, grid)
    fac1 = facies[0]

    x = 0.0
    y = 0.0
    z = -c+2*0.01
    trough.maybe_assign_facies_azim_dip(facies, angles, ids, x, y, z, 0, 0, grid)
    fac2 = facies[0]

    assert fac2 != fac1


def test_lag_trough():

    grid = Grid(0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 1.0, 1.0, 1.0, 0)

    dip = 0
    azim = 0
    type_params = {
            "structure":'flat',
            "facies":[5],
            "dip":[dip,dip],
            "azimuth":[azim,azim],
            'lag_height':0.05,
            'lag_facies':2,
    }

    a = 1.0
    b = 0.5
    c = 0.1
    alpha = 0.0
    trough = Trough(type_params, x=0.0, y=0.0, z=0.0, a=a, b=b, c=c, alpha=alpha)
    trough.num_ha = 0
    facies = np.zeros(1, dtype=np.int32)
    angles = np.zeros(2)
    ids = np.zeros(3, dtype=np.int32)

    trough.maybe_assign_facies_azim_dip(facies, angles, ids, 0.0, 0.0, -0.06, 0, 0, grid)
    assert facies[0] == 2
    trough.maybe_assign_facies_azim_dip(facies, angles, ids, 0.0, 0.0, -0.04, 0, 0, grid)
    assert facies[0] == 5

