from hyvr.optimized import assign_between
import numpy as np
import math


def test_assign_between():

    nx, ny, nz = 100, 100, 10
    data = np.zeros((nx, ny, nz), dtype=np.int16)

    z0 = 0
    dz = 0.1

    cs_z_bottom = np.zeros((nx, ny))
    cs_z_top = np.ones((nx, ny)) * 0.2

    assign_between(data, cs_z_bottom, cs_z_top, 1, z0, dz)

    assert np.all(data[:,:,0:2] == 1)
    assert np.all(data[:,:,2:] == 0)

    cs_z_bottom = cs_z_top
    cs_z_top = np.ones((nx, ny)) * 0.4

    assign_between(data, cs_z_bottom, cs_z_top, 2, z0, dz)

    assert np.all(data[:,:,0:2] == 1)
    assert np.all(data[:,:,2:4] == 2)
    assert np.all(data[:,:,4:] == 0)

if __name__ == '__main__':
    test_assign_between()
