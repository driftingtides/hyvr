import numpy as np
from hyvr.geo.grid cimport Grid
from hyvr.geo.contact_surface cimport ContactSurface
import hyvr.utils as hu
from libc.math cimport sin, cos, ceil, acos, sqrt
cimport cython
cimport numpy as np


cdef class Sheet:
    """
    This class holds parameters for single sheet objects.
    """

    cdef public:
        np.int32_t dipsets, num_facies
        np.float_t zmin, zmax
        np.float_t dip, azim
        np.float_t shift, layer_dist
        np.float_t normvec_x, normvec_y, normvec_z
        np.int32_t facies, num_ha
        np.int32_t [:] facies_array
        ContactSurface bottom_surface, top_surface


    def __init__(self, type_params, bottom_surface, top_surface, grid):
        cdef np.float_t sin_dip, cos_dip, sin_azim, cos_azim, xc, yc, zc

        self.bottom_surface = bottom_surface
        self.top_surface = top_surface
        self.zmin = self.bottom_surface.zmin
        self.zmax = self.top_surface.zmax
        self.dipsets = type_params['structure'] == 'dip'
        if self.dipsets:
            # the dipsets are layers of planes with distance dipset_dist such that
            # the 'zero plane' goes through the center of the bottom surface of
            # the sheet. The respective facies will be stored in an array such
            # that the zero-plane is in the center.
            self.layer_dist = type_params['dipset_dist']

            # The plane is defined by a normal vector:
            self.azim = np.random.uniform(*type_params['azimuth'])
            self.dip = np.random.uniform(*type_params['dip'])
            sin_dip = sin(self.dip*np.pi/180)
            cos_dip = cos(self.dip*np.pi/180)
            sin_azim = sin(self.azim*np.pi/180)
            cos_azim = cos(self.azim*np.pi/180)
            self.normvec_x = -sin_dip*cos_azim
            self.normvec_y = sin_dip*sin_azim
            self.normvec_z = cos_dip


            # The distance of a point to the zero-plane is:
            #
            # d(x, y, z) = nx*x + ny*y + nz*z - d
            #
            # and the distance of the center point is zero:
            #
            # d(xc, yc, zc) = nx*xc + ny*yc + nz*zc - d = 0
            #
            # or
            #
            # d = nx*xc + ny*y + nz*z
            xc = grid.x0 + grid.lx/2
            yc = grid.y0 + grid.ly/2
            # zc = grid.z0 + grid.lz/2
            zc = self.bottom_surface.zmean
            self.shift = self.normvec_x * xc + self.normvec_y * yc + self.normvec_z * zc
            self.shift += np.random.rand() * self.layer_dist


            # The maximum amount of necessary layers is max(lx, ly, lz)/layer_dist
            # 10 is added just to be sure
            self.num_facies = int(np.ceil((max(grid.ly, grid.ly, grid.lz) + self.shift)/self.layer_dist)) + 10
            self.facies_array = hu.get_alternating_facies(self.num_facies, type_params)
        else:
            # massive internal structure
            self.azim = np.random.uniform(*type_params['azimuth'])
            self.dip = np.random.uniform(*type_params['dip'])
            self.facies = np.random.choice(type_params['facies'])
