# cython: profile=True
import numpy as np
from hyvr.classes.grid cimport Grid
from hyvr.classes.contact_surface cimport ContactSurface
from hyvr.classes.geometrical_object cimport GeometricalObject
import hyvr.utils as hu
from libc.math cimport sin, cos, ceil, acos, sqrt
cimport cython
cimport numpy as np


cdef class Sheet(GeometricalObject):
    """
    This class holds parameters for single sheet objects.
    """

    cdef public:
        int dipsets, num_facies
        # double zmin, zmax
        double dip, azim
        double shift, layer_dist
        double normvec_x, normvec_y, normvec_z
        int facies, num_ha
        int [:] facies_array
        ContactSurface bottom_surface, top_surface


    def __init__(self, type_params, bottom_surface, top_surface, grid):
        cdef double sin_dip, cos_dip, sin_azim, cos_azim, xc, yc, zc

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

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    cpdef maybe_assign_facies_azim_dip(self, np.int32_t [:] facies, np.float_t [:] angles, np.int32_t [:] ids,
                                       double x, double y, double z,
                                       int x_idx, int y_idx, Grid grid):
        """
        This function checks whether the current grid cell with given
        coordinates is inside the sheet and assigns facies, azimuth and dip by
        altering the passed arrays.

        Parameters
        ----------
        facies : int array
            array of size 1 that holds the facies. This will be altered during
            the course of this function to either be -1 if the cell is not in
            the trough, or it is set to the facies number of the cell.
        angles : double array
            array of size 2 that holds the azimuth and dip of the cell. This
            will also be altered.
        ids : np.ndarray[np.int32, dim=1]
            Array of length three that holds AERealization-ID, ha, and hat
        x, y, z : double
            cell coordinates
        x_idx, y_idx : int
            indices of x and y position in grid.
        grid : Grid object
        """
        cdef double z_above, z_below, d
        cdef int n
        # if we reach this functions, we now that we are within [zmin, zmax]
        # the only way how this could be outside the domain is when the top or
        # bottom surfaces are not flat
        z_above = self.top_surface.surface[x_idx, y_idx]
        z_below = self.bottom_surface.surface[x_idx, y_idx]
        if z < z_below or z > z_above:
            facies[0] = -1
            return

        # at this point we now that the point is inside the sheet
        angles[0] = self.azim
        angles[1] = self.dip
        ids[1] = self.num_ha
        if self.dipsets:
            d = self.normvec_x*x + self.normvec_y*y + self.normvec_z*z - self.shift
            n = int(d/self.layer_dist) + self.num_facies//2
            facies[0] = self.facies_array[n]
            return
        else:
            facies[0] = self.facies
            return
