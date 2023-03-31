import numpy as np
import numpy.typing as npt
from hyvr.geo.grid import Grid
#from hyvr.geo.contact_surface import ContactSurface
import hyvr.utils as hu
#from libc.math cimport sin, cos, ceil, acos, sqrt


class Sheet:
    """
    This class holds parameters for single sheet objects.
    """


    def __init__(self, type_params, bottom_surface, ztop, grid, num_ha):

        self.bottom_surface = bottom_surface
        #self.top_surface = top_surface
        self.zmin = np.min(self.bottom_surface)
        self.zmax = ztop
        self.num_ha = num_ha
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
            sin_dip = np.sin(self.dip*np.pi/180)
            cos_dip = np.cos(self.dip*np.pi/180)
            sin_azim = np.sin(self.azim*np.pi/180)
            cos_azim = np.cos(self.azim*np.pi/180)
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
            self.num_facies = 1
            self.facies_array = hu.get_alternating_facies(self.num_facies, type_params)
