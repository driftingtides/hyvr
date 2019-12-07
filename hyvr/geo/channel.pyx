import numpy as np
cimport numpy as np
cimport cython
from hyvr.geo.grid cimport Grid
from libc.math cimport sqrt, sin, cos, atan2
import hyvr.utils as hu

cdef class Channel:
    """
    Channel implementation.

    This is work in progress, as we currently lack a centerline construction method.
    """

    cdef public:
        np.float_t[:] x_center, y_center, vx, vy
        np.int32_t[:,:] dont_check
        np.float_t a, width, depth, min_dx_dy
        np.int32_t len_centerline
        np.float_t zmin, zmax
        np.float_t ztop
        np.float_t dip, azim, sin_dip, cos_dip
        np.float_t shift, layer_dist
        np.float_t lag_height
        np.int32_t facies, num_ha, lag_facies, num_facies
        np.int32_t [:] facies_array
        np.int32_t dipsets


    def __init__(self, type_params, x_center, y_center, vx, vy, z, width, depth, grid):
        # x_center and y_center should be coordinates of points on the centerline.
        cdef int i
        cdef np.float_t dist_along_curve
        self.x_center = x_center
        self.y_center = y_center
        self.vx = vx
        self.vy = vy
        self.len_centerline = len(x_center)

        # we need this to reduce the amount of distance checks
        self.dont_check = np.zeros((grid.nx, grid.ny), dtype=np.int32)
        self.min_dx_dy = min(grid.dx, grid.dy)

        # for better searching
        self.ztop = z
        self.zmax = z
        self.zmin = z - depth


        # The channel shape is a parabolic function that maps the distance in
        # the x-y-plane onto a z-value:
        #
        #     z(d) = a*d**2 + b*d + c
        #
        # The parameters a, b, and c must be chosen such that ``z(0) = -depth``
        # and ``z(+-width/2) = 0``:
        #
        #     z(0) = c = -depth                                          (1)
        #
        #     z(width/2) = a*width**2/4 + b*width/2 - depth = 0
        #     z(-width/2) = a*width**2/4 - b*width/2 - depth = 0
        #
        # Therefore ``b = 0`` and ``a = 4*depth/width**2``
        self.depth = depth
        self.width = width
        self.a = 4*self.depth/self.width**2

        # lag
        self.lag_height = type_params['lag_height']
        self.lag_facies = type_params['lag_facies']


        # Get facies, dip, and azimuth
        self.dipsets = type_params['structure'] == 'dip'
        if self.dipsets:
            # the dipsets are layers of planes along the channel centerline
            # with distance 'dipset_dist'
            # To find which layer a point belongs to we will calculate its
            # distance along the centerline.
            # Since the layers are dipping, the real distance is longer than
            # the distance that we would go if we would go orthogonal to the
            # planes. Therefore the distance along the curve has to be
            # multiplied by the sin of the dip angle
            self.layer_dist = type_params['dipset_dist']
            self.dip = np.random.uniform(*type_params['dip'])
            self.sin_dip = sin(self.dip*np.pi/180)
            self.cos_dip = cos(self.dip*np.pi/180)
            # a random shift to be added to the distance
            self.shift = np.random.uniform(0, self.layer_dist)

            # The respective facies will be stored in an array such
            # that the zero-plane is in the center.
            # To find out how long the list has to be we have to find the length of the curve
            dist_along_curve = 0
            for i in range(1, self.len_centerline):
                dist_along_curve += sqrt((self.x_center[i] - self.x_center[i-1])**2
                                            +(self.y_center[i] - self.y_center[i-1])**2
                )
            self.num_facies = int(
                np.ceil(dist_along_curve*self.sin_dip/self.layer_dist + self.shift)
                + np.ceil(self.depth/self.layer_dist)
            )
            self.num_facies += 10  # just to be sure
            self.facies_array = hu.get_alternating_facies(self.num_facies, type_params)
        else:
            # massive internal structure
            self.dip = np.random.uniform(*type_params['dip'])
            self.facies = np.random.choice(type_params['facies'])
