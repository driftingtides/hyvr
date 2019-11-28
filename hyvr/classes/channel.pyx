# cython: profile=True
import numpy as np
cimport numpy as np
cimport cython
from hyvr.classes.grid cimport Grid
from hyvr.classes.geometrical_object cimport GeometricalObject
from libc.math cimport sqrt, sin, cos, atan2
import hyvr.utils as hu

cdef class Channel(GeometricalObject):
    """
    Channel implementation.

    This is work in progress, as we currently lack a centerline construction method.
    """

    cdef public:
        double[:] x_center, y_center, vx, vy
        int[:,:] dont_check
        double a, width, depth, min_dx_dy
        int len_centerline
        # double zmin, zmax
        double ztop
        double dip, azim, sin_dip, cos_dip
        double shift, layer_dist
        double lag_height
        int facies, num_ha, lag_facies, num_facies
        int [:] facies_array
        int dipsets


    def __init__(self, type_params, x_center, y_center, vx, vy, z, width, depth, grid):
        # x_center and y_center should be coordinates of points on the centerline.
        cdef int i
        cdef double dist_along_curve
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


    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    cpdef maybe_assign_facies_azim_dip(self, np.int32_t [:] facies, np.float_t [:] angles, np.int32_t [:] ids,
                                     double x, double y, double z,
                                     int x_idx, int y_idx, Grid grid):
        """
        This function checks whether the current grid cell with given
        coordinates is inside the trough and assigns facies, azimuth and dip by
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
        cdef double xy_dist, dz, radius, dist, weight, sum_weights, vx_now, vy_now
        cdef double dist_along_curve, dist_along_curve_tmp
        cdef int nx, ny, i, j, n, closest_idx


        # if the point is above the channel top, don't consider it
        dz = z - self.ztop
        if dz > 0:
            facies[0] = -1
            return

        # To check whether a point is inside the channel, we have to calculate
        # its distance to the centerline. This is quite costly, as it means we
        # have to iterate over the list of centerline points.
        # To avoid checks, after each distance calculation of a cell that does
        # belong to the channel, we try to find the indices of all points that
        # are closer to the current point than ``distance - channel width``.
        # These will then be set to 1 in ``dont_check``, so we can skip the
        # check for these.

        if self.dont_check[x_idx, y_idx]:
            facies[0] = -1
            return


        # Otherwise, get the distance to the centerline in xy-plane and z-direction
        # If we have dipsets, we will also directly calculate the velocity
        # (inverse distance weighting) and the distance along the curve
        sum_weights = 0
        vx_now = 0
        vy_now = 0
        dist_along_curve_tmp = 0
        dist_along_curve = 0
        xy_dist = 1e100
        if self.dipsets:
            for i in range(self.len_centerline):
                if i > 0:
                    dist_along_curve_tmp += sqrt((self.x_center[i] - self.x_center[i-1])**2
                                                +(self.y_center[i] - self.y_center[i-1])**2
                    )
                dist = sqrt((x - self.x_center[i])**2 + (y - self.y_center[i])**2)
                weight = 1/(dist + 1e-20)
                vx_now += self.vx[i] * weight
                vy_now += self.vy[i] * weight
                sum_weights += weight
                if dist < xy_dist:
                    dist_along_curve = dist_along_curve_tmp
                    xy_dist = dist
                    closest_idx = i
            vx_now /= sum_weights
            vy_now /= sum_weights
        else:
            # we only need the distance
            for i in range(self.len_centerline):
                dist = (x - self.x_center[i])**2 + (y - self.y_center[i])**2
                if dist < xy_dist:
                    xy_dist = dist
            xy_dist = sqrt(xy_dist)


        # Find other cells that are also not close enough
        radius = xy_dist - self.width
        if radius > self.min_dx_dy:
            # Idea: find number of cells in x-direction, then loop over these
            # and find y-directions for those.
            nx = int(radius/grid.dx)
            for i in range(nx):
                if x_idx-i >= 0 and x_idx+i < grid.nx:
                    # 'height' of the circle at distance x = i*dx from the center:
                    #
                    #     y(x) = sqrt(radius**2 - x**2)
                    #
                    ny = int(sqrt(radius**2 - (i*grid.dx)**2))
                    for j in range(ny):
                        if y_idx-j >= 0 and y_idx+j < grid.ny:
                            self.dont_check[x_idx+i, y_idx+j] = 1
                            self.dont_check[x_idx+i, y_idx-j] = 1
                            self.dont_check[x_idx-i, y_idx+j] = 1
                            self.dont_check[x_idx-i, y_idx-j] = 1

        # Now we calculate the value of the shape curve and check whether the
        # point is really inside (i.e. is above the parabola)
        if dz >= self.a*xy_dist**2 - self.depth:
            # it's inside: assign stuff
            ids[1] = self.num_ha
            if z < self.ztop - self.depth + self.lag_height:
                angles[0] = 0
                angles[1] = 0
                facies[0] = self.lag_facies
                return

            if self.dipsets:
                # azimuth from inverse distance weighted velocity
                angles[0] = atan2(vy_now, vx_now)/np.pi*180
                # dip as assigned
                angles[1] = self.dip

                # TODO: the distance along the curve is actually a bit
                # different, as the point is not directly orthogonal to the closest point.
                # However, if the point density on the centerline is high
                # enough, it won't really matter.
                # To correct for the distance in z-direction we subtract |dz| * cos(dip)
                # note that dz is negative here
                d = dist_along_curve * self.sin_dip + dz*self.cos_dip + self.shift
                n = int(d/self.layer_dist)
                facies[0] = self.facies_array[n]
                return
            else:
                angles[0] = self.azim
                angles[1] = self.dip
                facies[0] = self.facies
                return

        else:
            # it's not inside. This means lower points will also not be inside
            self.dont_check[x_idx, y_idx] = 1
            facies[0] = -1
            return



