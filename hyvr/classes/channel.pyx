import numpy as np
cimport numpy as np
cimport cython
from hyvr.classes.grid cimport Grid
from libc.math cimport sqrt, sin, cos
import hyvr.utils as hu

cdef class Channel:
    """
    Channel implementation.

    This is work in progress, as we currently lack a centerline construction method.
    """

    cdef public:
        double[:] x_center, y_center
        int[:,:] dont_check
        double a, width, depth, min_dx_dy
        int len_centerline
        double zmin, zmax, ztop
        double dip, azim
        double shift, layer_dist
        double normvec_x, normvec_y, normvec_z
        int facies, num_ha
        int [:] facies_array


    def __init__(self, type_params, x_center, y_center, z, width, depth, grid):
        # x_center and y_center should be coordinates of points on the centerline.
        self.x_center = x_center
        self.y_center = y_center
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


        # Get facies, dip, and azimuth
        self.dipsets = type_params['structure'] == 'dip'
        if self.dipsets:
            # the dipsets are layers of planes with distance dipset_dist such that
            # the 'zero plane' goes through the center of the domain (at the
            # top of the channel).
            # The respective facies will be stored in an array such
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
            zc = self.ztop
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
    cpdef maybe_assign_facies_azim_dip(self, int [:] facies, double [:] angles, int[:] ids,
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
        cdef double xy_dist, dz, radius, dist_sq
        cdef int nx, ny, i, j


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
        xy_dist = 1e100
        for i in range(self.len_centerline):
            dist_sq = (x - self.x_center[i])**2 + (y - self.y_center[i])**2
            if dist_sq < xy_dist:
                xy_dist = dist_sq
        xy_dist = sqrt(xy_dist)

        # Find other cells that are also not close enough
        radius = xy_dist - self.width
        if radius > self.min_dx_dy:
            # Idea: find number of cells in x-direction, then loop over these
            # and find y-directions for those.
            nx = int(radius/grid.dx)
            for i in range(nx):
                # 'height' of the circle at distance x = i*dx from the center:
                #
                #     y(x) = sqrt(radius**2 - x**2)
                #
                ny = int(sqrt(radius**2 - (i*grid.dx)**2))
                for j in range(ny):
                    self.dont_check[x_idx+i, y_idx+j] = 1
                    self.dont_check[x_idx+i, y_idx-j] = 1
                    self.dont_check[x_idx-i, y_idx+j] = 1
                    self.dont_check[x_idx-i, y_idx-j] = 1

        # Now we calculate the value of the shape curve and check whether the
        # point is really inside (i.e. is above the parabola)
        if dz >= self.a*xy_dist**2 - self.depth:
            # it's inside: assign stuff
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

        else:
            # it's not inside. This means lower points will also not be inside
            self.dont_check[x_idx, y_idx] = 1
            facies[0] = -1
            return



