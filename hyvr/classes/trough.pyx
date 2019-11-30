import numpy as np
cimport hyvr.optimized as ho
import hyvr.utils as hu
from libc.math cimport sin, cos, ceil, acos, sqrt
cimport numpy as np
from hyvr.classes.grid cimport Grid
cimport cython

cdef class Trough:
    """
    This class holds all parameters of a geometrical trough object.

    It has two methods: a constructor, that sets all parameters, e.g. dip set
    facies, etc., and a ``maybe_assign_facies_azim_dip`` method, that checks
    whether a given cell is inside the trough and returns its properties.
    """
    cdef public:
        double x, y, z
        double a, b, c
        double max_ab
        double zmin, zmax
        double alpha
        double sinalpha, cosalpha
        double lag, shift
        int structure
        double layer_dist
        double dip, azim
        int facies, lag_facies, num_ha
        double cosdip
        double normvec_x, normvec_y, normvec_z
        int num_facies
        int [:] facies_array

    def __init__(self, type_params, *, x, y, z, a, b, c, alpha):
        cdef double sin_dip, cos_dip, sin_azim, cos_azim

        self.x = x
        self.y = y
        self.z = z
        self.a = a
        self.b = b
        self.c = c
        self.max_ab = max(a, b)
        self.zmin = z - c
        self.zmax = z
        self.alpha = alpha

        # To calculate the 'normalized' distance (i.e. the distance normalized
        # by the parameters a, b, and c such that distance > 1 means outside of
        # the ellipsoid) of a rotated trough or to find the tangential planes
        # to calculate dip and azimuth of bulb sets, we could view our
        # ellipsoid as a deformed and rotated sphere. At first, a deformation
        # along the x-, y- and z-axes it performed, and then a rotation along
        # the z-axis with angle alpha.
        # In Matrix form, this operation looks as follows:
        #
        #     / cos -sin  0 \   / a  0  0 \   / a*cos -b*sin  0 \
        # A = | sin  cos  0 | * | 0  b  0 | = | a*sin  b*cos  0 |
        #     \  0    0   1 /   \ 0  0  c /   \   0      0    c /
        #
        # sin and cos are here sin and cos of alpha, the rotation angle
        # This matrix can be used to transform from the normalized space to the
        # real space. To go from the real space to the normalized space, the
        # inverse has to be applied.
        # Since this is a block diagonal matrix, we can just invert both blocks:
        #
        #        /  cos/a  sin/a   0  \
        # Ainv = | -sin/b  cos/b   0  |
        #        \    0       0   1/c /
        #
        #
        # Multiplication with a vector (x0, x1, x2) gives:
        #
        #          /  x0*cos/a + x1*sin/a \
        # Ainv*x = | -x0*sin/b + x1*cos/b |
        #          \         x2/c         /
        #
        #
        self.sinalpha = np.sin(alpha*np.pi/180)
        self.cosalpha = np.cos(alpha*np.pi/180)
        # self.A = np.asarray([[a*self.cosalpha, -b*self.sinalpha, 0],
        #                      [a*self.sinalpha, b*self.cosalpha, 0],
        #                      [0, 0, c]])
        # self.normalize_x = lambda x, y, z: np.array(
        #     [(x*self.cosalpha + y*self.sinalpha)/self.a,
        #      (-x*self.sinalpha + y*self.cosalpha)/self.b,
        #      z/c]
        # )
        # self.normalized_dist_square = lambda x, y, z: (
        #     ((x*self.cosalpha + y*self.sinalpha)/a)**2 +\
        #     ((-x*self.sinalpha + y*self.cosalpha)/b)**2 +\
        #     (z/c)**2
        #     )

        self.lag = type_params['lag_height']
        self.lag_facies = type_params['lag_facies']

        # Assign structure
        #
        # To avoid string comparisons we assign each structure a number as
        # follows:
        #
        # 0 : 'flat'
        # 1 : 'dip'
        # 2 : 'bulb'
        # 3 : 'bulb_sets'
        # -1: 'random'
        self.structure = ['flat', 'dip', 'bulb', 'bulb_sets', 'random'].index(type_params['structure'])
        if self.structure == 4: # random
            self.structure = np.random.randint(4)
        if self.structure == 0:  # flat
            self.facies = np.random.choice(type_params['facies'])
            self.azim = np.random.uniform(*type_params['azimuth'])
            self.dip = np.random.uniform(*type_params['dip'])
        elif self.structure == 1: # dip
            # In this case the internal structure consists of parallel dipping layers.
            # To describe these, we need their thickness, a point where we
            # could place the first layer, and a normal vector defined by two
            # angles, dip and azimuth.
            self.layer_dist = type_params['dipset_dist']

            # The normal vector is determined by the dip angle of the dipping
            # sets and their azimuth. We can again obtain the normal vector by
            # applying rotation matrices to the vector (0, 0, 1) (the upwards
            # pointing unit vector).
            # First, we apply the dip rotation by rotating around by the dip
            # angle along the y-axis to the right (anti-clockwise).
            # Then we rotate around the z-axis by the azimuth angle along the
            # z-axis to the left (clockwise)
            # In matrix notation:
            #
            #     / cos(azim)  sin(azim)  0 \   / cos(dip)  0  -sin(dip) \   / 0 \
            # n = |-sin(azim)  cos(azim)  0 | * |    0      1      0     | * | 0 |
            #     \    0            0     1 /   \ sin(dip)  0   cos(dip) /   \ 1 /
            #
            #
            #     / cos(azim)  sin(azim)  0 \   / -sin(dip) \
            # n = |-sin(azim)  cos(azim)  0 | * |    0      |
            #     \    0            0     1 /   \  cos(dip) /
            #
            #
            #     / -sin(dip)*cos(azim) \
            # n = |  sin(dip)*sin(azim) |
            #     \       cos(dip)      /
            #
            #
            # The azimuth is the given azimuth plus the rotation of the trough,
            # all angles are in degree
            self.azim = self.alpha + np.random.uniform(*type_params['azimuth'])
            self.dip = np.random.uniform(*type_params['dip'])
            sin_dip = sin(self.dip*np.pi/180)
            cos_dip = cos(self.dip*np.pi/180)
            sin_azim = sin(self.azim*np.pi/180)
            cos_azim = cos(self.azim*np.pi/180)
            self.normvec_x = -sin_dip*cos_azim
            self.normvec_y = sin_dip*sin_azim
            self.normvec_z = cos_dip

            # The 'starting point' will be the center of the trough, i.e.
            # (self.x, self.y, self.z).
            # To move the planes a bit we add a random component
            self.shift = np.random.uniform(0, self.layer_dist)

            # The facies changes with every distance `layer_dist` to the zero
            # plane. The maximum number of different facies is:
            self.num_facies = int(np.ceil((2*max(a, b, c) + self.shift)/self.layer_dist)) + 2
            self.facies_array = hu.get_alternating_facies(self.num_facies, type_params)
        elif self.structure == 2: # bulb
            self.facies = np.random.choice(type_params['facies'])
            # maximum possible dip
            self.dip = np.random.uniform(*type_params['dip'])
        elif self.structure == 3: # bulb_sets
            # The facies thickness is 'bulbset_dist', and is measured along the
            # z-axis. If a cell is not on the z-axis, this distance get's
            # deformed a bit. We add again a random shift, but this time in the
            # normalized domain.
            self.layer_dist = type_params['bulbset_dist']
            self.shift = np.random.uniform()
            self.num_facies = int(ceil((1 + self.shift)*c/self.layer_dist)) + 1
            self.facies_array = hu.get_alternating_facies(self.num_facies, type_params)
            # maximum possible dip
            self.dip = np.random.uniform(*type_params['dip'])
        else:
            print("Structure number:", self.structure)
            raise NotImplementedError('This structure is not implemented yet')




    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    cpdef maybe_assign_points(self,
                              np.int32_t [:] geo_ids,
                              np.float_t [:] angles,
                              double x, double y, double z,
                              int x_idx, int y_idx,
                              Grid grid):
        """
        This function checks whether the current grid cell with given
        coordinates is inside the trough and assigns facies, azimuth and dip by
        altering the passed arrays.

        Parameters
        ----------
        geo_ids : np.int32 array of size 4
            This holds the geological indices, i.e. facies number, architectural
            element number, hydrofacies assemblage (ha) and hydrofacies assemblage
            type (hat).
        angles : double array
            array of size 2 that holds the azimuth and dip of the cell. This
            will also be altered.
        x, y, z : double
            cell coordinates
        x_idx, y_idx : int
            indices of x and y position in grid.
        grid : Grid object
        """
        cdef double dx, dy, dz
        cdef double normalized_dx, normalized_dy, normalized_dz
        cdef double l2, plane_dist, aaa, bbb, ccc, len_tanvec
        cdef double tanvec_x, tanvec_y, tanvec_z
        cdef int n
        cdef dip, cos_azim, azim

        # Beware: This function is called a lot, so it should really be fast

        ####################################################################
        # CHECK IF INSIDE
        ####################################################################

        # At first, we check whether the cell is inside the domain. To keep the
        # number of calculations short, we first use very broad criterions
        # before calculating whether it's inside.
        dx = x - self.x
        dy = y - self.y
        dz = z - self.z
        if grid.periodic:
            if dx > grid.lx/2:
                dx -= grid.lx
            elif dx < -grid.lx/2:
                dx += grid.lx
            if dy > grid.ly/2:
                dy -= grid.ly
            elif dy < -grid.ly/2:
                dy += grid.ly
        # it seems to me that this is a bit faster than abs, and speed here is critical
        if dx > self.max_ab or -dx > self.max_ab or dy > self.max_ab or -dy > self.max_ab or dz > 0:
            geo_ids[0] = -1
            return

        # To decide whether the point is inside, we can just use the normalized
        # distance. Therefore we first calculate the normalized distance vector
        normalized_dx = (dx*self.cosalpha + dy*self.sinalpha)/self.a
        normalized_dy = (-dx*self.sinalpha + dy*self.cosalpha)/self.b
        normalized_dz = dz/self.c
        l2 = normalized_dx**2 + normalized_dy**2 + normalized_dz**2
        if l2 > 1:
            geo_ids[0] = -1
            return

        # At this point we know that the point is in the domain, so we have to
        # assign values

        # check whether it's below the lag surface
        if z < self.z - self.c + self.lag:
            # is inside lag surface
            geo_ids[0] = self.lag_facies
            angles[0] = 0 # no azimuth
            angles[1] = 0 # no dip
            geo_ids[2] = self.num_ha
            return

        if self.structure == 0: # flat
            geo_ids[0] = self.facies
            angles[0] = self.azim
            angles[1] = self.dip
            geo_ids[2] = self.num_ha
            return
        elif self.structure == 1: # dip
            # get distance from zero-plane
            plane_dist = dx*self.normvec_x +\
                         dy*self.normvec_y +\
                         dz*self.normvec_z + self.shift
            # get index in facies vector, where n = len(facies)//2 at d=0
            n = int(plane_dist/self.layer_dist) + self.num_facies//2
            geo_ids[0] = self.facies_array[n]
            angles[0] = self.azim
            angles[1] = self.dip
            geo_ids[2] = self.num_ha
            return
        elif self.structure == 2 or self.structure == 3: # bulb or bulb_sets
            # Since ellipsoids are isosurfaces of quadratic functions, we can
            # get the tangential vector by taking the (negative) gradient of
            # the quadratic function that has our ellipsoid as iso-surface.
            # This funciton can be written as:
            #
            # f(d(x)) = (Ainv * d(x))**2
            #
            # Since the gradient points outwards for this function, we have to
            # use the negative gradient as unit normal vector describing the
            # tangential plane.
            # The full derivation is a bit lengthy.
            aaa = normalized_dx/self.a
            bbb = normalized_dy/self.b
            ccc = normalized_dz/self.c
            len_tanvec = sqrt(aaa**2 + bbb**2 + ccc**2)
            if len_tanvec != 0:
                tanvec_x = (aaa*self.cosalpha + bbb*self.sinalpha)/len_tanvec
                tanvec_y = (-aaa*self.sinalpha + bbb*self.cosalpha)/len_tanvec
                tanvec_z = ccc/len_tanvec
                # The dip is the angle between the normal vector of the tangential
                # plane and the normal vector of the x-y-plane.
                # this means the cosine of the dip is dot(tanvec, [0, 0, 1]),
                # i.e. the third component of norm_vec
                dip = acos(tanvec_z)/np.pi*180
                # The angle is positive if the x-component of the normal vector is
                # positive
                # the normal vector projected and normalized to the x-y-plane,
                # multiplied with [0,1,0] is:
                len_tanvec_xy = sqrt(tanvec_x**2 + tanvec_y**2)
                if len_tanvec_xy == 0:
                    azim = 0
                else:
                    cos_azim = tanvec_y/len_tanvec_xy
                    azim = ho.sign(tanvec_x)*acos(cos_azim)/np.pi*180
                if dip > 90:
                    # if the dip is bigger than 90° we can also just rotate the azimuth
                    dip = 180 - dip
                    # add 180° to azimuth and shift such that it's between -180 and 180
                    azim = (azim + 180 + 180) % 360 - 180
            else:
                # the point is exactly at the center of the trough, this means
                # there is no azim and dip
                azim = 0
                dip = 0

            # use self.dip as maximum dip
            if dip > self.dip:
                dip = self.dip

            if self.structure == 2: # bulb
                geo_ids[0] = self.facies
                angles[0] = azim
                angles[1] = dip
                geo_ids[2] = self.num_ha
                return
            else:
                # create bulb sets
                # The facies can be found by dividing the normalized distance
                # by the normalized bublset_d:
                n = int(ceil((sqrt(l2) + self.shift)*self.c/self.layer_dist))
                geo_ids[0] = self.facies_array[n]
                angles[0] = azim
                angles[1] = dip
                geo_ids[2] = self.num_ha
                return
        else:
            print('Structure:', self.structure)
            raise NotImplementedError('This structure is not implemented yet')
