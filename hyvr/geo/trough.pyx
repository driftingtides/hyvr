import numpy as np
cimport hyvr.optimized as ho
import hyvr.utils as hu
from libc.math cimport sin, cos, ceil, acos, sqrt
cimport numpy as np
from hyvr.geo.grid cimport Grid
cimport cython

cdef class Trough:
    """
    This class holds all parameters of a geometrical trough object.
    
    The only method of this class is __init__, as after creating all trough
    objects, the list of troughs in TroughAE will be converted to multiple lists
    containing of attribute values for easier access from Cython.
    
    In principle this could also be implemented in pure Python, anyone
    interested in rewriting this class (and Channel and Sheet) can do so.
    """
    cdef public:
        np.float_t x, y, z
        np.float_t a, b, c
        np.float_t max_ab
        np.float_t zmin, zmax
        np.float_t alpha
        np.float_t sinalpha, cosalpha
        np.float_t lag, shift
        np.int32_t structure
        np.float_t layer_dist
        np.float_t dip, azim
        np.int32_t facies, lag_facies, num_ha
        np.float_t cosdip
        np.float_t normvec_x, normvec_y, normvec_z
        np.int32_t num_facies
        np.int32_t [:] facies_array

    def __init__(self, type_params, *, x, y, z, a, b, c, alpha):
        cdef np.float_t sin_dip, cos_dip, sin_azim, cos_azim

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
        # along the x-, y- and z-axes is performed, and then a rotation along
        # the z-axis with angle alpha.
        # In matrix form, this operation looks as follows:
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
        # 0 : 'massive'
        # 1 : 'dip'
        # 2 : 'bulb'
        # 3 : 'bulb_sets'
        # 4 : 'random'
        self.structure = ['massive', 'dip', 'bulb', 'bulb_sets', 'random'].index(type_params['structure'])
        if self.structure == 4: # random
            self.structure = np.random.randint(4)
        if self.structure == 0:  # massive
            # everything has the same facies, azim, and dip
            # alpha (rotation angle of trough) is added to azimuth
            self.facies = np.random.choice(type_params['facies'])
            self.azim = np.random.uniform(*type_params['azimuth']) + self.alpha
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
            # First, we apply the dip rotation by rotating by the dip
            # angle along the y-axis to the right (anti-clockwise).
            # Then we rotate around the z-axis by the azimuth angle along the
            # z-axis to the left (clockwise)
            # In matrix notation:
            #
            #     / cos(azim)  sin(azim)  0 \   /  cos(dip)  0  sin(dip) \   / 0 \
            # n = |-sin(azim)  cos(azim)  0 | * |     0      1     0     | * | 0 |
            #     \    0            0     1 /   \ -sin(dip)  0  cos(dip) /   \ 1 /
            #
            #
            #     / cos(azim)  sin(azim)  0 \   /  sin(dip) \
            # n = |-sin(azim)  cos(azim)  0 | * |    0      |
            #     \    0            0     1 /   \  cos(dip) /
            #
            #
            #     /  sin(dip)*cos(azim) \
            # n = | -sin(dip)*sin(azim) |
            #     \       cos(dip)      /
            #
            #
            # The azimuth is the given azimuth **with** the rotation of the
            # trough, all angles are in degree
            self.azim = np.random.uniform(*type_params['azimuth']) + self.alpha
            self.dip = np.random.uniform(*type_params['dip'])
            sin_dip = sin(self.dip*np.pi/180)
            cos_dip = cos(self.dip*np.pi/180)
            sin_azim = sin(self.azim*np.pi/180)
            cos_azim = cos(self.azim*np.pi/180)
            self.normvec_x = sin_dip*cos_azim
            self.normvec_y = -sin_dip*sin_azim
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
            # only used as maximum possible dip
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
