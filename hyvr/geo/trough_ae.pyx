
import numpy as np
from hyvr.geo.trough_utils import *
from hyvr.geo.trough import Trough

cimport cython
cimport numpy as np
from libc.math cimport sqrt, ceil, acos, pi
from hyvr.geo.ae_realization cimport AERealization
cimport hyvr.optimized as ho


cdef class TroughAE(AERealization):

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef create_object_arrays(self):
        # This is super ugly :(
        cdef int i
        self.object_x = np.zeros(self.n_objects, dtype=np.float)
        self.object_y = np.zeros(self.n_objects, dtype=np.float)
        self.object_z = np.zeros(self.n_objects, dtype=np.float)
        self.object_a = np.zeros(self.n_objects, dtype=np.float)
        self.object_b = np.zeros(self.n_objects, dtype=np.float)
        self.object_c = np.zeros(self.n_objects, dtype=np.float)
        self.object_max_ab = np.zeros(self.n_objects, dtype=np.float)
        self.object_alpha = np.zeros(self.n_objects, dtype=np.float)
        self.object_sinalpha = np.zeros(self.n_objects, dtype=np.float)
        self.object_cosalpha = np.zeros(self.n_objects, dtype=np.float)
        self.object_lag = np.zeros(self.n_objects, dtype=np.float)
        self.object_shift = np.zeros(self.n_objects, dtype=np.float)
        self.object_cosdip = np.zeros(self.n_objects, dtype=np.float)
        self.object_layer_dist = np.zeros(self.n_objects, dtype=np.float)
        self.object_normvec_x = np.zeros(self.n_objects, dtype=np.float)
        self.object_normvec_y = np.zeros(self.n_objects, dtype=np.float)
        self.object_normvec_z = np.zeros(self.n_objects, dtype=np.float)
        self.object_lag_facies = np.zeros(self.n_objects, dtype=np.int32)
        self.object_structure = np.zeros(self.n_objects, dtype=np.int32)

        for i, obj in enumerate(self.object_list):
            self.object_x[i] = obj.x
            self.object_y[i] = obj.y
            self.object_z[i] = obj.z
            self.object_a[i] = obj.a
            self.object_b[i] = obj.b
            self.object_c[i] = obj.c
            self.object_max_ab[i] = obj.max_ab
            self.object_alpha[i] = obj.alpha
            self.object_sinalpha[i] = obj.sinalpha
            self.object_cosalpha[i] = obj.cosalpha
            self.object_lag[i] = obj.lag
            self.object_shift[i] = obj.shift
            self.object_cosdip[i] = obj.cosdip
            self.object_layer_dist[i] = obj.layer_dist
            self.object_normvec_x[i] = obj.normvec_x
            self.object_normvec_y[i] = obj.normvec_y
            self.object_normvec_z[i] = obj.normvec_z
            self.object_lag_facies[i] = obj.lag_facies
            self.object_structure[i] = obj.structure
            

    def generate_objects(self, grid):
        """
        Generate trough objects and place them in the domain.
        """

        # get trough positions
        if self.type_params['te_xyz'] is None:
            te_xyz = generate_trough_positions(
                self.bottom_surface, self.top_surface, self.type_params, grid
            )
        else:
            te_xyz = self.type_params['te_xyz']

        # Generate troughs
        for xc, yc, zc in te_xyz:
            trough_params = rand_trough_params(self.type_params, zc, grid)
            trough_params['x'] = xc
            trough_params['y'] = yc
            trough_params['z'] = zc
            trough = Trough(self.type_params, **trough_params)
            self._add_object(trough)




    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    @cython.cdivision(True)
    cpdef maybe_assign_points_to_object(self, int oi,
                                        np.int32_t [:] geo_ids,
                                        np.float_t [:] angles,
                                        np.float_t x, np.float_t y, np.float_t z,
                                        int x_idx, int y_idx,
                                        Grid grid):
        """
        This function checks whether the current grid cell with given
        coordinates is inside the trough and assigns facies, azimuth and dip by
        altering the passed arrays.

        Parameters
        ----------
        oi : int
            object index
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
        cdef np.float_t dx, dy, dz
        cdef np.float_t normalized_dx, normalized_dy, normalized_dz
        cdef np.float_t l2, plane_dist, aaa, bbb, ccc, len_tanvec
        cdef np.float_t tanvec_x, tanvec_y, tanvec_z
        cdef np.int32_t n, structure
        cdef np.float_t dip, cos_azim, azim

        cdef np.float_t object_x = self.object_x[oi]
        cdef np.float_t object_y = self.object_y[oi]
        cdef np.float_t object_z = self.object_z[oi]


        # Beware: This function is called a lot, so it should really be fast

        ####################################################################
        # CHECK IF INSIDE
        ####################################################################

        # At first, we check whether the cell is inside the domain. To keep the
        # number of calculations short, we first use very broad criterions
        # before calculating whether it's inside.
        dx = x - object_x
        dy = y - object_y
        dz = z - object_z
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
        cdef np.float_t max_ab = self.object_max_ab[oi]
        if dx > max_ab or -dx > max_ab or dy > max_ab or -dy > max_ab or dz > 0:
            geo_ids[0] = -1
            return

        # To decide whether the point is inside, we can just use the normalized
        # distance. Therefore we first calculate the normalized distance vector
        cdef np.float_t cosalpha = self.object_cosalpha[oi]
        cdef np.float_t sinalpha = self.object_sinalpha[oi]
        normalized_dx = (dx*cosalpha + dy*sinalpha)/self.object_a[oi]
        normalized_dy = (-dx*sinalpha + dy*cosalpha)/self.object_b[oi]
        normalized_dz = dz/self.object_c[oi]
        l2 = normalized_dx**2 + normalized_dy**2 + normalized_dz**2
        if l2 > 1:
            geo_ids[0] = -1
            return

        # At this point we know that the point is in the domain, so we have to
        # assign values

        # check whether it's below the lag surface
        if z < object_z - self.object_c[oi] + self.object_lag[oi]:
            # is inside lag surface
            geo_ids[0] = self.object_lag_facies[oi]
            angles[0] = 0 # no azimuth
            angles[1] = 0 # no dip
            geo_ids[2] = self.object_num_ha[oi]
            return

        structure = self.object_structure[oi]
        if structure == 0: # flat
            geo_ids[0] = self.object_facies[oi]
            angles[0] = self.object_azim[oi]
            angles[1] = self.object_dip[oi]
            geo_ids[2] = self.object_num_ha[oi]
            return
        elif structure == 1: # dip
            # get distance from zero-plane
            plane_dist = dx*self.object_normvec_x[oi] +\
                         dy*self.object_normvec_y[oi] +\
                         dz*self.object_normvec_z[oi] + self.object_shift[oi]
            # get index in facies vector, where n = len(facies)//2 at d=0
            n = int(plane_dist/self.object_layer_dist[oi]) + self.object_num_facies[oi]//2
            geo_ids[0] = self.object_facies_array[oi,n]
            angles[0] = self.object_azim[oi]
            angles[1] = self.object_dip[oi]
            geo_ids[2] = self.object_num_ha[oi]
            return
        elif structure == 2 or structure == 3: # bulb or bulb_sets
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
            aaa = normalized_dx/self.object_a[oi]
            bbb = normalized_dy/self.object_b[oi]
            ccc = normalized_dz/self.object_c[oi]
            len_tanvec = sqrt(aaa**2 + bbb**2 + ccc**2)
            if len_tanvec != 0:
                tanvec_x = (aaa*self.object_cosalpha[oi] + bbb*self.object_sinalpha[oi])/len_tanvec
                tanvec_y = (-aaa*self.object_sinalpha[oi] + bbb*self.object_cosalpha[oi])/len_tanvec
                tanvec_z = ccc/len_tanvec
                # The dip is the angle between the normal vector of the tangential
                # plane and the normal vector of the x-y-plane.
                # this means the cosine of the dip is dot(tanvec, [0, 0, 1]),
                # i.e. the third component of norm_vec
                dip = acos(tanvec_z)/pi*180
                # The angle is positive if the x-component of the normal vector is
                # positive
                # the normal vector projected and normalized to the x-y-plane,
                # multiplied with [0,1,0] is:
                len_tanvec_xy = sqrt(tanvec_x**2 + tanvec_y**2)
                if len_tanvec_xy == 0:
                    azim = 0
                else:
                    cos_azim = tanvec_y/len_tanvec_xy
                    azim = ho.sign(tanvec_x)*acos(cos_azim)/pi*180
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

            # use self.object_dip[oi] as maximum dip
            if dip > self.object_dip[oi]:
                dip = self.object_dip[oi]

            if structure == 2: # bulb
                geo_ids[0] = self.object_facies[oi]
                angles[0] = azim
                angles[1] = dip
                geo_ids[2] = self.object_num_ha[oi]
                return
            else:
                # create bulb sets
                # The facies can be found by dividing the normalized distance
                # by the normalized bublset_d:
                n = int(ceil((sqrt(l2) + self.object_shift[oi])*self.object_c[oi]/self.object_layer_dist[oi]))
                geo_ids[0] = self.object_facies_array[oi,n]
                angles[0] = azim
                angles[1] = dip
                geo_ids[2] = self.object_num_ha[oi]
                return
        else:
            print('Structure:', structure)
            raise NotImplementedError('This structure is not implemented yet')