
import numpy as np
import numpy.typing as npt
from hyvr.geo.trough_utils import generate_trough_positions, rand_trough_params
from hyvr.geo.trough import Trough



#from libc.math cimport sqrt, ceil, acos, pi, fabs, atan
from hyvr.geo.ae_realization import AERealization
from hyvr.geo.grid import Grid
import hyvr.optimized as ho
from numba import jit,njit


class TroughAE(AERealization):


    def create_object_arrays(self):
        # This is super ugly :(
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
        self.object_xmax = np.zeros(self.n_objects, dtype=np.float)
        self.object_xmin = np.zeros(self.n_objects, dtype=np.float)

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
            self.object_xmin[i] = obj.xmin
            self.object_xmax[i] = obj.xmax
            try:
                self.object_shift[i] = obj.shift
            except AttributeError:
                continue
            try:
                self.object_cosdip[i] = obj.cosdip
            except AttributeError:
                continue
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
                self.zmin, self.zmax, self.type_params, grid
            )
        else:
            te_xyz = self.type_params['te_xyz']

        # Generate troughs
        # sorting troughs:
        te_xyz = sorted(te_xyz, key=lambda x: x[2], reverse= False) # Z ORDER MATTERS!!!
        num_ha = 0
        for xc, yc, zc in te_xyz:
            trough_params = rand_trough_params(self.type_params, zc, grid)
            trough_params['x'] = xc
            trough_params['y'] = yc
            trough_params['z'] = zc
            trough_params['num_ha'] = num_ha
            trough = Trough(self.type_params, **trough_params)
            self._add_object(trough)
            num_ha+=1





    def maybe_assign_points_to_object(self, oi: int,
                                        x,
                                        y,
                                        z,
                                        grid: Grid):
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

        object_x = self.object_x[oi]
        object_y = self.object_y[oi]
        object_z = self.object_z[oi]


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
        # if grid.periodic:
        #     if dx > grid.lx/2:
        #         dx -= grid.lx
        #     elif dx < -grid.lx/2:
        #         dx += grid.lx
        #     if dy > grid.ly/2:
        #         dy -= grid.ly
        #     elif dy < -grid.ly/2:
        #         dy += grid.ly
        # # it seems to me that this is a bit faster than abs, and speed here is critical
        # max_ab = self.object_max_ab[oi]
        # if dx > max_ab or -dx > max_ab or dy > max_ab or -dy > max_ab or dz > 0:
        #     geo_ids[0] = -1
        #     return

        # To decide whether the point is inside, we can just use the normalized
        # distance. Therefore we first calculate the normalized distance vector
        cosalpha = self.object_cosalpha[oi]
        sinalpha = self.object_sinalpha[oi]
        normalized_dx = (dx*cosalpha + dy*sinalpha)/self.object_a[oi]
        normalized_dy = (-dx*sinalpha + dy*cosalpha)/self.object_b[oi]
        normalized_dz = dz/self.object_c[oi]
        l2 = normalized_dx**2 + normalized_dy**2 + normalized_dz**2
        logic = l2 < 1
        # To ensure the semi-ellipse shape, we cut the ellipse in half by assuming a negative distance:
        logic = logic & (dz <= 0)
        # At this point we know that the point is in the domain, so we have to
        # assign values


        structure = self.object_structure[oi]
        if structure == 0: # massive
            facies = np.where(logic, self.object_facies[oi],-1)
            azim = np.where(logic,self.object_azim[oi],-1)
            dip = np.where(logic,self.object_dip[oi],-1)
            num_ha = np.where(logic, self.object_num_ha[oi], -1)
            
        elif structure == 1: # dip
            # get distance from zero-plane
            plane_dist = dx*self.object_normvec_x[oi] +\
                         dy*self.object_normvec_y[oi] +\
                         dz*self.object_normvec_z[oi] + self.object_shift[oi]
            # get index in facies vector, where n = len(facies)//2 at d=0
            ns = np.int(plane_dist/self.object_layer_dist[oi]) + np.int(self.object_num_facies[oi]//2)
            facies = np.concatenate([self.object_facies_array[oi,n] for n in ns])
            facies = np.where(logic, facies,-1)
            azim = np.where(logic,self.object_azim[oi],-1)
            dip = np.where(logic,self.object_dip[oi],-1)
            num_ha = np.where(logic, self.object_num_ha[oi], -1)
            
        elif structure == 2 or structure == 3: # bulb or bulb_sets
            # Since ellipsoids are isosurfaces of quadratic functions, we can
            # get the normal vector by taking the gradient of the quadratic
            # function that has our ellipsoid as iso-surface.
            # This function can be written as:
            #
            # f(d(x)) = (x*cos+y*sin)**2/a**2 + (-x*sin+y*cos)**2/b** + z**2/c**2
            #
            # The gradient is (up to a scalar)
            #
            #             /  nx*cos/a + ny*sin/b  \
            # grad f(x) = | -nx*sin/a + ny*cos/b  |
            #             \          nz/c         /
            #
            # where nx, ny, nz are the normalized distances.
            # The gradient points outwards.
            # The length of the vector is
            #
            # |grad f(x)| = (nx/a)**2 + (ny/b)**2 + (nz/c)**2
            #
            # The dip angle is the the angle between the normal
            # vector and the unit z-vector.
            # The azimuth is the angle between the projection of the normal
            # vector onto the x-y-plane and the unit x-vector.
            aaa = normalized_dx/self.object_a[oi]
            bbb = normalized_dy/self.object_b[oi]
            ccc = normalized_dz/self.object_c[oi]
            len_normvec = np.sqrt(aaa**2 + bbb**2 + ccc**2)
            normvec_x = (aaa*self.object_cosalpha[oi]
                         + bbb*self.object_sinalpha[oi])/len_normvec
            normvec_y = (-aaa*self.object_sinalpha[oi]
                         + bbb*self.object_cosalpha[oi])/len_normvec
            normvec_z = ccc/len_normvec
           
            len_normvec_xy = np.sqrt(normvec_x**2 + normvec_y**2)
            
            
            # The dip angle can be found as acos of the normalized
            # scalar product of unit z-vector and normal vector
            # It should be negative, if the gradient points in positive
            # x direction and in positive z direction, or if it points
            # in negative x and z direction.
            # direction.
            dip = -np.acos(np.abs(normvec_z))\
                * ho.sign(normvec_x*normvec_z)/ np.pi*180

            # The azimuth angle should also be between -90 and 90. It is
            # the angle between the projection of the normal vector into
            # the xy-plane and the x unit vector.
            # It can be calculated as negative atan of the y component
            # over the x component. To make sure the angle is between
            # -90 and 90, we change the sign of both components, if the
            # x component is negative. This means the x component is
            # always positive, and the y component is multiplied by the
            # sign of the x component
            azim = -np.atan(ho.sign(normvec_x)*normvec_y/np.abs(normvec_x))\
                   / np.pi*180
            
            #if len_normvec != 0:
            # the point is exactly at the center of the trough, this means
            # there is no azim and dip
            azim = np.where(len_normvec == 0, 0, azim)
            dip = np.where(len_normvec == 0, 0, dip)
            

            
            # normal vector points in z-direction -> dip = 0, azim = 0
            azim = np.where(len_normvec_xy == 0, 0, azim)
            dip = np.where(len_normvec_xy == 0, 0, dip)



            # use self.object_dip[oi] as maximum dip
            dip = np.where(dip > self.object_dip[oi],self.object_dip[oi], dip)
            dip = np.where(dip < - self.object_dip[oi], -self.object_dip[oi], dip)
            
            azim = np.where(logic, azim,-1)
            dip = np.where(logic, dip,-1)
            num_ha = np.where(logic, self.object_num_ha[oi],-1)

            if structure == 2: # bulb
                facies = np.where(logic, self.object_facies[oi],-1)
               

            else:
                # create bulb sets
                # The facies can be found by dividing the normalized distance
                # by the normalized bublset_d:
                ns = np.int(np.ceil((np.sqrt(l2) + self.object_shift[oi])*self.object_c[oi]/self.object_layer_dist[oi]))
                facies = np.concatenate([self.object_facies_array[oi,n] for n in ns])
                facies = np.where(logic, facies,-1)

        else:
            print('Structure:', structure)
            raise NotImplementedError('This structure is not implemented yet')
        
        # check whether it's below the lag surface
        lag_logic = z < (object_z - self.object_c[oi] + self.object_lag[oi])
        # is inside lag surface
        facies = np.where(lag_logic, self.object_lag_facies[oi], facies)
        azim = np.where(logic, 0, azim) # no azimuth
        dip = np.where(lag_logic,0, dip) # no dip
        
        return facies, azim, dip, num_ha
