
import numpy as np
from hyvr.geo.sheet import Sheet

cimport cython
cimport numpy as np
from hyvr.geo.ae_realization cimport AERealization

cdef class SheetAE(AERealization):

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef create_object_arrays(self):
        # This is super ugly :(
        cdef int i, nx, ny, j, k
        cdef np.float_t [:,:] top_surface, bottom_surface

        self.object_shift = np.zeros(self.n_objects, dtype=np.float)
        self.object_layer_dist = np.zeros(self.n_objects, dtype=np.float)
        self.object_normvec_x = np.zeros(self.n_objects, dtype=np.float)
        self.object_normvec_y = np.zeros(self.n_objects, dtype=np.float)
        self.object_normvec_z = np.zeros(self.n_objects, dtype=np.float)
        self.object_bottom_surface_zmean = np.zeros(self.n_objects, dtype=np.float)
        self.object_dipsets = np.zeros(self.n_objects, dtype=np.int32)

        if self.n_objects > 0:
            nx, ny = self.object_list[0].bottom_surface.surface.shape
        else:
            nx, ny = 0, 0
        self.object_bottom_surface = np.zeros((self.n_objects, nx, ny), dtype=np.float)
        self.object_top_surface = np.zeros((self.n_objects, nx, ny), dtype=np.float)

        for i, obj in enumerate(self.object_list):
            self.object_shift[i] = obj.shift
            self.object_layer_dist[i] = obj.layer_dist
            self.object_normvec_x[i] = obj.normvec_x
            self.object_normvec_y[i] = obj.normvec_y
            self.object_normvec_z[i] = obj.normvec_z
            self.object_bottom_surface_zmean[i] = obj.bottom_surface.zmean
            self.object_dipsets[i] = obj.dipsets
            top_surface = obj.top_surface.surface
            bottom_surface = obj.bottom_surface.surface
            for j in range(nx):
                for k in range(ny):
                    self.object_top_surface[i,j,k] = top_surface[j,k]
                    self.object_bottom_surface[i,j,k] = bottom_surface[j,k]


    def generate_objects(self, grid):
        """
        Generate sheet objects and place them in the domain.
        """

        # A sheet AE can consist of either one massive sheet (i.e. the sheet is
        # just the full AE) or of multiple sheets.
        # Internally, the sheets can again be homogeneous or consist of dipping
        # sets.
        #
        # Sheets have a bottom and top surface that is generally a flat contact
        # surface. Only the topmost and the bottom-most sheets have potentially
        # non-flat contact surfaces as they use the contact surface of the AE.

        if self.type_params['lens_thickness'] == -1:
            # massive bedding
            thickness = self.zmax - self.zmin
        elif self.type_params['size_ztrend'] is not None:
            zfactor = np.interp(np.mean([self.zmin, self.zmax]),
                                [grid.z0, grid.zmax],
                                self.type_params['size_ztrend'])
            thickness = self.type_params['lens_thickness'] * zfactor
        else:
            thickness = self.type_params['lens_thickness']

        zbottom = self.zmax - thickness
        top_surface = self.top_surface
        # create all sheets except the lowest one
        last_sheet = False
        while not last_sheet:
            # generate bottom surface
            if zbottom > self.zmin:
                # normal sheet
                bottom_surface = ContactSurface(grid, mode='flat', z=zbottom)
                # it's possible that this bottom surface is above the top
                # surface if we're close to the AE top. Then we will use the
                # lower value (i.e. the top surface value)
                bottom_surface.use_lower_surface_value(top_surface)

                # close to the bottom it might also be possible that the AE
                # bottom is higher than the current bottom, so we have to use
                # the higher value
                top_surface.use_higher_surface_value(self.bottom_surface)

            else:
                last_sheet = True
                bottom_surface = self.bottom_surface

            # generate sheet object
            sheet = Sheet(self.type_params, bottom_surface, top_surface, grid)
            self._add_object(sheet)

            # use current bottom as new top, zbottom decreases
            zbottom -= thickness
            top_surface = bottom_surface



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
        cdef np.float_t z_above, z_below, d
        cdef np.float_t normvec_x, normvec_y, normvec_z
        cdef np.float_t shift, azim, dip
        cdef np.int32_t n, num_facies

        normvec_x = self.object_normvec_x[oi]
        normvec_y = self.object_normvec_y[oi]
        normvec_z = self.object_normvec_z[oi]

        # if we reach this functions, we now that we are within [zmin, zmax]
        # the only way how this could be outside the domain is when the top or
        # bottom surfaces are not flat
        z_above = self.object_top_surface[oi, x_idx, y_idx]
        z_below = self.object_bottom_surface[oi, x_idx, y_idx]
        if z < z_below or z > z_above:
            geo_ids[0] = -1
            return

        # at this point we now that the point is inside the sheet
        angles[0] = self.object_azim[oi]
        angles[1] = self.object_dip[oi]
        geo_ids[2] = self.object_num_ha[oi]
        if self.object_dipsets[oi]:
            d = self.object_normvec_x[oi]*x + self.object_normvec_y[oi]*y + self.object_normvec_z[oi]*z - self.object_shift[oi]
            n = int(d/self.object_layer_dist[oi]) + self.object_num_facies[oi]//2
            geo_ids[0] = self.object_facies_array[oi,n]
            return
        else:
            geo_ids[0] = self.object_facies[oi]
            return
