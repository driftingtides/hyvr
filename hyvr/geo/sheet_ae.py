
import numpy as np
import numpy.typing as npt

from hyvr.geo.sheet import Sheet

from hyvr.geo.grid import Grid
from hyvr.geo.contact_surface import contact_surface


from hyvr.geo.ae_realization import AERealization

class SheetAE(AERealization):

    def create_object_arrays(self):

        self.object_shift = np.zeros(self.n_objects, dtype=np.float)
        self.object_layer_dist = np.zeros(self.n_objects, dtype=np.float)
        self.object_normvec_x = np.zeros(self.n_objects, dtype=np.float)
        self.object_normvec_y = np.zeros(self.n_objects, dtype=np.float)
        self.object_normvec_z = np.zeros(self.n_objects, dtype=np.float)
        self.object_bottom_surface_zmean = np.zeros(self.n_objects, dtype=np.float)
        self.object_dipsets = np.zeros(self.n_objects, dtype=np.int32)

        if self.n_objects > 0:
            nx, ny = self.object_list[0].bottom_surface.shape
        else:
            nx, ny = 0, 0
        self.object_bottom_surface = np.zeros((self.n_objects, nx, ny), dtype=np.float)
        self.object_top_surface = np.zeros((self.n_objects, nx, ny), dtype=np.float)
        #self.object_num_ha = np.zeros(self.n_objects, dtype=np.int32)

        for i, obj in enumerate(self.object_list):
            if obj.dipsets:
                self.object_dipsets[i] = obj.dipsets
                self.object_shift[i] = obj.shift
                self.object_layer_dist[i] = obj.layer_dist
                self.object_normvec_x[i] = obj.normvec_x
                self.object_normvec_y[i] = obj.normvec_y
                self.object_normvec_z[i] = obj.normvec_z
            self.object_bottom_surface_zmean[i] = np.mean(obj.bottom_surface)
            
            #top_surface = obj.top_surface
            bottom_surface = obj.bottom_surface
            self.object_num_ha[i] = obj.num_ha

            #self.object_top_surface[i,:,:] = top_surface
            self.object_bottom_surface[i,:,:] = bottom_surface


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

        bottom_surface = self.bottom_surface
        # create all sheets except the lowest one
        num_ha = 0
        z_iter = np.mean(bottom_surface)
        while z_iter < self.zmax:
            # normal sheet
            ztop = np.mean(bottom_surface) + thickness
            
            
            # generate sheet object
            sheet = Sheet(self.type_params, bottom_surface, ztop, grid, num_ha)
            self._add_object(sheet)
            num_ha+=1

            # use current top as new bottom
            bottom_surface = contact_surface(grid, mode='flat', z=ztop)
            z_iter = ztop
            




    def maybe_assign_points_to_object(self, oi: int,                                     
                                        x, y, z,
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

        normvec_x = self.object_normvec_x[oi]
        normvec_y = self.object_normvec_y[oi]
        normvec_z = self.object_normvec_z[oi]

        # if we reach this functions, we now that we are within [zmin, zmax]
        # the only way how this could be outside the domain is when the top or
        # bottom surfaces are not flat
        #z_above = self.object_zmaxs[oi]
        #z_below = self.object_bottom_surface[oi, :, :]

        # at this point we know that the point is inside the sheet
        azim = self.object_azim[oi]
        dip = self.object_dip[oi]
        num_ha = self.object_num_ha[oi]
        if self.object_dipsets[oi]:
            d = normvec_x*x + normvec_y*y + normvec_z*z - self.object_shift[oi]
            ns = np.int(d/self.object_layer_dist[oi] + self.object_num_facies[oi]//2)
            facies = np.concatenate([self.object_facies_array[oi,n] for n in ns])
        else:
            facies = self.object_facies[oi]
        
        return facies, azim, dip, num_ha
