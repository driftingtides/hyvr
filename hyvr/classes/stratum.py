import numpy as np
import hyvr.optimized as ho
import hyvr.utils as hu
from hyvr.classes.contact_surface import ContactSurface
from hyvr.classes.contact_surface_utils import *

def prob_choose(types, probs):
    # TODO: At the moment the given probabilities are only relative
    return np.random.choice(types, p=np.asarray(probs)/np.sum(probs))

class Stratum:

    def __init__(self, bottom_surface, top_surface, params, ae_types, grid):
        """
        Stratum constructor. This generates the stratum by generating the AEs
        within.

        Parameters
        ----------
        bottom_surface : ContactSurface object
            The bottom surface of the stratum
        top_surface : ContactSurface object
            The top surface of the stratum
        params : dict
            Dictionary with other necessary parameters
        ae_types : list of AEType objects
            The possible AE types within this stratum
        grid : Grid object
        """
        self.bottom_surface = bottom_surface
        self.top_surface = top_surface
        self.params = params
        self.name = params['strata']
        self.ae_types = ae_types
        self.ae_type_names = [ae_type.name for ae_type in ae_types]
        self.bg_facies = params['bg_facies']
        self.bg_azim = params['bg_azim']
        self.bg_dip = params['bg_dip']

        # generate AEs
        self.zmin = self.bottom_surface.zmin
        self.zmax = self.top_surface.zmax
        self.ae_zmins = []
        self.ae_zmaxs = []
        self.aes = []

        # first top is top of stratum
        ae_top_surface = self.top_surface
        # choose type of first AE
        curr_ae_type = prob_choose(self.ae_types, self.params['ae_prob'])
        znow = self.top_surface.zmean
        while znow > self.bottom_surface.zmean:

            num_type = self.ae_type_names.index(curr_ae_type.name)

            hu.print_to_stdout('Generating', curr_ae_type.name, 'from {:2.3f} m'.format(znow))

            # find height of AE
            mean_height = self.params['ae_z_mean'][num_type]
            # TODO: It might be nice to specify the height variance in the inifile
            height = np.random.normal(mean_height, mean_height*0.1)
            # potential avulsion
            if np.random.uniform() <= self.params['avul_prob'][0]:
                dz = -np.random.uniform(*self.params['avul'])
            else:
                dz = 0
            znow -= height + dz

            # create bottom contact surface
            # -----------------------------
            # the type of the bottom surface is determined by the type of the
            # AE below, if it is not the lowest AE
            next_ae_type = prob_choose(self.ae_types, self.params['ae_prob'])
            if znow <= self.bottom_surface.zmean:
                # limit top surface of uppermost AE to top surface of stratum
                znow = self.bottom_surface.zmean
                ae_bottom_surface = self.bottom_surface
            else:
                bottom_surface_params = next_ae_type.params['contact_model']
                bottom_surface_params['z'] = znow
                ae_bottom_surface = ContactSurface(grid, **bottom_surface_params)

            # if the bottom surface is above the top surface, we will assign
            # the (lower) top surface values instead
            use_lower_surface_value(ae_bottom_surface, ae_top_surface)

            # close to the bottom it might also be possible that the ae_bottom
            # is below the stratum bottom if the stratum bottom is not flat
            # in this case we will assign the (higher) bottom surface value instead
            use_higher_surface_value(ae_bottom_surface, self.bottom_surface)

            # generate element: this will place all the objects within an AE
            element = curr_ae_type.generate_realization(ae_bottom_surface, ae_top_surface, self, grid)
            self.aes.append(element)

            # if the minimum/maximum depth changed, update it
            ae_zmin = element.zmin
            if ae_zmin < self.zmin:
                self.zmin = ae_zmin
            self.ae_zmins.append(ae_zmin)
            ae_zmax = element.zmax
            if ae_zmax > self.zmax:
                self.zmax = ae_zmax
            self.ae_zmaxs.append(ae_zmax)

            # current bottom is top of next, next ae is new current ae
            ae_top_surface = ae_bottom_surface
            curr_ae_type = next_ae_type

        self.n_ae = len(self.aes)

    # def maybe_assign_facies_azim_dip(self, facies, angles, x, y, z, x_idx, y_idx, aei, oi, grid):
        # """
        # Assigns facies, azim and dip if the cell is inside this stratum or in
        # one of it's elements.

        # Since the objects inside the elements might reach into the domain
        # below, we can't just check the top and bottom surface but must check
        # every object.

        # Parameters
        # ----------
        # facies : np.ndarray[np.int32, dim=1]
            # Array of length one that hold the facies and will be altered by
            # this function. If this is -1 afterwards, the cell is not in the
            # stratum.
        # angles : np.ndarray[np.float, dim=1]
            # Array of length 2 that holds the azimuth and dip angles, will also
            # be altered.
        # x, y, z : float
            # x, y, and z positions of the cell in question
        # x_idx, y_idx : int
            # indices of x and y cells
        # aei : int
            # Index of where to start in the architectural element list. This is
            # typically the index of the last found element.
        # oi : int
            # Geometrical object index, similar to aei
        # grid : Grid object

        # Returns
        # -------
        # aei, oi : int
            # Updated indices
        # """
        # oi_orig = int(oi)
        # aei_orig = int(aei)
        # while aei < self.n_ae:
            # # if z > self.ae_zmaxs[aei]:
                # # if the current z is above the maximum z of the current AE, all other AEs will also
                # # be below the current cell, so we can stop searching for an AE
                # # This case should not happen
                # # break
            # if z < self.ae_zmins[aei]:
                # # if we're below the current AE, we can go search in the next one
                # oi = 0
                # aei += 1
                # continue
            # else:
                # self.aes[aei].maybe_assign_facies_azim_dip(
                    # facies, angles, x, y, z, x_idx, y_idx, oi, grid
                # )
                # if facies[0] != -1:
                    # # we found something
                    # return aei, oi
                # else:
                    # # if we didn't find something here, we check the next AE
                    # aei += 1
                    # oi = 0
                    # continue
