import numpy as np
import hyvr.optimized as ho
import hyvr.utils as hu
from hyvr.geo.contact_surface import ContactSurface

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
            ae_bottom_surface.use_lower_surface_value(ae_top_surface)

            # close to the bottom it might also be possible that the ae_bottom
            # is below the stratum bottom if the stratum bottom is not flat
            # in this case we will assign the (higher) bottom surface value instead
            ae_bottom_surface.use_higher_surface_value(self.bottom_surface)

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

        # convert lists to arrays
        self.ae_zmaxs = np.array(self.ae_zmaxs)
        self.ae_zmins = np.array(self.ae_zmins)
