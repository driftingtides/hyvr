"""
This class holds the abstraction of a HyVR model.
It can be broadly structured in two parts:

* **Model description**: This consists of all info about the model that is read
  from the ini-file.
* **Model realization**: This consists of several 3-d arrays which hold the
    data for the different variables on the grid. This is what get's put out in
    the end. The data is stored in a dictionary of arrays.
"""

import numpy as np
import hyvr.utils as hu
import hyvr.optimized as ho
from hyvr.classes.contact_surface import ContactSurface
from hyvr.classes.contact_surface_utils import *
from hyvr.classes.stratum import Stratum
from hyvr.classes.ae_types import AEType
from hyvr.classes.grid import Grid


class Model:

    def __init__(self, model_dict, strata_dict, elements, flowtrans_section):
        """
        Creates the model object based on the parsed parameter file.

        Parameters
        ----------
        model_dict : dictionary
            The parsed 'model' section of the parameter-file.
        strata_dict : dictionary
            The parsed 'strata' section of the parameter-file.
        elements : dictionary of dictionaries
            Dictionary containing the different element sections.
        flowtrans_section : dictionary
            The parsed 'flowtrans' section of the parameter-file.
        self.data = {}
        """
        self.grid = Grid(
                model_dict['x0'],
                model_dict['y0'],
                model_dict['z0'],
                model_dict['dx'],
                model_dict['dy'],
                model_dict['dz'],
                model_dict['lx'],
                model_dict['ly'],
                model_dict['lz'],
                model_dict['periodic'],
                )
        self.generate_ae_types(elements)

        self.strata_dict = strata_dict

        # anisotropy and heterogeneity
        self.generate_anisotropy = model_dict['anisotropy']
        self.generate_hydraulics = model_dict['hydraulics']
        self.generate_heterogeneity = model_dict['heterogeneity']
        self.heterogeneity_level = model_dict['heterogeneity_level']

        # flowtrans is necessary for model output
        self.flowtrans = flowtrans_section

        self.data = {}



    # STRATA AND ARCHITECTURAL ELEMENT GENERATION
    # ===================================================================================

    def generate_ae_types(self, elements):
        """
        Read in elements sections and generate type objects
        """
        # TODO: This could maybe be moved to model setup
        self.ae_types = {}
        for elem in elements:
            self.ae_types[elem] = AEType(elements[elem], elem)


    def generate_model(self):
        """
        Generates the model.

        This function generates the full model by creating all strata, the AEs
        within the strata and the geometrical objects within these.

        The method starts from the domain top and calls the strata generation
        function for each stratum, which handles the generation of
        architectural elements.
        """
        hu.print_to_stdout('Generating model...')
        strata_dict = self.strata_dict
        self.strata = []
        self.n_strata = len(strata_dict['strata'])

        # These are the additional parameters that a stratum needs
        param_names = [
            'bg_facies',
            'bg_azim',
            'bg_dip',
            'ae_prob',
            'ae_z_mean',
            'avul',
            'avul_prob',
            'strata',
        ]

        top_surface = ContactSurface(self.grid, mode='flat', z=self.grid.lz)
        for si in range(self.n_strata):
            # create top surface
            if si == self.n_strata - 1:
                bottom_surface = ContactSurface(self.grid, mode='flat', z=self.grid.z0)
            else:
                bottom_surface = ContactSurface(self.grid, **strata_dict['contact_models'][-1-si])
            # if the bottom surface is above the top surface, we will assign
            # the top surface values instead
            use_lower_surface_value(bottom_surface, top_surface)

            stratum_params = {name:strata_dict[name][-1-si] for name in param_names}
            ae_types = [self.ae_types[ae_type] for ae_type in strata_dict['ae_in_strata'][-1-si]]


            stratum = Stratum(bottom_surface,
                              top_surface,
                              stratum_params,
                              ae_types,
                              self.grid,
            )

            self.strata.append(stratum)

            # use old top as bottom
            top_surface = bottom_surface


        # enumerate all objects and AEs
        num_ae = 0
        num_ha = 0
        for stratum in self.strata:
            for ae in stratum.aes:
                ae.num = num_ae
                num_ae += 1
                for obj in ae.objects:
                    obj.num_ha = num_ha
                    num_ha += 1




    def generate_hydraulic_parameters(self, hydraulics):
        """
        Generates (heterogeneous) hydraulic parameters

        Parameters
        ----------
        hydraulics : dictionary
            parsed hydraulics section of the ini-file
        """

        hu.print_to_stdout('Generating hydraulic parameters')

        het = self.generate_heterogeneity
        ani = self.generate_anisotropy

        ha_arr = self.data['ha']
        fac = self.data['facies']
        ae_arr = self.data['ae']
        azim = self.data['azim']
        dip = self.data['dip']

        # Initialise storage arrays:
        # hydraulic conductivity
        self.data['k_iso'] = np.zeros((self.grid.nx, self.grid.ny, self.grid.nz),
                                      dtype=np.float)
        # porosity
        self.data['poros' ]= np.zeros((self.grid.nx, self.grid.ny, self.grid.nz),
                                      dtype=np.float)
        # anisotropy ratio
        self.data['anirat' ]= np.ones((self.grid.nx, self.grid.ny, self.grid.nz),
                                      dtype=np.float)

        k_iso = self.data['k_iso']
        poros = self.data['poros']
        anirat = self.data['anirat']

        if het:
            # Heterogeneous case
            for mi in np.unique(ha_arr):
                for fi in np.unique(fac[ha_arr == mi]):
                    mifi = (ha_arr == mi) & (fac == fi)    # Get mask for relevant values

                    if self.heterogeneity_level == 'internal':

                        # Generate internal heterogeneity - hydraulic conductivity
                        k_flat = hu.specsim(self.grid,
                                            hydraulics['sig_y'][fi],
                                            hydraulics['ycorlengths'][fi],
                                            covmod='exp',
                                            selection_mask=mifi)

                        k_flat = np.exp(k_flat) * hydraulics['k_h'][fi]
                        k_iso[mifi] = k_flat

                        # Generate internal heterogeneity - porosity
                        n_flat = hu.specsim(self.grid, 
                                            hydraulics['sig_n'][fi],
                                            hydraulics['ncorlengths'][fi],
                                            covmod='exp',
                                            selection_mask=mifi) + hydraulics['n'][fi]
                        poros[mifi] = n_flat

                    elif self.heterogeneity_level == 'facies':
                        # Assign heterogeneity at facies level only
                        k_iso[mifi] = hydraulics['k_h'][fi]
                        poros[mifi] = hydraulics['n'][fi]

                    # Assign anisotropy ratio
                    anirat[mifi] = hydraulics['k_ratio'][fi]


            """ Assign background heterogeneity per architectural element """
            if self.heterogeneity_level == 'internal':
                # for si, aei in enumerate(ae_lu):
                for stratum in self.strata:
                    for ae in stratum.aes:
                        m0 = ha_arr == 0
                        ms = ae_arr == ae.num
                        ae_mask = m0 & ms     # Get material that equals zero within in architectural element
                        if np.all(ae_mask == False):
                            continue

                        ae_bg_facies = ae.bg_facies

                        # Assign background material
                        k_ae_flat = hu.specsim(self.grid,
                                            hydraulics['sig_y'][ae_bg_facies],
                                            hydraulics['ycorlengths'][ae_bg_facies], 
                                            covmod='exp',
                                            selection_mask=ae_mask)
                        k_ae_flat = np.exp(k_ae_flat) * hydraulics['k_h'][ae_bg_facies]          # back-transform from log space
                        k_iso[ae_mask] = k_ae_flat

                        # Generate internal heterogeneity - porosity
                        n_flat = hu.specsim(self.grid,
                                            hydraulics['sig_n'][ae_bg_facies],
                                            hydraulics['ncorlengths'][ae_bg_facies],
                                            covmod='exp',
                                            selection_mask=ae_mask)
                        n_flat = n_flat + hydraulics['n'][ae_bg_facies]
                        poros[ae_mask] = n_flat

                        # Assign background anisotropy ratio
                        anirat[ae_mask] = hydraulics['k_ratio'][ae_bg_facies]

                        # Assign architectural element trends
                        if ae.type_params['k_ztrend'] is not None:
                            # Vertical trend
                            zf_vec = np.linspace(*ae.type_params['k_ztrend'], self.gridnz)  # Z factor at each elevation
                            k_iso *= zf_vec
                        if ae.type_params['k_xtrend'] is not None:
                            xf_vec = np.linspace(*ae.type_params['k_xtrend'], self.grid.nx)
                            k_iso = np.transpose(k_iso.transpose() * xf_vec)

                """ Assign trends to hydraulic parameters globally """
                if hydraulics['k_ztrend'] is not None:
                    zf_vec = np.linspace(*hydraulics['k_ztrend'], self.grid.nz)  # Z factor at each elevation
                    k_iso *= zf_vec
                if hydraulics['k_xtrend'] is not None:
                    xf_vec = np.linspace(*hydraulics['k_xtrend'], self.grid.nx)
                    k_iso = np.transpose(k_iso.transpose() * xf_vec)

                if hydraulics['n_ztrend'] is not None:
                    zf_vec = np.linspace(*hydraulics['n_ztrend'], self.grid.nz)  # Z factor at each elevation
                    poros *= zf_vec
                if hydraulics['n_xtrend'] is not None:
                   xf_vec = np.linspace(*hydraulics['n_xtrend'], self.grid.nx)
                   poros = np.transpose(poros.transpose() * xf_vec)

        else:
            # Homogeneous case
            for hy_idx, hyi in enumerate(hydraulics['hydro']):
                hyi = hy_idx
                k_iso[fac == hyi] = hydraulics['k_h'][hy_idx]
                poros[fac == hyi] = hydraulics['n'][hy_idx]
                anirat[fac == hyi] = hydraulics['k_ratio'][hy_idx]

        # Assignment of anisotropy
        self.data['ktensors'] = np.zeros((self.grid.nx, self.grid.ny, self.grid.nz, 3, 3), dtype=np.float)
        ho.set_anisotropic_ktensor(self.data['ktensors'], k_iso, azim*np.pi/180, dip*np.pi/180, anirat)



    def print_model(self):
        """
        Prints a long description of the model. This is mainly for debugging.
        """

        print("Number of strata:", self.n_strata)
        print()
        print('====================================')
        print()
        for si in range(self.n_strata):
            print("Stratum:", self.strata[si].name)
            print("number of AEs:", self.strata[si].n_ae)
            print("zmin:", self.strata[si].zmin)
            print("zmax:", self.strata[si].zmax)
            print("---------------------------")
            for aei in range(self.strata[si].n_ae):
                print("AE:", self.strata[si].aes[aei].type_name)
                print("number of objects:", self.strata[si].aes[aei].n_objects)
                print("zmin", self.strata[si].aes[aei].zmin)
                print("zmax", self.strata[si].aes[aei].zmax)
                print("object_zmins:")
                print(self.strata[si].aes[aei].object_zmins)
                print("object_zmaxs:")
                print(self.strata[si].aes[aei].object_zmaxs)
                print()
            print()
            print('====================================')
            print()
