# -*- coding: utf-8 -*-
""" Hydrogeological Virtual Reality simulation package.

    Hydrogeological virtual reality (HYVR) simulator for object-based modelling
    of sedimentary structures

    :Author: Jeremy P. Bennett

    :Notes:
         Grid nodes are cell-centred!


"""

import os
import sys
import numpy as np
import pickle
import math
import time
import random
import scipy.io as sio
import matplotlib.pyplot as plt
import h5py

import hyvr.grid as gr
import hyvr.utils as hu
import hyvr.optimized as ho
import hyvr.parameters as hp


def run(param_file):
    """
    Main function for HYVR generation

    Parameters:
        param_file (str): 	Parameter file location

    Returns:

        Save data outputs as parameter file

    """

    if param_file == 0:
        param_file = os.path.abspath(os.path.join(os.path.dirname( __file__ ), os.pardir, 'testcases', 'made.ini'))
        print('param file:', param_file)

    # Load parameter file
    run, model, strata, hydraulics, flowtrans, elements, mg = hp.model_setup(param_file)

    for sim in range(1, int(run['numsim'])+1):
        # Generate facies
        props, params = facies(run, model, strata, hydraulics, flowtrans, elements, mg)

        if hydraulics['gen'] is True:
            # Generate internal heterogeneity
            props, params = heterogeneity(props, params)

        # Save data
        if run['numsim'] > 1:
            realname = 'real_{:03d}'.format(sim)
            realdir = os.path.join(run['rundir'], realname)
        else:
            realname = run['runname']
            realdir = run['rundir']

        hu.try_makefolder(realdir)

        if run['dataoutputs'] is not None:
            if hydraulics['gen'] is True:
                outdict = {'fac': props['fac'], 'azim': props['azim'], 'dip': props['dip'],
                           'k_iso': props['k_iso'], 'poros': props['poros'], 'ae': props['ae_arr'],
                           'ha': props['ha_arr'], 'hat': props['hat_arr'], 'ssm': props['ssm_arr'], 'anirat': props['anirat'],
                           'ktensors': props['ktensors'], }
            else:
                outdict = {'fac': props['fac'], 'azim': props['azim'], 'dip': props['dip'],
                           'ae': props['ae_arr'], 'ha': props['ha_arr'], 'hat': props['hat_arr'],
                           'ssm': props['ssm_arr']}
            save_outputs(realdir, realname, run['dataoutputs'], mg, outdict)

        if run['modeloutputs'] is not None:
            if hydraulics['gen'] is False:
                print('No hydraulic parameters generated. No model outputs saved')
            else:
                save_models(realdir, realname, mg, run['modeloutputs'], flowtrans,
                            props['k_iso'],
                            props['ktensors'],
                            props['poros'],
                            props['anirat'],
                            props['dip'],
                            props['azim'])


def facies(run, model, strata, hydraulics, flowtrans, elements, mg):
    """  Generate hydrofacies fields

    Parameters:
        run:					Model run parameters like ``runname``, ``rundir``, ``l_dataoutputs``, ``l_modeloutputs``, etc.
        model:					Model domain parameters
        strata:				Details about the strata
        hydraulics:				Details about the hydraulics
        flowtrans:				Flow & transport simulation parameters
        elements:				Architectural elements and parameters
        mg:						Mesh grid object class

    Returns:
        (tuple): Tuple containing:

         - probs *(dict)*: Contains data of architectural element units and associated hydrofacies

        - params *(tuple)* - Contains parameters of model domain, strata, hydraulics, etc.
            - run (dict) - Model run parameters
            - model (dict) - Model domain parameters
            - strata (dict) - Strata parameters
            - hydraulics (dict) - Hydraulic properties parameters
            - flowtrans (dict) - Flow & transport simulation parameters
            - elements (dict) - Architectural elements and parameters
            - mg - Mesh grid object class
            - ae_lu - Architectural element lookup table

    """

    """--------------------------------------------------------------------------------------------------------------
    Simulate system contacts
    --------------------------------------------------------------------------------------------------------------"""
    num_ssm = len(strata['ssm'])
    if num_ssm > 1:
        """ Create contact surfaces """
        z_bot = np.zeros((mg.nx, mg.ny))
        ssm_arr = np.zeros((mg.nx, mg.ny, mg.nz), dtype=np.int32)  # Initialise system storage array
        ssm_top_z = np.zeros((mg.nx, mg.ny, num_ssm))
        strata['ssm_bot'] = [0.0] * num_ssm
        _, _, zzz = mg.meshup()

        for si in range(num_ssm):
            if strata['ssm_contact'] == 'random' and si != num_ssm - 1:
                sp = strata['ssm_contact_model'][si]       # geostatistical parameters of system

                # Generate random top contact
                z_top = hu.specsim(mg, sp[0], [sp[1], sp[2]], twod=True, covmod='gau') + strata['ssm_top'][si]
                z_top = hu.round_x(z_top, base=model['dz'])           # round the values to grid resolution
            else:
                # Flat top contact
                z_top = np.ones((mg.nx, mg.ny)) * strata['ssm_top'][si]

            # Update lowest and highest values due to randomness
            strata['ssm_bot'][si] = np.max([np.mean(z_bot), mg.oz])
            strata['ssm_top'][si] = np.min([np.mean(z_top), mg.oz + mg.lz])

            # Assign z_bot and z_top values to entire array
            z_bot_arr = np.tile(z_bot[..., np.newaxis], [1, 1, mg.nz])
            z_top_arr = np.tile(z_top[..., np.newaxis], [1, 1, mg.nz])
            zae = np.logical_and(zzz >= z_bot_arr, zzz < z_top_arr)
            ssm_arr[zae] = si

            z_bot = z_top                   # Update lower contact surface elevation
            ssm_top_z[:, :, si] = z_top     # Assign system top to storage array

    else:
        # Only one system present
        ssm_arr = np.zeros((mg.nx, mg.ny, mg.nz), dtype=np.int32)  # Initialise system storage array
        strata['ssm_bot'] = [mg.oz]
        strata['ssm_top'] = [mg.oz + mg.lz]

    """--------------------------------------------------------------------------------------------------------------
    Simulate architectural element units
    --------------------------------------------------------------------------------------------------------------"""
    if strata['ae_table'] is not None:
        """ Load architectural element lookup table (same directory as ini parameter input file """
        ae_lu = hu.read_lu(os.path.join(run['modeldir'],strata['ae_table']))

    elif len(strata['ae_z_mean']) < 0:
        """ Uniform Model """
        ae_lu = [[0, 0, mg.lz, strata['ssm_ae'][si]]]
        ae_arr = np.ones((mg.nx, mg.ny, mg.nz), dtype=np.int32)  # Initialise system storage array

    else:
        """ Assign architectural element units """
        print(time.strftime("%d-%m %H:%M:%S", time.localtime(time.time())) + ': Generating architectural element unit contacts')

        # Initialise architectural element unit lookup table
        # [architectural element unit #, z_bottom, z_top, architectural element type, system #
        ae_lu = []
        count = 0

        if max(strata['ssm_top']) < mg.lz:
            raise ValueError('Final major strata elevation and uppermost model grid elevation do not match.\n'
                             'Please update model.lz or strata.r_ssm_top.')

        for si in range(len(strata['ssm'])):
            # Randomly assign strata / architectural element contact surfaces
            znow = strata['ssm_bot'][si]
            while znow < np.min([mg.lz, strata['ssm_top'][si]]):
                # Loop over all depths in system
                aelu_z = [count, znow, 0, 0, si]       # Initialise AE entry in lookup table (and assign identifier)

                # Assign architectural element
                aelu_z[3] = prob_choose(strata['ssm_ae'][si],
                                        strata['ae_prob'][si])

                # Assign unit thickness
                ae_z_mean = strata['ae_z_mean'][si][strata['ssm_ae'][si].index(aelu_z[3])]
                ae_z = hu.round_x(np.random.normal(ae_z_mean, ae_z_mean * 0.1), base=model['dz'])
                aelu_z[2] = min(ae_z + znow, mg.lz)

                # Assign avulsion
                avul_prob = np.array(strata['avul_prob'][si])
                yn = prob_choose([-1, 0], [avul_prob, 1 - avul_prob])         # Avulsion yes/no
                avudr = strata['avul'][si]      # Avulsion depth range for system
                dz = np.random.uniform(avudr[0], avudr[1]) * yn
                znow += ae_z + dz

                # Append to lookup table
                ae_lu.append(aelu_z)
                count += 1

    """ Create contact surfaces """
    z_bot = np.zeros((mg.nx, mg.ny))
    ae_arr = np.zeros((mg.nx, mg.ny, mg.nz), dtype=np.int32)  # Initialise system storage array
    _, _, zzz = mg.meshup()

    for ae_i, ae_z in enumerate(ae_lu):
        ae_dict = elements[ae_z[3]]                             # Get architectural element dict
        if ae_i == len(ae_lu)-1:
            # If AE  is the upper-most in the domain
            z_top = np.ones((mg.nx, mg.ny)) * mg.lz             # Assign domain top as unit top
        elif ae_lu[ae_i+1][-1] != ae_z[-1]:
            # Use the system top contact if the AE unit is the top-most in the system
            z_top = ssm_top_z[:, :, ae_z[-1]]
        elif ae_dict['contact'] is not None and ae_dict['contact'] == 'random':
            # Generate random top contact
            sp = ae_dict['contact_model']
            z_top = hu.specsim(mg, sp[0], [sp[1], sp[2]], twod=True, covmod='gau') + ae_z[2]
            z_top = hu.round_x(z_top, base=model['dz'])           # round the values to grid resolution
            ae_z[2] = np.mean(z_top)

        else:
            # Flat top contact
            z_top = np.ones((mg.nx, mg.ny)) * ae_z[2]

        # Assign z_bot and z_top values to entire array
        z_bot_arr = np.tile(z_bot[..., np.newaxis], [1, 1, mg.nz])
        z_top_arr = np.tile(z_top[..., np.newaxis], [1, 1, mg.nz])
        zae = np.logical_and(zzz >= z_bot_arr, zzz < z_top_arr)
        zae = np.logical_and(zae, ssm_arr == ae_z[-1])
        ae_arr[zae] = ae_z[0]

        # Hack to make sure erosive elements aren't simulated in strata below
        if ae_dict['geometry'] in ['trunc_ellip', 'channel', 'ext_par']:
            ae_lu[ae_i][1] = np.max(z_bot.flatten())            # Update AE lookup table with highest value
        else:
            ae_lu[ae_i][1] = np.min(z_bot.flatten())            # Update AE lookup table with lowest value

        z_bot = z_top           # Update lower contact surface elevation

    # Save system lookup table
    # if 'ae_table' not in strata:
    #     lu_savetxt = rundir + '/ae_lu_' + time.strftime('%d-%m-%Y_%H.%M.%S.txt')
    #     with open(lu_savetxt, 'w') as fwr:
    #         print('strata summary')
    #         for i in ae_lu:
    #             fwr.write('%s\n' % str()[1:-1])
    #             print(i)

    """--------------------------------------------------------------------------------------------------------------
    Hydrofacies simulation
    --------------------------------------------------------------------------------------------------------------"""
    # Initialise storage arrays
    count = 1
    hat_arr, ha_arr, fac, azim, dip = save_arrays((mg.nx, mg.ny, mg.nz), bg=strata['bg'], mat_count=count)

    """ Create architectural elements and associated hydrofacies fields """
    # Loop over AE units rather than elevations
    for ae_i in ae_lu:
        print(time.strftime("%d-%m %H:%M:%S", time.localtime(time.time())) + ': generating ' + ae_i[3] + ' from ' + str(ae_i[1]) + 'm')
        ae_dict = elements[ae_i[3]]

        if ae_dict['geometry'] == 'trunc_ellip':
            # Generate truncated ellipsoid
            props_n, count = gen_trough(ae_dict, mg, model, ae_i, ae_arr, count)
            ae_mask = props_n['ae_arr_i'] == ae_i[0]

        elif ae_dict['geometry'] in ['ext_par', 'channel']:
            # Generate extruded parabola
            props_n, count = gen_extpar(ae_dict, mg, model, ae_i, ae_arr, count)
            ae_mask = props_n['ae_arr_i'] == ae_i[0]

        elif ae_dict['geometry'] == 'sheet':
            # Generate sheet
            props_n, count = gen_sheet(ae_dict, mg, ae_i, ae_arr, count)
            ae_mask = ae_arr == ae_i[0]

        """ ADD NEW GEOMETRIES HERE """

        # Assign simulated values to storage arrays
        ae_arr[ae_mask] = ae_i[0]
        ssm_arr[ae_mask] = ae_i[4]
        hat_arr[ae_mask] = props_n['hat_arr'][ae_mask]
        ha_arr[ae_mask] = props_n['ha_arr'][ae_mask]
        fac[ae_mask] = props_n['fac'][ae_mask]
        azim[ae_mask] = props_n['azim'][ae_mask]
        dip[ae_mask] = props_n['dip'][ae_mask]

    # Renumber material values from zero to remove eroded values
    ha_arr = ho.reindex(ha_arr)

    # Wrap storage arrays in a dictionary
    if run['anisotropy']:
        props = {'azim': azim, 'ha_arr': ha_arr, 'hat_arr': hat_arr, 'dip': dip, 'fac': fac, 'ae_arr': ae_arr, 'ssm_arr': ssm_arr}
    else:
        props = {'ha_arr': ha_arr, 'fac': fac, 'hat_arr': hat_arr, 'ae_arr': ae_arr, 'ssm_arr': ssm_arr}
    params = [run, model, strata, hydraulics, flowtrans, elements, mg, ae_lu]

    return props, params


def heterogeneity(props, params):
    """
    Generate internal heterogeneity

    Parameters:
        probs (list):			Data of architectural element units and associated hydrofacies (e.g. values of azimuth, material, dipping, etc.)
        params (list):			Parameters of model domain, system, hydraulics, etc.

    Returns:
        - probs *(list)* - Data of architectural element units and associated hydrofacies
        - params *(list)* - Input parameters, assigned with heterogenity

    """

    print(time.strftime("%d-%m %H:%M:%S", time.localtime(time.time())) + ': generating hydraulic parameters')
    run, model, strata, hydraulics, flowtrans, elements, mg, ae_lu = params


    ha_arr = props['ha_arr']
    fac = props['fac']
    ae_arr = props['ae_arr']
    hat_arr = props['hat_arr']
    ssm_arr = props['ssm_arr']
    if run['anisotropy']:
        azim = props['azim']
        dip = props['dip']
    else:
        azim = np.zeros((mg.nx, mg.ny, mg.nz), dtype=np.float)
        dip = np.zeros((mg.nx, mg.ny, mg.nz), dtype=np.float)

    # Initialise storage arrays
    k_iso = np.zeros((mg.nx, mg.ny, mg.nz), dtype=np.float)       # Horizontal hydraulic conductivity array
    poros = np.zeros((mg.nx, mg.ny, mg.nz), dtype=np.float)       # Porosity array
    anirat = np.ones((mg.nx, mg.ny, mg.nz), dtype=np.float)       # K_h/K_v anisotropy ratio

    if run['het'] is True:
        # Heterogeneous case
        for mi in np.unique(ha_arr):
            for fi in np.unique(fac[ha_arr == mi]):
                mifi = (ha_arr == mi) & (fac == fi)    # Get mask for relevant values

                if model['hetlev'] == 'internal':
                    # Generate internal heterogeneity
                    # Find outer limit of facies
                    fac_idx = np.where(mifi)                                # Get indices of facies
                    fac_nx = fac_idx[0].max() - fac_idx[0].min() + 1        # Get number of grid cells in x-direction
                    fac_ny = fac_idx[1].max() - fac_idx[1].min() + 1        # Get number of grid cells in y-direction
                    fac_nz = fac_idx[2].max() - fac_idx[2].min() + 1        # Get number of grid cells in z-direction

                    # Generate field with matching size
                    # Should include a condition that considers the characteristic lengths of the features
                    temp_gr = gr.Grid(dx=mg.dx, dy=mg.dy, dz=mg.dz, nx=fac_nx, ny=fac_ny, nz=fac_nz, gtype='cells')

                    # Generate internal heterogeneity - hydraulic conductivity
                    temp_k_small = hu.specsim(temp_gr, hydraulics['sig_y'][fi], hydraulics['ycorlengths'][fi], covmod='exp')
                    temp_k_small = np.exp(temp_k_small) * hydraulics['k_h'][fi]          # back-transform from log space
                    temp_k = np.zeros((mg.nx, mg.ny, mg.nz), dtype=np.float)

                    # Nest smaller array into larger array
                    # Get coordinates for 'nesting'
                    ix1 = fac_idx[0].min()
                    ix2 = ix1 + np.shape(temp_k_small)[0]
                    iy1 = fac_idx[1].min()
                    iy2 = iy1 + np.shape(temp_k_small)[1]
                    iz1 = fac_idx[2].min()
                    iz2 = iz1 + np.shape(temp_k_small)[2]

                    if np.shape(temp_k[ix1:ix2, iy1:iy2, iz1:iz2])[0] != np.shape(temp_k_small)[0]:
                        # QnD way to avoid indexing issues with nesting of the random field
                        ix1 -= 1
                        ix2 -= 1
                    elif np.shape(temp_k[ix1:ix2, iy1:iy2, iz1:iz2])[1] != np.shape(temp_k_small)[1]:
                        # QnD way to avoid indexing issues with nesting of the random field
                        iy1 -= 1
                        iy2 -= 1
                    elif np.shape(temp_k[ix1:ix2, iy1:iy2, iz1:iz2])[2] != np.shape(temp_k_small)[2]:
                        # QnD way to avoid indexing issues with nesting of the random field
                        iz1 -= 1
                        iz2 -= 1

                    # Insert into full-size array
                    temp_k[ix1:ix2, iy1:iy2, iz1:iz2] = temp_k_small
                    k_iso[mifi] = temp_k[mifi]

                    # Generate internal heterogeneity - porosity
                    temp_n_small = hu.specsim(temp_gr, hydraulics['sig_n'][fi], hydraulics['ncorlengths'][fi], covmod='exp') + hydraulics['n'][fi]
                    temp_n = np.zeros((mg.nx, mg.ny, mg.nz), dtype=np.float)
                    # Nest smaller array into larger array
                    temp_n[ix1:ix2, iy1:iy2, iz1:iz2] = temp_n_small
                    poros[mifi] = temp_n[mifi]

                elif model['hetlev'] == 'facies':
                    # Assign heterogeneity at facies level only
                    k_iso[mifi] = hydraulics['k_h'][fi]
                    poros[mifi] = hydraulics['n'][fi]

                # Assign anisotropy ratio
                anirat[mifi] = hydraulics['k_ratio'][fi]

        """ Assign background heterogeneity per architectural element """
        if model['hetlev'] == 'internal':
            for si, aei in enumerate(ae_lu):
                m0 = ha_arr == 0
                ms = ae_arr == int(aei[0])
                aemask = m0 & ms     # Get material that equals zero within in architectural element
                if elements[aei[3]]['bg'] is not None:
                    aebackfac = int(elements[aei[3]]['bg'][0])   # architectural element background facies
                else:
                    aebackfac = int(strata['bg'][0])

                # Assign background material
                temp_k = hu.specsim(mg, hydraulics['sig_y'][aebackfac], hydraulics['ycorlengths'][aebackfac], covmod='exp')
                temp_k = np.exp(temp_k) * hydraulics['k_h'][aebackfac]          # back-transform from log space
                k_iso[aemask] = temp_k[aemask]

                # Generate internal heterogeneity - porosity
                temp_n = hu.specsim(mg, hydraulics['sig_n'][aebackfac], hydraulics['ncorlengths'][aebackfac], covmod='exp')
                temp_n = temp_n + hydraulics['n'][aebackfac]
                poros[aemask] = temp_n[aemask]

                # Assign background anisotropy ratio
                anirat[aemask] = hydraulics['k_ratio'][aebackfac]

                # Assign architectural element trends
                if elements[aei[3]]['k_ztrend'] is not None:
                    # Vertical trend
                    zf_vec = np.linspace(elements[aei[3]]['k_ztrend'][0], elements[aei[3]]['k_ztrend'][1], mg.nz)  # Z factor at each elevation
                else:
                    # Longitudinal trend
                    zf_vec = np.ones((mg.nz,))
                if elements[aei[3]]['k_xtrend'] is not None:
                    xf_vec = np.linspace(elements[aei[3]]['k_xtrend'][0], elements[aei[3]]['k_xtrend'][1], mg.nx)
                else:
                    xf_vec = np.ones((mg.nx,))
                k_iso = np.transpose(k_iso.transpose() * xf_vec)
                k_iso *= zf_vec

            """ Assign trends to hydraulic parameters globally """
            if hydraulics['k_ztrend'] is not None:
                zf_vec = np.linspace(hydraulics['k_ztrend'][0], hydraulics['k_ztrend'][1], mg.nz)  # Z factor at each elevation
            else:
                # Longitudinal trend
                zf_vec = np.ones((mg.nz,))
            if hydraulics['k_xtrend'] is not None:
                xf_vec = np.linspace(hydraulics['k_xtrend'][0], hydraulics['k_xtrend'][1], mg.nx)
            else:
                xf_vec = np.ones((mg.nx,))
            k_iso = np.transpose(k_iso.transpose() * xf_vec)
            k_iso *= zf_vec

    else:
        # Homogeneous case
        for hy_idx, hyi in enumerate(hydraulics['hydro']):
            hyi = hy_idx
            k_iso[fac == hyi] = hydraulics['k_h'][hy_idx]
            poros[fac == hyi] = hydraulics['n'][hy_idx]
            anirat[fac == hyi] = hydraulics['k_ratio'][hy_idx]

    """ Assignment of anisotropy """
    # Initialise storage arrays
    ktensors = np.zeros((mg.nx, mg.ny, mg.nz, 3, 3), dtype=np.float)

    # convert angles to radians
    azim = azim * np.pi/180
    dip = dip * np.pi/180

    # Create hydraulic conductivity tensors
    # kplane = anirat ** 0.5
    # kperp = 1 / anirat ** 0.5

    # T =========================
    # R = np.array([[np.cos(azim), np.sin(azim), 0],
    #                   [-np.sin(azim), np.cos(azim), 0],
    #                   [0, 0, 1]], dtype=np.float) * \
    #         np.array([[np.cos(dip), 0, np.sin(dip)],
    #                   [0, 1, 0],
    #                   [-np.sin(dip), 0, np.cos(dip)]], dtype=np.float)
    # /T =========================

    # Iterate over all nodes
    ho.set_anisotropic_ktensor(ktensors, k_iso, azim, dip, anirat)

    # convert radians to angles
    azim = azim * 180/np.pi
    dip = dip * 180/np.pi

    props = {'azim': azim, 'ha_arr': ha_arr, 'dip': dip, 'fac': fac, 'ae_arr': ae_arr, 'ssm_arr': ssm_arr,
             'k_iso': k_iso, 'ktensors': ktensors, 'poros': poros, 'anirat': anirat, 'hat_arr': hat_arr}

    return props, params


def save_outputs(realdir, realname, outputs, mg, outdict):
    """
    Save data arrays to standard formats

    Parameters:
        realdir (str):  	File path to save to
        realname (str):  	File name
        outputs (str):		String codes for what type of outputs should be saved
        mg:					Mesh grid object class
        outdict:			Output directory

    Returns:
        Save data outputs as .vtk (Paraview), .mat (Matlab) or .dat (Python pickle output)

    """

    print('Saving files in {}'.format(realdir))
    realname = realname + '_hyvr'
    if 'vtk' in outputs:
        # VTK output for visualisation in ParaView
        hu.to_vtr({k: outdict[k] for k in outdict if k not in ['ktensors']},
                    os.path.join(realdir, realname), mg)
        print(time.strftime("%d-%m %H:%M:%S", time.localtime(time.time())) + ': VTR export complete')

    if 'mat' in outputs:
        # MATLAB output
        sio.savemat(os.path.join(realdir, realname), outdict)
        print(time.strftime("%d-%m %H:%M:%S", time.localtime(time.time())) + ': Matlab export complete')

    if 'py' in outputs:
        # Python pickle output
        with open(os.path.join(realdir, realname + '.dat'), 'wb') as outfile:
            pickle.dump(outdict, outfile, protocol=pickle.HIGHEST_PROTOCOL)
        print(time.strftime("%d-%m %H:%M:%S", time.localtime(time.time())) + ': Python pickle export complete')

    if 'npz' in outputs:
        # Numpy .npz
        np.savez_compressed(os.path.join(realdir, realname + '.npz'),
                            ae=outdict['ae'],
                            anirat=outdict['anirat'],
                            azim=outdict['azim'],
                            dip=outdict['dip'],
                            fac=outdict['fac'],
                            ha=outdict['ha'],
                            hat=outdict['hat'],
                            k_iso=outdict['k_iso'],
                            ssm=outdict['ssm'],
                            poros=outdict['poros'],
                            ktensors=outdict['ktensors'])
        print(time.strftime("%d-%m %H:%M:%S", time.localtime(time.time())) + ': Numpy export complete')

    if 'h5' in outputs:
        # HDF5 format
        with h5py.File(os.path.join(realdir, realname + '.h5'), 'w') as hf:
            print('G')
            for key in outdict.keys():
                hf.create_dataset(key, data=outdict[key])
        print(time.strftime("%d-%m %H:%M:%S", time.localtime(time.time())) + ': HDF5 export complete')


def save_models(realdir, realname, mg, outputs, flowtrans, k_iso, ktensors, poros, anirat, dip, azim):
    """
    Save HYVR outputs to standard modelling codes

    Parameters:
        run (dict):			Model run parameters
        mg:					Mesh grid object class
        flowtrans (dict):	Flow & transport simulation parameters
        k_iso:				Horizontal hydraulic conductivity array
        ktensors:			Array with tensor values of K
        poros:				Porosity array
        anirat:				Anistropic ratio ($K_h/K_v$)

    Returns:
        Save data outputs as .mf (MODFLOW) or .hgs (HydroGeoSphere)

    """
    realname = realname + '_hyvr'
    if 'mf' in outputs:
        # MODFLOW output

        # Create HGS output folder
        mfdir = os.path.join(realdir,'MODFLOW')
        mfname = os.path.join(mfdir, realname)
        hu.try_makefolder(mfdir)
        hu.to_modflow(mfname, mg, flowtrans, k_iso, anirat)
        print(time.strftime("%d-%m %H:%M:%S", time.localtime(time.time())) + ': MF export complete')

    if 'mf6' in outputs:
        # MODFLOW 6 output

        # Create HGS output folder
        mf6dir = os.path.join(realdir, 'mf6/hyvr')
        mf6name = realname
        if not os.path.exists(mf6dir):
            os.makedirs(mf6dir)

        hu.to_mf6(mf6dir, mf6name, mg, flowtrans, k_iso, anirat, dip, azim)
        print(time.strftime("%d-%m %H:%M:%S", time.localtime(time.time())) + ': MF6 export complete')

    if 'hgs' in outputs:
        # HydroGeoSphere output
        # Create HGS output folder
        hgsdir = os.path.join(realdir, 'HGS')
        hu.try_makefolder(hgsdir)

        # Write to HGS files
        hu.to_hgs(hgsdir, mg, flowtrans, ktensors, poros)
        print(time.strftime("%d-%m %H:%M:%S", time.localtime(time.time())) + ': HGS export complete')


"""-------------------------------------------------------------------------------------------------------------- """


"""--------------------------------------------------------------------------------------------------------------
Trough generators and utilities
--------------------------------------------------------------------------------------------------------------"""


def gen_trough(tr, mg, model, ae, ae_arr, count, ani=True):
    """
    Create trough shapes

    Parameters:
        tr (dict): 				Trough parameters
        mg (grid class): 		Model grid
        ae (list): 				Architectural element unit details
        ae_arr (ndarray): 		3D array of system numbers
        count (int): 			Material number and/or identifier
        ani (bool):         	Boolean if anisotropy is to be generated

    Returns:
        - probs *(numpy array)* - Grid properties
        - count *(int)* - Material number and/or identifier

    """

    x3, y3, z3 = mg.meshup()    # 3-D grid

    if ani:
        hat_arr, ha_arr, fac, azim, dip = save_arrays((mg.nx, mg.ny, mg.nz), mat_count=count, bg=tr['bg'])
    else:
        hat_arr, ha_arr, fac = save_arrays((mg.nx, mg.ny, mg.nz), mat_count=count, bg=tr['bg'], ani=False)
    count += 1

    ae_arr_i = np.zeros((mg.nx, mg.ny, mg.nz), dtype=int)

    # Assign background values
    ae_arr_i[ae_arr == ae[0]] = ae[0]

    # loop over trough top depths
    if tr['buffer'] is not None:
        tr_bot = ae[1] + tr['depth'] * tr['buffer']
    else:
        tr_bot = ae[1]
    tr_top = max(ae[2], tr_bot) + tr['agg']

    if tr['te_xyz'] is not None:
        gen_elevations = [xyz[2] for xyz in tr['te_xyz']]
        trough_gen_range = [1]
    else:
        trough_gen_range = np.arange(0, tr['el_z'] * mg.lx * mg.ly)         # Range of truncated ellipsoid to generate at each elevation
        gen_elevations = np.arange(tr_bot, tr_top, tr['agg'])               # aggradation elevations

    tr_xy_t = np.zeros((len(trough_gen_range), 2))                      # X,Y coordinates of truncated ellipsoid

    for zidx, znow in enumerate(gen_elevations):                                         # Loop over aggradation elevations
        for eli, elno in enumerate(trough_gen_range):                   # Loop over troughs
            # Reneration of trough parameters
            a, b, c = rand_trough(tr, mg=mg, ztr=znow)

            xnow = tr_xy_t[eli][0]
            ynow = tr_xy_t[eli][1]

            # center of trough
            if model['display']:
                # Place troughs in domain centre for display features
                xnow = mg.lx / 2
                ynow = 0
            elif tr['migrate'] is not None and znow > tr_bot:
                # Migration of troughs
                xnow += np.random.normal(tr['migrate'][0], tr['migrate'][1])
                ynow += np.random.normal(tr['migrate'][2], tr['migrate'][3])
            elif tr['te_xyz'] is not None:
                xnow = tr['te_xyz'][zidx][0]
                ynow = tr['te_xyz'][zidx][1]
            else:
                xnow = np.random.uniform(0, mg.lx)
                ynow = np.random.uniform(mg.ly/-2, mg.ly/2)

            tr_xy_t[eli][0] = xnow
            tr_xy_t[eli][1] = ynow

            alpha = np.random.uniform(tr['paleoflow'][0], tr['paleoflow'][1])   # orientation angle of trough ('paleoflow')
            angnow = np.random.uniform(tr['azimuth'][0], tr['azimuth'][1])      # orientation of material

            # Distances to ellipsoid centre
            xd = x3 - xnow
            yd = y3 - ynow
            zd = z3 - znow

            # Periodic boundary
            if model['periodic'] is True:
                xd[xd > mg.lx / 2] -= mg.lx
                xd[xd < -mg.lx / 2] += mg.lx
                yd[yd > mg.ly / 2] -= mg.ly
                yd[yd < -mg.ly / 2] += mg.ly
                zd[zd > mg.lz / 2] -= mg.lz
                zd[zd < -mg.lz / 2] += mg.lz

            # scaled and rotated distance squared
            select, R2 = ho.scale_rotate(xd, yd, zd, alpha, a, b, c)
            select = np.logical_and(select, ae_arr <= ae[0])                # Restrict selection to AE units equal or below current

            """" Assign internal structure """
            tr_struct = tr['structure']
            if tr_struct == 'random':
                tr_struct = random.choice(['dip', 'bulb_l'])

            if np.all(np.isnan(select)):
                # Skip section if no grid cells selected
                pass
            if ~np.any(select):
                # Skip section if no grid cells selected
                pass

            if model['hetlev'] == 'ha':
                # Add 'dip layers' into trough
                fac_now = random.choice(tr['facies'])
                fac[select] = fac_now
                ha_arr[select] = count
                hat_arr[select] = tr['ae_id']
                if ani:
                    azim[select] = angnow    # Save angle
                    dip[select] = np.random.uniform(tr['dip'][0], tr['dip'][1])                    # Assignment of architectural elements only
            elif tr_struct == 'bulb':
                """
                Add 'bulb' layers into trough
                    - Dip is derived from the gradient of the truncated ellipsoid boundary
                    - Azimuth is the angle of the ellipsoid
                """
                # Generate gradient information
                dip_tr, azim_tr = ellipsoid_gradient(xd, yd, zd, a, b, c, alpha, select, tr)

                # Assign generated values to grid cells
                fac_now = random.choice(tr['facies'])
                fac[select] = fac_now
                ha_arr[select] = count
                hat_arr[select] = tr['ae_id']
                if ani:
                    dip[select] = dip_tr[select]
                    azim[select] = azim_tr[select]

            elif tr_struct == 'bulb_l':

                # Ellipsoid 'c' radii
                c_range = c - np.arange(0, c, tr['bulbset_d'])

                # Iterate over internal truncated ellipsoids
                for c_idx, c_now in enumerate(c_range):
                    # Get scale factor for ellipsoids -
                    te_scale = c_now / c
                    a_now, b_now = np.array([a, b]) * te_scale

                    # Internal scaled and rotated distance squared
                    bulb_select, bulb_R2 = ho.scale_rotate(xd, yd, zd, alpha, a_now, b_now, c_now)

                    # Generate gradient information
                    dip_bulb, azim_bulb = ellipsoid_gradient(xd, yd, zd, a_now, b_now, c_now, alpha, bulb_select, tr)

                    # Assign generated values to grid cells
                    if c_idx == 0:
                        # Randomly choose facies
                        fac_now = random.choice(tr['facies'])
                    else:
                        # Choose next hydrofacies from alternating sets
                        pf_i = [i for i, x in enumerate(tr['facies']) if x == fac_now][0]    # Get facies index
                        fac_now = random.choice(tr['altfacies'][pf_i])                   # Get next alternating facies

                    fac[bulb_select] = fac_now                                      # Alternating facies
                    ha_arr[bulb_select] = count
                    hat_arr[bulb_select] = tr['ae_id']

                    if ani:
                        dip[bulb_select] = dip_bulb[bulb_select]
                        azim[bulb_select] = azim_bulb[bulb_select]

            elif tr_struct == 'dip':
                # Add 'dip layers' into trough, with alternating facies
                do, fd, dv, av = dip_sets(mg, tr, znow, select=select, azimuth_z=alpha)

                # Assign generated values to grid cells
                fac[select] = fd[select]
                ha_arr[select] = count
                hat_arr[select] = tr['ae_id']
                if ani:
                    dip[select] = dv
                    azim[select] = av

            else:
                # Add 'dip layers' into trough
                fac[select] = random.choice(tr['facies'])
                ha_arr[select] = count
                hat_arr[select] = tr['ae_id']
                if ani:
                    azim[select] = angnow    # Save angle
                    dip[select] = np.random.uniform(tr['dip'][0], tr['dip'][1])

            if tr['lag'] is not None:
                in_lag = (znow - tr['depth'] + float(tr['lag'][0])) > z3  # Is grid cell below top of lag surface
                fac[np.logical_and(select, in_lag)] = int(tr['lag'][1])

            count += 1
            ae_arr_i[select] = ae[0]

    if ani:
        props = {'ha_arr': ha_arr, 'azim': azim, 'dip': dip, 'fac': fac, 'ae_arr_i': ae_arr_i, 'hat_arr': hat_arr}
    else:
        props = {'ha_arr': ha_arr, 'fac': fac, 'ae_arr_i': ae_arr_i, 'hat_arr': hat_arr}

    return props, count


def scale_rotate(x, y, z, alpha=0, a=1, b=1, c=1):
    """
    Scale and rotate three-dimensional trough

    Parameters:
        x, y, z (float):	Spatial coordinates
        alpha (float):		Rotation angle about the z-axis
        a, b, c (float):	Axis lengths in x, y, z directions (ellipsoid length, width, depth)

    Returns:
        - select - Grid cells within ellipsoid
        - R2 - Grid of scaled and rotated values

    """




    alpha = np.radians(alpha)
    R2 = (x ** 2 * np.cos(alpha) ** 2 +
          2 * x * y * np.cos(alpha) * np.sin(alpha) +
          y ** 2 * np.sin(alpha) ** 2) / a ** 2 + \
         (x ** 2 * np.sin(alpha) ** 2 -
          2 * np.multiply(x, y) * np.cos(alpha) * np.sin(alpha) +
          np.power(y, 2) * np.cos(alpha) ** 2) / b ** 2 + \
          z ** 2 / c ** 2

    #  selection of cells
    mask1 = R2 <= 1
    mask2 = z <= 0
    select = mask1 & mask2

    return select, R2


def ellipsoid_gradient(x, y, z, a, b, c, alpha, select, tr):
    """
    Calculate dip and strike values in rotated ellipsoids

    Parameters:
        x, y, z:    	Distances to centre of ellipsoid
        a, b, c:    	Majox/minor axes of ellipsoid
        alpha:      	Rotation of ellipsoid from mean flow direction

    Returns:
        - dip_g - Dipping in 3D
        - azimuth_g - Azimuth in 3D

    """

    alpha = np.radians(-alpha)  # Convert alpha to radians

    # initialize arrays
    dip_g = np.zeros(np.shape(x))
    azimuth_g = np.zeros(np.shape(x))

    try:
        idx_z = np.where(select)[2].max()     # Find the 'surface cells' of the ellipsoid
    except ValueError:
        # Return if no values selected
        return dip_g, azimuth_g

    # Calcuate dip and strike for onion
    select_z_idx = np.where(select[:, :, idx_z])                    # Indices of grid cells in unit at top of unit
    ix = x[select_z_idx[0], select_z_idx[1], idx_z].flatten()
    iy = y[select_z_idx[0], select_z_idx[1], idx_z].flatten()
    iz = (1 - ((ix * np.cos(alpha) + iy * np.sin(alpha)) ** 2 / a ** 2 + (ix * np.sin(alpha) + iy * np.cos(alpha)) ** 2 / b ** 2) ** 0.5) * c

    # Get tangent plane coefficients
    fx = (2 * np.cos(alpha) ** 2 * ix + 2 * iy * np.cos(alpha) * np.sin(alpha)) / a ** 2 \
        + (2 * np.sin(alpha) ** 2 * ix + 2 * iy * np.cos(alpha) * np.sin(alpha)) / b ** 2
    fy = (2 * np.sin(alpha) ** 2 * iy + 2 * ix * np.cos(alpha) * np.sin(alpha)) / a ** 2 \
        + (2 * np.cos(alpha) ** 2 * iy + 2 * ix * np.cos(alpha) * np.sin(alpha)) / b ** 2
    fz = 2 * iz / c

    # Normal vectors of tangent plane, horizontal plane, vertical plane
    n_tan = np.array([fx, fy, fz]).T
    n_horizontal = np.array([0., 0., 1.])

    # Calculate the dip at each point
    dip_vec = np.minimum(tr['dip'][1], (angle(n_tan, n_horizontal) * 180/np.pi))

    # Insert into 2D array
    dip2d = np.zeros(np.shape(x)[0:2])
    dip2d[select_z_idx[0], select_z_idx[1]] = dip_vec
    dip2d = dip2d[:, :, None] * np.ones(np.shape(x))

    # Apply to 3D arrays
    dip_g[select] = dip2d[select]
    azimuth_g[select] = alpha * 180/np.pi

    return dip_g, azimuth_g


def rand_trough(tr, mg=False, ztr=[]):
    """
    Randomly generate ellipsoid geometry parameters:

    Parameters:
        tr:     	Ellipsoid parameters
        mg:     	Meshgrid object
        ztr:    	Elevation of ellipsoid centre point

    Returns:
        a, b, c - Length, width and depth of ellipsoid

    """
    if ztr and tr['geo_ztrend'] is not None:
        # zfactor = np.interp(ztr, [mg.oz, mg.oz + mg.lz], [tr['geo_ztrend'][0], tr['geo_ztrend'][1]])
        zfactor = tr['geo_ztrend'][0] + (ztr - mg.oz)*(tr['geo_ztrend'][1]-tr['geo_ztrend'][0])/(mg.lz)
    else:
        zfactor = 1

    a = tr['length'] * zfactor / 2
    b = tr['width'] * zfactor / 2
    c = tr['depth'] * zfactor

    return a, b, c

"""--------------------------------------------------------------------------------------------------------------
Extruded parabola generators and utilities
--------------------------------------------------------------------------------------------------------------"""


def gen_extpar(ch_par, mg, model, ssm, ae_array, count, ani=True):
    """
    Generate extruded parabola geometries:
        - Flow regime is assumed to be reasonably constant so the major geometry of the extruded parabolas doesn't change so much
        - 'Migration' of the extruded parabolas according to a shift vector

    Parameters:
        ch_par:         Extruded parabola parameters
        mg:         	Mesh grid object class
        model (dict):	Model domain parameters
        ssm (dict):		Strata parameters
        ae_array:		Array with architectural element unit details
        count (int): 	Material number and/or identifier
        ani (bool): 	Boolean if anisotropy is to be generated
        (z_in:       	starting depth)
        (thickness:  	Thickness of architectural element)

    Returns:
        (tuple): Tuple containing:

        - probs *(dict)*: Contains data of architectural element units and associated hydrofacies
            - mat - Material values
            - azim - Azimuth angles (ani == True)
            - dip - Dipping angles (ani == True)
            - fac - Facies values
            - ae_arr_i - Array with architectural element unit details

        - count *(int)* - Material number and/or identifier

    """

    # Vectors of spatial coordinates
    xvec, yvec, zvec = mg.vec()
    x2, y2 = np.meshgrid(xvec, yvec, indexing='ij')          # 2-D grid
    _, _, z3 = np.meshgrid(range(0, mg.nx), range(0, mg.ny), range(0, mg.nz), indexing='ij')          # 2-D grid

    # Initialize storage arrays
    if ani:
        hat_arr, ha_arr, fac, azim, dip = save_arrays((mg.nx, mg.ny, mg.nz), bg=ch_par['bg'], mat_count=count)
    else:
        hat_arr, ha_arr, fac = save_arrays((mg.nx, mg.ny, mg.nz), bg=ch_par['bg'], mat_count=count, ani=False)

    ae_arr_i = np.zeros(np.shape(ae_array), dtype=int)
    ae_arr_i[ae_array == ssm[0]] = ssm[0]

    # start location
    total_extpar = int(ch_par['channel_no'])
    if model['display']:
        # Place troughs in domain centre for display features
        xstart = np.random.uniform(0, 0, total_extpar)
        ystart = np.random.uniform(yvec[0], yvec[-1], total_extpar)
    else:
        # Randomly place curve starting points
        xstart = np.random.uniform(-10, 0, total_extpar)
        ystart = np.random.uniform(yvec[0], yvec[-1], total_extpar)

    # loop over extruded parabola top depths
    ch_bot = ssm[1]
    if ch_par['buffer'] is not None:
        ch_bot += ch_par['depth'] * ch_par['buffer']

    ch_top = ssm[2]
    ch_dz = ch_par['agg']
    ch_top += ch_dz

    for znow in np.arange(ch_bot, ch_top, ch_dz):
        print(time.strftime("%d-%m %H:%M:%S", time.localtime(time.time())) + ' z = ' + str(znow))
        # Assign linear trend to extruded parabola sizes
        if ch_par['geo_ztrend'] is not None:
            zfactor = np.interp(znow, [mg.oz, mg.oz + mg.lz], [ch_par['geo_ztrend'][0], ch_par['geo_ztrend'][1]])
        else:
            zfactor = 1
        z_ch_width = ch_par['width'] * zfactor
        z_ch_depth = ch_par['depth'] * zfactor

        if 'ch_start' in ch_par:
            cstart = ch_par['ch_start']
        else:
            cstart = []
        # Loop over total extruded parabolas per system
        for chan in range(0, total_extpar):
            """ Loop over multiple extruded parabolas at 'timestep' """
            aha = ferguson_curve(mg, ch_par['h'], ch_par['k'],  ch_par['ds'],
                                 ch_par['eps_factor'],
                                 disp=model['display'],
                                 ch_start=cstart)

            """ Get flow direction in extruded parabolas for azimuth """
            # For periodicity shift trajectories into model unit cell
            if model['periodic'] is True:
                aha[aha[:, 1] < yvec[0], 1] += mg.ly
                aha[aha[:, 1] > yvec[-1], 1] -= mg.ly

            # initialize 2-D distance matrix
            D = 1e20 * np.ones_like(x2)

            # initialize sum of inverse-distance weights
            sumW = np.zeros(np.shape(x2), dtype=float)
            W = np.zeros(np.shape(x2), dtype=float)

            # initialize velocity orientation at this level
            vx_znow = np.zeros(np.shape(x2), dtype=float)
            vy_znow = np.zeros(np.shape(x2), dtype=float)

            # loop over all points of trajectory
            for ii in range(0, len(aha)):
                # distance to current point
                R = np.sqrt((x2 - aha[ii][0]) ** 2 + (y2 - aha[ii][1]) ** 2)

                # smallest distance of entire grid to all points so far
                D[R < D] = R[R < D]

                # inverse-distance weight for velocity interpolation
                W[:] = 1e-20
                W[R < z_ch_width / 2] = 1 / (R[R < z_ch_width / 2] + 1e-20)

                # velocity interpolation in 2-D
                vx_znow += aha[ii][2] * W
                vy_znow += aha[ii][3] * W
                sumW += W

            vx_znow /= sumW
            vy_znow /= sumW

            # Assign facies sets with dip values
            if sum(ch_par['dip']) > 0:
                do, fd, dv, av = dip_sets(mg, ch_par, znow, curve=[aha[:, 0], aha[:, 1], vx_znow, vy_znow])
            else:
                do = np.ones((mg.nx, mg.ny, mg.nz)) + count
                fd = np.ones((mg.nx, mg.ny, mg.nz), dtype=int) * int(random.choice(ch_par['facies']))
                dv = 0.0
                av = np.zeros((mg.nx, mg.ny, mg.nz))

            """ Copy results into 3-D field """
            # Iterate over all nodes below current top elevation
            d_range = np.arange(max(0, mg.idx_z(znow - z_ch_depth)), min(mg.nz, mg.idx_z(znow)))      # Depth range
            if len(d_range) > 0:        # Only compute if extruded parabola depth range is finite

                # Get mask arrays for each condition
                in_extpar = D[:, :, None]**2 <= z_ch_width**2 / 4 - ((mg.idx_z(znow) - z3) * mg.dz * z_ch_width / (z_ch_depth*2)) ** 2     # is grid cell in extruded parabola
                finite_v = ~np.isnan(vx_znow)            # Only assign if velocity is finite
                below_top = ae_array <= ssm[0]          # Don't assign values to locations higher than top contact surface
                chan_mask = in_extpar * finite_v[:, :, None] * below_top
                if mg.idx_z(znow) <= z3[:, :, -1].max():
                    # Set mask above top of extruded parabola to False
                    chan_mask[:, :, mg.idx_z(znow):-1] = False

                # Assign properties
                fac[chan_mask] = fd[chan_mask]
                ha_arr[chan_mask] = count
                hat_arr[chan_mask] = ch_par['ae_id']
                ae_arr_i[chan_mask] = ssm[0]

                if ch_par['lag'] is not None:
                    in_lag = (znow - z_ch_depth + float(ch_par['lag'][0])) > z3 * mg.dz   # Is grid cell in extruded parabola
                    fac[np.logical_and(in_extpar, in_lag)] = int(ch_par['lag'][1])

                if ani:
                    # calcuate azimuth, to 1 degree
                    azim2d = np.round((np.arctan2(vx_znow, vy_znow) - np.pi/2) * 180/np.pi)
                    azim3d = azim2d[:, :, None] * np.ones((mg.nx, mg.ny, mg.nz))
                    azim[chan_mask] = azim3d[chan_mask]
                    dip[chan_mask] = dv

                count += 1

        # Shift starting values with migration vector from parameter file
        if 'mig' in ch_par:
            xstart += np.random.uniform(-ch_par['mig'][0], ch_par['mig'][0])
            ystart += np.random.uniform(-ch_par['mig'][1], ch_par['mig'][1])

    if ani:
        props = {'ha_arr': ha_arr, 'azim': azim, 'dip': dip, 'fac': fac, 'ae_arr_i': ae_arr_i, 'hat_arr': hat_arr}
    else:
        props = {'ha_arr': ha_arr, 'fac': fac, 'ae_arr_i': ae_arr_i, 'hat_arr': hat_arr}
    return props, count


def ferguson_curve(mg, h, k, ds, eps_factor, dist=0, disp=False, ch_start=[]):
    """
    Simulate extruded parabola centrelines using the Ferguson (1976) disturbed meander model
    Implementation of AR2 autoregressive model
    http://onlinelibrary.wiley.com/doi/10.1002/esp.3290010403/full

    Parameters:
        mg (object class):			Mesh grid object class
        h (float):					Height
        k (float):					Wave number
        ds(float):					Curve distance for calculations
        eps_factor (float):			Random background noise
        dist (float): 		        Distance to generate curves - defaults to mg.lx
        disp (bool): 		        Creating display extruded parabola - extruded parabola begins at (0,0)
        ch_start (tuple):           Starting location of channel (x,y coordinates)

    Returns:
        outputs (float array):		Simulated extruded parabola centerlines: storage array containing values for x coordinate, y coordinate, vx and vy

    """
    # Parameters
    ds += 1e-10
    if dist > 0:
        ns = dist
    else:
        ns = mg.lx*100
    s = np.arange(0, ns, ds)

    # Centreline starting point
    xp = 0
    yp = 0

    # Calculate curve directions
    theta = ferguson_theta(s, eps_factor, k, h)

    # Interpolate curve direction over interval of interest
    s_interp, th_interp = curve_interp(s, theta, 0.1)

    # Storage array
    outputs = np.zeros((len(th_interp), 4))

    for th_idx, th_i in enumerate(th_interp):
        vx = ds*np.cos(th_i)
        vy = ds*np.sin(th_i)
        xp += vx
        yp += vy

        # Assign to storage array
        outputs[th_idx, 0] = xp       # x coordinate
        outputs[th_idx, 1] = yp       # y coordinate
        outputs[th_idx, 2] = vx       # vx
        outputs[th_idx, 3] = vy       # vy

    # Rotate meanders into mean flow direction
    mean_th = -np.mean(th_interp)
    rotMatrix = np.array([[np.cos(mean_th), -np.sin(mean_th)],
                          [np.sin(mean_th),  np.cos(mean_th)]])
    roro = np.dot(rotMatrix, outputs[:, 0:2].transpose())

    outputs[:, 2:] = np.dot(rotMatrix, outputs[:, 2:].transpose()).transpose()

    # Move starting location in x-direction
    outputs[:, 0] = roro[0, :].transpose() - np.random.uniform(mg.lx/5, mg.lx/10)
    outputs[:, 1] = roro[1, :].transpose()

    # Remove values before model domain
    if dist > 0:
        indomain = outputs[:, 0] >= mg.ox
    else:
        indomain = np.logical_and(outputs[:, 0] >= mg.ox, outputs[:, 0] <= mg.lx)
    outputs = outputs[indomain, :]

    if len(ch_start) > 0:
        outputs[:, 0] += ch_start[0]
        outputs[:, 1] += ch_start[1]

    # Make sure streamlines begin within domain with respect to y
    yout = outputs[0, 1] > mg.ly/4 or outputs[0, 1] > mg.ly/4
    if disp is True:
        starty = outputs[0, 1]
    elif yout is True:
        starty = np.random.uniform(-mg.ly/4, mg.ly/4)
    else:
        starty = 0
    outputs[:, 1] = outputs[:, 1] - starty

    return outputs


def ferguson_theta(s, eps_factor, k, h):
    """
    Calculate curve direction angle

    Parameters:
        s:				    Steps within generated curve distance
        eps_factor:		    Random background noise
        k:				    Wave number
        h:				    Height

    Returns:
        th_store *(array)* - Curve direction angle

    """
    # Storage arrays
    th_store = np.zeros(len(s))

    for idex, si in enumerate(s):
        if idex == 0:
            t1 = 0
            t2 = 0
            eps = 0
        elif idex == 1:
            t1 = th_store[idex-1]
            t2 = 0
            eps = np.random.normal()*eps_factor
        else:
            t1 = th_store[idex-1]
            t2 = th_store[idex-2]
            eps = np.random.normal(1)*eps_factor

        th_store[idex] = thetaAR2(t1, t2, k, h, eps)

    return th_store


def thetaAR2(t1, t2, k, h, eps):
    """
    Implementation of AR2 autoregressive model (Ferguson, 1976, Eq.15)
    http://onlinelibrary.wiley.com/doi/10.1002/esp.3290010403/full

    Parameters:
        t1: 	theta(i-1)
        t2: 	theta(i-2)
        k:		Wavenumber
        h:		Height
        eps:	Random background noise

    Returns:
        2nd-order autoregression (AR2)

    """
    b1 = 2*np.exp(-k*h)*np.cos(k*np.arcsin(h))
    b2 = -np.exp(-2*k*h)
    return eps + b1*t1 + b2*t2


"""--------------------------------------------------------------------------------------------------------------
Sheet generators and utilities
--------------------------------------------------------------------------------------------------------------"""


def gen_sheet(sh, mg, ae_i, ae_array, count, ani=True):
    """
    Generate gravel sheet with internal heterogeneity

    Parameters:
        sh:         	Sheet parameters
        mg:         	Model grid class
        ae_i:       	Architectural element lookup details [system number, z_bottom, z_top, architectural element, geometry]
        ae_array:   	Architectural element array
        count (int): 	Material number and/or identifier
        ani (bool):		Boolean if anisotropy is to be generated

    Returns:
        - probs *(dict)* - Contains data of architectural element units and associated hydrofacies (e.g. values of azimuth, material, dipping, etc.)
        - count *(int)* - Material number and/or identifier

    """
    # Initialize storage arrays
    if ani:
        hat_arr, ha_arr, fac, azim, dip = save_arrays((mg.nx, mg.ny, mg.nz))
    else:
        hat_arr, ha_arr, fac = save_arrays((mg.nx, mg.ny, mg.nz), ani=False)
    hat_arr[ae_array == ae_i[0]] = sh['ae_id']

    # Massive bedding -----------------------------------
    if sh['lens_thickness'] == -1:
        count += 1
        # Generate dip
        if sh['dip'] is not None and max(sh['dip']) != 0:
            if len(sh['facies']) > 1:
                do, fd, dv, av = dip_sets(mg, sh, ae_i[1])               # Generate facies sets
                fac[ae_array == ae_i[0]] = fd[ae_array == ae_i[0]]
                ha_arr[ae_array == ae_i[0]] = do[ae_array == ae_i[0]] + count
                count += len(do)
                if ani:
                    azim[ae_array == ae_i[0]] = 0
                    dip[ae_array == ae_i[0]] = dv
            else:
                fac[ae_array == ae_i[0]] = sh['facies'][0]
                ha_arr[ae_array == ae_i[0]] =  count
                count += 1
                if ani:
                    azim[ae_array == ae_i[0]] = np.random.uniform(sh['azimuth'][0], sh['azimuth'][1])
                    dip[ae_array == ae_i[0]] = np.random.uniform(sh['dip'][0], sh['dip'][1])

        elif sh['azimuth'] is not None and max(sh['azimuth']) != 0:
            fac[ae_array == ae_i[0]] = sh['facies'][0]
            ha_arr[ae_array == ae_i[0]] =  count
            count += 1
            if ani:
                azim[ae_array == ae_i[0]] = np.random.uniform(sh['azimuth'][0], sh['azimuth'][1])
                dip[ae_array == ae_i[0]] = np.random.uniform(sh['dip'][0], sh['dip'][1])

        else:
            # No dip
            ha_arr[ae_array == ae_i[0]] = count
            fac[ae_array == ae_i[0]] = random.choice(sh['facies'])
            if ani:
                azim[ae_array == ae_i[0]] = 0
                dip[ae_array == ae_i[0]] = 0

    # Create lenses over depths ------------------------------
    else:
        # Assign lens thickness for system
        if sh['geo_ztrend'] is not None:
            zfactor = np.interp(np.mean(ae_i[1:3]), [mg.oz, mg.oz + mg.lz], [sh['geo_ztrend'][0], sh['geo_ztrend'][1]])
            z_lens_thick = sh['lens_thickness'] * zfactor
        else:
            z_lens_thick = sh['lens_thickness']
        z_lens = np.arange(ae_i[1], ae_i[2]*1.1, z_lens_thick)          # Buffer added to top elevation to avoid non-assignment

        # Loop over lenses
        for znow in z_lens:
            count += 1
            z_bottom = mg.idx_z(znow)
            z_top = min(mg.idx_z(znow + sh['lens_thickness']), mg.nz)
            z_range = range(z_bottom, z_top)

            # Generate dip
            if np.array(sh['dip']).ptp() > 0:
                do, fd, dv, av = dip_sets(mg, sh, znow)               # Generate facies sets

                # Iterate over all nodes - Brute force approach :(
                it = np.nditer(ae_array, flags=['multi_index'])
                while not it.finished:
                    if it.multi_index[2] in z_range and ae_array[it.multi_index] == ae_i[0]:
                        fac[it.multi_index] = fd[it.multi_index]
                        ha_arr[it.multi_index] = do[it.multi_index] + count

                    it.iternext()

                count += len(np.unique(do))
                if ani:
                    azim[:, :, z_range] = av
                    dip[:, :, z_range] = dv   # Assign facies sets to storage arrays
            else:
                fac[:, :, z_range] = random.choice(sh['facies'])
                ha_arr[:, :, z_range] = count
                if ani:
                    azim[:, :, z_range] = 0
                    dip[:, :, z_range] = 0

    if ani:
        props = {'ha_arr': ha_arr, 'azim': azim, 'dip': dip, 'fac': fac, 'hat_arr': hat_arr}
    else:
        props = {'ha_arr': ha_arr, 'fac': fac, 'hat_arr': hat_arr}
    return props, count


def dip_sets(mg, aep, znow, curve=[], select=[], azimuth_z=0):
    """
    Generate dip angles and assign to the dip matrix

    Parameters:
        mg:         Mesh grid object class
        aep:        Architectural element parameters (dict)
        cruve:    Tuple of x,y coordinates of extruded parabola (omitted for linear flows)
                        - x, y coordinates of extruded parabola
                        - vx, vy of extruded parabola flow
        select:     Model grid nodes to assign

    Returns:
        - dip_out - Array of assigned dip values
        - fac_out - Array of assigned hydrofacies

    """

    # Vectors of spatial coordinates of grid
    # xgvec, ygvec, zgvec = mg.vec()                              # Grid vectors
    xtemp, ytemp, ztemp = mg.meshup()      # 3-D grid

    # Define series of points for plane equations
    if len(curve) > 1:
        # Interpolate points along the extruded parabola trajectory
        x_dip, y_dip = curve_interp(curve[0], curve[1], aep['dipset_d'])
    else:
        if aep['azimuth'] is not None:
            azimuth_z += np.random.uniform(aep['azimuth'][0], aep['azimuth'][1])
        xst = -aep['dipset_d']*20 + np.random.uniform(0, aep['dipset_d'])        # Starting x-coordinate of plane points
        xend = mg.lx * 1.5 - xst                                                  # Final x-coordinate of plane points

        # Get coordinate differences
        lamb_dip = np.arange(xst, xend, aep['dipset_d'])
        xpvec = lamb_dip * np.cos(np.deg2rad(azimuth_z))
        ypvec = lamb_dip * np.sin(np.deg2rad(azimuth_z))

        # Calculate coordinates of dip points
        x_dip = xst + xpvec
        y_dip = 0 + ypvec

    # Calculate normal vector components in x/y by getting the difference between points
    p_setlamb = (xpvec**2 + ypvec**2) ** 0.5

    # Define normal vector (This might change if the plane is angled (i.e. extruded parabola settings)
    dip_z = np.random.uniform(aep['dip'][0], aep['dip'][1])
    dip_set = np.ones(np.shape(xpvec)) * dip_z
    dip_norm = np.array([xpvec, ypvec, p_setlamb * np.tan(np.deg2rad(90 - dip_set))]) #np.array((1, 0, np.tan(np.deg2rad(90 - aep['dip'][1]))))

    set_no = planepoint(dip_norm, x_dip, y_dip, znow, xtemp, ytemp, ztemp, select)
    # Re-index set_no, starting from 1 to work with 'count'
    set_no = ho.reindex(set_no) + 1

    """ Assign hydrofacies """
    # Initialise hydrofacies array
    fac_set = np.zeros((mg.nx, mg.ny, mg.nz), dtype=int)

    if aep['altfacies'] is not None:
        # Alternating hydrofacies
        ae_fac = np.asarray(aep['facies'], dtype=int)
        fac_now = random.choice(ae_fac)       # Initialise previous facies
        for idi in np.unique(set_no):
            pf_i = int(np.where(fac_now == ae_fac)[0])                                         # Get previous facies index
            fac_now = random.choice(aep['altfacies'][pf_i])  # Get next alternating facies
            fac_set[set_no == idi] = fac_now                                                    # Set previous facies

        for idx, facies in enumerate(aep['facies']):      # Cycle over hydrofacies in element
            fac_set[np.mod(set_no, idx + 1) == 0] = facies
    else:
        # Random assignment of hydrofacies
        for idi in np.unique(set_no):
            fac_set[set_no == idi] = np.random.choice(aep['facies'])

    return set_no, fac_set, dip_z, azimuth_z


def curve_interp(xc, yc, spacing):
    """
    Interpolate evenly spaced points along a curve. This code is based on code in an answer posted by 'Unutbu' on
    http://stackoverflow.com/questions/19117660/how-to-generate-equispaced-interpolating-values (retrieved 17/04/2017)

    Parameters:
        xc:			x coordinates of curve
        yc:			y coordinates of curve
        spacing:	Spacing between points

    Returns:
        - xn - x coordinates of interpolated points
        - yn - y coordinates of interpolated points

    """

    t = np.arange(xc[0], len(xc), spacing * 0.1)
    xc = np.interp(t, np.arange(len(xc)), xc)
    yc = np.interp(t, np.arange(len(yc)), yc)
    tol = spacing
    ic, idx = 0, [0]
    while ic < len(xc):
        total_dist = 0
        for j in range(ic+1, len(xc)):
            total_dist += math.sqrt((xc[j] - xc[j-1]) ** 2 + (yc[j] - yc[j-1]) ** 2)
            if total_dist > tol:
                idx.append(j)
                break
        ic = j + 1

    xn = xc[idx]
    yn = yc[idx]
    # fig, ax = plt.subplots()
    # ax.plot(xc, yc, '-')
    # ax.scatter(xn, yn)
    # ax.set_aspect('equal')
    # plt.show()

    return xn, yn


def dip_rotate(azimuth_in, dip_in):
    """
    Rotate dip angle based on azimuth
    Note that inputs and outputs are in degrees

    Parameters:
        azimuth_in:		Azimuth input angle
        dip_in:			Dipping input angle

    Returns:
        dip_out - Azimuth output angle
    """
    azimuth_in = azimuth_in * np.pi / 180
    dip_in = dip_in * np.pi / 180
    dip_out = np.arctan((np.sin(azimuth_in) + np.cos(azimuth_in) * np.tan(dip_in)) /
                        (np.cos(azimuth_in) - np.sin(azimuth_in) * np.tan(dip_in))) * 180 / np.pi
    return dip_out



"""--------------------------------------------------------------------------------------------------------------
Assignment of hydraulic properties
--------------------------------------------------------------------------------------------------------------"""

"""--------------------------------------------------------------------------------------------------------------
General functions
--------------------------------------------------------------------------------------------------------------"""


def save_arrays(arr_size, bg=None, mat_count=0, ani=True):
    """
    Generate arrays for material properties storage

    Parameters:
        arr_size:       Size of array
        bg:             List of background values for each array
        ani (bool):	    Boolean if anisotropy is to be generated

    Returns:
        - ha_arr - Material values
        - fac - Facies values

    """
    if not bg:
        bg = np.zeros((3,))

    hat_arr = np.ones(arr_size, dtype=np.int32) * -1                 # initialize architectural elements
    ha_arr = np.ones(arr_size, dtype=np.int32) * mat_count             # initialize material
    fac = np.ones(arr_size, dtype=np.int16) * int(bg[0])            # initialize hydrofacies

    if ani:
        azim = np.ones(arr_size, dtype=np.float) * bg[1]     # initialize azimuth angle
        dip = np.ones(arr_size, dtype=np.float) * bg[2]       # initialize dip angle
        return hat_arr, ha_arr, fac, azim, dip
    else:
        return hat_arr, ha_arr, fac


def prob_choose(choices, probs):
    """
    Get random values of an architectural element

    Parameters:
        choices:	Fixed number of choices
        probs: 		Contains data of architectural element units and associated hydrofacies

    Returns:
        choice *(int)* - Random value of architectural elements
    """

    ae_list = []
    for chi in range(len(choices)):
        ae_list += [choices[chi]] * int(probs[chi] * 1000)
    choice = random.choice(ae_list)

    return choice


def angle(v1, v2):
    """
    Return angle between two vectors in [°]

    Parameters:
        v1:	Vector 1
        v2:	Vector 2

    Returns:
        angle value *(float)* - Angle between v1 and v2
    """
    return np.arccos(np.abs(np.dot(v1, v2)) / (np.sqrt(np.sum(v1 ** 2, axis=1)) * np.sqrt(np.sum(v2 ** 2))))


def channel_checker(param_file, ae_name, no_extpar=1, dist=0):
    """
    extruded parabola checker function for quickly assessing the shape of extruded parabola inputs

    Parameters:
        param_file (str): 		Parameter file location
        ae_name (str):			Name of architectural element
        no_extpar (float):	Number of extruded parabolas
        dist (float): 			Distance to generate extruded parabolas - defaults to mg.lx

    Returns:
        Plots showing shape of Ferguson extruded parabolas

    """
    run, model, strata, hydraulics, flowtrans, elements, mg = hu.model_setup(param_file)
    ch_par = elements[ae_name]

    plt.figure()
    for i in range(no_extpar):
        chs = ferguson_curve(mg, ch_par['h'], ch_par['k'],  ch_par['ds'], ch_par['eps_factor'], dist=dist, disp=True)
        plt.plot(chs[:, 0], chs[:, 1])
        plt.axes().set_aspect('equal', 'datalim')

    plt.show()


def planepoint(dip_norm, x_dip, y_dip, znow, xtemp, ytemp, ztemp, select=[]):
    """
    Compute number of planes

    Parameters:
        dip_norm:
        x_dip:			X coordinates of points on dip planes
        y_dip:			Y coordinates of points on dip planes
        znow:			Current coordinates of Z, needed to compute Z coordinates of points on dip planes
        xtemp:			X dimension of model grid nodes
        ytemp:			Y dimension of model grid nodes
        ztemp:			Z dimension of model grid nodes
        select:         Model grid nodes to consider

    Returns:
        set_no - Number of planes with selected model grid nodes

    """
    # Get closest plane to points
    n_sets = dip_norm.shape[1]                   # Number of planes
    nx, ny, nz = xtemp.shape                       # Get number of model cells
    set_no = np.zeros((nx, ny, nz), dtype=np.int)  # Initialise set number array
    z_dip = np.ones(x_dip.shape) * znow

    points = np.array((xtemp[select].flatten(), ytemp[select].flatten(), ztemp[select].flatten()))      # Cartesian coordinates of model grid nodes
    plp = np.array((x_dip, y_dip, z_dip)).T                                     # Cartesian coordinates of points on dip planes
    pd = plp[:, None] - points.T                                                # subtract grid nodes from plane points

    select_idx = np.where(select)                                               # Get indices of selected model nodes

    # Loop over set planes
    for iset in range(n_sets-1):
        if iset > 1:
            pd_1 = pd_2
        else:
            abc_1 = dip_norm[:, iset]                                                           # Plane normal equation
            pd_1 = abc_1.dot(pd[iset, :, :].squeeze().T) / np.sqrt(sum(abc_1 * abc_1))          # Distance to plane
        pd1_c1 = pd_1 <= 0                                                                  # pd_1 meeting condition 1
        pd1_c1_idx = np.where(pd1_c1)

        if iset == 0:
            set_no[select_idx[0][pd1_c1_idx], select_idx[1][pd1_c1_idx], select_idx[2][pd1_c1_idx]] = iset+1
        # elif iset == n_sets: # this never happens
        #     pd1_c2_idx = np.where(pd_1 > 0)                     # index of pd_2 meeting condition 1
        #     set_no[select_idx[0][pd1_c2_idx], select_idx[1][pd1_c2_idx], select_idx[2][pd1_c2_idx]] = iset+1
        else:
            abc_2 = dip_norm[:, iset+1]
            # Points on plane
            pd_2 = abc_2.dot(pd[iset+1, :, :].squeeze().T) / np.sqrt(sum(abc_2 * abc_2))  # Distance to plane
            inset = np.logical_and(pd_1 <= 0, pd_2 > 0)                                   # grid cell between planes
            set_no[select_idx[0][inset], select_idx[1][inset], select_idx[2][inset]] = iset+1

    return set_no


"""--------------------------------------------------------------------------------------------------------------
Testing functions
--------------------------------------------------------------------------------------------------------------"""
if __name__ == '__main__':
    if len(sys.argv) > 1:
        param_file = str(sys.argv[1])
        print(sys.argv[1])
    else:
        param_file = 0
        #channel_checker(param_file, 'meander_channel', no_extpar=5, dist=10000)
        # param_file = 'E:\\Repositories\\WP3_effects\\case_studies\\trough\\tr1\\tr1dip.ini'
        # param_file = '..\\testcases\\made.ini'
    run(param_file)
