# -*- coding: utf-8 -*-
"""Some utility functions for HFM modelling

    :Authors: Jeremy P. Bennett, with help from Alessandro Comunian and Samuel Scherrer

    :Notes:



"""

import sys
import pickle
import numpy as np
import time
import configparser as cp
import pandas as pd
import linecache
import scipy.io
import scipy.spatial as sp
import os
import shutil
import errno
import scipy.stats as st
import hyvr.optimized as ho




def mf6_vtr(fhead, mg, fout):
    """
    Convert a MODFLOW 6 Binary head file into vtr suitable for visualisation in ParaView

    Parameters
    ----------



    """
    try:
        import flopy
    except ImportError:
        print('mf output not possible: Flopy not installed.')
        return
    hfile = flopy.utils.binaryfile.HeadFile(fhead)          #
    hdata = hfile.get_alldata()                             # Create numpy array with all data
    head_dict = dict()                                      # Initialise dict with the data
    for i in range(0, hdata.shape[0]):
        hf_i = np.squeeze(hdata[i, :, :, :])                # Get heads at individual time steps
        hf_i = np.transpose(hf_i, (2, 1, 0))                # Permute to be consistent with HyVR grids
        head_dict['head_timestep{}'.format(i)] = hf_i

    to_vtr(head_dict, fout, mg, points=True)


def dem_load(fn):
    """
    Load data from ESRI-style ASCII-file.

    Parameters:
        fn (str): 				Directory and file name for save

    Returns:
        - data *(numpy array)* - Data from ERSI-style ASCII-file
        - meta *(dict)* - Dict with grid metadata

    """

    # Extract header using linecache
    meta = {}
    meta['ncols'] = int(linecache.getline(fn, 1).split()[1])
    meta['nrows'] = int(linecache.getline(fn, 2).split()[1])
    meta['ox'] = linecache.getline(fn, 3).split()[1]
    meta['oy'] = linecache.getline(fn, 4).split()[1]
    meta['cell_size'] = linecache.getline(fn, 5).split()[1]
    meta['no_Data'] = linecache.getline(fn, 6).split()[1]

    # Extract data using pandas
    df = pd.read_csv(fn, header=None, delimiter=' ', skiprows=6, dtype=np.float)
    data = df.as_matrix()

    return data, meta


def dem_save(fn, data, gro):
    """
    Save DEM data to ESRI-style ASCII-file

    Parameters:
        fn (str):               Directory and file name for save
        data (numpy array):     DEM data
        gr (object class):      grid.Grid() object class

    Returns:
        Save DEM data to ESRI-style ASCII-file

    """

    header = ("ncols            {0.nx}\n"
              "nrows            {0.ny}\n"
              "xllcorner        {0.ox}\n"
              "yllcorner        {0.oy}\n"
              "cellsize         {0.cs2}\n"
              "NOODATA_value    -9999"
              ).format(gro)

    with open(fn, mode='wb') as out_file:
        np.savetxt(out_file,
                   data,
                   header=header,
                   fmt='%.4f',
                   comments='')



def load_gslib(fn):
    """
    Load .gslib files. This has been appropriated from the HPGL library
    https://github.com/hpgl/hpgl/blob/master/src/geo_bsd/routines.py
    commit b980e15ad9b1f7107fd4fa56ab117f45553be3aa

    Parameters:
        fn (str): 			.gslib file path and name

    Returns:
        gslib_dict *(dict)* - properties

    """
    gslib_dict = {}
    list_prop = []
    points = []

    f = open(fn)
    head = f.readline().split('\t')
    num_p = int(f.readline())
    #print num_p

    lx, ly, lz = [int(x) for x in head[0].split(' ')]
    nx, ny, nz = [float(x) for x in head[1].split(' ')]
    ox, oy, oz = [float(x) for x in head[2].split(' ')]

    for i in range(num_p):
        list_prop.append(str(f.readline().strip()))
    #print list_prop

    for i in range(len(list_prop)):
        gslib_dict[list_prop[i]] = np.zeros((lx * ly * lz))

    index = np.zeros(len(list_prop))

    for line in f:
        points = line.split()
        for j in range(len(points)):
            gslib_dict[list_prop[j]][index[j]] = float(points[j])
            index[j] += 1

    for dkey in gslib_dict.keys():
        gslib_dict[dkey] = gslib_dict[dkey].reshape((ly, lx, lz))

    f.close()

    return gslib_dict


''' HYVR-specific utilities'''


def read_lu(sq_fp):
    """
    Load user-defined strata (architectural element lookup table),
    split the data based on a delimiter and return it as a new list

    Parameters:
        sq_fp:			Load user-defined strata (architectural element lookup table)

    Returns:
        ssm_lu *(list)*: -Values of architectural element lookup table

    """
    # Load user-defined systems / architectural element lookup table
    print(time.strftime("%d-%m %H:%M:%S", time.localtime(time.time())) + ': Reading strata data from ' + sq_fp)
    with open(sq_fp) as f:
        lines = f.read().splitlines()

    ssm_lu = []
    for li in lines[1:]:
        temp = li.split(',')
        ssm_lu.append([int(temp[0]), float(temp[1]), float(temp[2]), str(temp[3]), int(temp[4])])

    return ssm_lu



def virtual_boreholes(data_dict, d, l, file_out=None, vals=[], opts=[]):
    """ Perform 'virtual' borehole sampling of parameter field

    Arguments:
        data_dict (dict):
            Data to sample
        d (list):
            3-tuple of model grid cell dimensions
        l (list):
            3-tuple of total model dimensions/lengths
        file_out (list: basefile path, list of file types):
            Output filename and path
        vals (list):
            Parameter fields to include
        opts (dict):
            Sampling options
            opts['noBH'] (int):             Random sampling
            opts['grid_spacing']:           float, or list of floats [x spacing, y spacing] Grid sample spacing
            opts['grid_n']:                 int, or list of ints [n in x, n in y] Number of grid nodes per x,y dimensions
            opt['lnK'] (bool):              Natural logarithm of isotropic hydraulic conductivity

    Returns:
        bh_df : Pandas DataFrame class

    """

    nx, ny, nz = np.shape(data_dict['fac'])

    # Set up column names
    cols = ['x', 'y', 'z']
    if len(vals) == 0:
        vals = data_dict.keys()
    cols.extend(vals)

    # Create dataframe
    bh_df = pd.DataFrame(columns=cols)

    # Sampling of grid
    xy_grid = []
    xv = np.arange(0.5 * d[0], l[0], d[0])
    yv = np.arange(0.5 * d[1], l[1], d[1])
    zv = np.arange(0.5 * d[2], l[2], d[2])

    if 'grid_spacing' in opts.keys():
        """ Uniform grid with set spacing """

        # Get cartesian coordinates in 2D (x,y)
        if len(opts['grid_spacing']) == 2:
            range_x = np.arange((opts['grid_spacing'][0] * 0.5), l[0], opts['grid_spacing'][0])
            range_y = np.arange((opts['grid_spacing'][1] * 0.5), l[1], opts['grid_spacing'][1])
        else:
            range_x = np.arange((opts['grid_spacing'] * 0.5), l[0], opts['grid_spacing'])
            range_y = np.arange((opts['grid_spacing'] * 0.5), l[1], opts['grid_spacing'])

        x_locs, y_locs = np.meshgrid(range_x, range_y)

        # Convert to array indices
        x_locs = np.floor(x_locs.flatten()/d[0]).astype(int)
        y_locs = np.floor(y_locs.flatten()/d[1]).astype(int)

    elif 'grid_n' in opts.keys():
        """ Sample over uniform grid """
        # Get cartesian coordinates in 2D (x,y)
        range_x = np.linspace(0, l[0]/d[0] - 1, opts['grid_n'][0])
        range_y = np.linspace(0, l[1]/d[1] - 1, opts['grid_n'][1]+2)[1:-1]
        x_locs, y_locs = np.meshgrid(range_x, range_y)

        # Convert to array indices
        x_locs = np.floor(x_locs.flatten()).astype(int)
        y_locs = np.floor(y_locs.flatten()).astype(int)

    elif 'noBH' in opts.keys():
        """ Randomly sample the xy plane """
        x_locs = np.random.choice(range(0, nx), opts['noBH'])    # Borehole location indices
        y_locs = np.random.choice(range(0, ny), opts['noBH'])    # Borehole location indices

    # Put data into dataframe
    for idx in range(len(x_locs)):
        # Get indices of location
        i = x_locs[idx]
        j = y_locs[idx]

        # Get vectors of Cartesian coordinates
        ibh = np.zeros((nz, 3 + len(vals)))
        ibh[:, 0] = np.ones((nz,)) * xv[i]          # x coordinates
        ibh[:, 1] = np.ones((nz,)) * yv[j]          # y coordinates
        ibh[:, 2] = zv                              # z coordinates

        for iv, v in enumerate(vals):
            # Append to list to be appended to dataframe
            ibh[:, iv+3] = data_dict[v][i, j, 0:nz]
        bh_df = bh_df.append(pd.DataFrame(ibh, columns=cols), ignore_index=True)

    if 'lnK' in opts and opts['lnK'] is True:
        vals.extend(['lnK'])
        bh_df['lnK'] = pd.Series(np.log(bh_df['k_iso']), index=bh_df.index)

    # Save borehole data
    fmtd = {'k_iso': '%.5e',
            'lnK': '%.5f',
            'poros': '%.5f',
            'fac': '%u',
            'dip': '%.2f',
            'azim': '%.2f'}
    if file_out is not None:
        if 'csv' in file_out[1]:
            bh_df.to_csv(file_out[0]+'.csv', index=False)
        if 'gslib' in file_out[1]:
            n_conddata = bh_df.shape[0]
            colsout = ['x', 'y', 'z']
            colsout.extend(vals)
            to_write = bh_df.as_matrix(columns=colsout)
            header = [str(n_conddata), str(3 + len(vals)), 'x', 'y', 'z']
            header.extend(vals)
            header = '\n'.join(header)

            fmts = '%.3f %.3f %.3f {}'.format(' '.join([fmtd[i] for i in vals]))
            np.savetxt(file_out[0] + '.gslib', to_write, delimiter=' ', header=header, comments='', fmt=fmts)

    return bh_df


def calc_norm(x):
    """
    Calculate norm (compute the complex conjugate from 'x')

    Parameters:
        x:	Input parameter

    Returns:
        Complex conjugate of x

    """
    s = (x.conj() * x).real
    return np.sqrt(np.add.reduce(s, axis=0))


def try_makefolder(makedir):
    """
    Create modflow output folder

    """
    # Create modflow output folder

    try:
        os.makedirs(makedir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def print_to_stdout(*args):
    """
    Prints a message to stdout with timestamp.
    """
    print(time.strftime("%d-%m %H:%M:%S", time.localtime(time.time())) + ':', *args)


####################################################################################################################
# Some utilities for generating shapes
####################################################################################################################

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

    if len(select) > 0:
        select_idx = np.where(select)                                               # Get indices of selected model nodes
    else:
        select = np.ones_like(xtemp, dtype=bool)
        select_idx = np.where(select)

    points = np.array((xtemp[select].flatten(), ytemp[select].flatten(), ztemp[select].flatten()))      # Cartesian coordinates of model grid nodes
    plp = np.array((x_dip, y_dip, z_dip)).T                                     # Cartesian coordinates of points on dip planes
    pd = plp[:, None] - points.T                                                # subtract grid nodes from plane points

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


def angle(v1, v2):
    """
    Return angle between two vectors in [°] between 0° and 180°

    Parameters:
        v1:	Vector 1
        v2:	Vector 2

    Returns:
        angle value *(float)* - Angle between v1 and v2
    """
    cos = np.dot(v1, v2) / np.sqrt(np.dot(v1, v1) * np.dot(v2, v2))
    angle = np.arccos(cos)/np.pi*180
    return angle

def get_alternating_facies(num_facies, type_params):
    """
    Returns a vector of alternating facies numbers based on the 'altfacies'
    setting in the inifile.
    """
    facies = np.zeros(num_facies, dtype=np.int32)
    facies[0] = np.random.choice(type_params['facies'])
    # The facies are changing according to the given 'altfacies'
    if type_params['altfacies'] is not None:
        for i in range(1, num_facies):
            fac_idx = type_params['facies'].index(facies[i-1])
            facies[i] = np.random.choice(type_params['altfacies'][fac_idx])
    else:
        for i in range(1, num_facies):
            facies[i] = np.random.choice(type_params['facies'])
    return facies

def norm(v):
    return np.sqrt(np.dot(v,v))


def specsim(grid, var, corl, selection_mask=None, two_dim=False, covmod='gaussian'):
    """
    Generate random variables with stationary covariance function using spectral
    techniques of Dietrich & Newsam (1993)

    Parameters
    ----------
    grid : Grid instance
    var : float
        Variance
    corl : tuple of floats
        Tuple of correlation lengths. 2-tuple for 2-d (x and y), and 3-tuple
        for 3-d (x, y, z).
    selection_mask : boolean numpy array
        This array should have the same size as the model grid and can be used to only generate
        random variables for a subset of the model grid. If ``two_dim=True``, it should have
        dimensions ``(nx, ny)``, otherwise ``(nx, ny, nz)``.
    two_dim : bool, optional (default: False)
        Whether to return a two-dimensional or a 3 dimensional field
    covmod : str, optional (default: "gaussian")
        Which covariance model to use ("gaussian" or "exp").

    Returns
    -------
    Y : 1d, 2d, or 3d numpy array
        Numpy array of random field. If no selection mask was given, this is either a 2d or 3d numpy
        array, depending on ``two_dim`` and of the size of the model grid.
        If a selection_mask was given, this is a flat array of values.
    """

    # THINK: If this is called multiple times in a row with the same parameters (i.e. same var,
    # corl, covmod), then it is not necessary to recalculate syy, and we even get two fields out.
    # Is there a way to use this?
    mask_given = selection_mask is not None
    if mask_given:
        # This method only works for rectangular grids, this means we have to find the smallest
        # rectangular selection that contains all selected cells
        if two_dim:
            x_idx, y_idx = np.where(selection_mask)
            min_x_idx = min(x_idx)
            max_x_idx = max(x_idx)
            min_y_idx = min(y_idx)
            max_y_idx = max(y_idx)
            rectangular_mask = np.zeros_like(selection_mask)
            rectangular_mask[min_x_idx:max_x_idx+1,min_y_idx:max_y_idx+1] = True
            selected_X = grid.X[min_x_idx:max_x_idx+1,min_y_idx:max_y_idx+1,0]
            selected_Y = grid.Y[min_x_idx:max_x_idx+1,min_y_idx:max_y_idx+1,0]
            selected_X_centered = selected_X - (np.min(selected_X) + np.max(selected_X))/2
            selected_Y_centered = selected_Y - (np.min(selected_Y) + np.max(selected_Y))/2
            selection_mask_small = selection_mask[min_x_idx:max_x_idx+1,min_y_idx:max_y_idx+1]

            h_square = (selected_X_centered/corl[0])**2 \
                    + (selected_Y_centered/corl[1])**2

        else:
            x_idx, y_idx, z_idx = np.where(selection_mask)
            try:
                min_x_idx = min(x_idx)
            except:
                import pdb; pdb.set_trace()
            max_x_idx = max(x_idx)
            min_y_idx = min(y_idx)
            max_y_idx = max(y_idx)
            max_z_idx = max(z_idx)
            min_z_idx = min(z_idx)
            rectangular_mask = np.zeros_like(selection_mask)
            rectangular_mask[min_x_idx:max_x_idx+1,min_y_idx:max_y_idx+1,min_z_idx:max_z_idx+1] = True
            selected_X = grid.X[min_x_idx:max_x_idx+1,min_y_idx:max_y_idx+1,min_z_idx:max_z_idx+1]
            selected_Y = grid.Y[min_x_idx:max_x_idx+1,min_y_idx:max_y_idx+1,min_z_idx:max_z_idx+1]
            selected_Z = grid.Z[min_x_idx:max_x_idx+1,min_y_idx:max_y_idx+1,min_z_idx:max_z_idx+1]
            selected_X_centered = selected_X - (np.min(selected_X) + np.max(selected_X))/2
            selected_Y_centered = selected_Y - (np.min(selected_Y) + np.max(selected_Y))/2
            selected_Z_centered = selected_Z - (np.min(selected_Z) + np.max(selected_Z))/2
            selection_mask_small = selection_mask[min_x_idx:max_x_idx+1,min_y_idx:max_y_idx+1,min_z_idx:max_z_idx+1]

            h_square = (selected_X_centered/corl[0])**2 \
                     + (selected_Y_centered/corl[1])**2 \
                     + (selected_Z_centered/corl[2])**2

    else:
        # full grid
        if two_dim:
            h_square = (np.asarray(grid.X_centered[:,:,0])/corl[0])**2 \
                    + (np.asarray(grid.Y_centered[:,:,0])/corl[1])**2
        else:
            h_square = (np.asarray(grid.X_centered)/corl[0])**2 \
                     + (np.asarray(grid.Y_centered)/corl[1])**2 \
                     + (np.asarray(grid.Z_centered)/corl[2])**2

    ntot = h_square.size

    # Covariance matrix of variables
    if covmod == 'gaussian':
        # Gaussian covariance model
        ryy = np.exp(-h_square) * var
    elif covmod == 'exp':
        # Exponential covariance model
        ryy = np.exp(-np.sqrt(h_square)) * var
    else:
        raise ValueError('Invalid covariance model')

    # Power spectrum of variable
    syy = np.fft.fftn(np.fft.fftshift(ryy)) / ntot
    syy = np.abs(syy)       # Remove imaginary artifacts
    syy[0] = 0

    real = np.random.randn(*syy.shape)
    imag = np.random.randn(*syy.shape)
    epsilon = real + 1j*imag
    rand = epsilon * np.sqrt(syy)
    Y = np.real(np.fft.ifftn(rand * ntot))

    if mask_given:
        return Y[selection_mask_small]

    return Y
