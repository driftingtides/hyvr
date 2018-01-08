# -*- coding: utf-8 -*-
"""Some utility functions for HFM modelling

    :Authors: Jeremy P. Bennett, with help from Alessandro Comunian

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
import errno
from pyevtk.hl import gridToVTK
import flopy
import scipy.stats as st
import hyvr.hyvr.grid as gr


''' File I/O and wrangling'''


def to_vtk(data, file_name, grid=None, sc_name=None):
    """ Save a numpy array into a VTK STRUCTURED_POINTS file.

    Parameters:
        data (numpy array):		Numpy array containing the data, `int` or `float` or `uint`. 
                                The dimensions should be between 1 and 3
        file_name (string):		Name of the file for the output, optional (None)
        grid (class Grid):		Information about the grid can be also provided as a Grid object            
        sc_name (string):		Name of the scalar quantities

    Returns:
        A VTK 'STRUCTURED_POINTS' dataset file containing the input numpy data.

        .. note:
            * Only 'STRUCTURED_POINTS' output allowed.
            * The type of the data is extracted from the input array.
            * Only 3D or 2D input data allowed.
            * The default output is set to 'POINT_DATA'.

    """

    if grid:
        nx = grid.nx
        ny = grid.ny
        nz = grid.nz
    else:
        # Create an internal default grid
        grid = gr.Grid()
        try:
            grid.nx, grid.ny, grid.nz = np.shape(data)
        except ValueError:
            print('    Warning (numpy2vtk): input data considered as 2D.')
            grid.nx, grid.ny = np.shape(data)
            grid.nz = 1

    # Set the correct format for the output
    if 'int' in data.dtype.name:
        fmt = '%i'
        fmt_head = 'int'
    elif 'float' in data.dtype.name:
        fmt = '%.4e'
        fmt_head = 'float'
    elif sc_name == 'facies':
        fmt = '%i'
        fmt_head = 'int'
    else:
        print(('    Error in "numpy2vtk", wrong data type "%s"' % data.dtype.name))

#    print("GRID:", grid)
    if grid.gtype == 'points':
        header = (
            "# vtk DataFile Version 3.4\n"
            "{0.gname}\n"
            "ASCII\n"
            "DATASET STRUCTURED_POINTS\n"
            "DIMENSIONS {0.nx:d} {0.ny:d} {0.nz:d}\n"
            "ORIGIN {0.ox:f} {0.oy:f} {0.oz:f}\n"
            "SPACING {0.dx:f} {0.dy:f} {0.dz:f}\n"
            "POINT_DATA {0.points:d}\n"
            "SCALARS type {1} 1\n"
            "LOOKUP_TABLE default"
            ).format(grid, fmt_head)
    elif grid.gtype == 'cells':
        header = (
            "# vtk DataFile Version 3.4\n"
            "{0.gname}\n"
            "ASCII\n"
            "DATASET RECTILINEAR_GRID\n"
            "DIMENSIONS {0.nx:d} {0.ny:d} {0.nz:d}\n"
            "X_COORDINATES {0.nx:d} float\n"
            "{2}\n"
            "Y_COORDINATES {0.ny:d} float\n"
            "{3}\n"
            "Z_COORDINATES {0.nz:d} float\n"
            "{4}\n"
            "POINT_DATA {0.points:d}\n"
            "SCALARS type {1} 1\n"
            "LOOKUP_TABLE default"
            ).format(grid, fmt_head, grid.vec_x(), grid.vec_y(), grid.vec_z())


    if sc_name:
        header = header.replace('type',sc_name)

    with open(file_name, mode='wb') as out_file:
        np.savetxt(out_file,
                   np.ravel(data, order='F'),
                   fmt=fmt, header=header, comments='')

    print(time.strftime('%X'), ': VTK export complete')


def to_vtr(data_dict, file_name, grid):
    """
    Save a numpy array into a ``.vtr`` rectilinear grid of voxels using pyevtk

    Parameters:
        data_dict (numpy array):		e.g. {'fac': fac, 'mat': mat}
        file_name (string):		        Name of the file for the output.
        grid (class Grid):		        The information about the grid can be also provided as a grid object.

    Returns:
        A VTK STRUCTURED_POINTS dataset file containing the input numpy data.

        .. note:
            * Only 'STRUCTURED_POINTS' output allowed.
            * The type of the data is extracted from the input array.
            * Only 3D of 2D input data allowed.
            * The default output is set to 'POINT_DATA'.

    """
    gvec = grid.vec_node()
    gridToVTK(file_name, gvec[0], gvec[1], gvec[2], cellData=data_dict)
    print(time.strftime('%X'), ': VTR export complete')


def vtk_mask(data, out_name, mask_val=-15, grid=None):
    """
    Prepare a mask file from a ``vtk`` input. All values should be 0 or 1 after this operation

    Parameters:
        data (numpy array):		Data to replace
        out_name (string):		Name of the file for the output.
        mask_val (int):			Value of the facies which should be masked
        grid:					Grid class

    Returns:
        A VTK 'STRUCTURED_POINTS' dataset file containing the mask data

    """

    # Get grid details
    if grid:
        nx = grid.nx
        ny = grid.ny
        nz = grid.nz
    else:
        # Create an internal default grid
        grid = gr.Grid(gname='IMPALA Mask input')
        try:
            grid.nx, grid.ny, grid.nz = np.shape(data)
        except ValueError:
            print('    Warning (numpy2vtk): input data considered as 2D.')
            grid.nx, grid.ny = np.shape(data)
            grid.nz = 1

    uni = np.unique(data)                             # Get unique values
    uni = [x for x in uni if x not in mask_val]     # Remove masked values from unique array
    mask = np.zeros_like(data)
    mask[np.in1d(data, uni)] = 1

    grid.gname = 'IMPALA Mask input'
    to_vtk(mask, out_name, grid, 'mask_code')
    return uni


def vtk_read(file_in):
    """
    Reads a ``.vtk`` file into a numpy array

    Parameters:
        file_in (str): 				Name and filepath to read

    Returns:
        gegrid (hyvr.grid class): 	Grid class
        props (numpy array): 		Grid properties


    """

    with open(file_in) as vtkfile:
        for line in vtkfile:
            if line.startswith('#') or not line.strip():
                # skip comments and blank lines
                continue
            elif line.startswith('DIMENSIONS'):
                nx, ny, nz = line.split()[1:4]
                nx = int(nx)
                ny = int(ny)
                nz = int(nz)
            elif line.startswith('SPACING'):
                dx, dy, dz = line.split()[1:4]
                dx = np.float16(dx)
                dy = np.float16(dy)
                dz = np.float16(dz)
            elif line.startswith('LOOKUP_TABLE'):
                props = []
                try:
                    while True:
                        val = int(next(vtkfile))
                        props.extend([val])
                except StopIteration:
                    break

    props = np.asarray(props, dtype=np.int8).reshape((nx, ny, nz), order='F')
    gegrid = gr.Grid(nx=nx, ny=ny, nz=nz,
                     dx=dx, dy=dy, dz=dz)

    print(time.strftime('%X'), ': VTK read complete')
    return gegrid, props


def vtk_trim(file_in, dims, file_out=None):
    """
    Trims a vtk file to the desired dimensions.
    Removes the effort of working out the indexing in a .grdecl file
    Saves as a .vtk file with everything the same except the dimensions

    Parameters:
        file_in: 	.vtk file to trim
        dims: 		3-tuple of dimensions
        file_out: 	Output file

    Returns:
        gegrid: 	Grid class of the data
        props: 		The data as a nx x ny x nz array

        ...notes:
            * Number of cells (nx, ny, nz) is changed
            * Spacing (dx, dy, dz) is NOT changed

    """

    # Read in VTK file
    v_in = vtk_read(file_in)
    gegrid = v_in[0]
    v_props = v_in[1]

    # Slice up property data
    props = v_props[0:dims[0], 0:dims[1], 0:dims[2]]

    # Create corresponding grid class
    gegrid.nx = dims[0]
    gegrid.ny = dims[1]
    gegrid.nz = dims[2]

    # Create output file name
    if not file_out:
        file_out = file_in[:file_in.find('.')] + '_trim.vtk'

    to_vtk(props, file_out, grid=gegrid)

    return gegrid, props


def dem_load(fn):
    """ 
    Load data from ESRI-style ASCII-file.

    Parameters:
        fn (str): 				Directory and file name for save

    Returns:
        data (numpy array):     Data from ERSI-style ASCII-file
        meta (dict):            Dict with grid metadata

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


def matlab_save(fn, data):
    """ 
    Save numpy array to .mat file for use in matlab.

    Parameters:
        fn (str):               File name (ending with .mat)
        data (numpy array):     Data to save

    Returns:
        Save a dictionary of names and arrays into a MATLAB-style .mat file.
        This saves the array objects in the given dictionary to a MATLAB- style .mat file.

    """

    scipy.io.savemat(fn, dict(data=data))


def load_gslib(fn):
    """ 
    Load .gslib files. This has been appropriated from the HPGL library
    https://github.com/hpgl/hpgl/blob/master/src/geo_bsd/routines.py
    commit b980e15ad9b1f7107fd4fa56ab117f45553be3aa

    Parameters:
        fn (str): 			.gslib file path and name

    Returns:
        gslib_dict (dict):  properties

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


def load_pickle(pickfile):
    """
    Pickle input file

    Parameters:
        pickfile:		Input file

    Return:
        data (dict):	Pickled data of input file 

    """
    with open(pickfile, 'rb') as f:
        data = pickle.load(f)

    return data


''' HYVR-specific utilities'''


def parameters(file_in):
    """ 
    Get parameters for hierarchical facies modelling

    Parameters:
        file_in (str):  		Parameter file path

    Returns:
        run (dict):             Model run parameters
        model (dict):           Model domain parameters
        sequences (dict):       Sequence parameters
        hydraulics (dict):      Hydraulic properties parameters
        flowtrans (dict):       Flow & transport simulation parameters
        elements (dict):        Architectural elements and parameters

    """

    p = cp.ConfigParser()
    p.read(file_in)
    if len(p.sections()) == 0:
        sys.exit('Parameter file not found')

    elements = {}
    str_values = 'geometry', 'structure', 'runname', 'contact', 'contact_file', 'modeldir', 'hetlev', 'ae_table',\
                 'k_trend', 'seq_contact'
    dtype_dict = {'ll_seq_ae': str, 'll_altfacies': int, 'll_contact_model': float,  'll_ae_prob': float,
                  'll_ae_z_mean': float, 'll_ycorlengths': float, 'll_ncorlengths': float, 'll_seq_contact_model': float,
                  'll_avul': float, 'll_avul_prob': float}



    for section in p.sections():
        ndict = dict(p[section])
        for key in ndict:
            if key.startswith('flag'):
                ndict[key] = p.getboolean(section, key)
            elif key.startswith(str_values):
                ndict[key] = str(ndict[key])
            elif key.startswith('r_'):
                ndict[key] = [float(i) for i in ndict[key].replace(' ', '').replace('[', '').replace(']', '').split(',')]
            elif key.startswith('hin') or key.startswith('hout'):
                ndict[key] = [float(i) for i in ndict[key].replace('[', '').replace(']', '').replace(' ', '').split(',')]
            elif key.startswith('l_'):
                ndict[key] = [i for i in ndict[key].replace('[', '').replace(']', '').replace(' ', '').split(',')]
            elif key.startswith('ll_'):
                ndict[key] = [[i for i in x.strip(" []").replace(' ', '').split(",")] for x in ndict[key].strip('[]').split("],")]
                if dtype_dict[key] is float:
                    ndict[key] = [[float(j) for j in i] for i in ndict[key]]
                elif dtype_dict[key] is int:
                    ndict[key] = [[int(j) for j in i] for i in ndict[key]]
            else:
                ndict[key] = float(ndict[key])

        if section == 'model':
            model = ndict
        elif section == 'sequences':
            sequences = ndict
        elif section == 'hydraulics':
            hydraulics = ndict
        elif section == 'run':
            run = ndict
        elif section == 'flowtrans':
            flowtrans = ndict
        else:
            elements[section] = ndict

    # Create runfile directory in main directory
    if 'modeldir' not in run:
        # Use same directory as parameter file
        run['modeldir'] = os.path.abspath('/'.join(file_in.split('/')[0:-1]))
    elif run['modeldir'] == 'select':
        # User input for save location
        run['modeldir'] = input('Please input the model directory save path, or press <enter> to save in default directory:\n')
        if len(run['modeldir']) == 0:
             run['modeldir'] = os.path.abspath(os.path.join(os.path.dirname(__file__)[0:-5], 'tests\\made_example\\'))
        try_makefolder(run['modeldir'])
    else:
        try_makefolder(run['modeldir'])

    run['rundir'] = run['modeldir'] + '\\' + run['runname']
    try_makefolder(run['rundir'])

    # Save parameter file
    if 'flag_ow' in run and run['flag_ow'] is False:
        file_save = run['rundir'] + '\\' + time.strftime('%d-%m-%Y_%H.%M.%S') + '_parameters.ini'
    else:
        file_save = run['rundir'] + '\\' + run['runname'] + '_parameters.ini'
    with open(file_save, 'w') as configfile:
        p.write(configfile)

    return run, model, sequences, hydraulics, flowtrans, elements,


def model_setup(pf):
    """ 
    Set up model using grid.Grid() class and assign parameters

    Parameters:
        pf (str):   					Parameter file path

    Returns:
        run (dict):                     Model run parameters
        mod (dict):                     Model domain parameters
        sequences (dict):               Sequence parameters
        hydraulics (dict):              Hydraulic properties parameters
        flowtrans (dict):               Flow & transport simulation parameters
        elements (dict):                Architectural elements and parameters
        model_grid (object class):      Grid object class

    """
    run, mod, sequences, hydraulics, flowtrans, elements = parameters(pf)
    model_grid = gr.Grid(dx=mod['dx'],
                         dy=mod['dy'],
                         dz=mod['dz'],
                         nx=int(mod['lx'] / mod['dx']),
                         ny=int(mod['ly'] / mod['dy']),
                         nz=int(mod['lz'] / mod['dz']),
                         gtype='cells', periodicity=mod['flag_periodic'])

    return run, mod, sequences, hydraulics, flowtrans, elements, model_grid


def read_lu(sq_fp):
    """
    Load user-defined sequences (architectural element lookup table),
    split the data based on a delimiter and return it as a new list

    Parameters:
        sq_fp:			Load user-defined sequences (architectural element lookup table)

    Returns:
        seq_lu (list):	Values of architectural element lookup table

    """
    # Load user-defined sequences / architectural element lookup table
    print(time.strftime("%d-%m %H:%M:%S", time.localtime(time.time())) + ': Reading sequence data from ' + sq_fp)
    with open(sq_fp) as f:
        lines = f.read().splitlines()

    seq_lu = []
    for li in lines:
        temp = li.split(', ')
        seq_lu.append([int(temp[0]), float(temp[1]), float(temp[2]), str(temp[3]), int(temp[4])])

    return seq_lu


def to_modflow(mfdir, mg, flowtrans, k_iso, anirat):
    """ 
    Convert HYVR outputs to MODFLOW inputs

    Parameters:
        mfdir:				Directory of MODFLOW model object
        mg:					Mesh grid object class
        flowtrans (dict):   Flow & transport simulation parameters
        k_iso:				Hydraulic conductivity of HYVR
        anirat:				Background anistropic ratio (K_h/K_v anisotropy ratio)

    Returns:
        mf:					MODFLOW model object
        dis:				Discretization of modflow object
        bas:				BAS package of modflow model
        lpf:				LPF package of modflow model
        oc:					OC package of modflow model
        pcg:				pcg package of modflow model

    """

    # Assign name and create modflow model object
    mf = flopy.modflow.Modflow(mfdir, exe_name='mf2005')

    # Create the discretization object
    ztop = mg.oz + mg.lz
    zbot = mg.oz
    botm = np.linspace(ztop, zbot, mg.nz + 1)
    dis = flopy.modflow.ModflowDis(mf, mg.nz, mg.nx, mg.ny, delr=mg.dx, delc=mg.dy, top=(mg.oz + mg.lz), botm=botm[1:])

    # Variables for the BAS package
    ibound = np.ones((mg.nz, mg.nx, mg.ny), dtype=np.int32)
    ibound[:, :, 0] = -1
    ibound[:, :, -1] = -1

    strt = np.ones((mg.nz, mg.nx, mg.ny), dtype=np.float32)
    strt[:, :, 0] = flowtrans['hin'][0]
    strt[:, :, -1] = flowtrans['hout'][0]

    bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt)

    # Assign hydraulic conductivity
    hyvr_hk = np.transpose(k_iso, (2, 0, 1))
    hyvr_layvka = 1                                           # VKA dataset is ratio of horizontal K
    hyvr_vka = np.transpose(anirat, (2, 0, 1))

    # Add LPF package to the MODFLOW model
    lpf = flopy.modflow.ModflowLpf(mf,                        # Modflow object
                                   hk=hyvr_hk,              # Horizontal hydraulic conductivity
                                   layvka=hyvr_layvka,      # Flag for each layer of anisotropic ratio
                                   vka=hyvr_vka)            # Anisotropy ratios.

    oc = flopy.modflow.ModflowOc(mf)        # Add OC package to the MODFLOW model
    pcg = flopy.modflow.ModflowPcg(mf)      # Add PCG package to the MODFLOW model
    mf.write_input()                        # Write the MODFLOW model input files

    return mf, dis, bas, lpf, oc, pcg


def to_mf6(mfdir, runname, mg, flowtrans, k_iso, anirat, dip, azim):
    """
    .. note:

        UNDER CONSTRUCTION
        https://github.com/modflowpy/flopy/blob/develop/examples/Notebooks/flopy3_mf6_tutorial.ipynb

    Convert HYVR outputs to MODFLOW6 inputs

    Parameters:
        mfdir:				Directory of MODFLOW model object
        runname:            Run name
        mg:					Mesh grid object class
        flowtrans (dict):   Flow & transport simulation parameters
        k_iso:				Hydraulic conductivity of HYVR
        anirat:				Background anistropic ratio (K_h/K_v anisotropy ratio)

    Returns:
        mf:					MODFLOW model object
        dis:				Discretization of modflow object
        bas:				BAS package of modflow model
        lpf:				LPF package of modflow model
        oc:					OC package of modflow model
        pcg:				pcg package of modflow model

    """

    # Transpose HyVR arrays for MF6 input
    transpose_order = (2, 1, 0)
    k_iso = np.transpose(k_iso, transpose_order)
    anirat = np.transpose(anirat, transpose_order)
    dip = np.transpose(dip, transpose_order)
    azim = np.transpose(azim, transpose_order)

    # Create the Flopy simulation object
    sim = flopy.mf6.MFSimulation(sim_name=runname, sim_ws=mfdir)

    # Create the Flopy temporal discretization object
    tdis = flopy.mf6.ModflowTdis(sim, time_units='DAYS', nper=1,
                                 tdisrecarray=[(1.0, 1, 1.0)])

    # Create the Flopy groundwater flow (gwf) model object
    gwf_name = runname
    gwf = flopy.mf6.MFModel(sim, model_type='gwf6',
                            model_name=gwf_name,
                            model_nam_file='{}.nam'.format(gwf_name),
                            sms_file_name='{}.sms'.format(gwf_name))

    # Create the Flopy iterative model solver Package object (taken from ex15-whirlsxt3d from MF6)
    # ims = flopy.mf6.ModflowIms(sim, print_option='SUMMARY', complexity='SIMPLE',
    #                            outer_hclose=1.e-8, outer_maximum=500,
    #                            under_relaxation='NONE', inner_maximum=100,
    #                            inner_hclose=1.e-8, rcloserecord=0.001,
    #                            linear_acceleration='BICGSTAB', scaling_method='NONE',
    #                            reordering_method='NONE', relaxation_factor=0.97)
    # sim.register_ims_package(ims, [gwf_name])

    # Create discretisation object
    ztop = mg.oz + mg.lz
    zbot = mg.oz
    botm = np.linspace(ztop, zbot, mg.nz + 1)
    dis = flopy.mf6.ModflowGwfdis(gwf, length_units='METERS',
                                  nlay=mg.nz, nrow=mg.ny, ncol=mg.nx,
                                  delr=mg.dy, delc=mg.dx,
                                  top=ztop, botm=botm[1:],
                                  fname='{}.dis'.format(gwf_name))

    # Create the initial conditions package
    # start = np.zeros((mg.nz, mg.ny, mg.nx))
    # ic = flopy.mf6.modflow.mfgwfic.ModflowGwfic(gwf, strt=start)

    # Create Node Property Flow package object
    npf_package = flopy.mf6.ModflowGwfnpf(gwf, save_flows=True, icelltype=0, xt3doptions='',
                                          k=k_iso,                                  # within-bedding hydraulic conductivity
                                          k33=k_iso/anirat,                         # across-bedding hydraulic conductivity
                                          angle1=azim,                              # azimuth
                                          angle2=dip,                               # dip
                                          angle3=np.zeros((mg.nz, mg.ny, mg.nx)))   # no rotation

    # Create the constant head package.
    # chd_rec = []
    # for layer in range(0, mg.nz):
    #     for row in range(0, mg.ny):
    #         chd_rec.append(((layer, row, 0),  flowtrans['hin'][0]))         # Apply at model inlet
    #         chd_rec.append(((layer, row, mg.nx-1), flowtrans['hout'][0]))   # Apply at model outlet
    # chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(gwf, maxbound=len(chd_rec),
    #                                                periodrecarray=chd_rec, save_flows=True)

    # Create the output control package
    # headfile = '{}.hds'.format(gwf_name)
    # head_filerecord = [headfile]
    # budgetfile = '{}.cbb'.format(gwf_name)
    # budget_filerecord = [budgetfile]
    # saverecord = [('HEAD', 'ALL'),
    #               ('BUDGET', 'ALL')]
    # printrecord = [('HEAD', 'LAST')]
    # oc = flopy.mf6.modflow.mfgwfoc.ModflowGwfoc(gwf, saverecord=saverecord,
    #                                             head_filerecord=head_filerecord,
    #                                             budget_filerecord=budget_filerecord,
    #                                             printrecord=printrecord)

    # write simulation to new location
    sim.write_simulation()
    # success, buff = sim.run_simulation()
    # print('\nSuccess is: ', success)


def to_hgs(hgspath, mg, flowtrans, ktensors, poros):
    """ 
    Convert HYVR outputs to HydroGeoSphere inputs

    Parameters:
        hgspath (str):		Path where to save HGS output file
        ktensors:			Array with tensor values of K
        poros:				Array with values of porosity

    Returns:
        val_fmts (dict):	Dictionary with values of K tensor and porosity
        val_filepath:		file name of HGS output file

    """

    uid = np.arange(1, len(ktensors[:, :, :, 1, 2].flatten()) + 1)                              # Create list of IDs
    vals_to_write = {'ktensors': np.column_stack((uid,
                                                  ktensors[:, :, :, 0, 0].flatten(),            # K_xx
                                                  ktensors[:, :, :, 1, 1].flatten(),            # K_yy
                                                  ktensors[:, :, :, 2, 2].flatten(),            # K_zz
                                                  ktensors[:, :, :, 0, 1].flatten(),            # K_xy
                                                  ktensors[:, :, :, 0, 2].flatten(),            # K_xz
                                                  ktensors[:, :, :, 1, 2].flatten())),          # K_yz
                     'porosity': np.column_stack((uid,
                                                  poros.flatten()))}
    val_fmts = {'ktensors': '%u %1.3e %1.3e %1.3e %1.3e %1.3e %1.3e',
                'porosity': '%u %1.3f'}

    # Loop over properties to write
    for val in vals_to_write:
        val_filepath = hgspath + val + '.txt'                      # File name of HGS output file
        np.savetxt(val_filepath, vals_to_write[val], fmt=val_fmts[val])


def round_x(x, base=1, prec=2):
    """
    Round to the nearest z-increment (Refer to http://stackoverflow.com/questions/2272149/round-to-5-or-other-number-in-python)

    Parameters:
        x (float):		Input parameter
        base (int):		Base parameter for avoiding floating-point values
        prec:			Precision of rounding

    Returns:	
        Rounded value of nearest z-increment

    """
    return np.round(base * np.round(x/base), prec)


def rotate_ktensor(count, aniso, azimuth, dip, k_in):
    """
    Create a rotated K tensor

    Parameters:
        count (int): 	Material number and/or identifier
        aniso:			Anisotropy
        azimuth:		Azimuth angles
        dip:			Dipping angles
        k_in:			-

    Returns:
        k_rotate:		Rotated K tensor
    
    """

    # convert dip and azimuth to radians
    dip = dip * np.pi / 180
    azimuth = azimuth * np.pi / 180

    kplane = np.ones(1, count) * np.sqrt(aniso)  # relative value
    kperp = np.ones(1, count) / np.sqrt(aniso)    # relative value

    k_rotate = np.empty((3,3,count), dtype=np.float16)
    for ii in np.arange(0, count):
        R = np.array([[np.cos(azimuth[ii]), np.sin(azimuth[ii]), 0],
                      [-np.sin(azimuth[ii]), np.cos(azimuth[ii]), 0],
                      [0, 0, 1]], dtype=np.float16) * ...
        np.array([[np.cos(dip[ii]), 0, np.sin(dip[ii])],
                  [0, 1, 0],
                  [-np.sin(dip[ii]), 0, np.cos(dip[ii])]], dtype=np.float16)
        k_rotate[:, :, ii] = R * np.diag(np.array([kplane[ii], kplane[ii], kperp[ii]])) * R.T

    return k_rotate


def get_boreholes(bh_loc, fin, fout=None):
    """
    Get virtual borehole data
    Returns values at centroids - this might not match the borehole inputs

    Parameters:
        bh_loc: 	Location of boreholes
        fin: 		Filepath of properties input
        fout: 		Filepath to save borehole information

    Returns:
        bhdf:	X, Y values at centroids

    .. note:

        UNDER CONSTRUCTION!

    """

    # Read in properties file
    df = pd.read_csv(fin, index_col='id', sep=' ')
    xy = df.as_matrix(['X', 'Y'])

    kdt = sp.KDTree(xy)             # Initialise search tree
    nn_idx = kdt.query(bh_loc)[1]   # Get indices of nearest centroids to boreholes
    bhdf = pd.DataFrame(columns=('bh', df.columns))

    for i in nn_idx:
        # Get all nn values with the same coordinates
        bhdf.append(df.ix[np.nonzero(xy[:, 0] == xy[i, 0])[0]])

    if fout:
        np.savetxt(fout, np.column_stack((np.arange(1, len(fe_zones)+1), centroids_pm, fe_facies, fe_K, fe_poros)),
                   fmt='%i %f %f %f %i %1.2e %.2f')

    return bhdf


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


def specsim(gr, var, corl, twod=False, covmod='gau'):
    """
    Generate random variables stationary covariance function using spectral techniques of Dietrich & Newsam (1993)

    Parameters:
        gr:     	Grid class object
        var:    	Variance
        corl:   	Tuple of correlation length of random variable
        twod:   	Flag for two-dimensional simulation
        covmod: 	Which covariance model to use
                        'gau': Gaussian
                        'exp': Exponential

    Returns:
        bigy:   	Random gaussian variable
                    Real part of a complex array, created via inverse discrete Fourier Transform

    """

    if twod is True:
        yy, xx = np.meshgrid(np.arange(-gr.ny*0.5*gr.dy, gr.ny*0.5*gr.dy, gr.dy),
                             np.arange(-gr.nx*0.5*gr.dx, gr.nx*0.5*gr.dx, gr.dx))
        h = ((xx / corl[0]) ** 2 + (yy / corl[1]) ** 2) ** 0.5      # Compute distance from origin

    else:
        yy, xx, zz = np.meshgrid(np.arange(-gr.ny*0.5*gr.dy, gr.ny*0.5*gr.dy, gr.dy),
                                 np.arange(-gr.nx*0.5*gr.dx, gr.nx*0.5*gr.dx, gr.dx),
                                 np.arange(-gr.nz*0.5*gr.dz, gr.nz*0.5*gr.dz, gr.dz))

        # Compute distance from origin
        h = ((xx / corl[0]) ** 2 + (yy / corl[1]) ** 2 + (zz / corl[2]) ** 2) ** 0.5

    ntot = np.size(xx)

    # Covariance matrix of variables
    if covmod == 'gau':
        # Gaussian covariance model
        ryy = np.exp(-h**2) * var
    elif covmod == 'exp':
        # Exponential covariance model
        ryy = np.exp(-np.abs(h)) * var
    else:
        sys.exit('Invalid covariance model')

    # Power spectrum of variable
    syy = np.fft.fftn(np.fft.fftshift(ryy)) / ntot
    syy = np.abs(syy)       # Remove imaginary artifacts
    if twod is True:
        syy[0, 0] = 0
    else:
        syy[0, 0, 0] = 0

    real = st.norm.rvs(size=syy.shape)
    imag = st.norm.rvs(size=syy.shape)
    epsilon = real + 1j*imag
    rand = epsilon * (syy**0.5)
    bigy = np.real(np.fft.ifftn(rand * ntot))

    return bigy


if __name__ == '__main__':
    """ Testing functions"""

