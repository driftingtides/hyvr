"""
This file contains functions to convert the internal representation of a HyVR
model output (a dictionary of arrays) to common data or model input formats.

The functions should be named `to_<format>` and take the following parameters:

model : Model instance (see model.py)
    The model object holding all the data arrays
fname : str
    Where to save the file

If you want to add an output format, add a function in the same style and link
the inifile-name of the format with the function and the file extension in the
``create_outputs`` function below (in the dictionary ``output_desc``).
"""

import os
import pickle
import scipy.io as sio
import numpy as np
from hyvr.utils import print_to_stdout, try_makefolder


def create_outputs(model, realization_dir, runname, formats):
    """
    This functions creates output files based on the model output in the
    current run directory.

    This calls mainly the functions defined below for the different datatypes.

    Parameters
    ----------
    model : Model instance (see model.py)
        The model object holding all the data arrays
    realization_dir : str
        Directory where the current realization results should be stored.
    runname : str
        Name of the model run
    formats : list of str
        List of output formats
    """

    # This dictionary links the output format names from the ini-file with output
    # functions and file endings.
    # If you want to add an output format, add a function and add the description here.
    output_desc = {
        'mat':to_mat,
        'py':to_pickle,
        'npz':to_npz,
        'h5':to_hdf5,
        'vtr':to_vtr,
        'mf':to_modflow,
        'mf6':to_mf6,
        'hgs':to_hgs,
    }

    for fmt in formats:
        if fmt not in output_desc:
            raise ValueError('No such output format: ' + fmt)
        fname = os.path.join(os.path.abspath(realization_dir), runname)
        output_desc[fmt](model, fname)
        print_to_stdout('Saved', fmt, 'output to', realization_dir)


def to_mat(model, fname):
    """
    Saves model output as .mat file.

    Parameters
    ----------
    model : Model instance (see model.py)
        The model object holding all the data arrays
    fname : str
        Where to save the file (without file format extension)
    """
    sio.savemat(fname+'.mat', model.data)


def to_pickle(model, fname):
    """
    Saves model output as .pickle file.

    Parameters
    ----------
    model : Model instance (see model.py)
        The model object holding all the data arrays
    fname : str
        Where to save the file (without file format extension)
    """
    with open(fname+'.pickle', 'wb') as outfile:
        pickle.dump(model.data, outfile, protocol=pickle.HIGHEST_PROTOCOL)

def to_npz(model, fname):
    """
    Saves model output as .npz file.

    Parameters
    ----------
    model : Model instance (see model.py)
        The model object holding all the data arrays
    fname : str
        Where to save the file (without file format extension)
    """
    np.savez_compressed(fname+'.npz', **model.data)


def to_hdf5(model, fname):
    """
    Saves model output as .h5 file. This requires h5py.

    Parameters
    ----------
    model : Model instance (see model.py)
        The model object holding all the data arrays
    fname : str
        Where to save the file (without file format extension)
    """
    try:
        import h5py
        with h5py.File(fname+'.h5', 'w') as hf:
            for key in model.data:
                hf.create_dataset(key, data=model.data[key], compression=True)
    except ImportError:
        raise ValueError('h5 output not possible: h5py not found.')


def to_vtr(model, fname):
    """
    Saves model output as .vtr file. This requires pyevtk

    Parameters
    ----------
    model : Model instance (see model.py)
        The model object holding all the data arrays
    fname : str
        Where to save the file (without file format extension)
    """
    # ktensor can not be saved as vtr
    try:
        from pyevtk.hl import gridToVTK
        data_dict = {key:model.data[key] for key in model.data if key != "ktensors"}
        xv = np.arange(model.grid.x0, model.grid.xmax+model.grid.dx, model.grid.dx)
        yv = np.arange(model.grid.y0, model.grid.ymax+model.grid.dy, model.grid.dy)
        zv = np.arange(model.grid.z0, model.grid.zmax+model.grid.dz, model.grid.dz)
        gridToVTK(fname, xv, yv, zv, cellData=data_dict)
    except ImportError:
        raise ValueError('vtr output not possible: pyevtk not found.')

def to_modflow(model, fname):
    """
    Saves model output in modflow format. This requires flopy.

    Parameters
    ----------
    model : Model instance (see model.py)
        The model object holding all the data arrays
    fname : str
        Where to save the file (without file format extension)
    """
    try:
        import flopy
    except ImportError:
        raise ValueError('modflow output not possible: flopy not found.')
    # For modflow we want to create a new folder instead of only a file. The folder name is the base
    # name of the passed filename
    realization_dir = os.path.dirname(fname)
    runname = os.path.basename(fname)
    mfdir = os.path.join(realization_dir, 'MODFLOW')
    mfname = os.path.join(mfdir, runname)
    try_makefolder(mfdir)

    # Assign name and create modflow model object
    mf = flopy.modflow.Modflow(mfname, exe_name='mf2005')

    # Create the discretization object
    ztop = model.grid.z0 + model.grid.lz
    zbot = model.grid.z0
    botm = np.linspace(ztop, zbot, model.grid.nz + 1)
    dis = flopy.modflow.ModflowDis(mf, model.grid.nz, model.grid.nx, model.grid.ny, 
                                   delr=model.grid.dx, delc=model.grid.dy,
                                   top=ztop, botm=botm[1:])

    # Variables for the BAS package
    ibound = np.ones((model.grid.nz, model.grid.nx, model.grid.ny), dtype=np.int32)
    ibound[:, :, 0] = -1
    ibound[:, :, -1] = -1

    strt = np.ones((model.grid.nz, model.grid.nx, model.grid.ny), dtype=np.float32)
    strt[:, :, 0] = model.flowtrans['hin'][0]
    strt[:, :, -1] = model.flowtrans['hout'][0]

    bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt)

    # Assign hydraulic conductivity
    hyvr_hk = np.transpose(model.data['k_iso'], (2, 0, 1))
    hyvr_layvka = 1                                           # VKA dataset is ratio of horizontal K
    if 'anirat' in model.data.keys():
        hyvr_vka = np.transpose(model.data['anirat'], (2, 0, 1))

        # Add LPF package to the MODFLOW model
        lpf = flopy.modflow.ModflowLpf(mf,                        # Modflow object
                                       hk=hyvr_hk,              # Horizontal hydraulic conductivity
                                       layvka=hyvr_layvka,      # Flag for each layer of anisotropic ratio
                                       vka=hyvr_vka)            # Anisotropy ratios.
    else:
        # Add LPF package to the MODFLOW model
        lpf = flopy.modflow.ModflowLpf(mf,                        # Modflow object
                                       hk=hyvr_hk)              # Horizontal hydraulic conductivity

    oc = flopy.modflow.ModflowOc(mf)        # Add OC package to the MODFLOW model
    pcg = flopy.modflow.ModflowPcg(mf)      # Add PCG package to the MODFLOW model
    mf.write_input()                        # Write the MODFLOW model input files
    



def to_mf6(model, fname):
    """
    Saves model output in mf6 format. This requires flopy.

    Parameters
    ----------
    model : Model instance (see model.py)
        The model object holding all the data arrays
    fname : str
        Where to save the file (without file format extension)
    """

    try:
        import flopy
    except ImportError:
        raise ValueError('mf6 output not possible: Flopy not installed.')

    # For modflow we want to create a new folder instead of only a file. The folder name is the base
    # name of the passed filename
    realization_dir = os.path.dirname(fname)
    runname = os.path.basename(fname)
    mfdir = os.path.join(realization_dir, 'mf6')
    mfname = os.path.join(mfdir, runname)
    try_makefolder(mfdir)

    # Transpose HyVR arrays for MF6 input
    transpose_order = (2, 1, 0)
    k_iso = np.transpose(model.data['k_iso'], transpose_order)
    if set(['anirat', 'dip', 'azim']).issubset(model.data.keys()):
        anirat = np.transpose(model.data['anirat'], transpose_order)
        dip = np.transpose(model.data['dip'], transpose_order)
        azim = np.transpose(model.data['azim'], transpose_order)
        xt3d = True
    else:
        xt3d = False

    """ create simulation """
    sim = flopy.mf6.MFSimulation(sim_name=runname,
                                 version='mf6',
                                 exe_name='mf6',
                                 sim_ws=mfdir)
                                 # sim_tdis_file='simulation.tdis')

    """ Create the Flopy temporal discretization object - STEADY-STATE """
    tdis = flopy.mf6.modflow.mftdis.ModflowTdis(sim,
                                                time_units='DAYS')

    """ create gwf model """
    gwf = flopy.mf6.MFModel(sim, modelname=runname)


    ims = flopy.mf6.ModflowIms(sim,
                               print_option='SUMMARY',
                               complexity='COMPLEX',
                               outer_hclose=1e-3,
                               outer_maximum=500,
                               under_relaxation='NONE',
                               inner_maximum=100,
                               inner_hclose=1e-4,
                               rcloserecord=0.001,
                               linear_acceleration='BICGSTAB',
                               scaling_method='NONE',
                               reordering_method='NONE',
                               relaxation_factor=0.97)
    sim.register_ims_package(ims, [gwf.name])

    """ Create discretization """
    ztop = model.grid.z0 + model.grid.lz
    zbot = model.grid.z0
    botm = np.around(np.arange(ztop, zbot-model.grid.dz, -model.grid.dz), decimals=3)
    dis = flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(gwf,
                                                   nlay=model.grid.nz,
                                                   nrow=model.grid.ny,
                                                   ncol=model.grid.nx,
                                                   delr=model.grid.dy,
                                                   delc=model.grid.dx,
                                                   top=ztop,
                                                   botm=botm[1:],
                                                   fname='{}.dis'.format(runname))

    """ Create Node Property Flow package object """
    if xt3d is True:
        npf_package = flopy.mf6.ModflowGwfnpf(gwf,
                                              save_flows=True, icelltype=0, xt3doptions='',
                                              k=k_iso,                                  # within-bedding hydraulic conductivity
                                              k33=k_iso/anirat,                         # across-bedding hydraulic conductivity
                                              angle1=azim,                              # azimuth
                                              angle2=dip,                               # dip
                                              angle3=np.zeros((model.grid.nz, model.grid.ny, model.grid.nx)))   # no rotation
    else:
        npf_package = flopy.mf6.ModflowGwfnpf(gwf,
                                              save_flows=True, icelltype=0,
                                              k=k_iso)                                # within-bedding hydraulic conductivity)

    """ Create constant head package """
    if model.flowtrans['hin'] is not None:
        hin = model.flowtrans['hin'][0]
        hout = model.flowtrans['hout'][0]

    elif flowtrans['gradh'] is not None:
        hout = 1
        hin = hout + model.grid.lx * flowtrans['gradh']

    if np.any([model.flowtrans['hin'] is not None, model.flowtrans['gradh'] is not None]):
        chd_rec = []
        for layer in range(0, model.grid.nz):
            for row in range(0, model.grid.ny):
                chd_rec.append(((layer, row, 0), hin))         # Apply at model inlet
                chd_rec.append(((layer, row, model.grid.nx-1), hout))   # Apply at model outlet
        # chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(gwf, maxbound=len(chd_rec),
        #                                                stress_period_data=chd_rec, save_flows=True)

        chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(gwf, maxbound=len(chd_rec),
                                                       stress_period_data=chd_rec, save_flows=True)
    elif model.flowtrans['q_in'] is not None:
        """ Apply fixed head at model outlet if fixed head at inlet"""
        hin = 1
        hout = 1
        chd_rec = []
        for layer in range(0, model.grid.nz):
            for row in range(0, model.grid.ny):
                chd_rec.append(((layer, row, model.grid.nx-1), hout))   # Apply at model outlet

        chd = flopy.mf6.modflow.mfgwfchd.ModflowGwfchd(gwf, maxbound=len(chd_rec),
                                                       stress_period_data=chd_rec, save_flows=True)
    else:
        hin = 1
        hout = 1

    """ Create the initial conditions package """
    # Create linear initial condition
    # hstart = np.ones_like(k_iso) *(hin - hout)/2
    hstart = np.ones_like(k_iso) * np.linspace(hin, hout, model.grid.nx)
    ic = flopy.mf6.modflow.mfgwfic.ModflowGwfic(gwf, strt=hstart)



    """ Create well package """
    # Apply constant discharges at model faces
    if 'q_in' in model.flowtrans:
        if model.flowtrans['q_in'] is not None:
            q_in = model.flowtrans['q_in']
        else:
            q_in = 0.01

        # if 'q_out' in flowtrans:
        #     q_out = flowtrans['q_out']
        # else:
        #     q_out = -q_in

        wel_rec = []
        for layer in range(0, model.grid.nz):
            for row in range(0, model.grid.ny):
                wel_rec.append(((layer, row, 0),  q_in, 'inlet'))          # Apply at model inlet
                # wel_rec.append(((layer, row, model.nx-1), q_out, 'outlet'))        # Apply at model outlet

        # Apply to model
        wel = flopy.mf6.ModflowGwfwel(gwf,
                                      print_input=True,
                                      print_flows=True,
                                      save_flows=True,
                                      boundnames=True,
                                      maxbound=len(wel_rec),
                                      stress_period_data=wel_rec)


    """ Create the output control package """
    headfile = '{}.hds'.format(runname)
    head_filerecord = [headfile]
    budgetfile = '{}.cbc'.format(runname)
    budget_filerecord = [budgetfile]
    saverecord = [('HEAD', 'ALL'),
                  ('BUDGET', 'ALL')]
    printrecord = [('HEAD', 'LAST')]
    oc = flopy.mf6.modflow.mfgwfoc.ModflowGwfoc(gwf,
                                                saverecord=saverecord,
                                                head_filerecord=head_filerecord,
                                                budget_filerecord=budget_filerecord,
                                                printrecord=printrecord)

    # write simulation
    sim.write_simulation()


def to_hgs(model, fname):
    """
    Saves model output in the HydroGeoSphere format.

    Parameters
    ----------
    model : Model instance (see model.py)
        The model object holding all the data arrays
    fname : str
        Where to save the file (without file format extension)
    """
    realization_dir = os.path.dirname(fname)
    runname = os.path.basename(fname)
    hgsdir = os.path.join(realization_dir, 'HGS')
    try_makefolder(hgsdir)

    uid = np.arange(1, len(model.data['ktensors'][:, :, :, 1, 2].flatten()) + 1)    # Create list of IDs
    vals_to_write = {'ktensors': np.column_stack((uid,
                                      model.data['ktensors'][:, :, :, 0, 0].flatten(),            # K_xx
                                      model.data['ktensors'][:, :, :, 1, 1].flatten(),            # K_yy
                                      model.data['ktensors'][:, :, :, 2, 2].flatten(),            # K_zz
                                      model.data['ktensors'][:, :, :, 0, 1].flatten(),            # K_xy
                                      model.data['ktensors'][:, :, :, 0, 2].flatten(),            # K_xz
                                      model.data['ktensors'][:, :, :, 1, 2].flatten())),          # K_yz
                     'porosity': np.column_stack((uid,
                                      model.data['poros'].flatten()))}
    val_fmts = {'ktensors': '%u %1.3e %1.3e %1.3e %1.3e %1.3e %1.3e',
                'porosity': '%u %1.3f'}

    # Loop over properties to write
    for val in vals_to_write:
        val_filepath = os.path.join(hgsdir, val + '.txt')       # File name of HGS output file
        np.savetxt(val_filepath, vals_to_write[val], fmt=val_fmts[val])
