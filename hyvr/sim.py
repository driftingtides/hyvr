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
import scipy.io as sio
import matplotlib.pyplot as plt

import hyvr.utils as hu
import hyvr.optimized as ho
import hyvr.input.parameters as hp
from hyvr.postprocess.output import create_outputs
from hyvr.classes.model import Model
from hyvr.assign_facies import assign_facies_azim_dip


def run(param_file, flag_ow=None):
    """
    Main function for HYVR generation

    Parameters:
        param_file (str): 	Parameter file location
        flag_ow : bool
            If True, an existing run directory is overwritten.

    Returns:

        Save data outputs as parameter file

    """
    # TODO: implement a seed option to get the same model
    # np.random.seed(84)

    # Load parameter file
    run, model, hydraulics = hp.model_setup(param_file)
    hp.set_up_directories(run, param_file, flag_ow)

    for sim in range(1, int(run['numsim'])+1):
        # Generate facies
        model.generate_model()
        # model.assign_facies_azim_dip()
        assign_facies_azim_dip(model)
        # return

        if model.generate_hydraulics is True:
            # Generate internal heterogeneity
            model.generate_hydraulic_parameters(hydraulics)
        # from hyvr.postprocess.plotting import cross_section_pcolor
        # cross_section_pcolor(model, 'k_iso', log=True, y=10, xlim=[15.5, 16.4], ylim=[7.1, 7.7]
        # cross_section_pcolor(model, 'facies', log=False, y=10, xlim=[15.5, 16.4], ylim=[7.1, 7.7])
        # cross_section_pcolor(model, 'ha', log=False, y=10, xlim=[15.5, 16.4], ylim=[7.1, 7.7])
        # cross_section_pcolor(model, 'hat', log=False, y=10, xlim=[15.5, 16.4], ylim=[7.1, 7.7])
        # cross_section_pcolor(model, 'ae', log=False, y=10, xlim=[15.5, 16.4], ylim=[7.1, 7.7])
        # cross_section_pcolor(model, 'strata', log=False, y=10, xlim=[15.5, 16.4], ylim=[7.1, 7.7])
        # cross_section_pcolor(model, 'dip')
        # cross_section_pcolor(model, 'azim')

        # Save data
        if run['numsim'] > 1:
            realname = 'real_{:03d}'.format(sim)
            realdir = os.path.join(run['rundir'], realname)
        else:
            realname = run['runname']
            realdir = run['rundir']
        hu.try_makefolder(realdir)
        create_outputs(model, realdir, realname, run['outputs'])


    if param_file == 0:
        # this is just the testcase, so we remove the output
        from pkg_resources import cleanup_resources
        cleanup_resources()
        import shutil
        runpath = os.path.abspath(run['rundir'])
        if os.path.exists(runpath):
            shutil.rmtree(runpath)
