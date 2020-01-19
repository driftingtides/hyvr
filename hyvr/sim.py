# -*- coding: utf-8 -*-
""" Hydrogeological Virtual Reality simulation package.

    Hydrogeological virtual reality (HYVR) simulator for object-based modelling
    of sedimentary structures

    :Author: Jeremy P. Bennett, Samuel Scherrer
"""

import os

import hyvr.utils as hu
import hyvr.optimized as ho
import hyvr.input.parameters as hp
from hyvr.postprocess.output import create_outputs
from hyvr.assign_points import assign_points


def run(inifile, flag_ow=None):
    """
    Runs HyVR based on a given ini-file.

    Parameters
    ----------
    inifile : path
        Path to inifile. If ``inifile`` is set to 0, the ini-file for the MADE
        test case is used.
    flag_ow : bool
        Whether to overwrite existing run directories.
    """

    run_settings, model, hydraulics = hp.setup_from_inifile(inifile, flag_ow=flag_ow)

    for sim in range(1, int(run_settings['numsim'])+1):
        # Generate facies
        model.generate_model()
        assign_points(model)

        # Generate internal heterogeneity
        if model.generate_hydraulics is True:
            model.generate_hydraulic_parameters(hydraulics)

        # plots, only for development purposes
        # from hyvr.postprocess.plotting import cross_section_pcolor
        # cross_section_pcolor(model, 'k_iso', log=True)
        # cross_section_pcolor(model, 'facies')
        # cross_section_pcolor(model, 'strata')
        # cross_section_pcolor(model, 'dip')
        # cross_section_pcolor(model, 'azim')

        # Save data
        if run_settings['numsim'] > 1:
            realdir = run_settings['rundir'] / 'real_{:03d}'.format(sim)
        else:
            realdir = run_settings['rundir']
        realdir.mkdir(parents=True, exist_ok=True)
        create_outputs(model, realdir, run_settings["runname"],
                       run_settings['outputs'])


    if inifile == 0:
        # this is just the testcase, so we remove the output
        import shutil
        runpath = run_settings['rundir']
        if runpath.exists():
            shutil.rmtree(runpath)
