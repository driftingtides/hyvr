"""
This is the main module for parsing the parameter file. There are two other
modules strongly linked ot this one:
* options: contains definitions of the possible options of a parameter file
* option_parsing: provides the classes Section and Option for simpler parsing and validation.

:Authors: Jeremy P. Bennett, Samuel Scherrer
"""


import configparser as cp
import os
import sys
import time
import shutil
from hyvr.utils import try_makefolder
import hyvr.grid as gr
from hyvr.options import options
from hyvr.option_parsing import *


def model_setup(pf):
    """
    Set up model using grid.Grid() class and assign parameters

    Parameters:
        pf (str):   					Parameter file path

    Returns:
        - run *(dict)* - Model run parameters
        - mod *(dict)* - Model domain parameters
        - sequences *(dict)* - Sequence parameters
        - hydraulics *(dict)* - Hydraulic properties parameters
        - flowtrans *(dict)* - Flow & transport simulation parameters
        - elements *(dict)* - Architectural elements and parameters
        - model_grid *(object class)* - Grid object class

    """
    run, mod, strata, hydraulics, flowtrans, elements = parameters(pf)
    model_grid = gr.Grid(dx=mod['dx'],
                         dy=mod['dy'],
                         dz=mod['dz'],
                         nx=int(mod['lx'] / mod['dx']),
                         ny=int(mod['ly'] / mod['dy']),
                         nz=int(mod['lz'] / mod['dz']),
                         gtype='cells', periodicity=mod['periodic'])

    # Assign architectural element identifiers
    for element in elements.keys():
        elements[element]['ae_id'] = strata['ae'].index(element)

    return run, mod, strata, hydraulics, flowtrans, elements, model_grid

def parameters(inifile):
    """
    Get parameters for hierarchical facies modelling from a .ini-file

    Parameters:
        inifile (str):  		Parameter file path

    Returns:
        - run *(dict)* - Model run parameters
        - model *(dict)* - Model domain parameters
        - sequences *(dict)* - Sequence parameters
        - hydraulics *(dict)* - Hydraulic properties parameters
        - flowtrans *(dict)* - Flow & transport simulation parameters
        - elements *(dict)* - Architectural elements and parameters
    """

    print(time.strftime("%d-%m %H:%M:%S", time.localtime(time.time())) + ': Reading parameter file')

    # read file
    p = cp.ConfigParser()
    try:
        p.read(inifile, encoding='utf-8')
    except cp.MissingSectionHeaderError:
        # this is probably caused by a wrong encoding
        p.read(inifile, encoding='utf-8-sig')
    if len(p.sections()) == 0:
        raise FileNotFoundError("Parameter file {:s} not found!".format(inifile))

    # TODO: these are still missing in the new parsing and validation
    # I don't know where I have to add them.
    str_values = 'contact_file', 'k_trend', 'linear_acceleration'

    sections = p.sections()
    section_parser = {}

    # The following code is not very nice, with lots of repetitions. This could
    # be much nicer if the strata section had only one possible name.
    must_haves = ['run', 'model', 'strata', 'hydraulics', 'flowtrans']
    for section in must_haves:
        if section not in sections:
            raise MissingSectionError(section)

    run = Section('run', options['run']).parse(dict(p['run']))
    del sections[sections.index('run')]

    model = Section('model', options['model']).parse(dict(p['model']))
    del sections[sections.index('model')]
    if model['dy'] == None:
        model['dy'] = model['dx']
    if model['dz'] == None:
        model['dz'] = model['dx']

    strata = Section('strata', options['strata']).parse(dict(p['strata']))
    del sections[sections.index('strata')]
    if strata['ae_table'] is not None:
        if not os.path.exists(os.path.join(modeldir), strata['ae_table']):
            raise FileNotFoundError('ae_table-file {:s} not found!'.format(strata['ae_table']))

    hydraulics = Section('hydraulics', options['hydraulics']).parse(dict(p['hydraulics']))
    del sections[sections.index('hydraulics')]

    flowtrans = Section('flowtrans', options['flowtrans']).parse(dict(p['flowtrans']))
    del sections[sections.index('flowtrans')]

    # remaining sections are architectural elements
    elements = {}
    for section in sections:
        dictionary = dict(p[section])
        assert_exists('geometry', dictionary, section)
        geometry = dictionary['geometry']
        elements[section] = Section(section, options[geometry]).parse(dict(p[section]))


    # Runname and Modeldir
    # ====================

    # if runname is not given, it's the part of the ini-file before
    # "_parameters.ini" or ".ini"
    ini_path = os.path.dirname(os.path.abspath(inifile))
    ini_name = os.path.basename(os.path.abspath(inifile))
    if 'runname' not in run:
        if ini_name[-15:] == "parameters.ini":
            run['runname'] = ini_name[0:-15]
        else:
            run['runname'] = ini_name[0:-4]



    # separate task: find modeldir

    # modeldir is either the given option or the path of the .ini file
    # rundir is modeldir/runname
    # But:
    # If the .ini-filename ends with "_parameters.ini", the directory of the
    # ini-file is the rundir and modeldir is the directory above.
    # This should overwrite all other settings (except flag_ow)
    if ini_name[-15:] == "_parameters.ini":
        run['modeldir'] = os.path.abspath(os.path.join(ini_path, '../'))
        run['rundir'] = os.path.abspath(ini_path)
    else:
        # either inipath or the chosen modeldir
        run['modeldir'] = run.get('modeldir', ini_path)
        if run['modeldir'] == 'select':
            run['modeldir'] = input(
                'Please input the model directory save path, or press <enter> to'
                'save in default directory:\n')
            if len(run['modeldir']) == 0:
                run['modeldir'] = inipath
        run['modeldir'] = os.path.abspath(run['modeldir'])
        run['rundir'] = os.path.join(run['modeldir'], run['runname'])


    # we're now done with parsing and can create the model directory.
    try_makefolder(run['modeldir'])

    # now test if the run directory already exists
    if os.path.exists(run['rundir']):
        if bool(run.get('flag_ow', True)):
            # If it exists and we want to overwrite we just delete it and create a new one.
            # Therefore we first move it and try to create a new directory. If this doesn't
            # work, we move it back to the old place and abort the program.
            dirname = os.path.dirname(run['rundir'])
            tmppath = os.path.join(dirname, 'hyvr_you_should_not_see_this_1234')
            shutil.move(run['rundir'], tmppath)
            try:
                os.makedirs(run['rundir'])
            except:
                # creating a new rundir somehow didn't work, therefore we abort the program
                shutil.move(tmppath, run['rundir'])
                print("Error when trying to create new run directory.", file=sys.stderr)
                raise
            shutil.rmtree(tmppath)
        else:
            raise FileExistsError("Run directory already exists and overwrite flag is"
                            "set to 'false'" ". Either change the runname or"
                            "change the overwrite flag ('flag_ow').")
    else:
        os.makedirs(run['rundir'])

    file_save = os.path.join(run['rundir'], run['runname'] + '_parameters.ini')
    with open(file_save, 'w') as configfile:
        p.write(configfile)

    return run, model, strata, hydraulics, flowtrans, elements,
