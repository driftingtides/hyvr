"""
This is the main module for parsing the parameter file. There are two other
modules strongly linked ot this one:

* options: contains definitions of the possible options of a parameter file
* option_parsing: provides the classes Section and Option for simpler parsing
  and validation.

:Authors: Jeremy P. Bennett, Samuel Scherrer
"""

import os
import pathlib
import sys
import shutil
import configparser as cp

from hyvr.utils import print_to_stdout
from hyvr.input.options import options
from hyvr.input.options_deprecated import options as old_options
from hyvr.input.option_parsing import *
from hyvr.geo.model import Model
from hyvr.geo.contact_surface_utils import parse_contact_model


def setup_from_inifile(inifile, flag_ow):
    """
    Parses the input file and sets up the directory structure for a HyVR run.

    Parameters
    ----------
    inifile : path
        Path to inifile. If ``inifile`` is set to 0, the ini-file for the MADE
        test case is used.
    flag_ow : bool
        Whether to overwrite existing run directories.
    
    Returns
    -------
    run : dict
        HyVR run parameters (number of simulations, where to store results)
    model : Model object
        A HyVR model object created from the ini-file.
    hydraulics : dict
        Parsed hydraulics section of the ini-file
    """
    run, model_dict, strata_dict, hydraulics, flowtrans, elements = parameters(inifile)

    # Assign architectural element identifiers
    for element in elements.keys():
        elements[element]['ae_id'] = strata_dict['ae'].index(element)

    # create model object
    model = Model(model_dict, strata_dict, elements, flowtrans)

    # set up directories
    set_up_directories(run, inifile, flag_ow)

    return run, model, hydraulics


def parameters(inifile):
    """
    Parses the inifile and returns all sections as dictionaries. Furthermore it
    sets to correct directory names for the run settings.

    Parameters
    ----------
    inifile : path
        Path to inifile. If ``inifile`` is set to 0, the ini-file for the MADE
        test case is used.
    
    Returns
    -------
    run : dict
        HyVR run settings (number of simulations, where to store results)
    model_dict : dict
        model setup parameters
    strata_dict : dict
        strata setup parameters
    hydraulics : dict
        Parsed hydraulics section of the ini-file
    flowtrans : dict
        Flow and transport settings for output
    elements : list of dicts
        List of dictionaries of settings for :class:`hyvr.model.ae_types.AEType`.
    """

    print_to_stdout("Reading parameter file")

    # test case
    if inifile == 0:
        from pkg_resources import resource_filename
        inifile = resource_filename(__name__, str(pathlib.Path('../made.ini')))

    # read file
    p = cp.ConfigParser()
    try:
        p.read(inifile, encoding='utf-8')
    except cp.MissingSectionHeaderError:
        # this is probably caused by a wrong encoding
        p.read(inifile, encoding='utf-8-sig')
    if len(p.sections()) == 0:
        raise FileNotFoundError("Parameter file {:s} not found!".format(inifile))

    old_format = 'dataoutputs' in dict(p['run'])
    if old_format:
        run, model, strata, hydraulics, flowtrans, elements = get_new_parameters_from_deprecated(*parse_deprecated_inifile(p))


    else:
        run, model, strata, hydraulics, flowtrans, elements = parse_inifile(p)

    # Runname and Modeldir
    # ====================
    # The following code sets
    # - runname: based on this, the output directory is named
    # - modeldir: the directory where to create the output directory
    # - rundir: the output directory

    # if runname is not given, it's the part of the ini-file before
    # "_parameters.ini" or ".ini"
    # in the test case runname is given
    ini_path = pathlib.Path(inifile).parent.resolve()
    ini_name = pathlib.Path(inifile).name
    if run['runname'] is None:
        if ini_name.endswith("autogenerated_backup.ini"):
            run['runname'] = ini_name[0:-25]
        else:
            run['runname'] = ".".join(ini_name.split('.')[0:-1])


    # separate task: find modeldir

    # modeldir is either the given option or the path of the .ini file
    # rundir is modeldir/runname
    # But:
    # If the .ini-filename ends with "autogenerated_backup.ini", the directory of the
    # ini-file is the rundir and modeldir is the directory above.
    # This should overwrite all other settings (except overwrite_old_output)
    if ini_name.endswith("_autogenerated_backup.ini"):
        run['modeldir'] = ini_path.parent.resolve()
        run['rundir'] = ini_path
    else:
        # either inipath or the chosen modeldir
        if run['modeldir'] is None:
            run['modeldir'] = ini_path
        if run['modeldir'] == 'select':
            directory = input(
                'Please input the model directory save path, or press <enter> to'
                'save in default directory:\n')
            if len(run['modeldir']) == 0:
                run['modeldir'] = inipath
            else:
                run['modeldir'] = pathlib.Path(directory).resolve()
        run['rundir'] = run['modeldir'] / run['runname']

    return run, model, strata, hydraulics, flowtrans, elements


def parse_inifile(p):
    """
    This function does the main work of parsing the input file.
    
    Parameters
    ----------
    p : ConfigParser
       A config parser that already read in the inifile.
    
    Returns
    -------
    run : dict
        HyVR run parameters (number of simulations, where to store results)
    model_dict : dict
        model setup parameters
    strata_dict : dict
        strata setup parameters
    hydraulics : dict
        Parsed hydraulics section of the ini-file
    flowtrans : dict
        Flow and transport settings for output
    elements : list of dicts
        List of dictionaries of settings for :class:`hyvr.model.ae_types.AEType`.
    """

    # TODO: these names were at some time keywords for hyvr, but I don't know
    # what they describe and they are not implemented anymore
    str_values = 'k_trend', 'linear_acceleration'

    sections = p.sections()
    section_parser = {}

    must_haves = ['run', 'model', 'strata', 'hydraulics']
    for section in must_haves:
        if section not in sections:
            raise MissingSectionError(section)

    # run section
    # -----------
    run = Section('run', options['run']).parse(dict(p['run']))
    del sections[sections.index('run')]

    # hydraulics section
    # ------------------
    hydraulics = Section('hydraulics', options['hydraulics']).parse(dict(p['hydraulics']))
    del sections[sections.index('hydraulics')]

    # model section
    # -------------
    model = Section('model', options['model']).parse(dict(p['model']))
    del sections[sections.index('model')]
    if model['dy'] == None:
        model['dy'] = model['dx']
    if model['dz'] == None:
        model['dz'] = model['dx']

    # strata section
    # --------------
    strata = Section('strata', options['strata']).parse(dict(p['strata']))
    del sections[sections.index('strata')]
    if len(strata['strata_contact_models']) != len(strata['strata']) - 1:
        raise ShapeError('strata_contact_models', 'strata')
    strata['contact_models'] = [
        parse_contact_model(model, depth=True) for model in strata['strata_contact_models']
    ]
    try:
        strata['bg_facies'] = parse_facies(strata['bg_facies'], hydraulics['hydrofacies'])
    except:
        raise ValueError('Invalid facies string in strata section in option bg_facies: ' + strata['bg_facies'])


    # flowtrans section
    # ------------------
    # this section is only necessary for model output
    for output in ['mf', 'mf6', 'hgs']:
        if output in run['outputs']:
            if 'flowtrans' not in sections:
                raise MissingSectionError('flowtrans')
            break
    if 'flowtrans' in sections:
        flowtrans = Section('flowtrans', options['flowtrans']).parse(dict(p['flowtrans']))
        del sections[sections.index('flowtrans')]
    else:
        flowtrans = {}

    # remaining sections are architectural elements
    elements = {}
    for section in sections:
        dictionary = dict(p[section])
        assert_exists('geometry', dictionary, section)
        geometry = dictionary['geometry']
        if geometry not in ['trough', 'channel', 'sheet']:
            raise ValueError('Invalid geometry: ' + geometry + ' in section ' + section)
        elements[section] = Section(section, options[geometry]).parse(dict(p[section]))
        elements[section]['contact_model'] = parse_contact_model(elements[section]['contact_model'], depth=False)
        # get facies number
        try:
            elements[section]['facies'] = parse_facies(elements[section]['facies'], hydraulics['hydrofacies'])
        except:
            raise ValueError('Invalid facies string in section ' + section + ' in option facies!')
        # get altfacies
        if elements[section]['altfacies'] is not None:
            try:
                elements[section]['altfacies'] = [
                    parse_facies(facies_list, hydraulics['hydrofacies']) for facies_list in elements[section]['altfacies']
                ]
            except:
                raise ValueError('Invalid facies string in section ' + section + ' in option altfacies!')
        # get bg_facies number
        if elements[section]['bg_facies'] == None:
            elements[section]['bg_facies'] = -1
        else:
            try:
                elements[section]['bg_facies'] = hydraulics['hydrofacies'].index(elements[section]['bg_facies'])
            except:
                raise ValueError('Invalid facies in section ' + section + ' in option bg_facies: ', elements[section]['bg_facies'])


        # if structure is 'dip', require 'dipset_dist'
        if elements[section]['structure'] == 'dip' or elements[section]['structure'] == 'random':
            assert_exists('dipset_dist', dictionary, section)
        # if structure is 'bulb_sets', require 'bulbset_dist'
        if elements[section]['structure'] == 'bulb' or elements[section]['structure'] == 'random':
            assert_exists('bulbset_dist', dictionary, section)

        # lag surface facies
        if geometry in ['trough', 'channel']:
            if elements[section]['lag_height'] != 0.0:
                try:
                    elements[section]['lag_facies'] = hydraulics['hydrofacies'].index(elements[section]['lag_facies'])
                except:
                    raise ValueError('lag_facies in section ' + section + ' is invalid!')
            else:
                elements[section]['lag_facies'] = -1



    return run, model, strata, hydraulics, flowtrans, elements

def get_new_parameters_from_deprecated(run, model, strata, hydraulics, flowtrans, elements):
    """
    This applies the necessary changes to parsed sections of an old-format inifile.
    """

    # run section
    run['outputs'] = run['modeloutputs'] + run['dataoutputs']
    run['overwrite_old_output'] = run['flag_ow']

    # hydraulics section
    hydraulics['hydrofacies'] = hydraulics['hydro']


    # model section
    model['anisotropy'] = run['anisotropy']
    model['hydraulics'] = hydraulics['gen']
    model['heterogeneity'] = run['het']
    model['heterogeneity_level'] = model['hetlev']
    if model['display']:
        raise ValueError('display option is not supported anymore')


    # strata section
    n_strata = len(strata['ssm'])
    strata['bg_facies'] = [int(strata['bg'][0])] * n_strata
    strata['bg_azim'] = [strata['bg'][1]] * n_strata
    strata['bg_dip'] = [strata['bg'][2]] * n_strata
    strata['strata'] = strata['ssm']
    if strata['ae_table'] is not None:
        raise ValueError('AE tables are not supported anymore!')
    if strata['save_aelu']:
        raise ValueError('save_aelu is not supported anymore!')
    contact_models = []
    for i in range(n_strata-1):
        cm = parse_old_strata_contact_model(strata, i)
        contact_models.append(cm)
    strata['contact_models'] = contact_models
    strata['ae_in_strata'] = strata['ssm_ae']

    # element sections
    for elem_name in elements:
        elem = elements[elem_name]
        elem['size_ztrend'] = elem['geo_ztrend']
        if elem['bg'] is not None:
            elem['bg_facies'] = elem['bg'][0]
            elem['bg_azim'] = elem['bg'][1]
            elem['bg_dip'] = elem['bg'][1]
        else:
            elem['bg_facies'] = -1
            elem['bg_azim'] = float('NaN')
            elem['bg_dip'] = float('NaN')
        elem['dipset_dist'] = elem['dipset_d']
        if 'bulbset_d' in elem:
            elem['bulbset_dist'] = elem['bulbset_d']
        if elem['geometry'] == 'trunc_ellip':
            elem['geometry'] = 'trough'
            elem['trough_density'] = elem['el_z']
        elem['contact_model'] = parse_old_elem_contact_model(elem)
        if elem['geometry'] in ['trough', 'channel']:
            if elem['lag'] is not None:
                elem['lag_height'] = elem['lag'][0]
                elem['lag_facies'] = elem['lag'][1]
            else:
                elem['lag_height'] = 0.0
                elem['lag_facies'] = -1



    return run, model, strata, hydraulics, flowtrans, elements


def parse_facies(facies_list, hydrofacies):
    """
    Reads a list of facies strings (facies_list) and returns the corresponding
    indices in the hydrofacies list.
    """
    return [
        hydrofacies.index(facies) for facies in facies_list
    ]

def get_new_facies_list(facies_list, hydrofacies):
    return [hydrofacies[i] for i in facies_list]


def parse_old_strata_contact_model(strata, i):
    mode = strata['ssm_contact']
    cm = {'mode':mode, 'z':strata['ssm_top'][i]}
    if mode.lower() == 'flat':
        return cm
    elif mode.lower() == 'random':
        cm['var'] = strata['ssm_contact_model'][i][0]
        cm['corlx'] = strata['ssm_contact_model'][i][1]
        cm['corly'] = strata['ssm_contact_model'][i][2]
        return cm
    else:
        raise ValueError('Unknown contact type in strata section: ' + mode)

def parse_old_elem_contact_model(elem):
    mode = elem['contact']
    cm = {'mode':mode}
    if mode.lower() == 'flat':
        return cm
    elif mode.lower() == 'random':
        cm['var'] = elem['contact_model'][0]
        cm['corlx'] = elem['contact_model'][1]
        cm['corly'] = elem['contact_model'][2]
        return cm
    else:
        raise ValueError('Unknown contact type in strata section: ' + mode)



def parse_deprecated_inifile(p):

    import warnings
    print()
    warnings.warn("You seem to be using the old ini-file format. We strongly recommend to use the new format as described in our documentation.", DeprecationWarning)
    print("------- Warning: You are using the old ini-file format ----------------------")
    print("We strongly recommend to use the new format as described in the documentation")
    print("-----------------------------------------------------------------------------")
    print()

    sections = p.sections()
    section_parser = {}

    # The following code is not very nice, with lots of repetitions. This could
    # be much nicer if the strata section had only one possible name.
    must_haves = ['run', 'model', 'strata', 'hydraulics', 'flowtrans']
    for section in must_haves:
        if section not in sections:
            raise MissingSectionError(section)

    run = Section('run', old_options['run']).parse(dict(p['run']))
    del sections[sections.index('run')]

    model = Section('model', old_options['model']).parse(dict(p['model']))
    del sections[sections.index('model')]
    if model['dy'] == None:
        model['dy'] = model['dx']
    if model['dz'] == None:
        model['dz'] = model['dx']

    strata = Section('strata', old_options['strata']).parse(dict(p['strata']))
    del sections[sections.index('strata')]
    if strata['ae_table'] is not None:
        ae_table = pathlib.Path(inifile).parent / strata['ae_table']
        if not ae_table.exists():
            raise FileNotFoundError('ae_table-file {:s} not found!'.format(ae_table))

    hydraulics = Section('hydraulics', old_options['hydraulics']).parse(dict(p['hydraulics']))
    del sections[sections.index('hydraulics')]

    flowtrans = Section('flowtrans', old_options['flowtrans']).parse(dict(p['flowtrans']))
    del sections[sections.index('flowtrans')]

    # remaining sections are architectural elements
    elements = {}
    for section in sections:
        dictionary = dict(p[section])
        assert_exists('geometry', dictionary, section)
        geometry = dictionary['geometry']
        elements[section] = Section(section, old_options[geometry]).parse(dict(p[section]))


    return run, model, strata, hydraulics, flowtrans, elements

def set_up_directories(run, inifile, overwrite_old_output=None):
    """
    This functions creates the necessary directories (modeldir, rundir). It also stores the used
    ini-file in the rundir.

    Parameters:
        run (dict):        parsed run-section of the config file
        inifile (str):  path to config file
        overwrite_old_output (bool):    Whether to overwrite the old run directory. If it is None, the option
                           from the config file will be chosen instead (default in inifile: False)
    """

    # for the test case we just create a temporary run directory/output
    # directory
    if inifile == 0:
        import tempfile
        run['rundir'] = pathlib.Path(tempfile.mkdtemp())
        return

    p = cp.ConfigParser()
    try:
        p.read(inifile, encoding='utf-8')
    except cp.MissingSectionHeaderError:
        # this is probably caused by a wrong encoding
        p.read(inifile, encoding='utf-8-sig')

    # we're now done with parsing and can create the model directory.
    run['modeldir'].mkdir(parents=True, exist_ok=True)

    # If the run directory already exists, it is overwritten if
    # overwrite_old_output is set to True.  If it doesn't exist, it will be
    # created.
    if overwrite_old_output is None:
        overwrite_old_output = run['overwrite_old_output']
    if run['rundir'].exists():
        if run['overwrite_old_output'] is True:
            # If it exists, we just delete everything in it
            for f in os.listdir(run['rundir']):
                path = os.path.join(run['rundir'], f)
                path = run['rundir'] / f
                if path.is_file():
                    path.unlink()
                elif path.is_dir():
                    shutil.rmtree(path)
        else:
            raise FileExistsError(
                "Run directory already exists and overwrite flag is "
                "set to 'false'" ". Either change the runname, or "
                "change the overwrite flag ('overwrite_old_output') "
                "in the config file or run hyvr with --overwrite.")
    else:
        run['rundir'].mkdir(parents=True, exist_ok=True)

    backup_file = run['rundir'] / (run['runname'] + '_autogenerated_backup.ini')
    with open(backup_file, 'w') as f:
        p.write(f)
