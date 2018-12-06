"""
This module contains all options for the ini-files. If you want to add an option, add it here and potentially do some post-processing in parameters.

:Authors: Samuel Scherrer, Jeremy P. Bennett
"""


from hyvr.input.option_parsing import Option

options = {}

# Run options
# ===========
options['run'] = [
    Option('runname', str, optional=True, default=None),
    Option('modeldir', str, optional=True, default=None),
    Option('numsim', int, optional=True, default=1),
    Option('outputs', list, optional=True, shape=[-1], default=[], datatype=str,
           validation_func=lambda x: x in ['vtr', 'py', 'mat', 'npz', 'h5', 'mf', 'mf6', 'hgs']),
    Option('overwrite_old_output', bool, optional=True, default=False),
]


# Model options
# =============
options['model'] = [
    Option('x0', float, optional=True, default=0.0),
    Option('y0', float, optional=True, default=0.0),
    Option('z0', float, optional=True, default=0.0),
    Option('dx', float, optional=False),
    Option('dy', float, optional=True), # use value of 'dx' if 'dy' is missing
    Option('dz', float, optional=True),
    Option('lx', float, optional=False),
    Option('ly', float, optional=False),
    Option('lz', float, optional=False),
    Option('periodic', bool, optional=True, default=False),
    Option('hydraulics', bool, optional=True, default=True),
    Option('anisotropy', bool, optional=True, default=True),
    Option('heterogeneity', bool, optional=True, default=True),
    Option('heterogeneity_level', str, optional=True, default='internal',
           validation_func=lambda x: x in ['facies', 'internal']),
]


# Strata options
# ==============
options['strata'] = [
    Option('strata', list, optional=False, shape=[-1], datatype=str),
    Option('strata_contact_models', list, optional=False, shape=[-1, -1], datatype=str),
    Option('ae', list, optional=False, shape=[-1], datatype=str),
    Option('ae_in_strata', list, optional=False, shape=['strata', -1], datatype=str),
    Option('ae_prob', list, optional=False, shape='ae_in_strata', datatype=float),
    Option('ae_z_mean', list, optional=False, shape='ae_in_strata', datatype=float),
    Option('avul_prob', list, optional=False, shape=['strata', 1], datatype=float),
    Option('avul', list, optional=False, shape=['strata', 2], datatype=float),
    Option('bg_facies', list, optional=False, shape='strata', datatype=str),
    Option('bg_azim', list, optional=False, shape='strata', datatype=float),
    Option('bg_dip', list, optional=False, shape='strata', datatype=float),
]


# Element options
# ===============
element_options = [
    Option('geometry', str, optional=False,
           validation_func=lambda x: x in ['trough', 'channel', 'sheet']),
    Option('contact_model', list, optional=True, shape=[-1], datatype=str),
    Option('facies', list, optional=False, shape=[-1], datatype=str),
    Option('altfacies', list, optional=True, shape=['facies', -1], datatype=str),
    Option('size_ztrend', list, optional=True, shape=2, datatype=float),
    Option('k_ztrend', list, optional=True, shape=2, datatype=float),
    Option('k_xtrend', list, optional=True, shape=2, datatype=float),
    Option('n_xtrend', list, optional=True, shape=2, datatype=float),
    Option('n_ztrend', list, optional=True, shape=2, datatype=float),
    Option('dip', list, optional=True, shape=2, default=[0.0, 0.0], datatype=float),
    Option('azimuth', list, optional=True, shape=2, datatype=float, default=[0., 0.]),
    # These are only necessary for the erosive elements, but they are not used for sheets anyways.
    Option('bg_facies', str, optional=True, default=None),
    Option('bg_azim', float, optional=True, default=float('NaN')),
    Option('bg_dip', float, optional=True, default=float('NaN')),
]

erosive_element_options = [
    Option('agg', float, optional=False),
    Option('buffer', float, optional=True, default=0.0),
    Option('lag', list, optional=True, shape=2, datatype=float),
    Option('lag_height', float, optional=True, default=0.0),
    Option('lag_facies', str, optional=True, default=None),
]

options['sheet'] = element_options + [
    Option('structure', str, optional=False,
           validation_func=lambda x: x in ['massive', 'dip']),
    Option('dipset_dist', float, optional=True),
    Option('lens_thickness', float, optional=False),
]


options['trough'] = element_options + erosive_element_options + [
    Option('structure', str, optional=False,
           validation_func=lambda x: x in ['massive', 'dip', 'bulb', 'bulb_sets', 'random', 'flat']),
    Option('dipset_dist', float, optional=True),
    Option('bulbset_dist', float, optional=True),
    Option('length', float, optional=False),
    Option('width', float, optional=False),
    Option('depth', float, optional=False),
    Option('trough_density', float, optional=False),
    Option('paleoflow', list, optional=False, shape=2, datatype=float),
    Option('te_xyz', list, optional=True, shape=[-1, 3], datatype=float),
    Option('migrate', list, optional=True, shape=4, datatype=float),
]


options['channel'] = element_options + erosive_element_options + [
    Option('structure', str, optional=False,
           validation_func=lambda x: x in ['massive', 'dip']),
    Option('dipset_dist', float, optional=True),
    Option('width', float, optional=False),
    Option('depth', float, optional=False),
    Option('h', float, optional=False),
    Option('k', float, optional=False),
    Option('ds', float, optional=False),
    Option('eps_factor', float, optional=False),
    Option('channel_no', int, optional=False),
    Option('mig', float, optional=True, datatype=float, default=0.),
    Option('flow_angle', float, optional=True, default=0.),
]


# Hydraulics options
# ==================
options['hydraulics'] = [
    Option('hydrofacies', list, optional=False, shape=[-1], datatype=str),
    Option('k_h', list, optional=False, shape='hydrofacies', datatype=float),
    Option('sig_y', list, optional=False, shape='hydrofacies', datatype=float),
    Option('k_ratio', list, optional=False, shape='hydrofacies', datatype=float),
    Option('n', list, optional=False, shape='hydrofacies', datatype=float),
    Option('sig_n', list, optional=False, shape='hydrofacies', datatype=float),
    Option('ycorlengths', list, optional=False, shape=['hydrofacies', 3], datatype=float),
    Option('ncorlengths', list, optional=False, shape=['hydrofacies', 3], datatype=float),
    Option('k_ztrend', list, optional=True, shape=2, datatype=float),
    Option('k_xtrend', list, optional=True, shape=2, datatype=float),
    Option('n_ztrend', list, optional=True, shape=2, datatype=float),
    Option('n_xtrend', list, optional=True, shape=2, datatype=float),
]


# Flowtrans options
# =================
options['flowtrans'] = [
    Option('hin', list, optional=True, shape=3, datatype=float),
    Option('hout', list, optional=True, shape=3, datatype=float),
    Option('gradh', float, optional=True),
    Option('q_in', float, optional=True),
    Option('kw_gradh', float, optional=True),
    Option('kw_q_in', float, optional=True)
]
