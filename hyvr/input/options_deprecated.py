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
    Option('anisotropy', bool, optional=True, default=True,
           alternatives='flag_anisotropy'),
    Option('het', bool, optional=True, default=True,
           alternatives='flag_het'),
    Option('dataoutputs', list, optional=True, shape=[-1], default=[], datatype=str,
           validation_func=lambda x: x in ['vtr', 'py', 'mat', 'npz', 'h5'],
           alternatives='l_dataoutputs'),
    Option('modeloutputs', list, optional=True, shape=[-1], default=[], datatype=str,
           validation_func=lambda x:  x in ['mf', 'mf6', 'hgs'],
           alternatives='l_modeloutputs'),
    Option('flag_ow', bool, optional=True, default=False),
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
    Option('periodic', bool, optional=True, default=False,
           alternatives='flag_periodic'),
    Option('display', bool, optional=True, default=False,
           alternatives='flag_display'),
    Option('hetlev', str, optional=False,
           validation_func=lambda x: x in ['ae', 'facies', 'internal']),
]


# Strata options
# ==============
options['strata'] = [
    Option('ssm', list, optional=False, shape=[-1], datatype=str,
           alternatives=['l_ssm', 'l_seq']),
    Option('ssm_top', list, optional=False, shape='ssm', datatype=float,
           alternatives='r_ssm_top'),
    Option('ssm_contact_model', list, optional=True, shape=['ssm', 3], datatype=float,
           alternatives='ll_ssm_contact_model'),
    Option('ssm_contact', str, optional=True, default='flat',
           validation_func=lambda x: x in ['flat', 'random', 'user']),
    Option('ae_table', str, optional=True, default=None),
    Option('save_aelu', bool, optional=True, default=False),
    Option('ae', list, optional=False, shape=[-1], datatype=str,
           alternatives='l_ae'),
    Option('ssm_ae', list, optional=False, shape=['ssm', -1], datatype=str,
           alternatives='ll_ssm_ae'),
    Option('ae_prob', list, optional=False, shape='ssm_ae', datatype=float,
           alternatives='ll_ae_prob'),
    Option('ae_z_mean', list, optional=False, shape='ssm_ae', datatype=float,
           alternatives='ll_ae_z_mean'),
    Option('avul_prob', list, optional=False, shape=['ssm', 1], datatype=float,
           alternatives='ll_avul_prob'),
    Option('avul', list, optional=False, shape=['ssm', 2], datatype=float,
           alternatives='ll_avul'),
    Option('bg', list, optional=False, shape=3, datatype=float),
]


# Element options
# ===============
element_options = [
    Option('geometry', str, optional=False,
           validation_func=lambda x: x in ['trunc_ellip', 'ext_par', 'sheet']),
    Option('structure', str, optional=False),
    Option('contact', str, optional=False,
           validation_func=lambda x: x in ['flat', 'random']),
    Option('contact_model', list, optional=True, shape=[3], datatype=float,
           alternatives='r_contact_model'),
    Option('facies', list, optional=False, shape=[-1], datatype=int,
           alternatives='l_facies'),
    Option('altfacies', list, optional=True, shape=['facies', -1], datatype=int,
           alternatives='ll_altfacies'),
    Option('geo_ztrend', list, optional=True, shape=2, datatype=float,
           alternatives='r_geo_ztrend'),
    Option('k_ztrend', list, optional=True, shape=2, datatype=float,
           alternatives='r_k_ztrend'),
    Option('k_xtrend', list, optional=True, shape=2, datatype=float,
           alternatives='r_k_xtrend'),
    Option('n_xtrend', list, optional=True, shape=2, datatype=float,
           alternatives='r_n_xtrend'),
    Option('n_ztrend', list, optional=True, shape=2, datatype=float,
           alternatives='r_n_ztrend'),
    Option('bg', list, optional=True, shape=3, datatype=float)
]

erosive_element_options = [
    Option('agg', float, optional=False),
    Option('buffer', float, optional=True, default=0.0),
    Option('dipset_d', float, optional=True),
    Option('migrate', list, optional=True, shape=4, datatype=float,
           alternatives='r_migrate'),
    Option('lag', list, optional=True, shape=2, datatype=float,
           alternatives='r_lag'),
    Option('dip', list, optional=False, shape=2, datatype=float,
           alternatives='r_dip'),

]

options['sheet'] = element_options + [
    Option('structure', str, optional=False,
           validation_func=lambda x: x in ['massive', 'dip']),
    Option('lens_thickness', float, optional=False),
    Option('dipset_d', float, optional=True),
    Option('dip', list, optional=True, shape=2, datatype=float, default=[0., 0.],
           alternatives='r_dip'),
    Option('azimuth', list, optional=True, shape=2, datatype=float, default=[0., 0.],
           alternatives='r_azimuth'),
]


options['trunc_ellip'] = element_options + erosive_element_options + [
    Option('structure', str, optional=False,
           validation_func=lambda x: x in ['massive', 'dip', 'bulb', 'bulb_l', 'random', 'flat']),
    Option('length', float, optional=False),
    Option('width', float, optional=False),
    Option('depth', float, optional=False),
    Option('el_z', float, optional=False),
    Option('paleoflow', list, optional=False, shape=2, datatype=float,
           alternatives='r_paleoflow'),
    Option('azimuth', list, optional=False, shape=2, datatype=float,
           alternatives='r_azimuth'),
    Option('bulbset_d', float, optional=True),
    Option('te_xyz', list, optional=True, shape=[-1, 3], datatype=float,
           alternatives='ll_te_xyz'),
]


options['ext_par'] = element_options + erosive_element_options + [
    Option('structure', str, optional=False,
           validation_func=lambda x: x in ['massive', 'dip']),
    Option('width', float, optional=False),
    Option('depth', float, optional=False),
    Option('h', float, optional=False),
    Option('k', float, optional=False),
    Option('ds', float, optional=False),
    Option('eps_factor', float, optional=False),
    Option('channel_no', float, optional=False),
]


# Hydraulics options
# ==================
options['hydraulics'] = [
    Option('gen', bool, optional=True, default=True, alternatives='flag_gen'),
    Option('hydro', list, optional=False, shape=[-1], datatype=str,
           alternatives='l_hydro'),
    Option('k_h', list, optional=False, shape='hydro', datatype=float,
           alternatives='r_k_h'),
    Option('sig_y', list, optional=False, shape='hydro', datatype=float,
           alternatives='r_sig_y'),
    Option('k_ratio', list, optional=False, shape='hydro', datatype=float,
           alternatives='r_k_ratio'),
    Option('n', list, optional=False, shape='hydro', datatype=float,
           alternatives='r_n'),
    Option('sig_n', list, optional=False, shape='hydro', datatype=float,
           alternatives='r_sig_n'),
    Option('ycorlengths', list, optional=False, shape=['hydro', 3], datatype=float,
           alternatives='ll_ycorlengths'),
    Option('ncorlengths', list, optional=False, shape=['hydro', 3], datatype=float,
           alternatives='ll_ncorlengths'),
    Option('k_ztrend', list, optional=True, shape=2, datatype=float,
           alternatives='r_k_ztrend'),
    Option('k_xtrend', list, optional=True, shape=2, datatype=float,
           alternatives='r_k_xtrend'),
    Option('n_ztrend', list, optional=True, shape=2, datatype=float,
           alternatives='r_n_ztrend'),
    Option('n_xtrend', list, optional=True, shape=2, datatype=float,
           alternatives='r_n_xtrend'),
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
