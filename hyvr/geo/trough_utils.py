import numpy as np
import hyvr.utils as hu

def generate_trough_positions(bottom_surface, top_surface, type_params, grid):
    # TODO:
    # Idea: Generating the objects on planes between the maximum of the
    # bottom surface and the top surface is not exactly perfect imho.
    # Instead, we could maybe calculate the volume of the AE and use
    # this for generating the number of troughs, e.g.:
    #
    #     trough_density = #troughs/area = (#throughs * height) / volume
    #
    # The height should be the 'agg' parameter. Therefore, the number
    # of troughs could be generated according to a poisson
    # distribution:
    #
    #     n_troughs = np.random.poisson(lam=trough_density * volume/agg)
    #
    # For migrating troughs, we could then first assign the first half
    # in the "upper part" (what should this be?) and the second half as
    # migrated troughs based on the previous ones.

    # TODO:
    # Also randomly place objects a bit outside the domain if the domain is not
    # periodic. Therefore we need a and b

    # Randomly generated trough positions
    bottom_z = bottom_surface.zmax + type_params['depth'] * type_params['buffer']
    top_z = top_surface.zmean
    # go from top to bottom
    planes_z = np.arange(top_z, bottom_z, -type_params['agg'])

    # construct list of random trough centers
    te_xyz = []
    for i, z in enumerate(planes_z):
        if i == 0 or type_params['migrate'] is None:
            # randomly draw positions according to spatial uniform distribution
            lam = type_params['trough_density'] * grid.lx * grid.ly
            n_troughs = np.random.poisson(lam=lam)
            for k in range(n_troughs):
                xc = np.random.uniform(grid.x0, grid.xmax)
                yc = np.random.uniform(grid.y0, grid.ymax)
                te_xyz.append([xc, yc, z])

        else:   # migrate previous troughs
            # n_troughs stays the same
            for k in range(n_troughs):
                above_idx = k + (i-1)*n_troughs
                xc = te_xyz[above_idx][0] + np.random.normal(type_params['migrate'][0],
                                                             type_params['migrate'][1])
                yc = te_xyz[above_idx][1] + np.random.normal(type_params['migrate'][2],
                                                             type_params['migrate'][3])
                te_xyz.append([xc, yc, z])

    return te_xyz


def rand_trough_params(type_params, z, grid):
    if type_params['size_ztrend'] is not None:
        delta_ztrend = type_params['size_ztrend'][1] - type_params['size_ztrend'][0]
        zfactor = type_params['size_ztrend'][0] +\
            (z - grid.z0)*(delta_ztrend)/grid.lz
    else:
        zfactor = 1

    params = {}
    params['a'] = type_params['length'] * zfactor / 2
    params['b'] = type_params['width'] * zfactor / 2
    params['c'] = type_params['depth'] * zfactor
    params['alpha'] = np.random.uniform(*type_params['paleoflow'])

    return params
