
import numpy as np
import hyvr.utils as hu

def contact_surface(grid, **kwargs):
    
    z = kwargs['z']
    if kwargs['mode'] == "random":
        for key in ["z", "var", "corlx", "corly"]:
            if not key in kwargs:
                raise ValueError(key + " is missing for Strata option 'random'!")
        surface = hu.specsim(grid,
                               kwargs['var'],
                               [kwargs['corlx'], kwargs['corly']],
                               two_dim=True,
                               covmod='gaussian')
        surface = surface + z

    elif kwargs['mode'] == "flat":
        surface = z * np.ones((grid.nx, grid.ny))

    elif kwargs['mode'] == "dem":
        raise NotImplementedError("This is not implemented yet!")
    
    return surface