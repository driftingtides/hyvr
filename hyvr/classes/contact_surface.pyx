
import numpy as np
cimport numpy as np
import hyvr.utils as hu


cdef class ContactSurface:

    def __init__(self, grid, **kwargs):

        self.z = kwargs["z"]

        if kwargs['mode'] == "random":
            for key in ["z", "var", "corlx", "corly"]:
                if not key in kwargs:
                    raise ValueError(key + " is missing for Strata option 'random'!")
            self.surface = hu.specsim(grid,
                                   kwargs['var'],
                                   [kwargs['corlx'], kwargs['corly']],
                                   two_dim=True,
                                   covmod='gaussian')
            self.surface += self.z*np.ones((grid.nx, grid.ny))

            self.zmean = np.mean(self.surface)
            self.zmax = np.max(self.surface)
            self.zmin = np.min(self.surface)

        elif kwargs['mode'] == "flat":
            self.surface = self.z * np.ones((grid.nx, grid.ny))
            self.zmean = self.z
            self.zmin = self.z
            self.zmax = self.z

        elif self.mode == "dem":
            raise NotImplementedError("This is not implemented yet!")

