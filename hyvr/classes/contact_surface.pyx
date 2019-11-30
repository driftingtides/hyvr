
import numpy as np
import hyvr.utils as hu

cimport cython
cimport numpy as np


cdef class ContactSurface:

    def __init__(self, grid, **kwargs):

        self.z = kwargs["z"]
        self.nx = grid.nx
        self.ny = grid.ny

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


    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef use_lower_surface_value(self, ContactSurface other_surface):
        cdef int i, j
        cdef np.float_t [:,:] other_surf = other_surface.surface
        for i in range(self.nx):
            for j in range(self.ny):
                if self.surface[i,j] > other_surf[i,j]:
                    self.surface[i,j] = other_surf[i,j]
        # using the higher of both values might have changed the minimum 
        self.zmin = np.min(self.surface)


    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef use_higher_surface_value(self, ContactSurface other_surface):
        cdef int i, j
        cdef np.float_t [:,:] other_surf = other_surface.surface
        for i in range(self.nx):
            for j in range(self.ny):
                if self.surface[i,j] < other_surf[i,j]:
                    self.surface[i,j] = other_surf[i,j]
        # using the higher of both values might have changed the maximum 
        self.zmax = np.max(self.surface)