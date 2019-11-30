cimport numpy as np

cdef class ContactSurface:

    cdef public int nx, ny
    cdef public np.float_t z, zmax, zmin, zmean
    cdef public np.float_t [:,:] surface

    cpdef use_lower_surface_value(self, ContactSurface other_surface)
    cpdef use_higher_surface_value(self, ContactSurface other_surface)
