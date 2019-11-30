cimport numpy as np

cdef class ContactSurface:

    cdef public np.float_t z, zmax, zmin, zmean
    cdef public np.float_t [:,:] surface
