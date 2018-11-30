
cdef class ContactSurface:

    cdef public double z, zmax, zmin, zmean
    cdef public double [:,:] surface
