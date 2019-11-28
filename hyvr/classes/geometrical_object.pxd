cimport numpy as np
from hyvr.classes.grid cimport Grid

cdef class GeometricalObject:

    cdef public double zmin, zmax

    cpdef maybe_assign_facies_azim_dip(self, np.int32_t [:] facies,
                                       np.float_t [:] angles,
                                       np.int32_t [:] ids,
                                       double x, double y, double z,
                                       int x_idx, int y_idx, Grid grid)