# "Header" file for Grid class
cdef class Grid:
    """Simple grid class for access from cython functions"""

    cdef public:
        double x0, y0, z0
        double lx, ly, lz
        double dx, dy, dz
        double xmax, ymax, zmax
        int nx, ny, nz
        int periodic
        double [:] x, y, z
        double [:,:,:] X, Y, Z, X_centered, Y_centered, Z_centered

    # cdef int get_x_index(self, double x_value)
    # cdef int get_y_index(self, double y_value)
    # cdef int get_z_index(self, double z_value)

