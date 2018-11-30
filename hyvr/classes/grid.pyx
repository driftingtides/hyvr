import numpy as np
cimport numpy as np

cdef class Grid:
    """Simple grid class for access from cython functions"""

    # cpdef int get_x_index(self, double x_value):
    #     return int(np.round((x_value - (self.x0 + self.dx/2))/self.dx))

    # cpdef int get_y_index(self, double y_value):
    #     return int(np.round((y_value - (self.y0 + self.dy/2))/self.dy))

    # cpdef int get_z_index(self, double z_value):
    #     return int(np.round((z_value - (self.z0 + self.dz/2))/self.dz))

    def __cinit__(self, double x0, double y0, double z0,
                  double dx, double dy, double dz,
                  double lx, double ly, double lz,
                  int periodic):
        """
        Create the model grid.

        The model grid is a 3-d rectangular grid. The x-axis is assumed to be
        eastwards, the y-axis northwards, and the z-axis upwards.

        TODO: implement rotation along z-axis

        Parameters
        ----------
        x0, y0, z0 : float
            Coordinates of grid origin.
        dx, dy, dz : float
            Grid spacing
        lx, ly, lz : float
            Grid size
        periodic : bool
            Whether the domain should be periodic.
        """
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.lx = lx
        self.ly = ly
        self.lz = lz
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.xmax = self.x0 + self.lx
        self.ymax = self.y0 + self.ly
        self.zmax = self.z0 + self.lz
        self.periodic = periodic

        # create grid vectors
        self.x = np.arange(x0+dx/2, x0+lx, dx)
        self.y = np.arange(y0+dy/2, y0+ly, dy)
        self.z = np.arange(z0+dz/2, z0+lz, dz)

        self.nx = len(self.x)
        self.ny = len(self.y)
        self.nz = len(self.z)

        # create meshgrid
        # The order here is such that the final matrices have shape (nx, ny, nz)
        self.Y, self.X, self.Z = np.meshgrid(self.y, self.x, self.z)
        # centered grids are necessary for specsim
        self.X_centered = self.X - (np.min(self.X) + np.max(self.X))/2
        self.Y_centered = self.Y - (np.min(self.Y) + np.max(self.Y))/2
        self.Z_centered = self.Z - (np.min(self.Z) + np.max(self.Z))/2



