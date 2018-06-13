""" Grid class

    A module containing some classes and function useful to work with structured grids.

:Usage:

    Import as a normal python module.

:Version:

    0.1 , 01-09-2016 : Forked from Alessandro Comunian

:Authors:

    Jeremy P. Bennett

    Notes:
        The grids for the moment are always considered as 3D.

"""

import numpy as np

class Grid:
    """ Grid class

    A simple class that contains the *Origin*, *Delta* and *Size* of a
    simulation. It can also be used as container for some information
    contained in a VTK structured grid file.

    Notes:
        * By default, the size of the grid is considered as points.

        .. seealso::

            :py:mod:`vtknumpy`

    """

    def __init__(self,
                 ox=0.0, oy=0.0, oz=0.0,
                 dx=1.0, dy=1.0, dz=1.0,
                 nx=200, ny=200, nz=10,
                 gtype='points', gname='image',
                 periodicity=False):

        """
        Define a structured grid, with some default values.

        Parameters:

            ox, oy, oz (float): 	Coordinates of the origin [optional]
            dx, dy, dz (float): 	Distance between points of the structured grid or side of the cells [optional]
            nx, ny, nz (int): 		Number of points. The number of cells is nx-1, ny-1 and nz-1
            gtype (string): 		'points' or 'cells'. Define the grid type, for grids made of points or grids made of cells
            gname (string): 		Give the grid a name. This is intended for vtk outputs

        Returns:
            Size of the grid with _lx, _ly and _lz

        """

        self.ox = ox
        self.oy = oy
        self.oz = oz
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.gtype = gtype
        self.gname = gname
        self.periodicity = periodicity

        if self.gtype == 'points':
            self.points = self.nx*self.ny*self.nz
        elif self.gtype == 'cells':
            self.cells = (self.nx-1 if self.nx> 1 else 1)*(self.ny-1 if self.ny> 1 else 1)*(self.nz-1 if self.nz> 1 else 1)

        # Compute the size of the grid
        self._lx = None #self.dx * self.nx - self.ox
        self._ly = None #self.dy * self.ny - self.oy
        self._lz = None #self.dz * self.nz - self.oz

    def __str__(self):
        """
        Print some of the informations provided in the class.

        Parameters:
            self:	An instance of the Grid class

        Returns:
            out:	Grid dimensions containing origin, delta, size, etc.

        """
        out = (
            '    *** Grid dimensions ***\n'
            '    Origin: ( {0.ox:f}, {0.oy:f}, {0.oz:f})\n'
            '    Delta:  ( {0.dx:f}, {0.dy:f}, {0.dz:f})\n'
            '    Size:   ( {0.lx:f}, {0.ly:f}, {0.lz:f})\n'
            '    N:      ( {0.nx:d}, {0.ny:d}, {0.nz:d})\n'
            '    type:     {0.gtype}\n'
            '    points:   {0.points}\n'
            '    cells:    {0.cells}\n'
            '    name:     {0.gname}\n'
            ).format(self)

        return out

    def print_intervals(self, axis='xyz'):
        """
        Print the intervals that constitute the simulation domain in a format like:

            [ ox, ox+nx*dx] [ oy, oy+ny*dy] [ oz, oz+nz*dz]

        where *ox* is the origin, *nx* is the number of points and *dx*
        is the delta between points (*idem* for *y* and *z*).

        Parameters:
            axis (string):  containing ['x','y','z'], [optional]
                            If the default value "xyz" is used, then all the intervals are printed.

        Returns:
            print intervals that constitute the simulation domain

        """

        ox = self.ox; nx = self.nx; dx = self.dx
        oy = self.oy; ny = self.ny; dy = self.dy
        oz = self.oz; nz = self.nz; dz = self.dz

        if 'x' in axis:
            print("        interval along *x* [ ", ox, ", ", nx*dx + ox, "] ", end=' ')
        if 'y' in axis:
            print("        interval along *y* [ ", oy, ", ", ny*dy + oy, "] ", end=' ')
        if 'z' in axis:
            print("        interval along *z* [ ", oz, ", ", nz*dz + oz, "] ", end=' ')

        print('')

    def spacing(self):
        """
        To print out the spacing of the grid as a tuple

        Parameters:
            self:	An instance of the Grid class

        Returns:
            A tuple containing the spacing defined in the grid

        """

        return self.dx, self.dy, self.dz

    def origin(self):
        """
        To print out the origin of the grid as a tuple

        Parameters:
            self:	An instance of the Grid class

        Returns:
            A tuple containing the origin defined in the grid

        """

        return self.ox, self.oy, self.oz

    def shape(self):
        """
        To print out the shape of the grid as a tuple

        Parameters:
            self:			An instance of the Grid class
            gtype(string):  'points' or 'cells'. String to decide to print the shape in terms of points or in terms of cells.

        Returns:
            A tuple containing the shape defined in the grid

        """

        if self.gtype == 'points':
            out = (self.nx, self.ny, self.nz)
        elif self.gtype == 'cells':
            # By default, nx, ny and nz contains the size of the grid
            # in terms of points. When one dimension is 1, then along
            # that direction number of cells == number of points.
            out = (self.nx-1 if self.nx > 1 else 1,
                   self.ny-1 if self.ny > 1 else 1,
                   self.nz-1 if self.nz > 1 else 1)
        return out

    def vec(self):
        """
        Create vectors of spatial coordinates

        Parameters:
            self:		An instance of the Grid class

        Returns:
            xv, yv, zv - Vectors of spatial coordinates

        """

        xv = np.arange(self.dx / 2, self.lx, self.dx)
        yv = np.arange(-self.ly / 2 + self.dy / 2, self.ly / 2, self.dy)
        zv = np.arange(self.oz, self.lz + self.oz, self.dz)

        if self.ox != 0:
            xv = np.arange(self.ox, self.lx + self.ox, self.dx)
            yv = np.arange(self.oy, self.ly + self.oy, self.dy)
            zv = np.arange(self.oz, self.lz + self.oz, self.dz)

        return xv, yv, zv

    def vec_node(self):
        """
        Create vectors of spatial coordinates of bounding nodes

        Parameters:
            self:	An instance of the Grid class

        Returns:
            xv, yv, zv - Vectors of spatial coordinates of bounding nodes

        """

        xv = np.arange(self.ox, self.lx + self.ox + self.dx, self.dx)
        yv = np.arange(self.oy, self.ly + self.oy + self.dy, self.dy)
        zv = np.arange(self.oz, self.lz + self.oz + self.dz, self.dz)

        return xv, yv, zv

    def vec_x(self):
        """
        Create vector of spatial x-coordinate

        Parameters:
            self:	An instance of the Grid class

        Returns:
            xv - Vector of spatial x-coordinate

        """	
        if self.ox != 0:
            ov = self.ox
            lv = self.self.lx + self.ox
        else:
            ov = self.dx / 2
            lv = self.lx

        xv = ""
        for num in np.arange(ov, lv, self.dx):
            xv += str(num) + " "

        return xv

    def vec_y(self):
        """
        Create vector of spatial y-coordinate

        Parameters:
            self:	An instance of the Grid class

        Returns:
            yv - Vector of spatial y-coordinate

        """		
        if self.oy != 0:
            ov = self.oy
            lv = self.self.ly + self.oy
        else:
            ov = self.dy / 2
            lv = self.ly

        yv = ""
        for num in np.arange(ov, lv, self.dy):
            yv += str(num) + " "

        return yv


    def vec_z(self):
        """
        Create vector of spatial z-coordinate

        Parameters:
            self:	An instance of the Grid class

        Returns:
            zv - Vector of spatial z-coordinate

        """			
        if self.oz != 0:
            ov = self.oz
            lv = self.self.lz + self.oz
        else:
            ov = self.dz / 2
            lv = self.lz

        zv = ""
        for num in np.arange(ov, lv, self.dz):
            zv += str(num) + " "

        return zv

    def meshup(self, ind='ij'):
        """
        Create a meshgrid representation of the grid

        Parameters:
            self: 	An instance of the Grid class

        Returns:
            A tuple containing the x,y,z-coordinates of the grid

        """
        xv, yv, zv = self.vec()
        x_reg, y_reg, z_reg = np.meshgrid(xv, yv, zv, indexing=ind)

        return x_reg, y_reg, z_reg

    def meshup2d(self, ind='ij'):
        """
        Create a 2D meshgrid representation of the grid

        Parameters:
            self:	An instance of the Grid class

        Returns:
            A tuple containing the x,y-coordinates of the grid

        """

        xv, yv, _ = self.vec()
        x_reg, y_reg = np.meshgrid(xv, yv, indexing=ind)

        return x_reg, y_reg

    def cart_coords(self):
        """
        Get x,y,z coordinates in a tuple

        Parameters:
            self:	An instance of the Grid class

        Returns:
            A Tuple containing the x,y and z-coordinates of the grid

        """

        mgx, mgy, mgz = self.meshup()
        coords = np.column_stack((np.ravel(mgx, order='C'), np.ravel(mgy, order='C'), np.ravel(mgz, order='C')))

        return coords

    def cart_coords2d(self):
        """
        Get x,y coordinates in a tuple

        Parameters:
            self:	An instance of the Grid class

        Returns:
            A tuple containing the x,y-coordinates of the grid

        """

        mgx, mgy = self.meshup2d()
        coords = np.column_stack((mgx.flatten(), mgy.flatten()))

        return coords

    def compute_max(self):
        """
        Compute the max values for *x*, *y* and *z* of the grid.

        Parameters:
            self:	An instance of the Grid class

        Returns:
            Max values of x,y and z

        """
        self.x_max = self.ox + self.dx*self.nx
        self.y_max = self.oy + self.dy*self.ny
        self.z_max = self.oz + self.dz*self.nz

    def get_lx(self):
        """
        Provide as output a tuple containing the x-size of a
        grid. Useful for the implementation of 'property'.

        Parameters:
            self:	An instance of the Grid class			

        Returns:
            Tuple containing the size of a grid in x-direction

        """
        return int(self.dx * self.nx - self.ox)

    def get_ly(self):
        """
        Provide as output a tuple containing the y-size of a
        grid. Useful for the implementation of 'property'.

        Parameters:
            self:	An instance of the Grid class

        Returns:
            Tuple containing the size of a grid in y-direction

        """
        return self.dy * self.ny - self.oy

    def get_lz(self):
        """
        Provide as output a tuple containing the z-size of a
        grid. Useful for the implementation of 'property'.

        Parameters:
            self:	An instance of the Grid class

        Returns:
            Tuple containing the size of a grid in z-direction

        """
        return self.dz * self.nz - self.oz

    def set_lx(self, val=None):
        """
        Set the x-size of a grid

        Parameters:
            self:	An instance of the Grid class

        Returns:
            Tuple containing the size of a grid in x-direction

        """
        self._lx = self.dx * self.nx - self.ox

    def set_ly(self):
        """
        Set the y-size of a grid

        Parameters:
            self:	An instance of the Grid class

        Returns:
            Tuple containing the size of a grid in y-direction

        """
        self._ly = self.dy * self.ny - self.oy

    def set_lz(self):
        """
        Set the z-size of a grid

        Parameters:
            self:	An instance of the Grid class

        Returns:
            Tuple containing the size of a grid in z-direction

        """
        self._lz = self.dz * self.nz - self.oz

    def del_lx(self):
        """
        Delete the x-dimension of a grid

        Parameters:
            self:	An instance of the Grid class

        Returns:
            Tuple containing the size of a grid in x-direction

        """
        del self._lx

    def del_ly(self):
        """
        Delete the y-dimension of a grid

        Parameters:
            self:	An instance of the Grid class

        Returns:
            Tuple containing the size of a grid in y-direction

        """
        del self._ly

    def del_lz(self):
        """
        Delete the z-dimension of a grid

        Parameters:
            self:	An instance of the Grid class

        Returns:
            Tuple containing the size of a grid in z-direction

        """
        del self._lz

    lx = property(get_lx, set_lx, del_lx, "'size' along *x* of the grid.")
    ly = property(get_ly, set_ly, del_ly, "'size' along *y* of the grid.")
    lz = property(get_lz, set_lz, del_lz, "'size' along *z* of the grid.")

    def get_points(self):
        """
        Update the number of points for a points grid

        Parameters:
            self:	An instance of the Grid class

        Returns:
            Number of points for a points grid

        """
        return self.nx*self.ny*self.nz

    def set_points(self, val=None):
        """
        Set the number of points for a points grid

        Parameters:
            self:	An instance of the Grid class

        Returns:
            Number of points for a points grid

        """
        self._points = self.nx*self.ny*self.nz

    def del_points(self):
        """
        Delete the points of a points grid

        Parameters:
            self:	An instance of the Grid class

        Returns:
            --
			
        """
        del self._points

    points = property(get_points, set_points, del_points, "Number of points")

    def get_cells(self):
        """
        Update the number of cells for a cells grid

        Parameters:
            self:	An instance of the Grid class

        Returns:
            Number of cells for x, y and z

        """
        return \
            (self.nx-1 if self.nx>1 else 1)* \
            (self.ny-1 if self.ny>1 else 1)* \
            (self.nz-1 if self.nz>1 else 1)

    def set_cells(self, val=None):
        """
        Set the number of cells for a cells grid

        Parameters:
            self:	An instance of the Grid class

        Returns:
            Number of cells for x, y and z

        """	
        self._cells = \
            (self.nx-1 if self.nx>1 else 1)* \
            (self.ny-1 if self.ny>1 else 1)* \
            (self.nz-1 if self.nz>1 else 1)
 
    def del_cells(self):
        """
        Delete the cells for a cells grid

        Parameters:
            self:	An instance of the Grid class

        Returns:
            --
        """	
        del self._cells

    cells = property(get_cells, set_cells, del_cells, "Number of cells")

    def cellsize_2d(self):
        """
        Compute the cell size in 2D

        Parameters:
            self:	An instance of the Grid class

        Returns:
            Cell size in 2D

        """	
        return self.dx * self.dy

    cs2 = property(cellsize_2d, "Cell size in 2D")

    """ Some helpful things for getting grid indices """

    def idx_z(self, zval):
        """
        Return the index of an elevation for the k-dimension of a 3D array

        Parameters:
            zval:	Elevation value	

        Returns:
            iz *(int)* - Index of elevation

        """
        iz = np.around((zval - self.oz) / self.dz)
        return int(iz)

if __name__ == '__main__':

    # Create a grid with the default parameters
    grid1 = Grid()

    print(grid1)

    grid1.print_intervals('xz')

    print("    Spacing:", grid1.spacing())
    print("    Origin: ", grid1.origin())
    print("    Shape:  ", grid1.shape())
    print("    Max:  ", grid1.compute_max())


    print("Before update")
    print(grid1)
    grid1.dx = 1.1
    print("After update")
    print(grid1)


