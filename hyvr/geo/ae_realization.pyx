"""
Base class for AE realizations
"""

import numpy as np

cimport cython
cimport numpy as np
from hyvr.geo.contact_surface cimport ContactSurface

cdef class AERealization:
    """
    AE realizations are the realization of a certain AE type in a HyVR
    simulation. For example, in the 'made.ini' test case, the lowest stratum
    consists of 'clay_sheet' architectural elements. While 'clay_sheet' is an AE
    type, each of these sheets is an AE realization.

    This is the base class for AE realizations. Any instantiated object should
    be an instance of one of the (currently) three subtypes: ``TroughAE``,
    ``ChannelAE``, and ``SheetAE``.

    Typically these objects should not be created directly, but via the
    ``generate_realization`` method of the respective ``AEType`` object.

    An AE realization contains a list of geometrical objects that are typical
    for this AE. For example, a TroughAE contains a list of Trough objects.
    These objects should be implemented in separate files as cython cdef
    classes, similar to the Trough class.

    AERealizations are designed to be used from Python and from Cython code.
    Therefore, after construction the object list is converted to multiple
    arrays of object properties.

    For example, a TroughAE has a list of Trough objects, and each Trough object
    has a property ``a`` of type float (lenth of semi-major axis of ellipsoid).
    Therefore, the TroughAE object also has a float array called ``object_a`` as
    property.
    This makes access from Cyton code faster, because we can then iterate over
    float arrays, instead of iterating over Python lists.
    
    The derived AERealization types must provide the methods
    ``generate_object``, ``maybe_assign_points_to_objects`` and
    ``maybe_assign_points_to_object``.

    """

    def __init__(self,
                 ContactSurface bottom_surface,
                 ContactSurface top_surface,
                 str type_name,
                 dict type_params,
                 stratum,
                 Grid grid):
        """
        Parameters
        ----------
        bottom_surface : ContactSurface object
            The bottom surface of the AE.
        top_surface : ContactSurface object
            The top surface of the AE.
        type_name : str
            Name of the AE type (e.g. 'clay_sheet').
        type_params : dict
            Parameters for this AE type. These are the parameters that are read
            from the ini-file.
        stratum : Stratum object
            The stratum to which the AE belongs
            NOTE/TODO: It seems like we only need the background values from
                 this, might be nicer to pass only those
        grid : Grid object
            The HyVR grid
        """
        self.bottom_surface = bottom_surface
        self.top_surface = top_surface
        self.type_params = type_params
        self.type_name = type_name
        self.type_id = type_params['ae_id']

        # containers for objects and their depths
        self.object_list = []
        self.object_zmin_list = []
        self.object_zmax_list = []

        # get initial minimum z
        self.zmin = self.bottom_surface.zmin
        self.zmax = self.top_surface.zmax

        # get background values
        if self.type_params['bg_facies'] != -1:
            self.bg_facies = self.type_params['bg_facies']
        else:
            self.bg_facies = stratum.bg_facies
        if not np.isnan(self.type_params['bg_azim']):
            self.bg_azim = self.type_params['bg_azim']
        else:
            self.bg_azim = stratum.bg_azim
        if not np.isnan(self.type_params['bg_dip']):
            self.bg_dip = self.type_params['bg_dip']
        else:
            self.bg_dip = stratum.bg_dip

        # generate objects
        self.generate_objects(grid)

        # create object arrays
        self.n_objects = len(self.object_list)
        self._create_common_object_arrays()
        self.create_object_arrays()


    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef _create_common_object_arrays(self):
        """
        Creates object property arrays for common element properties.
        """

        cdef int i, max_num_facies, num_facies, j
        self.object_zmins = np.array(self.object_zmin_list)
        self.object_zmaxs = np.array(self.object_zmax_list)
        # get the rest from the object list
        self.object_azim = np.zeros(self.n_objects, dtype=np.float)
        self.object_dip = np.zeros(self.n_objects, dtype=np.float)
        self.object_facies = np.zeros(self.n_objects, dtype=np.int32)
        self.object_num_ha = np.zeros(self.n_objects, dtype=np.int32)
        self.object_num_facies = np.zeros(self.n_objects, dtype=np.int32)
        for i, obj in enumerate(self.object_list):
            self.object_azim[i] = obj.azim
            self.object_dip[i] = obj.dip
            self.object_facies[i] = obj.facies
            self.object_num_ha[i] = obj.num_ha
            self.object_num_facies[i] = obj.num_facies

        if self.n_objects > 0:
            max_num_facies = np.max(self.object_num_facies)
        else:
            max_num_facies = 0
        self.object_facies_array = np.zeros((self.n_objects, max_num_facies), dtype=np.int32)
        for i, obj in enumerate(self.object_list):
            num_facies = self.object_num_facies[i]
            for j in range(num_facies):
                self.object_facies_array[i,j] = obj.facies_array[j]



    def generate_objects(self, grid):
        """
        This is the method that should be used for generating geometrical
        objects of the right type. Override this in all of your subclasses.

        It's very important that the objects are stored in the order such that
        the object with the highest zmax is the first one, i.e. such that we
        can go from top to bottom by iterating over the list of objects.

        To add the objects to the object list, use the _add_object function below
        """
        raise NotImplementedError("You must override the method 'generate_objects' in subclasses of AERealization!")

    def _add_object(self, object_):
        """
        Adds an object to the AE object list and updates the zmin and zmax
        values
        """
        self.object_list.append(object_)
        self.object_zmin_list.append(object_.zmin)
        # update zmins if necessary
        if object_.zmin < self.zmin:
            self.zmin = object_.zmin
        self.object_zmax_list.append(object_.zmax)
        if object_.zmax > self.zmax:
            self.zmax = object_.zmax

    cpdef create_object_arrays(self):
        """
        This method is called after the object list has been created, so after
        all objects have been generated.
        
        Its main purpose is to create fixed size arrays of object attributes as
        members of the AERealization subtype, so that they can be accessed
        without needing to access the python list of objects.
        
        All subclasses of AERealization must override this method.
        """
        raise NotImplementedError("You must override the method 'create_object_arrays' in subclasses of AERealization!")


    cpdef maybe_assign_points_to_object(self, int oi,
                                        np.int32_t [:] geo_ids,
                                        np.float_t [:] angles,
                                        double x, double y, double z,
                                        int x_idx, int y_idx,
                                        Grid grid):
        """
        This function checks whether the current grid cell with given
        coordinates is inside the trough and assigns facies, azimuth and dip by
        altering the passed arrays.

        Parameters
        ----------
        oi : int
            object index
        geo_ids : np.int32 array of size 4
            This holds the geological indices, i.e. facies number, architectural
            element number, hydrofacies assemblage (ha) and hydrofacies assemblage
            type (hat).
        angles : double array
            array of size 2 that holds the azimuth and dip of the cell. This
            will also be altered.
        x, y, z : double
            cell coordinates
        x_idx, y_idx : int
            indices of x and y position in grid.
        grid : Grid object
        """
        raise NotImplementedError("You must override the method 'maybe_assign_points_to_object' in subclasses of AERealization!")
    
