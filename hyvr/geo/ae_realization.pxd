"""
Header file for AERealization base class.

More info in ae_realization.pyx

:Author: Samuel Scherrer
"""

cimport numpy as np
from hyvr.geo.contact_surface cimport ContactSurface
from hyvr.geo.grid cimport Grid

cdef class AERealization:
    """
    Base class for AERealizations.
    """

    cdef public:
        np.int32_t type_id, num, bg_facies
        np.float_t zmin, zmax
        np.float_t bg_dip, bg_azim
        ContactSurface bottom_surface, top_surface
        dict type_params,
        str type_name
        # the number of objects is sometimes generated randomly, therefore
        # objects are first generated in a list. Object creation will happen in
        # the object classes. In the end, all object properties are stored in
        # fixed size arrays in the AERealization
        list object_list
        list object_zmax_list, object_zmin_list
        int n_objects
        np.float_t [:] object_zmins, object_zmaxs
        np.float_t [:] object_dip, object_azim
        np.int32_t [:] object_facies, object_num_ha

        # In case of alternating facies, these are
        # necessary. object_facies_array is a n_objects x max(num_facies)
        np.int32_t [:] object_num_facies
        np.int32_t [:,:] object_facies_array


    cpdef _create_common_object_arrays(self)
    cpdef create_object_arrays(self)
    cpdef maybe_assign_points_to_object(self, int oi,
                                        np.int32_t [:] geo_ids,
                                        np.float_t [:] angles,
                                        np.float_t x, np.float_t y, np.float_t z,
                                        int x_idx, int y_idx,
                                        Grid grid)
    