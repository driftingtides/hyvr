
cimport numpy as np
from hyvr.geo.contact_surface cimport ContactSurface
from hyvr.geo.grid cimport Grid
from hyvr.geo.ae_realization cimport AERealization

cdef class SheetAE(AERealization):

    cdef public:
        # These are arrays of object attributes for better access from cython
        np.float_t [:] object_shift, object_layer_dist
        np.float_t [:] object_normvec_x, object_normvec_y, object_normvec_z
        np.float_t [:,:,:] object_bottom_surface, object_top_surface
        np.float_t [:] object_bottom_surface_zmean
        np.int32_t [:] object_dipsets

    cpdef create_object_arrays(self)
    cpdef maybe_assign_points_to_object(self, int oi,
                                        np.int32_t [:] geo_ids,
                                        np.float_t [:] angles,
                                        np.float_t x, np.float_t y, np.float_t z,
                                        int x_idx, int y_idx,
                                        Grid grid)