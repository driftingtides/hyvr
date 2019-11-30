
cimport numpy as np
from hyvr.classes.contact_surface cimport ContactSurface
from hyvr.classes.grid cimport Grid
from hyvr.classes.ae_realization cimport AERealization


cdef class ChannelAE(AERealization):

    cdef public:
        # These are arrays of object attributes for better access from cython
        np.float_t [:,:] object_x_center, object_y_center
        np.float_t [:,:] object_vx, object_vy
        np.float_t [:] object_a, object_width, object_depth, object_min_dx_dy
        np.int32_t [:] object_len_centerline
        np.float_t [:] object_ztop
        np.float_t [:] object_sin_dip, object_cos_dip
        np.float_t [:] object_shift, object_layer_dist
        np.float_t [:] object_lag_height
        np.int32_t [:] object_lag_facies
        np.int32_t [:] dipsets
        np.int32_t [:,:,:] dont_check



    cpdef create_object_arrays(self)
    cpdef maybe_assign_points_to_object(self, int oi,
                                        np.int32_t [:] geo_ids,
                                        np.float_t [:] angles,
                                        double x, double y, double z,
                                        int x_idx, int y_idx,
                                        Grid grid)