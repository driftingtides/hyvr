cimport numpy as np

cdef struct FaciesAnglesIDs:
     np.int32_t facies
     np.float_t azim
     np.float_t dip
     np.int32_t aeID
     np.int32_t ha
     np.int32_t hat
