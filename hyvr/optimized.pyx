"""
This module contains optimized versions of routines that were once in sim or utils.

:Author: Samuel Scherrer
"""

import numpy as np
import math
cimport numpy as np
cimport cython


@cython.boundscheck(False)
@cython.wraparound(False)
def reindex(inray):
    """
    Reindex array from 0

    Parameters:
        inray:		Array with indices

    Returns:
        inray:      Array with altered indices. ATTENTION: This changes ``inray``!!

    """

    uniq = np.unique(inray)
    remat = dict(zip(uniq, range(len(uniq))))
    linear_inray = inray.reshape(inray.size)
    for i in range(len(linear_inray)):
        linear_inray[i] = remat[linear_inray[i]] + 1
    return inray


@cython.boundscheck(False)
@cython.wraparound(False)
def select_trough(np.ndarray[np.float_t, ndim=3] Xd,
                  np.ndarray[np.float_t, ndim=3] Yd,
                  np.ndarray[np.float_t, ndim=3] Zd,
                  np.float_t a=1.0, np.float_t b=1.0, np.float_t c=1.0,
                  np.float_t alpha=0.0):
    """
    Select cells that belong to the trough.

    Parameters
    ----------
    Xd, Yd, Zd : np.ndarray[np.float_t, ndim=3]
        x, y, and z distances of grid cells to trough center
    a, b, c : np.float
        Trough parameters (principal semi-axes lengths)
    alpha : np.float
        Rotation angle in degree

    Returns
    -------
    select : np.ndarray[np.bool, ndim=3]
        Mask that selects only points inside the trough
    """
    cdef np.float_t xi, yi, zi, cos, sin, alpha_rad, r2, max_ab
    cdef np.int_t nx, ny, nz, i, j, k
    nx = Xd.shape[0]
    ny = Xd.shape[1]
    nz = Xd.shape[2]

    # calculate cos and sin
    alpha_rad = alpha/180*np.pi
    cos = np.cos(alpha_rad)
    sin = np.sin(alpha_rad)

    # The trough is at most max(a, b) from the center
    max_ab = max(a, b)

    # assign mask
    cdef np.ndarray[np.uint8_t, ndim=3, cast=True] select = np.empty((nx, ny, nz), dtype=np.bool)
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                xi = Xd[i,j,k]
                yi = Yd[i,j,k]
                zi = Zd[i,j,k]
                if zi <= 0 and abs(xi) <= max_ab and abs(yi) <= max_ab:
                    r2 = (xi**2*cos**2 + 2*xi*yi*cos*sin + yi**2*sin**2)/a**2
                    r2 += (xi**2*sin**2 - 2*xi*yi*cos*sin + yi**2*cos**2)/b**2
                    r2 += zi**2/c**2
                    select[i,j,k] = (r2 <= 1)
                else:
                    select[i,j,k] = 0

    return select



@cython.boundscheck(False)
@cython.wraparound(False)
def set_anisotropic_ktensor(np.ndarray[np.float_t, ndim=5] ktensors,
            np.ndarray[np.float_t, ndim=3] k_iso,
            np.ndarray[np.float_t, ndim=3] azim,
            np.ndarray[np.float_t, ndim=3] dip,
            np.ndarray[np.float_t, ndim=3] anirat):
    """
    Calculate and set the anisotropic conductivity tensor for every grid point

    Parameters
    ----------
    ktensors : np.ndarray, dtype = np.float, shape = (nx, ny, nz, 3, 3)
        The matrix of conductivity tensors. This matrix will be altered in place.
    k_iso : np.ndarray, dtype = np.float, shape = (nx, ny, nz)
        The matrix of isotropic hydraulic conductivities
    azim : np.ndarray, dtype = np.float, shape = (nx, ny, nz)
        The matrix of azimuth angles (in radians)
    dip : np.ndarray, dtype = np.float, shape = (nx, ny, nz)
        The matrix of dip angles (in radians)
    anirat : np.ndarray, dtype = np.float, shape = (nx, ny, nz)
        The matrix of anisotropy ratios

    Returns
    -------
    None
    """
    cdef np.float_t kappa, psi, a, sink, cosk, sinp, cosp
    cdef np.float_t a11, a12, a13, a22, a23, a33
    cdef np.float_t kiso
    cdef np.int_t imax, jmax, kmax, i, j, k
    imax = k_iso.shape[0]
    jmax = k_iso.shape[1]
    kmax = k_iso.shape[2]
    for i in range(imax):
        for j in range(jmax):
            for k in range(kmax):
                kappa = azim[i,j,k]
                psi = dip[i,j,k]
                a = 1/anirat[i,j,k]
                sink = np.sin(kappa)
                cosk = np.cos(kappa)
                sinp = np.sin(psi)
                cosp = np.cos(psi)
                kiso = k_iso[i,j,k]
                a00 = kiso*(sink**2 + cosk**2*(cosp**2 + a*sinp**2))
                a01 = kiso*(1-a)*sink*cosk*sinp**2
                a02 = kiso*((a-1)*(np.sin(kappa+2*psi)-np.sin(kappa-2*psi)))/4
                a11 = kiso*((a-1)*sink**2*sinp**2+1)
                a12 = kiso*((1-a)*(np.cos(kappa-2*psi)-np.cos(kappa+2*psi)))/4
                a22 = kiso*(a*cosp**2+sinp**2)

                ktensors[i,j,k,0,0] = a00
                ktensors[i,j,k,0,1] = a01
                ktensors[i,j,k,0,2] = a02
                ktensors[i,j,k,1,0] = a01
                ktensors[i,j,k,1,1] = a11
                ktensors[i,j,k,1,2] = a12
                ktensors[i,j,k,2,0] = a02
                ktensors[i,j,k,2,1] = a12
                ktensors[i,j,k,2,2] = a22

# @cython.boundscheck(False)
# @cython.wraparound(False)
# def assign_between(np.ndarray[np.int16_t, ndim=3] strata_number,
                   # np.ndarray[np.float_t, ndim=2] cs_z_bottom,
                   # np.ndarray[np.float_t, ndim=2] cs_z_top,
                   # np.int_t n_strata,
                   # np.float_t z0, np.float_t dz):
    # """
    # Assigns all points below the contact surface to the strata by decreasing
    # its strata number by 1.

    # Parameters
    # ----------
    # strata_number : np.npdarray[np.int16_t, ndim=3]
        # Array of strata numbers. This will be altered.
    # cs_z_bottom : np.ndarray[np.float_t, ndim=2]
        # z-coordinate of the bottom contact surface at every point in the x-y-plane
    # cs_z : np.ndarray[np.float_t, ndim=2]
        # z-coordinate of the top contact surface at every point in the x-y-plane
    # n_strata : int
        # number that should be assigned
    # z0: np.float
        # z-coordinate of the lowest grid layer
    # dz : np.float
        # distance between grid layers
    # """
    # cdef np.int_t imax, jmax, kmax, below_idx, above_idx, i, j
    # cdef np.int16_t n_strata_16
    # imax = strata_number.shape[0]
    # jmax = strata_number.shape[1]
    # kmax = strata_number.shape[2]
    # n_strata_16 = np.int16(n_strata)
    # for i in range(imax):
        # for j in range(jmax):
            # z_bottom = cs_z_bottom[i,j]
            # z_top = cs_z_top[i,j]
            # # find index that is below/above contact surface
            # bottom_idx = int(math.ceil((z_bottom - (z0+0.5*dz))/dz))
            # top_idx = int(math.floor((z_top - (z0+0.5*dz))/dz))
            # for k in range(bottom_idx,top_idx+1):
                # strata_number[i,j,k] = n_strata_16



# @cython.boundscheck(False)
# @cython.wraparound(False)
# def trough_normalized_distance_squared(np.float_t x, np.float_t y, np.float_t z,
#                                        np.float_t a, np.float_t b, np.float_t c,
#                                        np.float_t sin, np.float_t cos):
#     return ((x*cos + y*sin)/a)**2 + ((-x*sin+y*cos)/b)**2 + (z/c)**2

# @cython.boundscheck(False)
# @cython.wraparound(False)
# def trough_normalized_vector(np.float_t x, np.float_t y, np.float_t z,
#                              np.float_t a, np.float_t b, np.float_t c,
#                              np.float_t sin, np.float_t cos):
#     cdef np.ndarray[np.float_t, ndim=1] d = np.empty(3, dtype=np.float)
#     d[0] = (x*cos + y*sin)/a
#     d[1] = (-x*sin + y*cos)/b
#     d[2] = z/c
#     return d


cdef int sign(double x):
    return (x > 0) - (x < 0)



@cython.boundscheck(False)
@cython.wraparound(False)
def curve_interp(double [:] xc, double [:] yc, double spacing):
    """
    Interpolate evenly spaced points along a curve. This code is based on code in an answer posted by 'Unutbu' on
    http://stackoverflow.com/questions/19117660/how-to-generate-equispaced-interpolating-values (retrieved 17/04/2017)

    Parameters:
        xc:			x coordinates of curve
        yc:			y coordinates of curve
        spacing:	Spacing between points

    Returns:
        - xn - x coordinates of interpolated points
        - yn - y coordinates of interpolated points

    """
    cdef int ic, j, i
    cdef double total_dist, tol
    cdef double [:] t
    cdef double [:] xn
    cdef double [:] yn
    cdef list idx

    t = np.arange(xc[0], len(xc), spacing * 0.1)
    xc = np.interp(t, np.arange(len(xc)), xc)
    yc = np.interp(t, np.arange(len(yc)), yc)
    tol = spacing
    ic, idx = 0, [0]
    while ic < len(xc):
        total_dist = 0
        j = ic + 1 # this is to make sure j is initialized
        for j in range(ic+1, len(xc)):
            total_dist += math.sqrt((xc[j] - xc[j-1]) ** 2 + (yc[j] - yc[j-1]) ** 2)
            if total_dist > tol:
                idx.append(j)
                break
        ic = j + 1

    xn = np.zeros(len(idx))
    yn = np.zeros(len(idx))
    for i, ic in enumerate(idx):
        xn[i] = xc[ic]
        yn[i] = yc[ic]

    return xn, yn
