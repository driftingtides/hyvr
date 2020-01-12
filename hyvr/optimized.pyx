"""
This module contains optimized versions of routines that were once in sim or utils.

:Author: Samuel Scherrer
"""

import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport sin, cos, sqrt


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef set_anisotropic_ktensor(np.ndarray[np.float_t, ndim=5] ktensors,
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
    cdef double kappa, psi, a, sink, cosk, sinp, cosp
    cdef double a11, a12, a13, a22, a23, a33
    cdef double kiso
    cdef int imax, jmax, kmax, i, j, k
    imax = k_iso.shape[0]
    jmax = k_iso.shape[1]
    kmax = k_iso.shape[2]
    for i in range(imax):
        for j in range(jmax):
            for k in range(kmax):
                kappa = azim[i,j,k]
                psi = dip[i,j,k]
                a = 1/anirat[i,j,k]
                sink = sin(kappa)
                cosk = cos(kappa)
                sinp = sin(psi)
                cosp = cos(psi)
                kiso = k_iso[i,j,k]
                a00 = kiso*(sink**2 + cosk**2*(cosp**2 + a*sinp**2))
                a01 = kiso*(1-a)*sink*cosk*sinp**2
                a02 = kiso*((a-1)*(sin(kappa+2*psi)-sin(kappa-2*psi)))/4
                a11 = kiso*((a-1)*sink**2*sinp**2+1)
                a12 = kiso*((1-a)*(cos(kappa-2*psi)-cos(kappa+2*psi)))/4
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


cdef int sign(double x):
    return (x >= 0) - (x < 0)




@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
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
    cdef int ic, j, i, nc
    cdef double total_dist, tol
    cdef double [:] t
    cdef double [:] xn, yn
    cdef list idx

    nc = len(xc)

    t = np.arange(xc[0], nc, spacing * 0.1)
    xc = np.interp(t, np.arange(nc), xc)
    yc = np.interp(t, np.arange(nc), yc)
    tol = spacing
    ic, idx = 0, [0]
    while ic < nc:
        total_dist = 0
        j = ic + 1 # this is to make sure j is initialized
        for j in range(ic+1, nc):
            total_dist += sqrt((xc[j] - xc[j-1]) ** 2 + (yc[j] - yc[j-1]) ** 2)
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
