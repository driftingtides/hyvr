"""
This module contains optimized versions of routines that were once in sim or utils.

:Author: Samuel Scherrer
"""


import numpy as np
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
def scale_rotate(np.ndarray[np.float_t, ndim=3] x,
                 np.ndarray[np.float_t, ndim=3] y,
                 np.ndarray[np.float_t, ndim=3] z,
                 np.float_t alpha=0.0, np.float_t a=1.0,
                 np.float_t b=1.0, np.float_t c=1.0):
    """
    Scale and rotate three-dimensional trough

    Parameters:
        x, y, z (float):	Spatial coordinates
        alpha (float):		Rotation angle about the z-axis
        a, b, c (float):	Axis lengths in x, y, z directions (ellipsoid length, width, depth)

    Returns:
        - select - Grid cells within ellipsoid
        - R2 - Grid of scaled and rotated values

    """

    cdef np.float_t xi, yi, zi, cos, sin, alpha_rad
    alpha_rad = alpha/180*np.pi
    cos = np.cos(alpha_rad)
    sin = np.sin(alpha_rad)
    cdef np.int_t imax, jmax, kmax
    imax = x.shape[0]
    jmax = x.shape[1]
    kmax = x.shape[2]
    cdef np.ndarray[np.float_t, ndim=3] R2 = np.empty((imax, jmax, kmax), dtype=np.float)
    cdef np.ndarray[np.uint8_t, ndim=3, cast=True] select = np.empty((imax, jmax, kmax), dtype=np.bool)
    for i in range(imax):
        for j in range(jmax):
            for k in range(kmax):
                xi = x[i,j,k]
                yi = y[i,j,k]
                zi = z[i,j,k]
                R2[i,j,k] = (xi**2*cos**2 + 2*xi*yi*cos*sin + yi**2*sin**2)/a**2
                R2[i,j,k] += (xi**2*sin**2 - 2*xi*yi*cos*sin + yi**2*cos**2)/b**2
                R2[i,j,k] += zi**2/c**2
                select[i,j,k] = (R2[i,j,k] <= 1 and zi <= 0)

    return select, R2


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
        The matrix of TODO: what exactly is this?

    Returns
    -------
    None
    """
    cdef np.float_t kappa, psi, a, sink, cosk, sinp, cosp
    cdef np.float_t a11, a12, a13, a22, a23, a33
    cdef np.float_t kiso
    cdef np.int_t imax, jmax, kmax
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

