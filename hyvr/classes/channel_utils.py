import numpy as np
from hyvr.optimized import curve_interp
import matplotlib.pyplot as plt

def ferguson_curve(grid, h, k, ds, eps_factor, flow_angle, domain_length, xstart, ystart, width):
    """
    Simulate extruded parabola centrelines using the Ferguson (1976) disturbed meander model
    Implementation of AR2 autoregressive model
    http://onlinelibrary.wiley.com/doi/10.1002/esp.3290010403/full

    Parameters
    ----------
    grid : Grid object
    h : float
        Height (Ferguson model parameter)
    k : float
        Wave number (Ferguson model parameter)
    ds : float
        Curve distance for calculations (Ferguson model parameter)
    eps_factor : float
        Random background noise (Ferguson model parameter)
    flow_angle : float
        Angle of mean flow direction, in radians
    domain_length : float
        Length of the domain in mean flow direction
    xstart, ystart : float
        Starting coordinates of the channel centerline
    width : float
        Width of the channel

    Returns
    -------
    outputs : float matrix
        Simulated extruded parabola centerlines: storage array containing values for x coordinate, y
        coordinate, vx and vy
    """
    # Parameters
    ds += 1e-10
    ns = domain_length*2
    s = np.arange(0, ns, ds)

    # Centreline starting point
    xp = 0
    yp = 0

    # Calculate curve directions
    theta = ferguson_theta(s, eps_factor, k, h)

    # Interpolate curve direction over interval of interest
    s_interp, th_interp = curve_interp(s, theta, 0.1)

    # Storage array
    outputs = np.zeros((len(th_interp), 4))

    for th_idx, th_i in enumerate(th_interp):
        vx = ds*np.cos(th_i)
        vy = ds*np.sin(th_i)

        # Assign to storage array
        outputs[th_idx, 0] = xp       # x coordinate
        outputs[th_idx, 1] = yp       # y coordinate
        outputs[th_idx, 2] = vx       # vx
        outputs[th_idx, 3] = vy       # vy

        # go to next point
        xp += vx
        yp += vy

    # Rotate meanders into mean flow direction
    rot_angle = -np.mean(th_interp) + flow_angle
    rotMatrix = np.array([[np.cos(rot_angle), -np.sin(rot_angle)],
                          [np.sin(rot_angle),  np.cos(rot_angle)]])
    roro = np.dot(rotMatrix, outputs[:, 0:2].transpose())

    outputs[:, 2:] = np.dot(rotMatrix, outputs[:, 2:].transpose()).transpose()
    outputs[:, 0] = roro[0, :].transpose()
    outputs[:, 1] = roro[1, :].transpose()

    # move to channel start
    outputs[:, 0] += xstart
    outputs[:, 1] += ystart

    # plt.plot(outputs[:,0], outputs[:,1])
    # plt.axhline(y=grid.y0, xmin=grid.x0, xmax=grid.xmax, color='r')
    # plt.axhline(y=grid.ymax, xmin=grid.x0, xmax=grid.xmax, color='r')
    # plt.axvline(x=grid.x0, ymin=grid.y0, ymax=grid.ymax, color='r')
    # plt.axvline(x=grid.lx, ymin=grid.y0, ymax=grid.ymax, color='r')
    # plt.show()

    # Remove values before model domain
    indomain = np.logical_and(outputs[:, 0] >= grid.x0-width, outputs[:, 0] <= grid.xmax+width)
    outputs = outputs[indomain, :]

    # plt.plot(outputs[:,0], outputs[:,1])
    # plt.axhline(y=grid.y0, xmin=grid.x0, xmax=grid.xmax, color='r')
    # plt.axhline(y=grid.ymax, xmin=grid.x0, xmax=grid.xmax, color='r')
    # plt.axvline(x=grid.x0, ymin=grid.y0, ymax=grid.ymax, color='r')
    # plt.axvline(x=grid.lx, ymin=grid.y0, ymax=grid.ymax, color='r')
    # plt.show()

    return outputs


def ferguson_theta(s, eps_factor, k, h):
    """
    Calculate curve direction angle

    Parameters:
        s:				    Steps within generated curve distance
        eps_factor:		    Random background noise
        k:				    Wave number
        h:				    Height

    Returns:
        th_store *(array)* - Curve direction angle

    """
    # Storage arrays
    th_store = np.zeros(len(s))

    for idex, si in enumerate(s):
        if idex == 0:
            t1 = 0
            t2 = 0
            eps = 0
        elif idex == 1:
            t1 = th_store[idex-1]
            t2 = 0
            eps = np.random.normal()*eps_factor
        else:
            t1 = th_store[idex-1]
            t2 = th_store[idex-2]
            eps = np.random.normal(1)*eps_factor

        th_store[idex] = thetaAR2(t1, t2, k, h, eps)

    return th_store


def thetaAR2(t1, t2, k, h, eps):
    """
    Implementation of AR2 autoregressive model (Ferguson, 1976, Eq.15)
    http://onlinelibrary.wiley.com/doi/10.1002/esp.3290010403/full

    Parameters:
        t1: 	theta(i-1)
        t2: 	theta(i-2)
        k:		Wavenumber
        h:		Height
        eps:	Random background noise

    Returns:
        2nd-order autoregression (AR2)

    """
    b1 = 2*np.exp(-k*h)*np.cos(k*np.arcsin(h))
    b2 = -np.exp(-2*k*h)
    return eps + b1*t1 + b2*t2
