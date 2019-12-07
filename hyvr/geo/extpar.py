
import numpy as np
import hyvr.utils as hu

def gen_extpar(channel_type, model, ae, count, ani=True):
    """
    Generate extruded parabola geometries:

    * Flow regime is assumed to be reasonably constant so the major geometry of the extruded
      parabolas doesn't change so much
    * 'Migration' of the extruded parabolas according to a shift vector


    Parameters
    ----------
    extpar_type : ExtrudedParabolaType
    model : Model object
    ae : AERealization object
    ae : list
        AE parameters
        TODO: use AERealization instead
    count : int
        Material number and/or identifier
    ani : bool, optional (default: True)
        Whether anisotropy should be generated
        
    Returns
    -------
    count : int
        Material number and/or identifier
    """
    # This is just for convenience
    ae_arr = model.data['ae']
    ha_arr = model.data['ha']
    hat_arr = model.data['hat']
    fac = model.data['fac']
    azim = model.data['azim']
    dip = model.data['dip']

    x2 = model.X[:,:,0]
    y2 = model.Y[:,:,0]
    # TODO: replace z3 with model.Z
    _, _, z3 = np.meshgrid(range(0, model.nx), range(0, model.ny), range(0, model.nz), indexing='ij')          # 2-D grid

    # start location
    total_extpar = int(channel_type.channel_no)
    if model.display:
        # Place troughs in domain centre for display features
        xstart = (model.x0 + model.lx)/2 * np.ones(total_extpar)
        ystart = np.random.uniform(model.y0, model.y0+model.ly, total_extpar)
    else:
        # Randomly place curve starting points
        # TODO: this was the original setting for xstart, why was this hard-coded?
        # xstart = np.random.uniform(-10, 0, total_extpar)
        xstart = np.random.uniform(model.x0, model.xmax, total_extpar)
        ystart = np.random.uniform(model.y0, model.ymax, total_extpar)

    # loop over extruded parabola top depths
    ch_bot = ae.z_bottom
    ch_bot += channel_type.depth * channel_type.buffer

    ch_top = ae.z_top
    ch_dz = channel_type.agg
    ch_top += ch_dz

    for znow in np.arange(ch_bot, ch_top, ch_dz):
        hu.print_to_stdout('z = {:2.4f}m'.format(znow))
        # Assign linear trend to extruded parabola sizes
        if channel_type.geo_ztrend is not None:
            zfactor = np.interp(znow, [model.z0, model.z0 + model.lz], [channel_type.geo_ztrend[0], channel_type.geo_ztrend[1]])
        else:
            zfactor = 1
        z_ch_width = channel_type.width * zfactor
        z_ch_depth = channel_type.depth * zfactor

        if channel_type.channel_start is not None:
            cstart = channel_type.channel_start
        else:
            cstart = []
        # Loop over total extruded parabolas per system
        for chan in range(0, total_extpar):
            """ Loop over multiple extruded parabolas at 'timestep' """
            aha = ferguson_curve(model, channel_type.h, channel_type.k,  channel_type.ds,
                                 channel_type.eps_factor,
                                 disp=model.display,
                                 ch_start=cstart)

            """ Get flow direction in extruded parabolas for azimuth """
            # For periodicity shift trajectories into model unit cell
            if model.periodic:
                aha[aha[:, 1] < yvec[0], 1] += model.ly
                aha[aha[:, 1] > yvec[-1], 1] -= model.ly

            # initialize 2-D distance matrix
            D = 1e20 * np.ones_like(x2)

            # initialize sum of inverse-distance weights
            sumW = np.zeros(np.shape(x2), dtype=float)
            W = np.zeros(np.shape(x2), dtype=float)

            # initialize velocity orientation at this level
            vx_znow = np.zeros(np.shape(x2), dtype=float)
            vy_znow = np.zeros(np.shape(x2), dtype=float)

            # loop over all points of trajectory
            for ii in range(0, len(aha)):
                # distance to current point
                R = np.sqrt((x2 - aha[ii][0]) ** 2 + (y2 - aha[ii][1]) ** 2)

                # smallest distance of entire grid to all points so far
                D[R < D] = R[R < D]

                # inverse-distance weight for velocity interpolation
                W[:] = 1e-20
                W[R < z_ch_width / 2] = 1 / (R[R < z_ch_width / 2] + 1e-20)

                # velocity interpolation in 2-D
                vx_znow += aha[ii][2] * W
                vy_znow += aha[ii][3] * W
                sumW += W

            vx_znow /= sumW
            vy_znow /= sumW

            # Assign facies sets with dip values
            if sum(channel_type.dip) > 0:
                do, fd, dv, av = hu.dip_sets(model, channel_type, znow, curve=[aha[:, 0], aha[:, 1], vx_znow, vy_znow])
            else:
                do = np.ones((model.nx, model.ny, model.nz)) + count
                fd = np.ones((model.nx, model.ny, model.nz), dtype=int) * int(np.random.choice(channel_type.facies))
                dv = 0.0
                av = np.zeros((model.nx, model.ny, model.nz))

            """ Copy results into 3-D field """
            # Iterate over all nodes below current top elevation
            znow_index = int(np.around((znow - model.z0)/model.dz)) # TODO: fix around stuff
            zbelow_index = int(np.around((znow - z_ch_depth - model.z0)/model.dz))
            d_range = np.arange(max(0, zbelow_index), min(model.nz, znow_index))      # Depth range
            if len(d_range) > 0:        # Only compute if extruded parabola depth range is finite

                # Get mask arrays for each condition
                in_extpar = D[:, :, None]**2 <= z_ch_width**2 / 4 - ((znow_index - z3) * model.dz * z_ch_width / (z_ch_depth*2)) ** 2     # is grid cell in extruded parabola
                finite_v = ~np.isnan(vx_znow)            # Only assign if velocity is finite
                below_top = ae_arr <= ae.num          # Don't assign values to locations higher than top contact surface
                chan_mask = in_extpar * finite_v[:, :, None] * below_top
                if znow_index <= z3[:, :, -1].max():
                    # Set mask above top of extruded parabola to False
                    chan_mask[:, :, znow_index:-1] = False

                # Assign properties
                fac[chan_mask] = fd[chan_mask]
                ha_arr[chan_mask] = count
                hat_arr[chan_mask] = channel_type.ae_id
                ae_arr[chan_mask] = ae.num

                if channel_type.lag is not None:
                    in_lag = (znow - z_ch_depth + float(channel_type.lag[0])) > z3 * model.dz   # Is grid cell in extruded parabola
                    fac[np.logical_and(in_extpar, in_lag)] = int(channel_type.lag[1])

                if ani:
                    # calcuate azimuth, to 1 degree
                    azim2d = np.round((np.arctan2(vx_znow, vy_znow) - np.pi/2) * 180/np.pi)
                    azim3d = azim2d[:, :, None] * np.ones((model.nx, model.ny, model.nz))
                    azim[chan_mask] = azim3d[chan_mask]
                    dip[chan_mask] = dv

                count += 1

        # Shift starting values with migration vector from parameter file
        if channel_type.mig: # TODO: fix this
            xstart += np.random.uniform(-channel_type.mig[0], channel_type.mig[0])
            ystart += np.random.uniform(-channel_type.mig[1], channel_type.mig[1])

    return count


def ferguson_curve(model, h, k, ds, eps_factor, dist=0, disp=False, ch_start=[]):
    """
    Simulate extruded parabola centrelines using the Ferguson (1976) disturbed meander model
    Implementation of AR2 autoregressive model
    http://onlinelibrary.wiley.com/doi/10.1002/esp.3290010403/full

    Parameters
    ----------
    model : Model object
    h : float
        Height
    k : float
        Wave number
    ds : float
        Curve distance for calculations
    eps_factor : float
        Random background noise
    dist : float, optional (default: 0)
        Distance to generate curves - defaults to model.lx (if dist=0)
    disp : bool, optional (default: False)
        Creating display extruded parabola - extruded parabola begins at (0,0)
    ch_start : tuple of floats
        Starting location of channel (x,y coordinates)

    Returns
    -------
    outputs : numpy array 
        Simulated extruded parabola centerlines: storage array containing values for x coordinate, y
        coordinate, vx and vy
    """
    # Parameters
    ds += 1e-10
    if dist > 0:
        ns = dist
    else:
        ns = model.lx*100
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
        xp += vx
        yp += vy

        # Assign to storage array
        outputs[th_idx, 0] = xp       # x coordinate
        outputs[th_idx, 1] = yp       # y coordinate
        outputs[th_idx, 2] = vx       # vx
        outputs[th_idx, 3] = vy       # vy

    # Rotate meanders into mean flow direction
    mean_th = -np.mean(th_interp)
    rotMatrix = np.array([[np.cos(mean_th), -np.sin(mean_th)],
                          [np.sin(mean_th),  np.cos(mean_th)]])
    roro = np.dot(rotMatrix, outputs[:, 0:2].transpose())

    outputs[:, 2:] = np.dot(rotMatrix, outputs[:, 2:].transpose()).transpose()

    # Move starting location in x-direction
    outputs[:, 0] = roro[0, :].transpose() - np.random.uniform(model.lx/5, model.lx/10)
    outputs[:, 1] = roro[1, :].transpose()

    # Remove values before model domain
    if dist > 0:
        indomain = outputs[:, 0] >= model.x0
    else:
        indomain = np.logical_and(outputs[:, 0] >= model.x0, outputs[:, 0] <= model.lx)
    outputs = outputs[indomain, :]

    if len(ch_start) > 0:
        outputs[:, 0] += ch_start[0]
        outputs[:, 1] += ch_start[1]

    # Make sure streamlines begin within domain with respect to y
    yout = outputs[0, 1] > model.ly/4 or outputs[0, 1] > model.ly/4
    if disp is True:
        starty = outputs[0, 1]
    elif yout is True:
        starty = np.random.uniform(-model.ly/4, model.ly/4)
    else:
        starty = 0
    outputs[:, 1] = outputs[:, 1] - starty

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


def curve_interp(xc, yc, spacing):
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

    # TODO: what does this function do exactly? Can we use np.interp instead?
    N = len(xc)
    xc_in = np.array(xc)
    yc_in = np.array(yc)

    t = np.arange(xc[0], N, spacing * 0.1)
    xc = np.interp(t, np.arange(N), xc)
    yc = np.interp(t, np.arange(N), yc)
    tol = spacing
    ic, idx = 0, [0]
    while ic < N:
        total_dist = 0
        for j in range(ic+1, N):
            total_dist += np.sqrt((xc[j] - xc[j-1]) ** 2 + (yc[j] - yc[j-1]) ** 2)
            if total_dist > tol:
                idx.append(j)
                break
        ic = j + 1

    xn = xc[idx]
    yn = yc[idx]
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot(xc_in, yc_in, 'rx-')
    ax.scatter(xn, yn)
    # ax.set_aspect('equal')
    plt.show()

    return xn, yn
