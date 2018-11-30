import sys
import numpy as np
cimport cython
cimport numpy as np
from hyvr.classes.grid cimport Grid
from hyvr.classes.ae_realizations import SheetAE
import hyvr.utils as hu


@cython.boundscheck(False)
@cython.wraparound(False)
def assign_facies_azim_dip(model):
    """
    This function does all the dirty work of assigning grid cells to their
    stratum, ae, or geometrical object.
    """
    cdef int x_idx, y_idx, z_idx, total_n, current_n
    cdef int si, aei, oi
    cdef double x, y, z
    cdef double [:] angles
    cdef int [:] facies, ids
    cdef int [:,:,:] facies_array, ha_array, hat_array, ae_array
    cdef double [:,:,:] azim_array, dip_array

    hu.print_to_stdout("Assigning parameters to grid cells:")

    facies_array = np.zeros((model.grid.nx, model.grid.ny, model.grid.nz),
                                    dtype=np.int32)
    # number of current stratum
    strata_array = np.zeros((model.grid.nx, model.grid.ny, model.grid.nz),
                            dtype=np.int32)
    # This holds the number of the current AE starting from the top. This is
    # mostly for visualising the model.
    ae_array = -1*np.ones((model.grid.nx, model.grid.ny, model.grid.nz),
                        dtype=np.int32)
    # ha = hydrofacies assemblage
    # this number is basically a counter for the different geometrical objects
    # within the architectural elements
    ha_array = -1*np.ones((model.grid.nx, model.grid.ny, model.grid.nz),
                            dtype=np.int32)
    # hat = hydrofacies assemblage type
    # This holds the corresponding number of AEType to which this belongs
    hat_array = -1*np.ones((model.grid.nx, model.grid.ny, model.grid.nz),
                         dtype=np.int32)

    # angles
    azim_array = np.zeros((model.grid.nx, model.grid.ny, model.grid.nz),
                                    dtype=np.float)
    dip_array = np.zeros((model.grid.nx, model.grid.ny, model.grid.nz),
                                 dtype=np.float)


    # containers for passing facies, azimuth and dip, and ae, ha, and hat
    facies = np.empty(1, dtype=np.int32)
    angles = np.empty(2, dtype=np.float)
    ids = np.empty(3, dtype=np.int32)

    total_n = model.grid.nx * model.grid.ny
    # Loop over all grid cells, starting from the top, "towerwise"
    # The advantage of the towerwise looping is that we can remember in
    # which stratum/AE/object the previous cell was and know that the
    # current can only be in the same or below


    # TODO: make this a parallel loop
    for x_idx in range(model.grid.nx):
        for y_idx in range(model.grid.ny):
            current_n = model.grid.ny * x_idx + y_idx
            percentage_done = int(np.round(100*current_n/total_n))
            sys.stdout.write("Currently done: {:>3d} %\r".format(percentage_done))

            # Here we start to go through one tower from top to bottom
            x = model.grid.X[x_idx, y_idx, 0]
            y = model.grid.Y[x_idx, y_idx, 0]
            si = 0    # stratum index
            aei = 0   # AE index
            oi = 0    # object index

            for z_idx in range(model.grid.nz)[::-1]:
                z = model.grid.z[z_idx]
                # set default values
                facies[0] = -1
                angles[0] = np.nan
                angles[1] = np.nan
                ids[0] = -1
                ids[1] = -1
                ids[2] = -1

                # go through all strata and check whether the current cell is inside by
                # delegating the task to the strata, which then calls the AE, which calls the
                # individual geometrical objects
                # This part here is probably the performance-wise most critical part of the
                # whole program.
                # Therefore, any checks made by the stratum, AE and geometrical object will
                # first consist of very simple tests that don't require much computations, since
                # most of the time these will be enough to make sure that a cell is not inside
                # an object
                while si < model.n_strata:
                    if z < model.strata[si].zmin:
                        # if the current cell is below the current stratum, we can go to the
                        # next stratum
                        si += 1
                        aei = 0
                        oi = 0
                        continue
                    else:
                        # otherwise we have to check really to get an answer
                        aei, oi = stratum_maybe_assign_facies_azim_dip(
                            model.strata[si], facies, angles, ids, x, y, z,
                            x_idx, y_idx, aei, oi, model.grid
                        )
                        if facies[0] != -1:
                            facies_array[x_idx, y_idx, z_idx] = facies[0]
                            azim_array[x_idx, y_idx, z_idx] = angles[0]
                            dip_array[x_idx, y_idx, z_idx] = angles[1]
                            ae_array[x_idx, y_idx, z_idx] = ids[0]
                            ha_array[x_idx, y_idx, z_idx] = ids[1]
                            hat_array[x_idx, y_idx, z_idx] = ids[2]
                            strata_array[x_idx, y_idx, z_idx] = si
                            break
                        else:
                            si += 1
                            aei = 0
                            oi = 0
                            continue

                if si == model.n_strata:
                    print("Warning: Something went wrong: Cell can not be assigned to any stratum")
    model.data['facies'] = np.asarray(facies_array)
    model.data['azim'] = np.asarray(azim_array)
    model.data['dip'] = np.asarray(dip_array)
    model.data['strata'] = np.asarray(strata_array)
    model.data['ae'] = np.asarray(ae_array)
    model.data['ha'] = np.asarray(ha_array)
    model.data['hat'] = np.asarray(hat_array)
    hu.print_to_stdout("Done  ")


@cython.boundscheck(False)
@cython.wraparound(False)
def stratum_maybe_assign_facies_azim_dip(
        stratum, int [:] facies, double [:] angles, int[:] ids,
        double x, double y, double z, int x_idx, int y_idx,
        int aei, int oi, Grid grid):

    cdef int oi_orig, aei_orig

    oi_orig = int(oi)
    aei_orig = int(aei)
    while aei < stratum.n_ae:
        # if z > stratum.ae_zmaxs[aei]:
            # if the current z is above the maximum z of the current AE, all other AEs will also
            # be below the current cell, so we can stop searching for an AE
            # This case should not happen
            # break
        if z < stratum.ae_zmins[aei]:
            # if we're below the current AE, we can go search in the next one
            oi = 0
            aei += 1
            continue
        else:
            oi = ae_maybe_assign_facies_azim_dip(
                stratum.aes[aei], facies, angles, ids, x, y, z, x_idx, y_idx, oi, grid
            )
            # stratum.aes[aei].maybe_assign_facies_azim_dip(
            #     facies, angles, x, y, z, x_idx, y_idx, oi, grid
            # )
            if facies[0] != -1:
                # we found something
                return aei, oi
            else:
                # if we didn't find something here, we check the next AE
                aei += 1
                oi = 0
                continue

    # At this point we know that the cell is not in any AE and therefore not in this stratum
    facies[0] = -1
    return 0, 0

@cython.boundscheck(False)
@cython.wraparound(False)
def ae_maybe_assign_facies_azim_dip(
        ae, int [:] facies, double [:] angles, int [:] ids,
        double x, double y, double z, int x_idx, int y_idx,
        int oi, Grid grid):
        """
        Return facies, azim and dip if the cell is inside this AE.

        Parameters
        ----------
        ae : AERealization object
        facies : np.ndarray[np.int32, dim=1]
            Array of length one that hold the facies and will be altered by
            this function. If this is -1 afterwards, the cell is not in the
            stratum.
        angles : np.ndarray[np.float, dim=1]
            Array of length 2 that holds the azimuth and dip angles, will also
            be altered.
        ids : np.ndarray[np.int32, dim=1]
            Array of length three that holds AERealization-ID, ha, and hat
        x, y, z : float
            cell coordinates
        x_idx, y_idx : int
            indices of x and y position in grid.
        oi : int
            Object index of the last found object. This can be used as a
            starting point for the next cell.
        grid : Grid object

        Returns
        -------
        oi : int
            Object index of the last found object. This can be used as a
            starting point for the next cell.
        """
        cdef int oi_orig

        # First we check if the current point is already below the current
        # object. If this is the case, we check the next object. We do this
        # until either no objects can be checked anymore, or the next objects
        # are all below.
        # Since the objects are sorted in a way that the top-most objects come
        # first, we can just check the zmax of the current object. If it's
        # below the current z, the zmax of all following objects will also be
        # below z and we can therefore abort our search.

        oi_orig = int(oi)
        while oi < ae.n_objects:
            if z > ae.object_zmaxs[oi]:
                # we won't find any object below
                break
            elif z < ae.object_zmins[oi]:
                oi += 1
                continue
            else:
                ae.objects[oi].maybe_assign_facies_azim_dip(
                    facies, angles, ids, x, y, z, x_idx, y_idx, grid
                )
                if facies[0] != -1:
                    ids[0] = ae.num
                    ids[2] = ae.type_id
                    return oi
                else:
                    oi += 1
                    continue
        # at this point we know that the point does not belong to any object,
        # but it might still be in the background matrix
        # If this is a SheetAE, this is not possible
        # if isinstance(ae, SheetAE):
        #     return 0
        if isinstance(ae, SheetAE):
            return 0
        z_above = ae.top_surface.surface[x_idx, y_idx]
        z_below = ae.bottom_surface.surface[x_idx, y_idx]
        if z_below < z and z < z_above:
            facies[0] = np.int32(ae.bg_facies)
            angles[0] = float(ae.bg_azim)
            angles[1] = float(ae.bg_dip)
            ids[0] = ae.num
            ids[2] = ae.type_id
            # If we didn't find an object, but it belongs to the current AE, we can start the next
            # time from the same object
            return oi_orig
        else:
            # The cell does not belong to this AE
            return 0


