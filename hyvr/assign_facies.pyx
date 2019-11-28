# cython: profile=True
import sys
import numpy as np
cimport cython
cimport numpy as np
from hyvr.classes.grid cimport Grid
from hyvr.classes.ae_realizations import AERealization, SheetAE
from hyvr.classes.geometrical_object cimport GeometricalObject
import hyvr.utils as hu

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef assign_facies_azim_dip(model):
    """
    This function does all the dirty work of assigning grid cells to their
    stratum, ae, or geometrical object.
    """
    cdef int x_idx, y_idx, z_idx, loop_idx, total_n
    cdef int nx, ny, nz, n_strata
    cdef int si, aei, oi
    cdef double x, y, z
    cdef np.float_t [:] angles
    cdef np.int32_t [:] facies, ids
    cdef np.int32_t [:,:,:] facies_array, ha_array, hat_array, ae_array, strata_array
    cdef np.float_t [:,:,:] azim_array, dip_array
    cdef np.float_t [:,:,:] gridX, gridY
    cdef np.float_t [:] gridz
    cdef np.float_t [:] strata_zmins

    # getting grid information
    nx = model.grid.nx
    ny = model.grid.ny
    nz = model.grid.nz
    gridX = model.grid.X
    gridY = model.grid.Y
    gridz = model.grid.z
    n_strata = model.n_strata
    strata_zmins = model.strata_zmins

    hu.print_to_stdout("Assigning parameters to grid cells:")

    facies_array = np.zeros((nx, ny, nz), dtype=np.int32)

    # number of current stratum
    strata_array = np.zeros((nx, ny, nz), dtype=np.int32)

    # This holds the number of the current AE starting from the top. This is
    # mostly for visualising the model.
    ae_array = -1*np.ones((nx, ny, nz), dtype=np.int32)

    # ha = hydrofacies assemblage
    # this number is basically a counter for the different geometrical objects
    # within the architectural elements
    ha_array = -1*np.ones((nx, ny, nz), dtype=np.int32)

    # hat = hydrofacies assemblage type
    # This holds the corresponding number of AEType to which this belongs
    hat_array = -1*np.ones((nx, ny, nz), dtype=np.int32)

    # angles
    azim_array = np.zeros((nx, ny, nz), dtype=np.float)
    dip_array = np.zeros((nx, ny, nz), dtype=np.float)


    # container for passing facies, azimuth and dip, and ae, ha, and hat
    facies = np.empty(1, dtype=np.int32)
    angles = np.empty(2, dtype=np.float)
    ids = np.empty(3, dtype=np.int32)


    total_n = nx * ny

    # Loop over all grid cells, starting from the top, "towerwise"
    # The advantage of the towerwise looping is that we can remember in
    # which stratum/AE/object the previous cell was and know that the
    # current can only be in the same or below

    # TODO: make this a parallel loop
    for loop_idx in range(total_n):
        x_idx = loop_idx // ny
        y_idx = loop_idx % ny
        
        percentage_done = int(np.round(100*loop_idx/total_n))
        sys.stdout.write("Currently done: {:>3d} %\r".format(percentage_done))

        # Here we start to go through one tower from top to bottom
        x = gridX[x_idx, y_idx, 0]
        y = gridY[x_idx, y_idx, 0]
        si = 0    # stratum index
        aei = 0   # AE index
        oi = 0    # object index

        for z_idx in range(nz-1, -1, -1):
            z = gridz[z_idx]
            # set default values
            #facies_angles_ids.facies = -1
            #facies_angles_ids.azim = np.nan
            #facies_angles_ids.dip = np.nan
            #facies_angles_ids.aeID = -1
            #facies_angles_ids.ha = -1
            #facies_angles_ids.hat = -1
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
            while si < n_strata:
                if z < strata_zmins[si]:
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
                        #facies_array[x_idx, y_idx, z_idx] = facies_angles_ids.facies
                        #azim_array[x_idx, y_idx, z_idx] = facies_angles_ids.azim
                        #dip_array[x_idx, y_idx, z_idx] = facies_angles_ids.dip
                        #ae_array[x_idx, y_idx, z_idx] = facies_angles_ids.aeID
                        #ha_array[x_idx, y_idx, z_idx] = facies_angles_ids.ha
                        #hat_array[x_idx, y_idx, z_idx] = facies_angles_ids.hat
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

            if si == n_strata:
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
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef stratum_maybe_assign_facies_azim_dip(
        stratum, np.int32_t [:] facies, np.float_t [:] angles, np.int32_t[:] ids,
        double x, double y, double z, int x_idx, int y_idx,
        int aei, int oi, Grid grid):

    cdef int oi_orig, aei_orig, n_ae
    cdef np.float_t [:] ae_zmaxs = stratum.ae_zmaxs
    cdef np.float_t [:] ae_zmins = stratum.ae_zmins

    oi_orig = int(oi)
    aei_orig = int(aei)
    n_ae = stratum.n_ae
    while aei < n_ae:
        # if z > ae_zmaxs[aei]:
            # if the current z is above the maximum z of the current AE, all other AEs will also
            # be below the current cell, so we can stop searching for an AE
            # This case should not happen
            # break
        if z < ae_zmins[aei]:
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
@cython.nonecheck(False)
@cython.cdivision(True)
cpdef int ae_maybe_assign_facies_azim_dip(
        ae, np.int32_t [:] facies, np.float_t [:] angles, np.int32_t[:] ids,
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
        cdef int oi_orig, n_objects
        oi_orig = int(oi)
        n_objects = ae.n_objects

        cdef np.float_t [:] object_zmaxs = ae.object_zmaxs
        cdef np.float_t [:] object_zmins = ae.object_zmins

        cdef np.float_t [:,:] top_surface = ae.top_surface.surface
        cdef np.float_t [:,:] bottom_surface = ae.bottom_surface.surface
        cdef np.float_t z_above, z_below

        # First we check if the current point is already below the current
        # object. If this is the case, we check the next object. We do this
        # until either no objects can be checked anymore, or the next objects
        # are all below.
        # Since the objects are sorted in a way that the top-most objects come
        # first, we can just check the zmax of the current object. If it's
        # below the current z, the zmax of all following objects will also be
        # below z and we can therefore abort our search.

        while oi < n_objects:
            if z > object_zmaxs[oi]:
                # we won't find any object below
                break
            elif z < object_zmins[oi]:
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
        z_above = top_surface[x_idx, y_idx]
        z_below = bottom_surface[x_idx, y_idx]
        if z_below < z and z < z_above:
            facies[0] = ae.bg_facies
            angles[0] = ae.bg_azim
            angles[1] = ae.bg_dip
            ids[0] = ae.num
            ids[2] = ae.type_id
            # If we didn't find an object, but it belongs to the current AE, we can start the next
            # time from the same object
            return oi_orig
        else:
            # The cell does not belong to this AE
            return 0


