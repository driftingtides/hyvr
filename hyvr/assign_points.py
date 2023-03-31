#import sys
import numpy as np
import numpy.typing as npt
import hyvr.utils as hu

#cimport cython
#cimport numpy as np
from hyvr.geo.grid import Grid
from hyvr.geo.ae_realization import AERealization
from hyvr.geo.sheet_ae import SheetAE

from numba import njit, jit, prange

#@cython.boundscheck(False)
#@cython.wraparound(False)
#@cython.nonecheck(False)
#@cython.cdivision(True)

def assign_points(model):
    """
    This function does all the dirty work of assigning grid cells to their
    stratum, ae, or geometrical object.
    """
    # cdef int x_idx, y_idx, z_idx, loop_idx, total_n
    # cdef int nx, ny, nz, n_strata
    # cdef int si, aei, oi
    # cdef double x, y, z
    # cdef np.float_t [:] angles
    # cdef np.int32_t [:] geo_ids
    # cdef np.int32_t [:,:,:] facies_array, ha_array, hat_array, ae_array, strata_array
    # cdef np.float_t [:,:,:] azim_array, dip_array
    # cdef np.float_t [:,:,:] gridX, gridY
    # cdef np.float_t [:] gridz
    # cdef np.float_t [:] strata_zmins

    # getting grid information
    nx = model.grid.nx
    ny = model.grid.ny
    nz = model.grid.nz
    # gridX = model.grid.X
    # gridY = model.grid.Y
    # gridz = model.grid.z
    # n_strata = model.n_strata
    # strata_zmins = model.strata_zmins

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
    azim_array = np.zeros((nx, ny, nz), dtype=np.float64)
    dip_array = np.zeros((nx, ny, nz), dtype=np.float64)


    # container for passing facies, ae, ha, and hat, and azimuth and dip
    # geo_ids = np.empty(4, dtype=np.int32)
    # angles = np.empty(2, dtype=np.float64)

    # total_n = nx * nyW

    #Refactoring from the original HyVR model, here we loop only in the domain of each element:
    for si in range(model.n_strata):
        strata = model.strata[si]
        print(np.mean(strata.bottom_surface))
        strata_array[(model.grid.Z >= strata.bottom_surface[:,:,np.newaxis])] = si
        print((model.grid.Z >= strata.bottom_surface[:,:,np.newaxis]))
        for aei, zmin, zmax in zip(strata.aes, strata.ae_zmins, strata.ae_zmaxs):
            facies_array[(model.grid.Z >= aei.bottom_surface[:,:,np.newaxis])] = aei.bg_facies
            azim_array[(model.grid.Z >= aei.bottom_surface[:,:,np.newaxis])] = aei.bg_azim
            dip_array[(model.grid.Z >= aei.bottom_surface[:,:,np.newaxis])] = aei.bg_dip
            ae_array[(model.grid.Z >= aei.bottom_surface[:,:,np.newaxis])] = aei.num
            hat_array[(model.grid.Z >= aei.bottom_surface[:,:,np.newaxis])] = aei.type_id
            print(aei.n_objects)
            for oi in range(aei.n_objects):
                try:
                    xmin = aei.object_xmin[oi]
                    xmax = aei.object_xmin[oi]
                    ymin = aei.object_ymin[oi]
                    ymax = aei.object_ymax[oi]
                except (NameError, AttributeError) as e:
                    xmin = np.min(model.grid.X)
                    ymin = np.min(model.grid.Y)
                    xmax = np.max(model.grid.X)
                    ymax = np.max(model.grid.Y)
                
                logic = (model.grid.Z >= aei.bottom_surface[:,:,np.newaxis]) & (model.grid.X <= xmax) & (model.grid.X >= xmin) &\
                    (model.grid.Y <= ymax) & (model.grid.Y >= ymin)
                #if np.sum(logic) < 1:
                #    continue
                grid_elementX = model.grid.X[logic]
                grid_elementY = model.grid.Y[logic]
                grid_elementZ = model.grid.Z[logic]
                
                xs = grid_elementX
                ys = grid_elementY
                zs = grid_elementZ
                
                facies_oi, azim_oi, dip_oi, num_ha_oi = aei.maybe_assign_points_to_object(oi, xs,ys,zs, model.grid)
                
                facies_array[logic] = np.where(facies_oi != -1, facies_oi, facies_array[logic])
                azim_array[logic] = np.where(azim_oi != -1, azim_oi, azim_array[logic])
                dip_array[logic] = np.where(dip_oi != -1, dip_oi, dip_array[logic])
                ha_array[logic] = np.where(num_ha_oi != -1, num_ha_oi, ha_array[logic])
    
    model.data['facies'] = facies_array
    model.data['azim'] = azim_array
    model.data['dip'] = dip_array
    model.data['strata'] = strata_array
    model.data['ae'] = ae_array
    model.data['ha'] = ha_array
    model.data['hat'] = hat_array
    hu.print_to_stdout("Done  ")
    

#     # TODO: make this a parallel loop
#     for loop_idx in prange(total_n):
#         x_idx = loop_idx // ny
#         y_idx = loop_idx % ny
        
#         # percentage_done = int(np.round(100*loop_idx/total_n))
#         # sys.stdout.write("Currently done: {:>3d} %\r".format(percentage_done))

#         # Here we start to go through one tower from top to bottom
#         x = gridX[x_idx, y_idx, 0]
#         y = gridY[x_idx, y_idx, 0]
#         si = 0    # stratum index
#         aei = 0   # AE index
#         oi = 0    # object index

#         for z_idx in range(nz-1, -1, -1):
#             z = gridz[z_idx]
#             geo_ids[0] = -1
#             angles[0] = np.nan
#             angles[1] = np.nan
#             geo_ids[1] = -1
#             geo_ids[2] = -1
#             geo_ids[3] = -1


#             # go through all strata and check whether the current cell is inside
#             # by delegating the task to the strata, which then calls the AE,
#             # which calls the individual geometrical objects This part here is
#             # probably the performance-wise most critical part of the whole
#             # program.
#             # Therefore, any checks made by the stratum, AE and geometrical
#             # object will first consist of very simple tests that don't require
#             # much computations, since most of the time these will be enough to
#             # make sure that a cell is not inside an object
#             while si < n_strata:
#                 if z < strata_zmins[si]:
#                     # if the current cell is below the current stratum, we can
#                     # go to the next stratum
#                     si += 1
#                     aei = 0
#                     oi = 0
#                     continue
#                 else:
#                     # otherwise we have to check really to get an answer
#                     stratum = model.strata[si]
#                     aei, oi = maybe_assign_points_to_stratum(
#                         stratum.aes, stratum.n_ae,
#                         stratum.ae_zmaxs, stratum.ae_zmins,
#                         geo_ids, angles, x, y, z,
#                         x_idx, y_idx, aei, oi, model.grid
#                     )
#                     if geo_ids[0] != -1:
#                         facies_array[x_idx, y_idx, z_idx] = geo_ids[0]
#                         azim_array[x_idx, y_idx, z_idx] = angles[0]
#                         dip_array[x_idx, y_idx, z_idx] = angles[1]
#                         ae_array[x_idx, y_idx, z_idx] = geo_ids[1]
#                         ha_array[x_idx, y_idx, z_idx] = geo_ids[2]
#                         hat_array[x_idx, y_idx, z_idx] = geo_ids[3]
#                         strata_array[x_idx, y_idx, z_idx] = si
#                         break
#                     else:
#                         si += 1
#                         aei = 0
#                         oi = 0
#                         continue

#             if si == n_strata:
#                 print("Warning: Something went wrong: Cell can not be assigned to any stratum")

#     # we're done now, and can write the results back to the model
#     model.data['facies'] = np.asarray(facies_array)
#     model.data['azim'] = np.asarray(azim_array)
#     model.data['dip'] = np.asarray(dip_array)
#     model.data['strata'] = np.asarray(strata_array)
#     model.data['ae'] = np.asarray(ae_array)
#     model.data['ha'] = np.asarray(ha_array)
#     model.data['hat'] = np.asarray(hat_array)
#     hu.print_to_stdout("Done  ")


# def maybe_assign_points_to_stratum(
#     stratum_aes: list,
#     n_ae: int,
#     ae_zmaxs: npt.NDArray[np.float64],
#     ae_zmins: npt.NDArray[np.float64],
#     geo_ids: npt.NDArray[np.int32],
#     angles: npt.NDArray[np.float64],
#     x: float, y: float, z: float,
#     x_idx: int, y_idx: int,
#     aei: int, oi: int,
#     grid: Grid):
#     """
#     Tries to assign the point in question to an architectural element in the
#     stratum, if it is inside the stratum.
    
#     Parameters
#     ----------
#     stratum_aes : list of architectural element objects
#     n_ae : int
#         Number of AEs in stratum
#     geo_ids : np.int32 array of size 4
#         This holds the geological indices, i.e. facies number, architectural
#         element number, hydrofacies assemblage (ha) and hydrofacies assemblage
#         type (hat).
#     angles : np.float array of size 2
#         Holds azimuth and dip angle
#     x, y, z : float
#         x, y, and z position of the point that should be assigned
#     x_idx, y_idx : int
#         x and y indices in the grid
#     aei : int
#         Architectural element index. This is passed, so that we only have to
#         search in elements above the previous one for the right architectural
#         element.
#     oi : int
#         Object index.
#     grid : Grid object
#     """

#     #oi_orig = int(oi)
#     #aei_orig = int(aei)

#     while aei < n_ae:
#         # if z > ae_zmaxs[aei]:
#             # if the current z is above the maximum z of the current AE, all other AEs will also
#             # be below the current cell, so we can stop searching for an AE
#             # This case should not happen
#             # break
#         if z < ae_zmins[aei]:
#             # if we're below the current AE, we can go search in the next one
#             oi = 0
#             aei += 1
#             continue
#         else:
#             oi = maybe_assign_points_to_ae(
#                 stratum_aes[aei], geo_ids, angles, x, y, z, x_idx, y_idx, oi, grid
#             )
#             # stratum.aes[aei].maybe_assign_facies_azim_dip(
#             #     facies, angles, x, y, z, x_idx, y_idx, oi, grid
#             # )
#             if geo_ids[0] != -1:
#                 # we found something
#                 return aei, oi
#             else:
#                 # if we didn't find something here, we check the next AE
#                 aei += 1
#                 oi = 0
#                 continue

#     # At this point we know that the cell is not in any AE and therefore not in this stratum
#     geo_ids[0] = -1
#     return 0, 0


# def maybe_assign_points_to_ae(
#     ae: AERealization,
#     geo_ids: npt.NDArray[np.int32],
#     angles: npt.NDArray[np.float64],
#     x: float, y: float, z: float,
#     x_idx: int, y_idx: int,
#     oi: int,
#     grid: Grid):
#     """
#     Return facies, azim and dip if the cell is inside this AE.
#     Assigns values to point if point is in current ae

#     Parameters
#     ----------
#     ae : AERealization object
#     geo_ids : np.int32 array of size 4
#         This holds the geological indices, i.e. facies number, architectural
#         element number, hydrofacies assemblage (ha) and hydrofacies assemblage
#         type (hat).
#     angles : np.ndarray[np.float, dim=1]
#         Array of length 2 that holds the azimuth and dip angles, will also
#         be altered.
#     x, y, z : float
#         cell coordinates
#     x_idx, y_idx : int
#         indices of x and y position in grid.
#     oi : int
#         Object index of the last found object. This can be used as a
#         starting point for the next cell.
#     grid : Grid object

#     Returns
#     -------
#     oi : int
#         Object index of the last found object. This can be used as a
#         starting point for the next cell.
#     """
#     oi_orig = int(oi)
#     n_objects = ae.n_objects

#     object_zmaxs = ae.object_zmaxs
#     object_zmins = ae.object_zmins

#     top_surface = ae.top_surface.surface
#     bottom_surface = ae.bottom_surface.surface
    

#     # First we check if the current point is already below the current
#     # object. If this is the case, we check the next object. We do this
#     # until either no objects can be checked anymore, or the next objects
#     # are all below.
#     # Since the objects are sorted in a way that the top-most objects come
#     # first, we can just check the zmax of the current object. If it's
#     # below the current z, the zmax of all following objects will also be
#     # below z and we can therefore abort our search.

#     while oi < n_objects:
#         if z > object_zmaxs[oi]:
#             # we won't find any object below
#             break
#         elif z < object_zmins[oi]:
#             oi += 1
#             continue
#         else:
#             ae.maybe_assign_points_to_object(
#                 oi, geo_ids, angles, x, y, z, x_idx, y_idx, grid
#             )
#             if geo_ids[0] != -1:
#                 geo_ids[1] = ae.num
#                 geo_ids[3] = ae.type_id
#                 return oi
#             else:
#                 oi += 1
#                 continue
#     # at this point we know that the point does not belong to any object,
#     # but it might still be in the background matrix
#     # If this is a SheetAE, this is not possible
#     # if isinstance(ae, SheetAE):
#     #     return 0
#     if isinstance(ae, SheetAE):
#         return 0
#     z_above = top_surface[x_idx, y_idx]
#     z_below = bottom_surface[x_idx, y_idx]
#     if z_below < z and z < z_above:
#         geo_ids[0] = ae.bg_facies
#         angles[0] = ae.bg_azim
#         angles[1] = ae.bg_dip
#         geo_ids[1] = ae.num
#         geo_ids[3] = ae.type_id
#         # If we didn't find an object, but it belongs to the current AE, we can start the next
#         # time from the same object
#         return oi_orig
#     else:
#         # The cell does not belong to this AE
#         return 0


