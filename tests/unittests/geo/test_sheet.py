# from hyvr.classes.grid import Grid
# from hyvr.classes.contact_surface import ContactSurface
# from hyvr.classes.sheet import Sheet
# import numpy as np

# def test_dipsets():

#     grid = Grid(0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 1.0, 1.0, 1.0, 0)
#     bottom_surface = ContactSurface(grid, mode='flat', z=0.0)
#     top_surface = ContactSurface(grid, mode='flat', z=1.0)


#     dip = 89.99
#     azim = 10
#     type_params = {
#             "structure":'dip',
#             "dipset_dist":0.1,
#             "facies":[5, 2],
#             "altfacies":[[2],[5]],
#             "dip":[dip,dip],
#             "azimuth":[azim,azim],
#     }

#     sheet = Sheet(type_params, bottom_surface, top_surface, grid)
#     sheet.num_ha = 0
#     facies = np.zeros(1, dtype=np.int32)
#     angles = np.zeros(2)
#     ids = np.zeros(3, dtype=np.int32)

#     x_idx = 2
#     y_idx = 2
#     sheet.maybe_assign_facies_azim_dip(facies, angles, ids,
#                                        grid.X[x_idx, y_idx, 0],
#                                        grid.Y[x_idx, y_idx, 0],
#                                        0.3, x_idx, y_idx, grid)
#     fac1 = facies[0]
#     assert angles[0] == azim
#     assert angles[1] == dip
#     assert facies[0] != -1


#     sheet.maybe_assign_facies_azim_dip(facies, angles, ids,
#                                        grid.X[x_idx, y_idx, 0],
#                                        grid.Y[x_idx, y_idx, 0],
#                                        0.3+0.1, x_idx, y_idx, grid)
#     assert facies[0] == fac1


#     x_idx += 1
#     sheet.maybe_assign_facies_azim_dip(facies, angles, ids,
#                                        grid.X[x_idx, y_idx, 0],
#                                        grid.Y[x_idx, y_idx, 0],
#                                        0.3, x_idx, y_idx, grid)
#     assert facies[0] != fac1


#     sheet.maybe_assign_facies_azim_dip(facies, angles, ids,
#                                        grid.X[x_idx, y_idx, 0],
#                                        grid.Y[x_idx, y_idx, 0],
#                                        1.3, x_idx, y_idx, grid)
#     assert facies[0] == -1

#     sheet.maybe_assign_facies_azim_dip(facies, angles, ids,
#                                        grid.X[x_idx, y_idx, 0],
#                                        grid.Y[x_idx, y_idx, 0],
#                                        -0.1, x_idx, y_idx, grid)
#     assert facies[0] == -1
