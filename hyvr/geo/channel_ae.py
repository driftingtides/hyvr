
import numpy as np
#import numpy.typing as npt
from hyvr.geo.channel_utils import ferguson_curve
from hyvr.geo.channel import Channel
from hyvr.geo.grid import Grid
from hyvr.geo.ae_realization import AERealization


class ChannelAE(AERealization):

    def create_object_arrays(self):
        # This is super ugly :(
        self.object_a = np.zeros(self.n_objects, dtype=np.float)
        self.object_width = np.zeros(self.n_objects, dtype=np.float)
        self.object_depth = np.zeros(self.n_objects, dtype=np.float)
        self.object_min_dx_dy = np.zeros(self.n_objects, dtype=np.float)
        self.object_len_centerline = np.zeros(self.n_objects, dtype=np.int32)
        self.object_ztop = np.zeros(self.n_objects, dtype=np.float)
        self.object_sin_dip = np.zeros(self.n_objects, dtype=np.float)
        self.object_cos_dip = np.zeros(self.n_objects, dtype=np.float)
        self.object_shift = np.zeros(self.n_objects, dtype=np.float)
        self.object_layer_dist = np.zeros(self.n_objects, dtype=np.float)
        self.object_lag_height = np.zeros(self.n_objects, dtype=np.float)
        self.object_lag_facies = np.zeros(self.n_objects, dtype=np.int32)
        self.object_dipsets = np.zeros(self.n_objects, dtype=np.int32)
        self.object_dip = np.zeros(self.n_objects, dtype=np.float)
        self.object_azim = np.zeros(self.n_objects, dtype=np.float)
        self.object_facies = np.zeros(self.n_objects, dtype=np.int32)

        if self.n_objects > 0:
            nx, ny = self.object_list[0].dont_check.shape
        else:
            nx, ny = 0, 0
        self.object_dont_check = np.zeros((self.n_objects, nx, ny), dtype=np.int32)

        for i, obj in enumerate(self.object_list):
            self.object_a[i] = obj.a
            self.object_width[i] = obj.width
            self.object_depth[i] = obj.depth
            self.object_min_dx_dy[i] = obj.min_dx_dy
            self.object_len_centerline[i] = obj.len_centerline
            self.object_ztop[i] = obj.ztop
            self.object_sin_dip[i] = obj.sin_dip
            self.object_cos_dip[i] = obj.cos_dip
            self.object_shift[i] = obj.shift
            self.object_layer_dist[i] = obj.layer_dist
            self.object_lag_height[i] = obj.lag_height
            self.object_lag_facies[i] = obj.lag_facies
            self.object_dipsets[i] = obj.dipsets
            for j in range(nx):
                for k in range(ny):
                    self.object_dont_check[i,j,k] = obj.dont_check[j,k]
            # if not self.object_dipsets[i]:
            #     self.object_azim[i] = obj.azim
            #     self.object_dip[i] = obj.dip
                #self.object_facies[i] = obj.facies
                
                

        # find length of centerline discretizations
        if self.n_objects > 0:
            max_len_centerline = np.max(self.object_len_centerline)
        else:
            max_len_centerline = 0
        self.object_x_center = np.zeros((self.n_objects, max_len_centerline), dtype=np.float)
        self.object_y_center = np.zeros((self.n_objects, max_len_centerline), dtype=np.float)
        self.object_vx = np.zeros((self.n_objects, max_len_centerline), dtype=np.float)
        self.object_vy = np.zeros((self.n_objects, max_len_centerline), dtype=np.float)

        for i, obj in enumerate(self.object_list):
            for j in range(self.object_len_centerline[i]):
                self.object_x_center[i,j] = obj.x_center[j]
                self.object_y_center[i,j] = obj.y_center[j]
                self.object_vx[i,j] = obj.vx[j]
                self.object_vy[i,j] = obj.vy[j]

    
    def generate_objects(self, grid):
        """
        Generate channel objects and place them in the domain.



        The following was just an idea I had, ignore it.

        The channels follow a mean flow angle that is given in the parameter
        file. From this angle, a flow axis is constructed as the "center axis"
        in direction of the flow angle. (This means the flow axis is chosen such
        that the projection of the (x-y-)domain onto the axis perpendicular to
        the flow axis is cut in half by the flow axis.)

        The channels are then created as random 1-d functions of the position
        on the flow axis, with random starting points (within the projection of
        the domain onto the axis perpendicular to the flow axis).
        This can then be rotated to get an approximation of the channel curve
        (a set of points which lie on the curve).

        On every level, 'channel_no' channels are constructed. They share the
        same flow angle but are shifted w.r.t. the flow axis.

        Its also possible to let the channels migrate in time, i.e. slightly
        change the flow angle and the distance between the channels (by
        changing their starting points).
        """

        n_channels = self.type_params['channel_no']

        flow_angle = self.type_params['flow_angle']*np.pi/180
        domain_width = np.abs(grid.lx * np.sin(flow_angle)) + np.abs(grid.ly * np.cos(flow_angle))
        domain_length = np.abs(grid.lx * np.cos(flow_angle)) + np.abs(grid.ly * np.sin(flow_angle))
        ystart = np.random.uniform(-domain_width/2, domain_width/2, n_channels)


        # get number of layers
        dz = self.type_params['agg']
        if self.type_params['size_ztrend'] is not None:
            z = 0.5*(self.zmax + np.mean(self.bottom_surface))
            zfactor = np.interp(z, [grid.z0, grid.zmax], *self.type_params['size_ztrend'])
        else:
            zfactor = 1
        depth = self.type_params['depth'] * zfactor
        zbottom = np.mean(self.bottom_surface) + depth*self.type_params['buffer']
        ztop = self.zmax

        z = zbottom
        ztop - dz
        while z < ztop - dz:


            # get current width and depth
            if self.type_params['size_ztrend'] is not None:
                zfactor = np.interp(z, [grid.z0, grid.zmax], self.type_params['size_ztrend'])
            else:
                zfactor = 1
            width = self.type_params['width'] * zfactor
            depth = self.type_params['depth'] * zfactor

            for i in range(n_channels):
                # get start of channel
                ystart_grid = (grid.y0 + grid.ly)/2 - (domain_length+width)/2 * np.sin(flow_angle) + ystart[i] * np.cos(flow_angle)
                xstart_grid = (grid.x0 + grid.lx)/2 - (domain_length+width)/2 * np.cos(flow_angle) - ystart[i] * np.sin(flow_angle)
                # generate centerline
                outputs = ferguson_curve(grid,
                                         self.type_params['h'],
                                         self.type_params['k'],
                                         self.type_params['ds'],
                                         self.type_params['eps_factor'],
                                         flow_angle,
                                         domain_length,
                                         xstart_grid,
                                         ystart_grid,
                                         width,
                )
                x_center = outputs[:,0]
                y_center = outputs[:,1]
                vx = outputs[:,2]
                vy = outputs[:,3]


                # generate channels
                azimuth = self.type_params['azimuth'][i]
                channel = Channel(self.type_params, x_center, y_center, vx, vy, z, width, depth, azimuth, grid)
                self._add_object(channel)

                if self.type_params['mig'] is not None:
                    ystart[i] -= np.random.uniform(-self.type_params['mig'], self.type_params['mig'])

            z += dz



    # @jit
    def maybe_assign_points_to_object(self, oi: int,
                                        #geo_ids,
                                        #angles: npt.NDArray[np.int32],
                                        xs,
                                        ys,
                                        zs,
                                        #x_idx: int,
                                        #y_idx: int,
                                        grid: Grid):
        """
        This function checks whether the current grid cell with given
        coordinates is inside the trough and assigns facies, azimuth and dip by
        altering the passed arrays.

        Parameters
        ----------
        oi : int
            object index
        geo_ids : np.int32 array of size 4
            This holds the geological indices, i.e. facies number, architectural
            element number, hydrofacies assemblage (ha) and hydrofacies assemblage
            type (hat).
        angles : double array
            array of size 2 that holds the azimuth and dip of the cell. This
            will also be altered.
        x, y, z : double
            cell coordinates
        x_idx, y_idx : int
            indices of x and y position in grid.
        grid : Grid object
        """


        # if the point is above the channel top, don't consider it
        dz = zs - self.object_ztop[oi]
        
        logic_z = dz <= 0
        # if dz > 0:
        #     geo_ids[0] = -1
        #     return

        # To check whether a point is inside the channel, we have to calculate
        # its distance to the centerline. This is quite costly, as it means we
        # have to iterate over the list of centerline points.
        # To avoid checks, after each distance calculation of a cell that does
        # belong to the channel, we try to find the indices of all points that
        # are closer to the current point than ``distance - channel width``.
        # These will then be set to 1 in ``dont_check``, so we can skip the
        # check for these.

        # if self.object_dont_check[oi,x_idx, y_idx]:
        #     geo_ids[0] = -1
        #     return


        # Otherwise, get the distance to the centerline in xy-plane and z-direction
        # If we have dipsets, we will also directly calculate the velocity
        # (inverse distance weighting) and the distance along the curve
        sum_weights = 0
        vx_now = 0
        vy_now = 0
        dist_along_curve_tmp = 0
        dist_along_curve = np.zeros(ys.shape)
        xy_dist = 1e100
        
        for i in range(self.object_len_centerline[oi]):
            #if i > 0:
            dist = np.sqrt((xs - self.object_x_center[oi,i])**2 + (ys - self.object_y_center[oi,i])**2)
            if self.object_dipsets[oi]:
                dist_along_curve_tmp += np.sqrt((self.object_x_center[oi,i] - self.object_x_center[oi,i-1])**2
                                            +(self.object_y_center[oi,i] - self.object_y_center[oi,i-1])**2
                )
                weight = 1/(dist + 1e-20)
                vx_now += self.object_vx[oi,i] * weight
                vy_now += self.object_vy[oi,i] * weight
                sum_weights += weight
                vx_now /= sum_weights
                vy_now /= sum_weights
                dist_along_curve = np.where(dist < xy_dist, dist_along_curve_tmp, dist_along_curve)
                #dist_along_curve = dist_along_curve_tmp
            xy_dist = np.where(dist < xy_dist, dist, xy_dist)
                
                    #closest_idx = i
            
        # else:
        #     # we only need the distance
        #     for i in range(self.object_len_centerline[oi]):
        #         dist = (xs - self.object_x_center[oi,i])**2 + (ys - self.object_y_center[oi,i])**2
        #         if dist < xy_dist:
        #             xy_dist = dist
        xy_dist = np.sqrt(xy_dist)


        # # Find other cells that are also not close enough
        # radius = xy_dist - self.object_width[oi]
        # if radius > self.object_min_dx_dy[oi]:
        #     # Idea: find number of cells in x-direction, then loop over these
        #     # and find y-directions for those.
        #     nx = int(radius/grid.dx)
        #     for i in range(nx):
        #         if x_idx-i >= 0 and x_idx+i < grid.nx:
        #             # 'height' of the circle at distance x = i*dx from the center:
        #             #
        #             #     y(x) = sqrt(radius**2 - x**2)
        #             #
        #             ny = int(np.sqrt(radius**2 - (i*grid.dx)**2))
        #             for j in range(ny):
        #                 if y_idx-j >= 0 and y_idx+j < grid.ny:
        #                     self.object_dont_check[oi,x_idx+i, y_idx+j] = 1
        #                     self.object_dont_check[oi,x_idx+i, y_idx-j] = 1
        #                     self.object_dont_check[oi,x_idx-i, y_idx+j] = 1
        #                     self.object_dont_check[oi,x_idx-i, y_idx-j] = 1

        # Now we calculate the value of the shape curve and check whether the
        # point is really inside (i.e. is above the parabola)
        logic_inside = dz >= self.object_a[oi]*xy_dist**2 - self.object_depth[oi]
        
        num_ha = np.where(logic_inside & logic_z, self.object_num_ha[oi], -1)
        
        
        # if dz >= self.object_a[oi]*xy_dist**2 - self.object_depth[oi]:
        #     # it's inside: assign stuff
        #     geo_ids[2] = self.object_num_ha[oi]


        if self.object_dipsets[oi]:
            # azimuth from inverse distance weighted velocity
            azim = np.where(logic_inside & logic_z, np.atan2(vy_now, vx_now)/np.pi*180, -1)
            # dip as assigned
            

            # TODO: the distance along the curve is actually a bit
            # different, as the point is not directly orthogonal to the closest point.
            # However, if the point density on the centerline is high
            # enough, it won't really matter.
            # To correct for the distance in z-direction we subtract |dz| * cos(dip)
            # note that dz is negative here
            d = dist_along_curve * self.object_sin_dip[oi] + dz*self.object_cos_dip[oi] + self.object_shift[oi]
            n = int(d/self.object_layer_dist[oi])
            facies = np.where(logic_inside & logic_z, self.object_facies_array[oi,n], -1)
            #return
        else:
            azim = np.where((logic_inside & logic_z), self.object_azim[oi], -1)
            #angles[1] = self.object_dip[oi]
            #geo_ids[0] = self.object_facies[oi]
            facies = np.where(logic_inside & logic_z, self.object_facies[oi], -1)
            #return
        
        dip = np.where(logic_inside & logic_z, self.object_dip[oi], -1)

        logic_lag = zs < (self.object_ztop[oi] - self.object_depth[oi] + self.object_lag_height[oi])
        azim = np.where(logic_lag & (azim > 0), 0, azim)
        dip = np.where(logic_lag & (dip > 0), 0, dip)
        facies = np.where(logic_lag & (facies > 0), self.object_lag_facies[oi], facies)
                
        return facies, azim, dip, num_ha