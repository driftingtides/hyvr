
import numpy as np
from hyvr.geo.channel_utils import ferguson_curve
from hyvr.geo.channel import Channel

cimport cython
cimport numpy as np
from libc.math cimport sqrt, atan2
from hyvr.geo.ae_realization cimport AERealization


cdef class ChannelAE(AERealization):

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef create_object_arrays(self):
        # This is super ugly :(
        cdef int i, nx, ny, j, k, max_len_centerline
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
            z = 0.5*(self.top_surface.zmean + self.bottom_surface.zmean)
            zfactor = np.interp(z, [grid.z0, grid.zmax], *self.type_params['size_ztrend'])
        else:
            zfactor = 1
        depth = self.type_params['depth'] * zfactor
        zbottom = self.bottom_surface.zmean + depth*self.type_params['buffer']
        ztop = self.top_surface.zmean

        z = ztop - dz
        while z > zbottom:


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
                channel = Channel(self.type_params, x_center, y_center, vx, vy, z, width, depth, grid)
                self._add_object(channel)

                if self.type_params['mig'] is not None:
                    ystart[i] -= np.random.uniform(-self.type_params['mig'], self.type_params['mig'])

            z -= dz



    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    @cython.cdivision(True)
    cpdef maybe_assign_points_to_object(self, int oi,
                                        np.int32_t [:] geo_ids,
                                        np.float_t [:] angles,
                                        np.float_t x, np.float_t y, np.float_t z,
                                        int x_idx, int y_idx,
                                        Grid grid):
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
        cdef double xy_dist, dz, radius, dist, weight, sum_weights, vx_now, vy_now
        cdef double dist_along_curve, dist_along_curve_tmp
        cdef int nx, ny, i, j, n, closest_idx


        # if the point is above the channel top, don't consider it
        dz = z - self.object_ztop[oi]
        if dz > 0:
            geo_ids[0] = -1
            return

        # To check whether a point is inside the channel, we have to calculate
        # its distance to the centerline. This is quite costly, as it means we
        # have to iterate over the list of centerline points.
        # To avoid checks, after each distance calculation of a cell that does
        # belong to the channel, we try to find the indices of all points that
        # are closer to the current point than ``distance - channel width``.
        # These will then be set to 1 in ``dont_check``, so we can skip the
        # check for these.

        if self.object_dont_check[oi,x_idx, y_idx]:
            geo_ids[0] = -1
            return


        # Otherwise, get the distance to the centerline in xy-plane and z-direction
        # If we have dipsets, we will also directly calculate the velocity
        # (inverse distance weighting) and the distance along the curve
        sum_weights = 0
        vx_now = 0
        vy_now = 0
        dist_along_curve_tmp = 0
        dist_along_curve = 0
        xy_dist = 1e100
        if self.object_dipsets[oi]:
            for i in range(self.object_len_centerline[oi]):
                if i > 0:
                    dist_along_curve_tmp += sqrt((self.object_x_center[oi,i] - self.object_x_center[oi,i-1])**2
                                                +(self.object_y_center[oi,i] - self.object_y_center[oi,i-1])**2
                    )
                dist = sqrt((x - self.object_x_center[oi,i])**2 + (y - self.object_y_center[oi,i])**2)
                weight = 1/(dist + 1e-20)
                vx_now += self.object_vx[oi,i] * weight
                vy_now += self.object_vy[oi,i] * weight
                sum_weights += weight
                if dist < xy_dist:
                    dist_along_curve = dist_along_curve_tmp
                    xy_dist = dist
                    closest_idx = i
            vx_now /= sum_weights
            vy_now /= sum_weights
        else:
            # we only need the distance
            for i in range(self.object_len_centerline[oi]):
                dist = (x - self.object_x_center[oi,i])**2 + (y - self.object_y_center[oi,i])**2
                if dist < xy_dist:
                    xy_dist = dist
            xy_dist = sqrt(xy_dist)


        # Find other cells that are also not close enough
        radius = xy_dist - self.object_width[oi]
        if radius > self.object_min_dx_dy[oi]:
            # Idea: find number of cells in x-direction, then loop over these
            # and find y-directions for those.
            nx = int(radius/grid.dx)
            for i in range(nx):
                if x_idx-i >= 0 and x_idx+i < grid.nx:
                    # 'height' of the circle at distance x = i*dx from the center:
                    #
                    #     y(x) = sqrt(radius**2 - x**2)
                    #
                    ny = int(sqrt(radius**2 - (i*grid.dx)**2))
                    for j in range(ny):
                        if y_idx-j >= 0 and y_idx+j < grid.ny:
                            self.object_dont_check[oi,x_idx+i, y_idx+j] = 1
                            self.object_dont_check[oi,x_idx+i, y_idx-j] = 1
                            self.object_dont_check[oi,x_idx-i, y_idx+j] = 1
                            self.object_dont_check[oi,x_idx-i, y_idx-j] = 1

        # Now we calculate the value of the shape curve and check whether the
        # point is really inside (i.e. is above the parabola)
        if dz >= self.object_a[oi]*xy_dist**2 - self.object_depth[oi]:
            # it's inside: assign stuff
            geo_ids[2] = self.object_num_ha[oi]
            if z < self.object_ztop[oi] - self.object_depth[oi] + self.object_lag_height[oi]:
                angles[0] = 0
                angles[1] = 0
                geo_ids[0] = self.object_lag_facies[oi]
                return

            if self.object_dipsets[oi]:
                # azimuth from inverse distance weighted velocity
                angles[0] = atan2(vy_now, vx_now)/np.pi*180
                # dip as assigned
                angles[1] = self.object_dip[oi]

                # TODO: the distance along the curve is actually a bit
                # different, as the point is not directly orthogonal to the closest point.
                # However, if the point density on the centerline is high
                # enough, it won't really matter.
                # To correct for the distance in z-direction we subtract |dz| * cos(dip)
                # note that dz is negative here
                d = dist_along_curve * self.object_sin_dip[oi] + dz*self.object_cos_dip[oi] + self.object_shift[oi]
                n = int(d/self.object_layer_dist[oi])
                geo_ids[0] = self.object_facies_array[oi,n]
                return
            else:
                angles[0] = self.object_azim[oi]
                angles[1] = self.object_dip[oi]
                geo_ids[0] = self.object_facies[oi]
                return

        else:
            # it's not inside. This means lower points will also not be inside
            self.object_dont_check[oi,x_idx, y_idx] = 1
            geo_ids[0] = -1
            return