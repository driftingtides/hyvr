"""
Contains the different AE realizations.
"""

import numpy as np
from hyvr.classes.trough import Trough
from hyvr.classes.trough_utils import *
from hyvr.classes.sheet import Sheet
from hyvr.classes.channel import Channel
from hyvr.classes.contact_surface import ContactSurface
from hyvr.classes.contact_surface_utils import *
from hyvr.classes.channel_utils import ferguson_curve

class AERealization:
    """
    AE realizations are the realization of a certain AE type in a HyVR
    simulation. For example, in the 'made.ini' test case, the lowest stratum
    consists of 'clay_sheet' architectural elements.
    Each of these sheets is an AE realization.

    This is the base class for AE realizations. Any instantiated object should
    be an instance of one of the (currently) three subtypes: ``TruncEllipAE``,
    ``ExtParAE``, and ``SheetAE``.

    Typically these objects should not be created directly, but via the
    ``generate_realization`` method of the respective ``AEType`` object.

    An AE realization contains a list of geometrical objects that are typical
    for this AE. For example, a TruncEllipAE contains a list of Trough objects.
    These objects should be implemented in separate files as cython cdef
    classes, similar to the Trough class.
    """

    def __init__(self, bottom_surface, top_surface, type_name, type_params, stratum, grid):
        """
        Parameters
        ----------
        bottom_surface : ContactSurface object
            The bottom surface of the AE.
        top_surface : ContactSurface object
            The top surface of the AE.
        type_name : str
            Name of the AE type (e.g. 'clay_sheet').
        type_params : dict
            Parameters for this AE type. These are the parameters that are read
            from the ini-file.
        """
        self.bottom_surface = bottom_surface
        self.top_surface = top_surface
        self.type_params = type_params
        self.type_name = type_name
        self.type_id = type_params['ae_id']

        # containers for objects and their depths
        self.objects = []
        self.object_zmins = []
        self.object_zmaxs = []

        # get initial minimum z
        self.zmin = self.bottom_surface.zmin
        self.zmax = self.top_surface.zmax

        # get background values
        if self.type_params['bg_facies'] != -1:
            self.bg_facies = self.type_params['bg_facies']
        else:
            self.bg_facies = stratum.bg_facies
        if not np.isnan(self.type_params['bg_azim']):
            self.bg_azim = self.type_params['bg_azim']
        else:
            self.bg_azim = stratum.bg_azim
        if not np.isnan(self.type_params['bg_dip']):
            self.bg_dip = self.type_params['bg_dip']
        else:
            self.bg_dip = stratum.bg_dip

        # generate objects
        self.generate_objects(grid)
        self.n_objects = len(self.objects)
        self.y_idx_cached = None


    def generate_objects(self, grid):
        """
        This is the method that should be used for generating geometrical
        objects of the right type. Override this in all of your subclasses.

        It's very important that the objects are stored in the order such that
        the object with the highest zmax is the first one, i.e. such that we
        can go from top to bottom by iterating over the list of objects.

        To add the objects to the object list, use the _add_object function below
        """
        pass

    def _add_object(self, object_):
        """
        Adds an object to the AE object list and updates the zmin and zmax
        values
        """
        self.objects.append(object_)
        self.object_zmins.append(object_.zmin)
        # update zmins if necessary
        if object_.zmin < self.zmin:
            self.zmin = object_.zmin
        self.object_zmaxs.append(object_.zmax)
        if object_.zmax > self.zmax:
            self.zmax = object_.zmax


class TruncEllipAE(AERealization):

    def generate_objects(self, grid):
        """
        Generate trough objects and place them in the domain.
        """

        # get trough positions
        if self.type_params['te_xyz'] is None:
            te_xyz = generate_trough_positions(
                self.bottom_surface, self.top_surface, self.type_params, grid
            )
        else:
            te_xyz = self.type_params['te_xyz']

        # Generate troughs
        for xc, yc, zc in te_xyz:
            trough_params = rand_trough_params(self.type_params, zc, grid)
            trough_params['x'] = xc
            trough_params['y'] = yc
            trough_params['z'] = zc
            trough = Trough(self.type_params, **trough_params)
            self._add_object(trough)





class SheetAE(AERealization):

    def generate_objects(self, grid):
        """
        Generate sheet objects and place them in the domain.
        """

        # A sheet AE can consist of either one massive sheet (i.e. the sheet is
        # just the full AE) or of multiple sheets.
        # Internally, the sheets can again be homogeneous or consist of dipping
        # sets.
        #
        # Sheets have a bottom and top surface that is generally a flat contact
        # surface. Only the topmost and the bottom-most sheets have potentially
        # non-flat contact surfaces as they use the contact surface of the AE.

        if self.type_params['lens_thickness'] == -1:
            # massive bedding
            thickness = self.zmax - self.zmin
        elif self.type_params['size_ztrend'] is not None:
            zfactor = np.interp(np.mean([self.zmin, self.zmax]),
                                [grid.z0, grid.zmax],
                                type_params['size_ztrend'])
            thickness = self.type_params['lens_thickness'] * zfactor
        else:
            thickness = self.type_params['lens_thickness']

        zbottom = self.zmax - thickness
        top_surface = self.top_surface
        # create all sheets except the lowest one
        last_sheet = False
        while not last_sheet:
            # generate bottom surface
            if zbottom > self.zmin:
                # normal sheet
                bottom_surface = ContactSurface(grid, mode='flat', z=zbottom)
                # it's possible that this bottom surface is above the top
                # surface if we're close to the AE top. Then we will use the
                # lower value (i.e. the top surface value)
                use_lower_surface_value(bottom_surface, top_surface)

                # close to the bottom it might also be possible that the AE
                # bottom is higher than the current bottom, so we have to use
                # the higher value
                use_higher_surface_value(bottom_surface, self.bottom_surface)

            else:
                last_sheet = True
                bottom_surface = self.bottom_surface

            # generate sheet object
            sheet = Sheet(self.type_params, bottom_surface, top_surface, grid)
            self._add_object(sheet)

            # use current bottom as new top, zbottom decreases
            zbottom -= thickness
            top_surface = bottom_surface


class ExtParAE(AERealization):


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
            zfactor = np.interp(z, [grid.z0, grid.zmax], *ch_par['size_ztrend'])
        else:
            zfactor = 1
        depth = self.type_params['depth'] * zfactor
        zbottom = self.bottom_surface.zmean + depth*self.type_params['buffer']
        ztop = self.top_surface.zmean

        z = ztop - dz
        while z > zbottom:


            # get current width and depth
            if self.type_params['size_ztrend'] is not None:
                zfactor = np.interp(z, [grid.z0, grid.zmax], *ch_par['size_ztrend'])
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
