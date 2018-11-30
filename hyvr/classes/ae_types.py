"""
This module contains the architectural element types.
"""

from hyvr.classes.ae_realizations import *

class AEType:
    """
    Architectural element types are the kind of things that are read from the
    supplied ini-file. For example, in the ``made.ini`` testcase, there are the
    architectural element types "clay_sheet", "sand_sheet", "clay_lens",
    "crossbedded_scour", ...
    In HyVR, these are instances of the ``AEType`` class.

    The actual AEs that are realized in a HyVR simulation are instances of
    ``AERealization`` and are implemented in 'ae_realization.py'.

    To generate a realization, you should normally call the
    ``generate_realization`` method of a type instance.
    """

    def __init__(self, ae_section, name):
        """
        AE type constructor.

        Parameters
        ----------
        ae_section : dict
            Parsed element section from the ini-file.
        name : str
            The name of the element type (e.g. 'crossbedded_scour').
        """
        self.params = ae_section
        self.geometry = ae_section['geometry']
        self.name = name

    def generate_realization(self, bottom_surface, top_surface, stratum, grid):
        """
        Generates a realization of the AE type.

        For example, if the AE geometry is 'trunc_ellip', this generates a
        ``TruncEllip`` object, i.e. a sub-layer containing trough features
        within a stratum.

        Parameters
        ----------
        bottom_surface : ContactSurface object
            The bottom surface of the AE.
        top_surface : ContactSurface object
            The top surface of the AE.
        model : Model object
        """
        if self.geometry == 'trough':
            return TruncEllipAE(bottom_surface, top_surface, self.name, self.params, stratum, grid)
        elif self.geometry == 'channel':
            return ExtParAE(bottom_surface, top_surface, self.name, self.params, stratum, grid)
        elif self.geometry == 'sheet':
            return SheetAE(bottom_surface, top_surface, self.name, self.params, stratum, grid)
        else:
            raise ValueError('No such geometry in ' + self.name + ': ' + self.geometry)






