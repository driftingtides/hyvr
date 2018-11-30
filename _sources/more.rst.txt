====================================================================
Extending HyVR
====================================================================

The HyVR package is a work in progress. It has been implemented in Python in order to make it
accessible for researchers, easily customisable, and hopefully somewhat easy to extend to suit the
needs of future research. In this section I have included some tips on ways to extend the HyVR
source code.


------------------------------------------------------------------------
Adding more geometries
------------------------------------------------------------------------

HyVR has been set up in such a way to facilitate the implementation of additional hydrofacies
assemblage geometries. 

In order to generate new types of geometries a new architectural element class has to be added
(probably in ``classes/ae_realizations.py``) and a new geometry class has to be generated (similar
to ``classes/trough.pyx``).

Each geometry class needs a constructor that sets all the random values, e.g. dip, azimuth, facies
etc. Additionally a assign-function is necessary that assigns facies, azimuth and dip to a grid cell
with given coordinates.

When adding geometries, we suggest reviewing the code for the existing geometries to see how it is
currently implemented. This should provide a reasonable idea of how to put a new geometry together.

------------------------------------------------------------------------
The HyVR wish list
------------------------------------------------------------------------

Any modelling project will have 'areas for growth' (as opposed to weaknesses). I have identified
some things that I would like HyVR to have, but that are outside of the scope of my PhD research
(and funds...). Perhaps you have the time and need to complete these yourself?

* Some level of conditioning, or improved interfacing with multiple-point geostatistical packages.
* Interaction of extruded parabolas, as well as more complex/realisitic configurations of channel deposits (e.g. point bars).
* Utilities for deriving HyVR simulation parameters from transitional probability geostatistics.
* Simulation of chemofacies.
* Better trough placement
* Horizontally rotated grid
