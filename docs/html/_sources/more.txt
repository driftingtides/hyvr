====================================================================
Extending HyVR
====================================================================

The HyVR package is a work in progress. It has been implemented in Python in order to make it accessible for researchers, easily customisable, and hopefully somewhat easy to extend to suit the needs of future research. In this section I have included some tips on ways to extend the HyVR source code.


------------------------------------------------------------------------
Adding more geometries
------------------------------------------------------------------------

HyVR has been set up in such a way to facilitate the implementation of additional architectural element geometries. 

In order to generate new types of geometries a new function needs to be written in the ``hyvr`` module that will be called from ``hyvr.hyvr_main()`` where individual architectural elements and hydrofacies are simulated (around line 188 of ``hyvr.sim.main()``). 

Any new geometry function needs to return the following properties:

* ``mat`` (numpy array)
* ``azim`` (numpy array)
* ``dip`` (numpy array)
* ``fac`` (numpy array)
* ``ae_arr_i`` (numpy array)


------------------------------------------------------------------------
The HyVR wish list
------------------------------------------------------------------------

Any modelling project will have 'areas for growth' (as opposed to weaknesses). I have identified some things that I would like HyVR to have, but that are outside of the scope of my PhD research (and funds...). Perhaps you have the time and need to complete these yourself?

* Extensions in ``C`` programming language to speed up bottlenecks in simulations, particularly in ``hyvr.scale_rotate``, ``hyvr.reindex``, and ``hyvr.planepoint``.
* Some level of conditioning, or improved interfacing with multiple-point geostatistical packages.
* Interaction of channels, as well as more complex/realisitic configurations of channel deposits (e.g. point bars).
* Utilities for deriving HyVR simulation parameters from transitional probability geostatistics.
* Simulation of chemofacies.