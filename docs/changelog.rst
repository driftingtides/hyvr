====================
Changelog/Bug Fixes
====================

HyVR 1.1.0
----------

Release Date: 

Contributors
""""""""""""

* Samuel Scherrer
* Jeremy Bennett
* Pablo Ortega Tong

Changes
"""""""

* more Cython -> code faster
* profiling setup
* coverage.py setup
* better tests
* better developers guide
* 'massive' is now a valid structure for troughs, as described in the documentation
* changed 'flat' structure in made.ini to 'massive'
* some other minor changes


HyVR 1.0.0
----------

Release Date: 06 December 2018

Contributors
""""""""""""

* Samuel Scherrer

Changes
"""""""

Channels ported to new code (introduced in version 0.9.9).


HyVR 0.9.9
----------

Release Date: 30 November 2018

Contributors
""""""""""""

* Samuel Scherrer
* Jeremy Bennett

Changes
"""""""
This release introduces major changes to how HyVR works behind the scenes.
This is a pre-release of the coming version 1.0.0. This will only be published
on github, since the channel generation code is not yet ported to the new
version. Version 1.0.0 will also contain a reimplementation of the channel
generation code and will also be published on PyPI.

These are the most important changes:

* Object oriented code: HyVR is now heavily object oriented. This will hopefully
  make it easier to add new geometries to HyVR, and will enable us to introduce
  features like contact surface conditioning in future versions.
* Better separation of grid and model: A HyVR model is now almost fully
  independent from the grid on which it should be evaluated. Instead of directly
  mapping objects onto the grid at the time of object creation, the model is
  first constructed as list of strata, which again hold lists to their
  architectural elements.
  Only after the full model is created we go over all grid points, decide to
  which object the point belongs, and assign values accordingly.
  This comes with some performance loss due to python slowness. We'll try to
  alleviate these in future releases.
* Changes to ``.ini``-file: Some changes have been made to the parameter file.
  These were partly necessary, e.g. the new contact models, and are partly only
  cosmetic changes to get more verbose option names (e.g. ``strata`` instead of
  ``ssm``, ``ae_in_strata`` instead of ``ssm_ae``, etc.).
  Old-format ini-files will probably still work, but support for them will be
  dropped soon, so we urge all users to make the necessary changes to their
  ini-files. (It's really not that much work.)
* **No channels**: This version does not yet support channels. This will be part
  of version 1.0.0.



HyVR 0.2.3
----------

Release Date: 19 July 2018

Contributors
""""""""""""

* Jeremy Bennett
* Samuel Scherrer


Changes
"""""""
* Improvements to file/directory handling, including separation of parameter file parsing and directory setup. 
* New name for back-up parameter files.
* Flopy is now an optional dependency.
* Virtual boreholes can now be sampled on grids defined by number of boreholes in x,y-directions, or by a specific grid spacing.
* Journal article references added.



HyVR 0.2.2
----------

Release Date: 12 June 2018

Contributors
""""""""""""

* Jeremy Bennett
* Samuel Scherrer


Changes
"""""""

* Removed ``hyvr.utils.to_vtk`` function.
* HyVR now uses Flopy 3.2.9, and incorporates more of that package's features.
* Some changes to MODFLOW 6 utilities.
* HyVR can now be installed from PyPI using pip
* Improvements to h5 I/O.
* Improvements to ``hyvr.utils.virtual_boreholes`` function


HyVR 0.2.1
----------

Release Date: 9 May 2018

Contributors
""""""""""""

* Jeremy Bennett
* Samuel Scherrer
* Emilio Sanchez


Changes
"""""""

* Fixed bug in parsing of boolean options: previously all existing boolean
  options were parsed as ``True``
* Outputs for ParaView .vtr files are now specified with ``vtr`` instead of ``vtk`` as in previous versions.
* Some small changes to the testcase parameter file examples.
* Trends in porosity microstructure are now working.
* Architectural element lookup tables can now be saved to text files following simulation.
* Addition of ``virtual_boreholes`` function to ``HyVR.utils`` module. This can be used for generating borehole data from HyVR simulations.
* Some improvements to creation of MODFLOW 6 input files, including linear hydraulic head initial condition.
* Added testing functions.



HyVR 0.2
--------

Release Date: April 2018

Contributors
""""""""""""

* Jeremy Bennett
* Samuel Scherrer

Changes
"""""""

* First Release
