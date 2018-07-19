====================
Changelog/Bug Fixes
====================

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
