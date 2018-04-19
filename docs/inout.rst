.. _inout:

==========================================================
HYVR inputs and outputs
==========================================================

-----------------------------------
Parameter inputs and model outputs
-----------------------------------

Model outputs include three-dimensional fields with values at each model grid cell for the following parameters:

- **Strata**, ``ssm`` - strata identifier
- **Architectural element** ``ae`` - architectural element identifier
- **Hydrofacies assemblage** ``ha`` - Unique identifiers for each hydrofacies assemblage generated
- **Hydrofacies assemblage type** ``hat`` - Type of hydrofacies assemblage generated
- **Hydrofacies** ``fac`` - Type of hydrofacies
- **Isotropic hydraulic conductivity** ``k_iso``
- **Dip** ``dip`` - bedding parameter
- **Azimuth** ``azim`` - bedding parameter
- **Anisotropy ratio** ``ani_rat`` 
- **Porosity** ``poros``
- **Full hydraulic-conductivity tensor** ``ktensors`` - based on isotropic hydraulic conductivity, dip, azimuth and anisotropy ratio (:ref:`methods <tensorgen>`)


------------------------------------------------------------------------
The ``.ini`` configuration file
------------------------------------------------------------------------

The key piece of information requried by HYVR is the parameter file. This is an ``.ini`` configuration file that contains all the parameters required to run a HYVR simulation. 

The parameter file is separated into sections denoted by headers surrounded by brackets (e.g. ``[section_name]``). Parameters (or keys) and their associated values are then stipulated using the equals sign (e.g. ``key = value``). Each key is associated with the section in which it is located. Section and variable names should be in lower case. String values do nore require quotation marks in ``.ini`` files.

In HYVR there following sections are necessary:

*   ``[run]`` - This contains parameters related to the model run.
*   ``[model]`` - This contains details about the model dimensions 
*   ``[strata]`` - In this section the strata parameters are set.
*   ``[*architectural_elements*]`` - Each architectural element to be simulated in HYVR is included in its own section. Please see the subsection below for more information.
*   ``[hydraulics]`` - This section includes information for setting the hydraulic properties of simulated features.

An additional section ``[flowtrans]`` is included in some parameter files - this section included parameters used in groundwater flow and solute transport simulation packages not implemented in HYVR. ``.ini`` files are readable by a range of programming languages so the user can also store and read flow and transport parameters from the configuration file.

------------------------------------------------------------------------
Model setup sections
------------------------------------------------------------------------

^^^^^^^^^^^^^^^^^^^^^^
``[run]`` section
^^^^^^^^^^^^^^^^^^^^^^

This section contains general sections how the program should be run. It
controls where outputs are stored and which kind of output files are generated,
what kind of parameters are generated and how many realisations are generated.

HYVR simulations are structured in the following way:

Model -> Run -> Realisation

Typically we assume that the ``.ini`` file is stored in a model directory::

    mymodel/
    |-- myconfig.ini

When HyVR is run with this parameter file, it will create a run directory.
The name of the run directory can be set with the ``runname`` option, e.g. to ``myrun``.
HyVR then stores a copy of the ``.ini`` file inside the run directory under
the name ``myrun_parameters.ini``. If ``runname`` is not given, the filename of
the ini-file without ``.ini`` will be used. If you are only running one
realization, the output of this realization is also stored directly in this
directory::

    mymodel/
    |-- myconfig.ini
    |-- myrun/
    |   |-- myrun_parameters.ini
    |   |-- myrun_hyvr.dat
    |   |-- myrun_hyvr.vtk
    |   |...

If you are running multiple realizations, HyVR creates a subdirectory for
each realization output::

    mymodel/
    |-- myconfig.ini
    |-- myrun/
    |   |-- myrun_parameters.ini
    |   |-- real_001/
    |   |   |-- myrun_real_001_hyvr.dat
    |   |   |-- myrun_real_001_hyvr.vtk
    |   |-- real_002/
    |   |   |-- myrun_real_002_hyvr.dat
    |   |   |-- myrun_real_002_hyvr.vtk
    |   |...

**ATTENTION**: If the name of the ini-file ends with ``_parameters.ini`` we
assume the file was created automatically and is already in the run directory.
In this case the model directory will default to the directory one above the
directory of the ini-file and the run directory will be the directory of the
ini-file. Also, if the ``runname`` is not given in the ini-file, the part before
``_parameters.ini`` will be chosen as runname.

Example: If you run the file ``myrun_parameters.ini`` in the example above, the
run directory will be ``myrun`` and the runname will default to ``myrun`` if it
is not given.

The following settings are possible in the run section:

- ``runname``: *(optional)* Name of the model simulation run. Defaults to the name of the
  ini-file without the last for characters (``.ini``)
- ``numsim``: *(optional, default: 1)* Number of realisations to be generated.
- ``dataoutputs``: *(optional)* List of simulation output data formats (see `Model outputs
  <modelout>`), e.g. ``[vtk, mat, py]``
- ``modeloutputs``: *(optional)* List of simulation output formats for model input (see `Model
  outputs <modelout>`), e.g. ``[mf, hgs]``
- ``flag_ow``: *(optional, default: true)* Whether to overwrite previous model
  results. If ``true`` *(default)* model outputs are stored in the current
  directory and previous results will be overwritten. If ``false``, HyVR
  will check if the directory already exists and will ask you to change the
  runname in case it exists.
- ``anisotropy``: *(optional, default: true)* Generate anisotropy parameters?
- ``het``: *(optional, default: true)* Generate heterogeneity?

^^^^^^^^^^^^^^^^^^^^^^
``[model]`` section
^^^^^^^^^^^^^^^^^^^^^^

- ``dx``, ``dy``, ``dz``: *(required/optional)* Model grid cell dimensions. If
  ``dy`` or ``dz`` are not given, ``dx`` will be used instead.
- ``lx``, ``ly``, ``lz``: *(required)* Model domain dimensions.
- ``periodic``: *(optional, default: false)* Periodic model domain?
  (Sheets/truncated ellipsoids only)
- ``display``: *(optional, default: false)* 'Display'-type simulation? If this
  flag is set to ``true``, the simulated architectural elements are centred in
  the model domain so they can be viewed easily.
- ``hetlev``: *(required)* Hierarchical level at which heterogeneity should be
  simulated. Can be ``ae``, ``facies`` or ``internal``

^^^^^^^^^^^^^^^^^^^^^^
``[strata]`` section
^^^^^^^^^^^^^^^^^^^^^^

- ``ssm``: *(required)* List of sequence names. This should be a list of
  strings.
- ``ssm_top``: *(required)* List of mean strata contact elevations. This should
  be a list of floats of the same length as ``ssm``.
- ``ssm_contact_model``: *(required)* Statistical parameters for strata contact
  model. This can either be a list of floats of length 3, e.g. ``[0.05, 6, 6]``,
  or a list of the same length as ``ssm`` of lists of floats of length 3, e.g.
  ``[[0.05, 6, 6], [0.05, 5, 4], ...]``
- ``ssm_contact``: *(optional, default: flat)* Contact surface type, either flat,
  random, or user
- ``ae_table``: *(optional)* Relative filepath (starting from the modeldir) for
  a architectural element lookup table.
- ``ae``: List of architectural elements. This is a list of strings, which are
  the names of the ``[*architectural_elements*]`` sections below.
- ``ssm_ae``: *(required)* Which architectural elements are in each stratum.
  This should be a list of lists of strings. The outer list must have the same
  length as ``ssm``, the inner list can be of variable length. The elements of
  the inner lists must be strings from ``ae``.
- ``ae_prob``: *(required)* Probability of an architectural element occuring.
  This must be a list of lists of floats with the same shape as ``ssm_ae``.
- ``ae_z_mean``: *(required)* Mean thickness of architectural element unit. This
  must be a list of lists of floats with the same shape as ``ssm_ae``.
- ``avul_prob``: *(required)* Probability of avulsion. List of floats with the
  same length as ``ssm``.
- ``avul``: *(required)* Avulsion depth range. List of lists of floats. The
  outer list must have the same length as ``ssm``, the inner lists must be of
  length 2 and are the start and end point of the depth range.
- ``bg``: *(optional)* Background parameters for unassigned cells in the
  architectural elements. This should be three float values: facies, azimuth, and
  dip background values.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
``[element]`` sections for architectural elements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sections that describe architectural elements are entitled with an identifying
name (e.g. ``[sparse_scour]``). Note that section names should not include
spaces. The first parameter to be set it the ``geometry``. The current
implementation of HYVR includes three geometries: truncated ellipsoids
(``trunc_ellip``), channels (``channel``), and sheets (``sheet``).

Unless otherwise noted, ranges (``r_``) represent the lower and upper limits of
uniform distributions from which values are randomly generated.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
General ``[*element]`` parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``geometry``: *(required)* Geometry of hydrofacies assemblages within
  architectural element, either ``trunc_ellip``, ``ext_par``, or ``sheet``
- ``structure``: *(required)* Internal structure of hydrofacies assemblages.
  This can be ``massive`` or ``dip`` and also ``bulb``, ``bulb_l``, or
  ``random`` for truncated ellipsoids.
- ``contact``: *(required)* Type of bedding contact between element units.
  Either ``flat`` or ``random``.
- ``contact_model`` *(required)* Statistical parameters for bedding contact
  model. This should be a list of floats of length 3.
- ``facies``: *(required)* Hydrofacies included in hydrofacies assemblage. These
  are indices referring to ``[hydraulics].hydro`` (starting from 0).
- ``altfacies``: *(optional)* Alternating facies specification. This is a list of
  lists where the outer list has the same length as ``facies``.
- ``bg``: *(optional)* Background parameters for unassigned cells in the architectural element. This should be three float values: facies, azimuth, and dip background values.
- ``geo_ztrend``: *(optional)* Linear trend in geometry sizes with elevation.
  Given as a percentage change mulitplier in mean value from bottom to top of
  domain, i.e. :math:`[\lambda_{bottom}, \lambda_{top}]`
- ``k_ztrend``: *(optional)* Linear trend in isotropic hydraulic conductivity
  from bottom to top of domain :math:`[\xi_{bottom},\xi_{top}]`
- ``k_xtrend``: *(optional)* Linear trend in isotropic hydraulic conductivity from model inlet to outlet :math:`[\xi_{inlet},\xi_{outlet}]`
- ``n_ztrend``: *(optional)* Linear trend in porosity from bottom to top of domain :math:`[\xi_{bottom},\xi_{top}]`
- ``n_xtrend``: *(optional)* Linear trend in porosity from model inlet to outlet :math:`[\xi_{inlet},\xi_{outlet}]`

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Erosive element-specific parameters (truncated_ellipsoid, extruded parabola)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
- ``agg``: *(required)* Aggradation thickness added between each generation elevation. 
- ``buffer``: *(optional)* Buffer to reduce erosion of underlying units (see :ref:`methods <temethod>`).
- ``dipset_d``: *(optional)* Thickness of dipping internal structures.
- ``migrate``: *(optional)* Lateral migration of ellipsoid centrepoints drawn from a random normal distribution, given as mean and variance in :math:`x` and :math:`y` directions :math:`[\overline{\Delta x}, \sigma^2_{\Delta x}, \overline{\Delta y}, \sigma^2_{\Delta y}]`. 
- ``lag``: *(optional)* Parameters for lag surface *[lag thickness, hydrofacies ID]*
- ``dip``: *(required)* Range of the uniform distribution from which the dip will be randomly drawn.


.. _teparams:

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Truncated ellipsoid parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- ``el_z``: *(required)* Number of elements to be simulated per simulation elevation and layer area
- ``length``, ``width``, ``depth``: *(required)* Mean geometry of truncated ellipsoids
- ``paleoflow``: *(required)* Range of the uniform distribution from which the paleoflow orientation will be randomly drawn. 
- ``azimuth``: *(required)* Range of the uniform distribution from which the azimuth will be randomly drawn.
- ``bulbset_d``: *(optional)* Thickness of nested-bulb structures at the maximum depth of the truncated ellipsoid.
- ``te_xyz``: *(optional)* List of 3D coordinated for manually setting the
  centrepoint of truncated ellipsoids. This should be a list of lists. The inner
  lists must have length 3.

.. _chparams:

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Extruded parabola parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
TODO: add required/optional
- ``width``, ``depth`` -  Mean geometry of channel
- ``h`` - Extruded parabola centreline curve shape parameter
- ``k`` - Extruded parabola centreline curve shape wave number
- ``ds`` - Distance between centreline points along trajectory
- ``eps_factor`` - Variance of random fluctuations of channel centreline.
- ``channel_no`` - Number of Extruded parabolas to generate at each elevation
- ``dipset_d`` - Thickness of dipping internal structures.

.. _shparams:

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Sheet parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- ``lens_thickness`` - Thickness of individual sheets. If set to ``-1`` then no individual sheets are generated within each sheet architectural element unit.


------------------------------------------------------------------------
``[hydraulics]`` section
------------------------------------------------------------------------
The input parameters in this section are associated with the simulation of hydraulic parameters. It is also possible to only simulate the geometries of architectural elements and hydrofacies if required.

- ``gen``: *(optional, default: true)* Generate hydraulic parameters (i.e. hydraulic conductivity)?
- ``hydro``: *(required)* List of hydrofacies codes
- ``k_h``: Mean horizontal hydraulic conductivity. This must be either a float
  if it is the same for all hydrofacies, or a list of the same length as
  ``hydro``.
- ``sig_y`` - Variance of log hydraulic conductivity. This must be either a
  float if it is the same for all hydrofacies, or a list of the same length as
  ``hydro``.
- ``ycorlengths``: *(required)* Default correlation lengths for
  :math:`\log(K_{iso})` in each hydrofacies in :math:`x,y,z`-directions. This
  can be either a single float, if it's the same in all directions for all
  hydrofacies, a list of floats of length 3 if it's the same for all
  hydrofacies, or a list of lists of floats, where the outer list has the same
  length as ``hydro`` and the inner lists have length 3
- ``k_ratio``: *(required)* List of perpendicular anisotropy ratios (i.e
  :math:`\frac{K_h}{K_v}`) or single value if it's the same for all hydrofacies.
- ``n``: *(required)* List of mean porosity values or single value if it's the
  same for all hydrofacies.
- ``sig_n``: *(required)* Variance of porosity values. List of floats or single
  float if it's the same for all hydrofacies.
- ``ncorlengths``: *(required)* Default correlation lengths for porosity in each
  hydrofacies in :math:`x,y,z`-directions This can be either a single float, if
  it's the same in all directions for all hydrofacies, a list of floats of
  length 3 if it's the same for all hydrofacies, or a list of lists of floats,
  where the outer list has the same length as ``hydro`` and the inner lists have
  length 3


------------------------------------------------------------------------
``[flowtrans]`` section
------------------------------------------------------------------------
This section contains parameters to be used for groundwater flow and solute
transport simulations. This allows all input parameters for field generation and
subsequent modelling to be stored in the same ``.ini`` file.

- ``hin``: *(required)* boundary condition (head in). List of 3 floats
- ``hout``: *(required)* boundary condition (head out). List of 3 floats


.. _modelout:

-----------------------------------
Model outputs
-----------------------------------
HyVR has a number model outputs that can be set in the input parameter file. A
copy of the ``ini`` model parameter file is saved in the model directory
automatically. The following data output files include model outputs as
three-dimensional arrays:

- ``dat`` Python 'pickle' file - this is a native Python format that can be loaded into Python using ``hyvr.utils.load_pickle()``.
- ``mat`` MATLAB file
- ``vtr`` VTK rectilinear grid file -  this can be opened in ParaView for improved three-dimensional visualisation.
- ``h5`` HDF5 format 
- ``npz`` Numpy compressed format

HyVR can also create files that can be used as model inputs for some flow and transport modelling packages These currently include:

- MODFLOW-2005 - ``bas``, ``dis``, ``lpf``, ``nam``, ``oc``, and ``pcg`` model input files. Provided suitable flow and transport parameters are set in the ``[flowtrans]`` section of the input parameter file, this simulation can be executed. 
- MODFLOW 6 - ``dis``, ``nam``, and ``npf`` model input files. A complete set of MODFLOW 6 input files cannot be generated in HyVR at this stage.
- HydroGeoSphere - *K* tensors and porosity at each grid node. 

Note that these model inputs can only have regular model grids. They have not been tested for use in the above-named packages. 
