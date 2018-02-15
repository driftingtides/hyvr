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

A number of key prefixes are used in the configuration files to assist identification when the HYVR package reads the configuration file:

* ``r_`` - This prefix denotes the lower and upper bounds of a range of float values.
* ``l_`` - This denotes a list of string values (no punctuation marks are required).
* ``ll_`` - This prefix denote a list of lists, each enclosed in brackets.
* ``flag_`` - keys with this prefix will be intepreted as true/false (boolean) values.

------------------------------------------------------------------------
Model setup sections
------------------------------------------------------------------------

^^^^^^^^^^^^^^^^^^^^^^
``[run]`` section
^^^^^^^^^^^^^^^^^^^^^^

- ``runname``				- Name of the model simulation run
- ``numsim`` 				- Number of realisations to be generated
- ``l_dataoutputs``			- Simulation output data formats [vtk, mat, py]
- ``l_modeloutputs``		- Simulation output model formats [mf, hgs]
- ``modeldir`` 				- filepath/directory for simulation outputs
- ``flag_ow``				- overwrite previous simulation with same name?
- ``flag_anisotropy``		- Generate anisotropy parameters?
- ``flag_het``				- Generate heterogeneity?

^^^^^^^^^^^^^^^^^^^^^^
``[model]`` section
^^^^^^^^^^^^^^^^^^^^^^

- ``dx``, ``dy``, ``dz``	- Model grid cell dimensions
- ``lx``, ``ly``, ``lz``	- Model domain dimensions
- ``flag_periodic``			- Periodic model domain? (Sheets/truncated ellipsoids only)
- ``flag_display``			- 'Display'-type simulation? If this flag is set to ``true``, the simulated architectural elements are centred in the model domain so they can be viewed easily.
- ``hetlev`` 				- Hierarchical level at which heterogeneity should be simulated *[ae, facies, internal]*

^^^^^^^^^^^^^^^^^^^^^^
``[strata]`` section
^^^^^^^^^^^^^^^^^^^^^^

- ``l_ssm``					- List of sequences
- ``r_ssm_top`` 			- List of mean strata contact elevations
- ``ll_ssm_contact_model``	- Statistical parameters for strata contact model
- ``ae_table`` 				- Filepath for architectural element lookup table
- ``ssm_contact``			- Contact surface type *[flat, random, user]*
- ``ll_ssm_ae``				- Which architectural elements are in each stratum
- ``ll_ae_prob``			- Probability of an architectural element occuring
- ``ll_ae_z_mean``			- Mean thickness of architectural element unit
- ``ll_avul_prob``			- Probability of avulsion
- ``ll_avul``				- Avulsion depth range
- ``r_bg``					- Background values for unassigned cells 

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
``[element]`` sections for architectural elements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sections that describe architectural elements are entitled with an identifying name (e.g. ``[sparse_scour]``). Note that section names should not include spaces. The first parameter to be set it the ``geometry``. The current implementation of HYVR includes three geometries: truncated ellipsoids (``trunc_ellip``), channels (``channel``), and sheets (``sheet``).

Unless otherwise noted, ranges (``r_``) represent the lower and upper limits of uniform distributions from which values are randomly generated.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
General ``[*element]`` parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- ``geometry``  		- Geometry of hydrofacies assemblages within architectural element *[trunc_ellip, ext_par, sheet]*
- ``structure``			- Internal structure of hydrofacies assemblages *[massive, dip]*, also , *[bulb, bulb_l, random]* for truncated ellipsoids.
- ``contact`` 			- Type of bedding contact between element units *[flat, random]*
- ``r_contact_model``	- Statistical parameters for bedding contact model
- ``l_facies`` 			- Hydrofacies included in hydrofacies assemblage. These refer to ``[hydraulics].l_hydro`` and are zero-indexed.
- ``ll_altfacies`` 		- Alternating facies specification. This is a list of lists and should have one entry for each value in ``[*element].l_facies``.
- ``r_bg`` 				- Background parameters for unassigned cells in the architectural element. This should be three values: facies, azimuth, and dip background values.
- ``r_geo_ztrend``		- Linear trend in geometry sizes with elevation. Given as a percentage change mulitplier in mean value from bottom to top of domain, i.e. :math:`[\lambda_{bottom}, \lambda_{top}]`
- ``r_k_ztrend``		- Linear trend in isotropic hydraulic conductivity from bottom to top of domain :math::math:`\xi_{bottom},\xi_{top}`
- ``r_k_xtrend``		- Linear trend in isotropic hydraulic conductivity from model inlet to outlet :math:`\xi_{inlet},\xi_{outlet}`
- ``r_n_ztrend``		- Linear trend in porosity from bottom to top of domain :math:`\xi_{bottom},\xi_{top}`
- ``r_n_xtrend``		- Linear trend in porosity from model inlet to outlet :math:`\xi_{inlet},\xi_{outlet}`
- ``r_dip`` 			- Range of dip

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Erosive element-specific parameters (truncated_ellipsoid, extruded parabola)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
- ``agg`` 		- Aggradation thickness added between each generation elevation. 
- ``buffer``	- Buffer to reduce erosion of underlying units (see :ref:`methods <temethod>`).
- ``r_migrate``	- Lateral migration of ellipsoid centrepoints drawn from a random normal distribution, given as mean and variance in :math:`x` and :math:`y` directions :math:`[\overline{\Delta x}, \sigma^2_{\Delta x}, \overline{\Delta y}, \sigma^2_{\Delta y}]`. 
- ``l_lag`` 	- Parameters for lag surface *[lag thickness, hydrofacies ID]*


.. _teparams:

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Truncated ellipsoid parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- ``el_z``		- Number of elements to be simulated per simulation elevation and layer area
- ``length``, ``width``, ``depth`` -  Mean geometry of truncated ellipsoids
- ``r_paleoflow`` 	- Range of the uniform distribution from which the paleoflow orientation will be randomly drawn. 
- ``r_dip``- Range of the uniform distribution from which the dip will be randomly drawn.
- ``r_azimuth`` - Range of the uniform distribution from which the azimuth will be randomly drawn.
- ``bulbset_d`` - Thickness of nested-bulb structures at the maximum depth of the truncated ellipsoid.
- ``dipset_d`` - Thickness of dipping internal structures.

.. _chparams:

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Extruded parabola parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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

- ``flag_gen`` 			- Generate hydraulic parameters (i.e. hydraulic conductivity)?
- ``l_hydro`` 			- List of hydrofacies codes
- ``r_k_h`` 			- Mean horizontal hydraulic conductivity 
- ``r_sig_y`` 			- Variance of log hydraulic conductivity
- ``ll_ycorlengths`` 	- Default correlation lengths for :math:`\log(K_{iso})` in each hydrofacies in :math:`x,y,z`-directions
- ``r_k_ratio`` 		- List of perpendicular anisotropy ratios (i.e :math:`\frac{K_h}{K_v}`)
- ``r_n`` 				- List of mean porosity values
- ``r_sig_n``			- Variance of porosity values
- ``ll_ncorlengths`` 	- Default correlation lengths for porosity in each hydrofacies in :math:`x,y,z`-directions


------------------------------------------------------------------------
``[flowtrans]`` section
------------------------------------------------------------------------
This section contains parameters to be used for groundwater flow and solute transport simulations. This allows all input parameters for field generation and subsequent modelling to be stored in the same ``.ini`` file. 


----------------------------------
File structure
----------------------------------

HYVR simulations are structured in the following way:

Model -> Run -> Realisation

The default operation of HYVR is to save simulation outputs in the same directory as the ``.ini`` parameter input file. However, it is possible to save the simulation outputs into another directory by specifying ``run.modeldir`` in the parameter file.


-----------------------------------
Model outputs
-----------------------------------
HyVR has a number model outputs that can be set in the input parameter file. A copy of the ``ini`` model parameter file is saved in the model directory automatically. The following data output files include model outputs as three-dimensional arrays:

- ``dat`` Python 'pickle' file - this is a native Python format that can be loaded into Python using ``hyvr.utils.load_pickle()``.
- ``mat`` MATLAB file
- ``vtr`` VTK rectilinear grid file -  this can be opened in ParaView for improved three-dimensional visualisation.

HyVR can also create files that can be used as model inputs for some flow and transport modelling packages These currently include:

- MODFLOW-2005 - ``bas``, ``dis``, ``lpf``, ``nam``, ``oc``, and ``pcg`` model input files. Provided suitable flow and transport parameters are set in the ``[flowtrans]`` section of the input parameter file, this simulation can be executed. 
- MODFLOW 6 - ``dis``, ``nam``, and ``npf`` model input files. A complete set of MODFLOW 6 input files cannot be generated in HyVR at this stage.
- HydroGeoSphere - *K* tensors and porosity at each grid node. 

Note that these model inputs can only have regular model grids. They have not been tested for use in the above-named packages. 