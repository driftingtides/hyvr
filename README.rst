======================================================
The Hydrogeological Virtual Reality simulation package
======================================================

**HyVR: Turning your geofantasy into reality!** 

Welcome to the Hydrogeological Virtual Reality simulation package (HYVR). This is an *alpha* release of the package and has not been thoroughly tested. Check out the `Technical documentation <https://driftingtides.github.io/hyvr/index.html>`_

HyVR is a Python module that helps researchers and practitioners generate subsurface models with multiple scales of heterogeneity that are based on geological concepts. The simulation outputs can then be used to explore groundwater flow and solute transport behaviour. This is facilitated by HyVR outputs in common flow-and-transport simulation input formats. As each site is unique, HyVR has been designed that users can take the code and extend it to suit their particular simulation needs.

The original motivation for HyVR was the lack of tools for modelling sedimentary deposits that include bedding structure model outputs (i.e. dip and azimuth). Such bedding parameters were required to approximate full hydraulic-conductivity tensors for groundwater flow and solute transport modelling. HyVR is able to simulate these bedding parameters and generate spatially distributed parameter fields, including full hydraulic-conductivity tensors.

I hope you enjoy using HyVR much more than I enjoyed putting it together! I look forward to seeing what kind of funky fields you created in the course of your work. 

Installing the HYVR package (Windows)
--------------------------------------

With conda
^^^^^^^^^^

To install HyVR we recommend first installing the `Anaconda distribution <https://www.anaconda.com/download/>`_ of Python 3. This distribution has the majority of dependencies that HyVR requires.

It is a good idea to install the HyVR package into a `virtual environment <https://conda.io/docs/user-guide/tasks/manage-environments.html>`_. Do this by opening a command prompt window and typing the following::    

    python -m virtualenv hyvr_env
	    
You need to then activate this environment::

    call hyvr_env/scripts/activate
	
Install the necessary python packages by downloading the ``requirements.txt`` file in the HyVR repository and then running::
	
	pip install -r https://raw.githubusercontent.com/driftingtides/hyvr/master/requirements.txt
    
Once this is completed install HyVR using pip::

    (hyvr_env) <working directory> pip install hyvr
	
You can test whether HyVR is working by running the default parameters in the Python console::
	
	>>> import hyvr
	>>> hyvr.sim.main(0)
	
A number of messages should then be displayed in the console, including the location where the simulation outputs have been saved. 
   
Source
------
The most current version of HyVR will be available at the `github repository <https://github.com/driftingtides/hyvr/>`_; a version will also be available on the `PyPI index <https://pypi.python.org/pypi/hyvr/>`_ which can be installed using ``pip``.


Requirements
------------

Python
^^^^^^
HyVR was developed for use with Python 3.4 or greater. It may be possible to use with earlier versions of Python 3, however this has not been tested.

Modules
^^^^^^^

* scipy
* pandas
* numpy
* matplotlib
* `flopy <https://github.com/modflowpy/flopy>`_
* `pyevtk <https://pypi.python.org/pypi/PyEVTK>`_


Development
-----------
You can contact the developer(s) of HyVR by `email <mailto:hyvr.sim@gmail.com>`_.  HyVR is currently being developed by Jeremy Bennett (`website <https://jeremypaulbennett.weebly.com>`_) as part of his doctoral research at the University of Tübingen.
