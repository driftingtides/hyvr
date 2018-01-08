======================================================
The Hydrogeological Virtual Reality simulation package
======================================================

Welcome to the Hydrogeological Virtual Reality simulation package (HYVR). This is an *alpha* release of the package and has not been thoroughly tested.

At this stage if you are interested in using HYVR please download the repository. The simulation package can be run from the command line by opening a command line window in the src/hyvr directory and then typing the following::

    python sim.py
    
This will run the simulation function with the default example parameters. It may be necessary to install the Python packages listed below if they aren't already on your system.If you would like to use an alternative input parameter file you will need to edit the parameter file path at the end of sim.py. This will be improved in future commits of the HYVR package.

-----------------
Install (Windows)
-----------------

With conda

^^^^^^^^^^

To install HYVR we recommend first installing the `Anaconda distribution <https://www.anaconda.com/download/>`_ of Python 3 This distribution has the majority of dependencies that HYVR requires.

It is a good idea to install the HYVR package into a `virtual environment <https://conda.io/docs/user-guide/tasks/manage-environments.html>`_. Do this by opening a command prompt window and typing the following::    

    python -m virtualenv hyvr_env
	    
You need to then activate this environment::

    call hyvr_env/scripts/activate
	
Install the necessary python packages by downloading the ``requirements.txt`` file in the HyVR repository and then running::
	
	pip install -r <path to>requirements.txt
    
Once this is completed install HYVR using pip::

    (hyvr_env) <working directory> pip install hyvr
   
    

Source
------
The most current version of HyVR will be available at the 'github repository <https://github.com/driftingtides/hyvr/>`_; a version will also be available on the 'PyPI index <https://pypi.python.org/pypi/hyvr/>'_ which can be installed using ``pip``.


Requirements
------------

Python
^^^^^^
Python 3.4 or greater

Modules
^^^^^^^

* scipy
* `pyevtk <https://pypi.python.org/pypi/PyEVTK>`_
* pandas
* numpy
* matplotlib
* `flopy <https://github.com/modflowpy/flopy>`_

