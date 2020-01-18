====================================================================
Introduction
====================================================================

**HyVR: Turning your geofantasy into reality!** 

The Hydrogeological Virtual Reality simulation package (HyVR) is a Python module
that helps researchers and practitioners generate subsurface models with
multiple scales of heterogeneity that are based on geological concepts. The
simulation outputs can then be used to explore groundwater flow and solute
transport behaviour. This is facilitated by HyVR outputs in common flow
simulation packages' input formats. As each site is unique, HyVR has been
designed that users can take the code and extend it to suit their particular
simulation needs.

The original motivation for HyVR was the lack of tools for modelling sedimentary
deposits that include bedding structure model outputs (i.e., dip and azimuth).
Such bedding parameters were required to approximate full hydraulic-conductivity
tensors for groundwater flow modelling. HyVR is able to simulate these bedding
parameters and generate spatially distributed parameter fields, including full
hydraulic-conductivity tensors. More information about HyVR is available in the
online `technical documentation <https://driftingtides.github.io/hyvr/index.html>`_.

I hope you enjoy using HyVR much more than I enjoyed putting it together! I look
forward to seeing what kind of funky fields you created in the course of your
work.

*HyVR can be attributed by citing the following journal article: Bennett, J. P.,
Haslauer, C. P., Ross, M., & Cirpka, O. A. (2018). An open, object-based
framework for generating anisotropy in sedimentary subsurface
models. Groundwater.
DOI:* `10.1111/gwat.12803 <https://onlinelibrary.wiley.com/doi/abs/10.1111/gwat.12803>`_.
*A preprint version of the article is available* `here <https://github.com/driftingtides/hyvr/blob/master/docs/Bennett_GW_2018.pdf>`_.

Installing the HyVR package
--------------------------------------

Installing Python
^^^^^^^^^^^^^^^^^


Windows
"""""""

If you are using Windows, we recommend installing the `Anaconda distribution
<https://www.anaconda.com/download/>`_ of Python 3. This distribution has the
majority of dependencies that HyVR requires.

It is also a good idea to install the HyVR package into a `virtual environment
<https://conda.io/docs/user-guide/tasks/manage-environments.html>`_. Do this by
opening a command prompt window and typing the following::

    conda create --name hyvr_env

You need to then activate this environment::

    conda activate hyvr_env
	

Linux
"""""

Depending on your preferences you can either use the Anaconda/Miniconda
distribution of python, or the version of your package manager. If you choose
the former, follow the same steps as for Windows.

If you choose the latter, you probably already have Python 3 installed. If not,
you can install it using your package manager (e.g. ``apt`` on Ubuntu/Debian).

In any way we recommend using a virtual environment. Non-conda users can use for
example `virtualenvwrapper <https://virtualenvwrapper.readthedocs.io/en/latest/>`_.


Installing HyVR
^^^^^^^^^^^^^^^

Once you have activated your virtual environment, you can install HyVR from PyPI
using ``pip`` (with Anaconda you might have to install pip first into your
environment using ``conda install pip``)::

    pip install hyvr

It might be necessary to install numpy separately before installing HyVR.
To check whether installation was successful, you can run::

    python -m hyvr

If this runs without error (warnings are ok), ``hyvr`` should be successfully
installed.

If installing from PyPI doesn't work for you, please let us know (see
below). You might want to try installing from source instead.

Installation from conda-forge will (hopefully) be coming soon.

Installing from source
""""""""""""""""""""""

If installation via ``pip`` fails, you can try to install from source by cloning
or downloading the github repository.
To install from source you need a C/C++ compiler. On Windows you can get one by
installing "Build Tools for Visual Studio".
If you have ``git`` installed on your machine, you can do::

    git clone https://github.com/driftingtides/hyvr.git
    pip install ./hyvr

Otherwise you can download the code as a zip file, unzip it, and then install it via::

    pip install <path/to/unzipped/hyvr>

If you are in the same directory as the unzipped ``hyvr`` directory make sure to
use ``./hyvr`` instead of ``hyvr``, otherwise pip tries to install from PyPI.


Dependencies
^^^^^^^^^^^^

``pip`` should normally install all required dependencies. Optional dependencies are:

- ``Cython`` to recreate the C-extensions
- ``h5py`` for HDF5 output
- ``flopy`` for output of MODFLOW files
- ``pyevtk`` for VTR output


Usage
-----

To use HyVR you have to create a configuration file with your settings.
You can then run HyVR the following way::

    (hyvr_env) $ python -m hyvr my_configfile.ini

HyVR will then run and store all results in a subdirectory. If no configfile is
given, it will run a test case instead::

    (hyvr_env) $ python -m hyvr

If you want to use HyVR in a script, you can import it and use the ``run`` function::

    import hyvr
    hyvr.run('my_configfile.ini')
    
Examples can be found in the ``tests/testcases`` directory of the `github
repository <https://github.com/driftingtides/hyvr/>`_, the general setup and
possible options of the config-file are described in the
documentation. Currently only ``tests/testcaes/made.ini`` is ported to version 1.0.0.


Development
-----------
HyVR has been developed by Jeremy Bennett (`website <https://jeremypaulbennett.weebly.com>`_)
as part of his doctoral research at the University of Tübingen and by Samuel
Scherrer as a student assistant.

You can contact the developer(s) of HyVR by `email <mailto:hyvr.sim@gmail.com>`_
or via github.

Problems, Bugs, Unclear Documentation
-------------------------------------

If you have problems with HyVR have a look at the `troubleshooting
<https://driftingtides.github.io/hyvr/troubleshooting.html>`_ section. If this
doesn't help, don't hesitate to contact us via email or at github.

If you find that the documentation is unclear, lacking, or wrong, please also
contact us.
