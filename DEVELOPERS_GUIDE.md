General Coding Guidelines
=========================

The following are some general guidelines that you should follow as much as
possible if you add code to HyVR. Unfortunately the existing code does not meet
all these criterions, but you are very welcome to change this, as long as it
doesn't affect how HyVR works (e.g. it doesn't change ini-file option names).

* 4 spaces for indentation
* Use underscores for word separations of variable names
* Variable names: Try to use variable names from which it is easy to infer what
  the variable is. Use abbreviations only if you think they are easy to
  understand for people that have no idea how HyVR works. Unfortunately the old
  code is still full of cryptic abbreviations like e.g. `ha_arr`.  Instead, use
  something like `hydrofacies_assemblages`, at least for member variables of the
  model object and keys for the model object's data dictionary.
* **Write docstrings!** If you add a function, make sure it has a docstring in
  the **numpy format**, i.e. a short summary of what the function does, followed
  by a description of input and output variables and their types.  Only if your
  function is very short and it's very obvious what it does (e.g. from the
  function and variable names), you can omit the docstring.
  The existing code still has some docstrings in google format, which should be
  converted to numpy format.
* Comments: You don't have to comment every line you write, but make sure to
  write a short note why you did things the way you did them if they are not too
  clear. You are also encouraged to use short descriptive comments/headers
  before blocks of code to add some structure to your code.
* Functions: Try to encapsulate the tasks that you are trying to solve into
  functions. A function should typically not be too long, it should fit on one
  or at most two pages in your text editor. If your function is very long,
  consider splitting it into several shorter functions.


Uploading a new version
=======================

Please run the tests before you push anything to github. Furthermore, you should
also publish the new version on PyPI. Therefore, you have to create a source
distribution or a wheel using
```
python setup.py sdist
```
or
```
python setup.py bdist_wheel
```
Make sure to set the versionnumber correctly. If you want to publish version
1.0.1, you should use the versionnumber 1.0.1.x for PyPI-test, where x is the
number of times you already tested whether installing from PyPI-test works.

Then you can upload the distributions first to PyPI-test using `twine`:
```
python -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*
```
To test if it worked, create a new virtual python environment, e.g. named test,
and then install the package from PyPI test:
```
pip install --index-url https://test.pypi.org/simple/
```

If everything worked, you can upload it using the correct versionnumber (i.e. in
this example 1.0.1) to PyPI:
```
python -m twine upload dist/*
```


How HyVR Works
==============

There are roughly 4 big tasks that HyVR performs:

1) Reading the parameter file and validating input
2) Generating the model
3) Assigning values to grid points
4) Storing the results in different formats

Below you will find short descriptions of all of these steps, but first some
general guidelines for working on HyVR.


Reading Parameters
------------------

The first part of HyVR is reading (and validating) the input from the `*.ini`
configuration file.  The code for this part can be found in
``hyvr/input/``. There are 3 important files in there:

* ``option_parsing.py``: This is our option parsing module, you probably
  shouldn't touch this.
* ``options.py``: This is the list of possible options. If you add
  functionality, you will need to add your option here so that it will be
  parsed.
* ``parameters.py``: This does the actual work of opening the file, parsing the
  different sections, and some additional validation of the parameters
  (everything that couldn't be done by the option parsing module).


Generating model
----------------

The model is generated from top to bottom, by first generating contact surfaces
of the uppermost strata, then generating contact surfaces of AEs within this
strata, and then generating the objects that are placed inside the AEs. This is
done for every stratum.
The model is stored in separate objects: The model objects holds a list of
strata, each strata holds a list of the AEs inside, and each AE holds a list of
the geometrical objects/hydrofacies assemblages it contains.


Assigning Values to Grid
------------------------

Assigning values to the grid is mostly done in the `assign_facies.pyx` file.
This basically just iterates over all grid points and assigns facies, dip,
azimuth, strata number, AE number, etc. to the grid point.
The code currently only works with regular grids. The grid points are iterated
over "towerwise", that is, the innermost loop is the z-loop. This is important,
as it reduces the number of strata, AEs and objects we have to check, since we
assume that a grid point can only be in the same object/AE/stratum as the grid
point before, or in a lower object/AE/stratum. Therefore we can always start our
search at the last found object/AE/stratum.

If a new geometry is added, it must also supply a
``maybe_assign_facies_azim_dip`` function as the other geometries. This function
should be implemented in Cython if possible, to get better speed. In the future,
the ``ae_realizations.py`` module should also be converted to Cython.


Output
------

The last step is to output the HyVR results to different formats. This is
implemented in `postprocess/output.py`. The central method is `create_outputs`,
which in turn calls the several `to_*` methods, that return the specified
format. New formats can be added by providing a `to_*` method, adding it in the
list in `create_outputs`, and as possible option in `options.py`.
