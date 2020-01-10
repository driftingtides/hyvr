<!-- markdown-toc start - Don't edit this section. Run M-x markdown-toc-refresh-toc -->
# Table of Contents

- [General Coding Guidelines](#general-coding-guidelines)
- [Version Control](#version-control)
    - [Working on an issue](#working-on-an-issue)
    - [Releasing a new version](#releasing-a-new-version)
        - [PyPI-test](#pypi-test)
        - [PyPI](#pypi)
        - [Github](#github)
        - [Documentation](#documentation)
- [How HyVR Works](#how-hyvr-works)
    - [Reading Parameters](#reading-parameters)
    - [Generating model](#generating-model)
    - [Assigning Values to Grid](#assigning-values-to-grid)
    - [Output](#output)

<!-- markdown-toc end -->


# General Coding Guidelines

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

# Version Control

The two major branches on [github](https://github.com/driftingtides/hyvr) are
`main` and `develop`.
The `main` branch should always contain the current released version, that
should also be available on PyPI.
The `develop` contains the current state of development.
For working on an issue, create a new branch based on the `develop`
branch. Only for fixing critical bugs you can directly base your work on
`master` and merge your changes to master afterwards.


## Working on an issue

If you don't have cloned the repository yet, you can do it with
```
git clone git@github.com:driftingtides/hyvr.git
```
Otherwise, do a `git pull` to get the current version. If you are not already on
the develop branch, change to it with `git checkout develop`.
Now you can create a new branch with
```
git checkout -b <type>/<short issue title>
```
where `<type>` should be the type of the issue and `<short issue title>` a short
description about what the issue is. This also works if you forgot to create a
new branch but didn't add/commit any changes yet.

`<type>` should be one of the following:
- `wip` (work in progress) for stuff that will take a while and require a lot of
  coding.
- `feat` for adding minor features or minor fixes (spelling etc.)
- `bug` for bugfixes
- `doc` for documentation changes
`<short issue title>` should normally be one word describing what the issue is about.

Commit messages should always have one header line and then be followed by a
block description. The header line should start with the issue type and
description in brackets `[<type>][<short issue title>]`. You can reference a
github issue by `#<issue number>`, and automatically close an issue with
`Resolves #<issue number>`. In case you do this, put this statement on the last
line of the commit message.

After fixing the issue, you can merge the branch into `develop` or
create a pull request. Therefore you have to change into the `develop` branch
again:
```
git checkout develop
git merge <type>/<short issue title>
```
Then you can push your changes to github with
```
git push -u origin develop
```
In case you are working on a bigger issue, you might want to push your feature
branch regularly to github, so you have a backup in case of disaster:
```
git push -u origin <branch name>
```


## Releasing a new version

To release a new version, you first have to merge the current changes into
`master` locally (either from `develop` for a real new version or from a bugfix
branch for critical bugs).

Then you should run the tests to make sure everything works as
expected. Documentation for the tests can be found [here](tests/README.md).

In case everything works as expected you should perform the following steps:
- set the correct version number for PyPI-test and create distributions
- upload the new version to PyPI-test and try installing it from there
- set the correct version number and rebuild distributions
- upload the new version to PyPI
- push to the github master branch
- update the online documentation

### PyPI-test

The versionnumber should be in the file `versionnumber`. It follows the format
`<major version number>.<minor version number>.<number of fixes>`.
As PyPI and PyPI-test don't allow reuploading the same version again, we append
a further number to the version number for PyPI-test, e.g. the <n>-th try on
PyPI-test should get the version number `<major>.<minor>.<fix>.<n>`.

After you have set the version number, create distributions of the python
package by running `setup.py`:
```
python setup.py sdist bdist_wheel
```
(Note: On linux you cannot create wheels.)

Then you can upload the distributions first to PyPI-test using `twine`:
```
python -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*
```
To test if it worked, create a new virtual python environment, e.g. named test,
and then install the package from PyPI test:
```
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple hyvr
```

If everything worked, move on, otherwise try to fix the issues first.

### PyPI

Now you can set the correct version number. Then remove the distributions you
created for PyPI-test and rebuild them:
```
rm dist/*
python setup.py sdist bdist_wheel
```
Then you can upload the new version to PyPI:
```
python -m twine upload dist/*
```

### Github

Running `setup.py` might have recreated the C-files. You can add them and commit
all other changes if you didn't already do so. The commit message header should
start with `[v<version number>]`.

Then tag the current commit with the version number:
```
git tag v<version number>
```
e.g. `git tag v1.1.0` for version 1.1.0.

Now you can push the changes and the version tag:
```
git push origin
git push origin <tag>
```

### Documentation

To rebuild the documentation, change to the `docs` directory and run the
following commands:
```
sphinx-apidoc -fM -H 'Module Reference' -o modules ../hyvr
make html
```

To push the new documentation to github pages, you have to check out the branch
`gh-pages` and recreate the documentation:
```
git checkout gh-pages
rm _sources _static
git checkout master --
git reset HEAD
cd docs
sphinx-apidoc -fM -H 'Module Reference' -o modules ../hyvr
make html
cd ..
cp docs/_build/html/* .
rm -rf hyvr docs README.rst testcases versionnumber
git add -A
git commit -m <commit message>
git push origin gh-pages
```


# How HyVR Works

There are roughly 4 big tasks that HyVR performs:

1) Reading the parameter file and validating input
2) Generating the model
3) Assigning values to grid points
4) Storing the results in different formats


## Reading Parameters

The first part of HyVR is reading (and validating) the input from the `*.ini`
configuration file.  The code for this part can be found in
``hyvr/input/``. There are 3 important files in there:

* ``option_parsing.py``: This is our option parsing module, you probably
  don't have to change this.
* ``options.py``: This is the list of possible options. If you add
  functionality, you will need to add your option here so that it will be
  parsed.
* ``parameters.py``: This does the actual work of opening the file, parsing the
  different sections, and some additional validation of the parameters
  (everything that couldn't be done by the option parsing module).


## Generating model

The model is generated from top to bottom, by first generating contact surfaces
of the uppermost strata, then generating contact surfaces of AEs within this
strata, and then generating the objects that are placed inside the AEs. This is
done for every stratum.
The model is stored in separate objects: The model objects holds a list of
strata, each strata holds a list of the AEs inside, and each AE holds a list of
the geometrical objects/hydrofacies assemblages it contains.

Implementation of architectural elements is split in several ways. First,
there's a distinction between AE types (``hyvr/geo/ae_types.py``), and AE
realizations (``hyvr/geo/ae_realizations.py``).
The AE types correspond to the sections in the ini-file, e.g. ``clay_sheet`` in
``made.ini``. From these templates, AE realizations are randomly created
according to the parameters.
AE realization implementation is also split. Firstly, there's the Python part,
where each AE realization has a list of objects, for access/creation from Python
code. However, to make the assignment of grid points to objects faster, each AE
realization also contains arrays of object properties, since we can iterate much
faster over arrays compared to Python lists.

Most of the code for model generation can be found in ``hyvr/geo/``.


## Assigning Values to Grid

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


## Output

The last step is to output the HyVR results to different formats. This is
implemented in `postprocess/output.py`. The central method is `create_outputs`,
which in turn calls the several `to_*` methods, that return the specified
format. New formats can be added by providing a `to_*` method, adding it in the
list in `create_outputs`, and as possible option in `options.py`.
