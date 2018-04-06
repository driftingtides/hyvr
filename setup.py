from setuptools import setup, find_packages, Extension
# To use a consistent encoding
from codecs import open
from os import path

# this here is a bit ugly, but we need to make sure that cython and numpy are installed
try:
    import numpy as np
except ModuleNotFoundError:
    print("Building hyvr requires numpy. Installing numpy.")
    import pip
    pip.main(['install', "numpy"])
    import numpy as np

try:
    from Cython.Build import cythonize
except ModuleNotFoundError:
    print("Building hyvr requires cython. Installing numpy.")
    import pip
    pip.main(['install', "cython"])
    from Cython.Build import cythonize

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='hyvr',
    version='0.2',
    description='A python package for simulating hydrogeological virtual realities',
    long_description=long_description,
    url='https://github.com/driftingtides/hyvr',
    author='Jeremy Bennett',
    author_email='hyvr.sim@gmail.com',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Education',
        'Topic :: Scientific/Engineering',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
    ],
    keywords=['hydrogeology', 'sediment', 'simulator'],
    packages=find_packages(),
    python_requires='>=3.4',
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'pandas',
        'pyevtk',
        'flopy==3.2.8',
        'cython',
        ],

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
    include_package_data = True,
    package_data={
        'tests': ['made.ini', 'test_lu.ini', 'test_aelu.txt'],
    },

    # Cython stuff
    ext_modules = cythonize((Extension("hyvr.optimized", sources=["hyvr/optimized.pyx"], include_dirs=[np.get_include()])))

)
