from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='hyvr',
    version='0.1.dev9',
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
        'flopy',
        ],

    # If there are data files included in your packages that need to be
    # installed, specify them here.  If using Python 2.6 or less, then these
    # have to be included in MANIFEST.in as well.
	include_package_data = True,
    package_data={
        '': ['*.ini'],
    },

)