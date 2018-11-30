from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.sdist import sdist as _sdist
# To use a consistent encoding
from codecs import open
from os import path



# Numpy stuff
# ----------
# To build the extensions, we need numpy headers. This makes sure numpy is correctly installed.
# See https://stackoverflow.com/questions/27021270/how-to-handle-dependency-on-scipy-in-setup-py
class build_ext(_build_ext):

    def finalize_options(self):
        _build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process:
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())




# Cython stuff
# ------------
# Here I'm trying to import cython to recreate the *.c files if possible. If not, I'm using the
# existing *.c files for creating the extensions. See here for more info:
# https://stackoverflow.com/questions/4505747/how-should-i-structure-a-python-package-that-contains-cython-code
try:
    from Cython.Build import cythonize
    use_cython = True
    ext = ".pyx"
except ModuleNotFoundError:
    use_cython = False
    ext = ".c"

# Always run cython before creating a source distribution
class sdist(_sdist):
    def run(self):
        # Make sure the compiled Cython files in the distribution are up-to-date
        cythonize(extensions)
        _sdist.run(self)

# Extensions
# ----------
extensions = [Extension("hyvr.optimized", sources=["hyvr/optimized"+ext]),
              Extension("hyvr.classes.grid", sources=["hyvr/classes/grid"+ext]),
              Extension("hyvr.classes.contact_surface", sources=["hyvr/classes/contact_surface"+ext]),
              Extension("hyvr.classes.trough", sources=["hyvr/classes/trough"+ext]),
              Extension("hyvr.classes.sheet", sources=["hyvr/classes/sheet"+ext]),
              Extension("hyvr.classes.channel", sources=["hyvr/classes/channel"+ext]),
              Extension("hyvr.assign_facies", sources=["hyvr/assign_facies"+ext]),
              ]
                         # include_dirs=[np.get_include()])]
if use_cython:
    ext_modules = cythonize(extensions)
else:
    ext_modules = extensions



# Get the long description from the README file
here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

# get version number from file
with open('versionnumber', 'r') as f:
    version = f.read().strip('\n')


# copy the made.ini test case to the package directory
from shutil import copyfile
copyfile(path.join(here, 'testcases', 'made.ini'), path.join(here, 'hyvr', 'made.ini'))

setup(
    name='hyvr',
    version=version,
    description='A python package for simulating hydrogeological virtual realities',
    long_description=long_description,
    long_description_content_type='text/x-rst',
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
    cmdclass={'build_ext':build_ext,
              'sdist':sdist},
    ext_modules=ext_modules,
    setup_requires=['numpy'],
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib',
        'pandas',
        ],
    extras_require = {
        'HDF5': ['h5py'],
        'MODFLOW': ['flopy'],
        'VTR': ['pyevtk'],
        },

    # include testcase config file
    package_data={
        'hyvr': ['made.ini'],
    },
    zip_safe=False, # Due to the .pxd files
)
