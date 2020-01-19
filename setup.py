# setup.py for HyVR
# -----------------
#
# The setup.py is a bit more complicated than in pure python packages, as we
# have to make sure the extensions are built.
#
# Normally, running this will use the already pre-generated *.c-files from
# running Cython during creating a source distribution.
# If you want to cythonize the extensions yourself, pass the "--cythonize"
# option to ``setup.py install`` or ``--install-option="--cythonize"`` to ``pip
# install``.
# Below, you will find also some options for cythonizing (profiling, annotating,
# etc.) in the function ``custom_cythonize``.

from os import path
from shutil import copyfile

from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.sdist import sdist as _sdist
from setuptools.command.install import install
from setuptools.command.develop import develop


# Numpy stuff
# ----------
# To build the extensions, we need numpy headers. This makes sure numpy is
# correctly installed. See
# https://stackoverflow.com/questions/27021270/how-to-handle-dependency-on-scipy-in-setup-py
class build_ext(_build_ext):

    def finalize_options(self):
        _build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process:
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())


# Cython stuff
# ------------
# custom settings for cythonize
def custom_cythonize():
    from Cython.Build import cythonize
    cythonize(get_extensions(".pyx"),
              language_level=3,
              annotate=False,
              compiler_directives={'profile':False})
    

# cythonize when ``--cythonize`` command line option is passed
class CythonizeMixin(object):

    user_options = [
        ("cythonize", None, "recreate the c extionsions with cython")
    ]

    def initialize_options(self):
        super().initialize_options()
        self.cythonize = False

    def run(self):
        if self.cythonize:
            custom_cythonize()
        super().run()

class develop_with_cythonize(CythonizeMixin, develop):
    user_options = getattr(develop, 'user_options', [])\
                   + CythonizeMixin.user_options
    
class install_with_cythonize(CythonizeMixin, install):
    user_options = getattr(install, 'user_options', [])\
                   + CythonizeMixin.user_options

# Always run cython before creating a source distribution
class sdist(_sdist):
    def run(self):
        # Make sure the compiled Cython files in the distribution are up-to-date
        custom_cythonize()
        # copy the made.ini test case to the package directory as default test
        copyfile(path.join(here, 'tests', 'full_testcases', 'made.ini'), path.join(here, 'hyvr', 'made.ini'))
        _sdist.run(self)

# Extensions
# ----------
here = path.abspath(path.dirname(__file__))
hyvr_path = path.join(here, "hyvr")



def get_extensions(ext):

    extensions = {
        "hyvr.optimized": path.join(hyvr_path, "optimized"),
        "hyvr.geo.grid": path.join(hyvr_path, "geo", "grid"),
        "hyvr.geo.contact_surface": path.join(hyvr_path, "geo", "contact_surface"),
        "hyvr.geo.ae_realization": path.join(hyvr_path, "geo", "ae_realization"),
        "hyvr.geo.trough": path.join(hyvr_path, "geo", "trough"),
        "hyvr.geo.trough_ae": path.join(hyvr_path, "geo", "trough_ae"),
        "hyvr.geo.sheet": path.join(hyvr_path, "geo", "sheet"),
        "hyvr.geo.sheet_ae": path.join(hyvr_path, "geo", "sheet_ae"),
        "hyvr.geo.channel": path.join(hyvr_path, "geo", "channel"),
        "hyvr.geo.channel_ae": path.join(hyvr_path, "geo", "channel_ae"),
        "hyvr.assign_points": path.join(hyvr_path, "assign_points"),
    }

    return [Extension(name, sources=[extensions[name]+ext]) for name in extensions]


# Get the long description from the README file
with open(path.join(here, 'README.rst'), 'r', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='hyvr',
    version="1.1.0",
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
    packages=find_packages(exclude=['tests', 'tests.*']),
    python_requires='>=3.4',
    cmdclass={'build_ext':build_ext,
              'sdist':sdist,
              'install':install_with_cythonize,
              'develop':develop_with_cythonize,
    },
    ext_modules=get_extensions(".c"),
    setup_requires=['numpy'],
    install_requires=[
        'numpy',
        'scipy',
        ],
    extras_require = {
        'HDF5': ['h5py'],
        'MODFLOW': ['flopy'],
        'VTR': ['pyevtk'],
        'develop': ['Cython',
                    'pytest',
                    'twine',
                    'Sphinx',
                    'sphinxcontrib-bibtex',
                    'sphinxcontrib-fulltoc',
                    'sphinxcontrib-fulltoc',
                    'h5py',
                    'flopy',
                    'pyevtk',
                    'restructuredtext-lint',
        ],
        },

    # include testcase config file
    package_data={
        "" : ["LICENCE", "*.c", "*.pyx", "*.pxd"],
        "hyvr": ["made.ini"],
    },
    exclude_package_data={"": ["tests"]},
    zip_safe=False, # Due to the .pxd files
)
