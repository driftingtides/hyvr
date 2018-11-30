.. _troubleshooting:

==========================================================
Troubleshooting
==========================================================

This section contains problems we encountered during testing and how we solved
them. It might or might not help you.
If you have problems that are not listed here, don't hesitate to contact us
either via opening an issue in the github repository or via email.
If you had problems that are not listed here but you solved them yourself, let
us know anyway so we can either fix it or add the solution here.

----------------------------
Installation
----------------------------

::

   ModuleNotFoundError: No module named 'numpy'

   --------------------------------------------
   Command "python setup.py egg_info" failed with error code 1 in <some path>/flopy/

This happens because in order to install flopy you need to have numpy installed first. Just install numpy and try again. In newer versions flopy will be an optional dependency, so this should not happen anymore.


::

   ImportError: libtk8.6.so: cannot open shared object file: No such file or directory


On Linux this can be fixed by installing the package ``tk`` using your
distribution's package manager. Depending on your distribution, this might be
called ``tk-devel``, ``tkinter``, ``tk-dev``, ``tk``, ``python-tk`` or maybe
another variation. A google search could be helpful.


::

    Failed building wheel for hyvr
    ...
    <very long command containing the words 'compile' and 'install'>

This happens if there's only a source distribution on PyPI and you don't have
access to a compiler. Contact us if this happens and we'll try to upload a
pre-built binary.
