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

   File "/tmp/pip-build-3mre9i5d/flopy/flopy/mbase.py", line 21, in <module>
   import numpy as np
   ModuleNotFoundError: No module named 'numpy'

In this case it helps to install numpy first separately.

::

   ImportError: libtk8.6.so: cannot open shared object file: No such file or directory


On Linux this can be fixed by installing the package ``tk`` using your
distribution's package manager. Depending on your distribution, this might be
called ``tk-devel``, ``tkinter``, ``tk-dev``, ``tk``, ``python-tk`` or maybe
another variation. A google search could be helpful.


-----------------------------
Runtime errors
-----------------------------

::

    TypeError: __init__() got an unexpected keyword argument 'tdisrecarray'

This happens because HyVR is currently incompatible with ``flopy 3.2.9``.
Use ``flopy 3.2.8`` instead. To avoid conflicts with dependencies of other
package consider using a virtual environment for python.
