************
Installation
************


Linux
=====

In almost any up-to-date Linux distribution the dependencies of **pypath** are
built-in, or provided by the distributors. You only need to install a couple
of things in your package manager (cairo, py(2)cairo, igraph,
python(2)-igraph, graphviz, pygraphviz), and after install **pypath** by *pip*
(see below). If any module still missing, you can install them the usual way
by *pip* or your package manager.

igraph C library, cairo and pycairo
-----------------------------------

*python(2)-igraph* is a Python interface to use the igraph C library. The
C library must be installed. The same goes for *cairo*, *py(2)cairo* and
*graphviz*.

Directly from git
-----------------

.. code:: bash

    pip install git+https://github.com/saezlab/pypath.git

With pip
--------

Download the package from /dist, and install with pip:

.. code:: bash

    pip install pypath-x.y.z.tar.gz

Build source distribution
-------------------------

Clone the git repo, and run setup.py:

.. code:: bash

    python setup.py sdist


Mac OS X
========

On OS X installation is not straightforward primarily because cairo needs to
be compiled from source. We provide 2 scripts here: the
**mac-install-brew.sh** installs everything with HomeBrew, and
**mac-install-conda.sh** installs from Anaconda distribution. With these
scripts installation of igraph, cairo and graphviz goes smoothly most of the
time, and options are available for omitting the 2 latter. To know more see
the description in the script header. There is a third script
**mac-install-source.sh** which compiles everything from source and presumes
only Python 2.7 and Xcode installed. We do not recommend this as it is time
consuming and troubleshooting requires expertise.

Troubleshooting
---------------

* ``no module named ...`` when you try to load a module in Python. Did
theinstallation of the module run without error? Try to run again the specific
part from the mac install shell script to see if any error comes up. Is the
path where the module has been installed in your ``$PYTHONPATH``? Try ``echo
$PYTHONPATH`` to see the current paths. Add your local install directories if
those are not there, e.g.
``export PYTHONPATH="/Users/me/local/python2.7/site-packages:$PYTHONPATH"``.
If it works afterwards, don't forget to append these export path statements to
your ``~/.bash_profile``, so these will be set every time you launch a new
shell.

* ``pkgconfig`` not found. Check if the ``$PKG_CONFIG_PATH`` variable is
set correctly, and pointing on a directory where pkgconfig really can be
found.

* Error while trying to install py(2)cairo by pip. py(2)cairo could not be
installed by pip, but only by waf. Please set the ``$PKG_CONFIG_PATH`` before.
See **mac-install-source.sh** on how to install with waf.

* Error at pygraphviz build: ``graphviz/cgraph.h file not found``. This is
because the directory of graphviz detected wrong by pkgconfig. See
**mac-install-source.sh** how to set include dirs and library dirs by
``--global-option`` parameters.

* Can not install bioservices, because installation of jurko-suds fails. Ok,
this fails because pip is not able to install the recent version of
setuptools, because a very old version present in the system path. The
development version of jurko-suds does not require setuptools, so you can
install it directly from git as it is done in **mac-install-source.sh**.

* In **Anaconda**, *pypath* can be imported, but the modules and classes are
missing. Apparently Anaconda has some built-in stuff called *pypath*. This
has nothing to do with this module. Please be aware that Anaconda installs a
completely separated Python distribution, and does not detect modules in the
main Python installation. You need to install all modules within Anaconda's
directory. **mac-install-conda.sh** does exactly this. If you still
experience issues, please contact us.


Microsoft Windows
=================

Not many people have used *pypath* on Microsoft computers so far. Please share
your experiences and contact us if you encounter any issue. We appreciate
your feedback, and it would be nice to have better support for other computer
systems.

With Anaconda
-------------

The same workflow like you see in ``mac-install-conda.sh`` should work for
Anaconda on Windows. The only problem you certainly will encounter is that not
all the channels have packages for all platforms. If certain channel provides
no package for Windows, or for your Python version, you just need to find an
other one. For this, do a search:

.. code:: bash

    anaconda search -t conda <package name>

For example, if you search for *pycairo*, you will find out that *vgauther*
provides it for osx-64, but only for Python 3.4, while *richlewis* provides
also for Python 3.5. And for win-64 platform, there is the channel of
*KristanAmstrong*. Go along all the commands in ``mac-install-conda.sh``, and
modify the channel if necessary, until all packages install successfully.

With other Python distributions
-------------------------------

Here the basic principles are the same as everywhere: first try to install all
external dependencies, after *pip* install should work. On Windows certain
packages can not be installed by compiled from source by *pip*, instead the
easiest to install them precompiled. These are in our case *fisher, lxml,
numpy (mkl version), pycairo, igraph, pygraphviz, scipy and statsmodels*. The
precompiled packages are available here:
http://www.lfd.uci.edu/~gohlke/pythonlibs/. We tested the setup with Python
3.4.3 and Python 2.7.11. The former should just work fine, while with the
latter we have issues to be resolved.

Known issues
------------

* *"No module fabric available."* -- or *pysftp* missing: this is not
important, only certain data download methods rely on these modules, but
likely you won't call those at all.
* Progress indicator floods terminal: sorry about that, will be fixed soon.
* Encoding related exceptions in Python2: these might occur at some points in
the module, please send the traceback if you encounter one, and we will fix
as soon as possible.

*Special thanks to Jorge Ferreira for testing pypath on Windows!*
