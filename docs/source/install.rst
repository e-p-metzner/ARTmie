**Date**: |today|, **Version**: |version|, **Author**: E. P. Metzner

Installation and tests
======================




Getting *ARTmie*
----------------

from the Python Package Index (PyPI)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

	pip install ARTmie


from source
~~~~~~~~~~~

1. Download the source code from this git by clicking the triangular button next to *Code* on the starting page and selecting e.g. `zip <https://gitlab.kit.edu/enrico.metzner/artmie/-/archive/main/artmie-main.zip>`_ underneath *Download source code*.
2. Extract to source code to a desired location.
3. Open a terminal with an appropriate python environment (Python 3.0 or higher) and Numpy (1.5 or higher, tested with 2.0)
4. Change to the directory of the ARTmie project (dir above src):
   
   .. code-block:: shell

       cd /path/to/ARTmie
   
5. Install it with the following command:
   
   .. code-block:: python

       python -m pip install .




Testing
-------

Make sure, you have the latest version of Pytest:

.. code-block:: python

    pip --update pytest

Then you can run the testsuite with the following commands:

.. code-block:: python

    import pytest
    test ARTmie


