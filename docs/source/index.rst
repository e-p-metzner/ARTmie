.. ARTmie documentation master file, created by
   sphinx-quickstart on Mon Apr 14 12:55:01 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


ARTmie documentation
====================


*ARTmie* is a python library to calculate optical properties of spherical particles.
Those can be simple, single coated or drawn from a particle size distribution.
*ARTmie* uses a lot of improved C++ in the backend for accelaration of various calculations.

The Code was mainly ported from Matlab and Fortran to C++.


.. Add your content using ``reStructuredText`` syntax. See the
   `reStructuredText <https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html>`_
   documentation for details.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Simple particles <simple.rst>
   Coated particles <coated.rst>
   Particle distributions <psd.rst>
   Mathematical accuracy <bessel.rst>
   Installation and tests <install.rst>
   API Reference <api.rst>
