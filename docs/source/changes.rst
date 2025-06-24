
Change-Log
==========

v0.1.2
******

New:

- This library can now be downloaded from `PyPI <https://pypi.org/project/ARTmie/>`_
- Documentation (hosted on `readthedocs/artmie <https://artmie.readthedocs.io/en/latest/>`_ )
- Testing with pytest (:math:`\sim`\ 90% of core functionality)

Changes:

- z attribute in bessel function can now be 1dimensional
- :func:`~ARTmie.MieQ`, :func:`~ARTmie.MieCoatedQ` accepts 1dimensional arrays for refractive indices and wavelength

Fixes:

- fix indexing of helper variables :math:`D_u, D_v, D_w` in :func:`~ARTmie.MieCoated_ab` to be comparable with Python library PyMieScatt




v0.1.1
******

First release