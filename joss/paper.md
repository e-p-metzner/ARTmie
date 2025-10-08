---
title: 'ARTmie:  Fast and Flexible Computation of Aerosol Optical Properties for Atmospheric Modeling and Remote Sensing'
tags:
  - Python
  - Mie
  - scattering
  - backscattering
  - extinction
  - absorption
  - size distribution
  - spheres
  - optics
authors:
  - name: Enrico P. Metzner
    orcid: 0000-0001-6523-9205
    corresponding: true
    affiliation: '&dagger;'
affiliations:
  - index: '&dagger;'
    name: Karlsruhe Institute of Technology (KIT), Germany
date: 10 April 2025
bibliography: paper.bib
---


# Summary

In atmospheric science, astronomy, and related fields, accurate computation of particle optical properties — such as mass extinction coefficient, single-scattering albedo, asymmetry factor, and phase function — is essential for modeling interactions with radiation.
For spherical particles, these properties can be calculated using Mie theory, a well-established analytical solution.
The foundational formulas and early computational algorithms for Mie scattering were first formalized and improved by @BH1983.

ARTmie implements these Mie formulas to compute extinction, scattering, absorption, and backscattering efficiencies, as well as the corresponding cross sections.
It ports key routines from [@Maetzler2002matlab] (originally in MATLAB) and Bessel function implementations from [@Amos1986] (Fortran 77) into optimized C++.

ARTmie was developed with a strong focus on computational speed, particularly for applications involving particle size distributions and angle-weighted backscattering efficiencies — important for lidar-based remote sensing.
Its performance gains stem not only from the C++ backend, but also from algorithmic optimizations that reduce redundant calculations, especially when computing ensemble-averaged properties (e.g., for log-normal distributions).

Compared to widely used Mie codes, ARTmie offers significantly faster execution, making it ideal for large-scale aerosol optical property calculations in atmospheric models, remote sensing applications, and the creation of high-resolution datasets for machine learning (ML).
It supports standard optical property outputs and is designed for seamless integration into modeling pipelines or data-driven workflows.


# Statement of need

Accurate calculation of aerosol optical properties is fundamental to research in physics, atmospheric science, and astronomy.
As the complexity of aerosol mixtures increases and the number of particles considered grows, the computational demand for Mie scattering calculations becomes significant.
Traditional Mie codes written in FORTRAN or Python often struggle with performance in such scenarios.

ARTmie addresses this challenge by offering a fast and efficient solution for computing Mie scattering properties.
It significantly outperforms pure Python implementations, making it suitable for large-scale applications, including high-resolution look-up table generation and ensemble simulations.

Beyond the standard Mie efficiencies for extinction, scattering, and absorption, ARTmie also computes the backscatter efficiency through a weighted average over the backward scattering angles — an important metric for interpreting lidar-based remote sensing measurements.

Despite its performance-oriented design, ARTmie remains user-friendly and requires minimal setup, enabling easy integration into research workflows.

ARTmie has already been successfully applied in the ICON-ART module of the ICON numerical weather and climate prediction model [@Rieger2015; @icon202504].
It was used to recalculate the optical properties of aerosols routinely used in ICON-ART and to extend the model with additional aerosol types and coating configurations [@art202504].


# State of the field

There are many Mie scattering libraries available in a variety of programming languages.
Python is probably the most accessible programming language for quickly and easily analyzing and testing data within the scientific community.
ARTmie is not the first Python module for Mie scattering; therefore, it will be compared to three other commonly used Python libraries:

1. miepython [@swMiepython] :

   miepython is a minimal, pure Python implementation of classical Mie theory that provides efficient scattering, absorption and extinction calculations for homogeneous spheres only.
   While it is ideal for simple, educational or small-scale use cases, it lacks support for coated particles or size distributions.
   Its performance is limited due to the absence of a compiled backend accelerator.

2. PyMieScatt [@articlePyMieScatt] :

   PyMieScatt is a feature-rich Python library that supports homogeneous and single-coated spheres. It also offers additional tools for data analysis and angular scattering, as well as inversion methods.
   However, like miepython, it is implemented entirely in Python, making it significantly slower than ARTmie for coated particles and polydisperse systems, where computational demands are higher.

3. PyMieSim [@articlePyMieSim; @swPyMieSim] :

   Designed for simulating optical systems, PyMieSim provides a modular interface for defining sources, detectors and particle configurations.
   It supports both homogeneous and single-coated spheres.
   Although it uses a C++ backend to optimize performance, the modular design introduces additional setup overhead.
   PyMieSim is best suited to integrating Mie scattering into broader optical modelling workflows but does not yet support polydisperse particles.


# Acknowledgement

This work was done during the funding by the Deutsche Forschungsgemeinschaft (DFG) under the project FE.5250.0010.7772.


# References
