---
title: 'ARTmie: a fast python package for mie scattering and backscattering calculations of single particles, coated particles and log-normal distributed particles'
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

In atmospheric science, astronomy and similar areas, the optical properties of particles like dust, ash and cloud droplets are necessary when calculating interactions with radiation.
The simplest assumption for these particles is that they are spherical, which leads to the mathematical fully derived theory of Mie.
The corresponding formulae and first computational algorithms were improved and published by @BH1983.

ARTmie uses these formulae to calculate extinction, scattering, absorption and backscattering efficiencies as well as cross sections by implementing and porting of the relevant functions [@Maetzler2002matlab] from Matlab and the bessel functions [@Amos1986] from Fortran 77 to C++.

ARTmie was developed with focus on speed, especially for particle size distribution and scattering angle weighted backscattering efficiency. ARTmie mainly gains its speed not only by using a C++ backend but also by optimizing the amount of necessary calculation when computing combined informations, e.g. for log-normal particle size distributions.


# Statement of need

A lot of research in physics and atmospheric and astronomical science depends on the optical properties of aerosols. The more variability in the constituents and the number of particles there is, the more calculations are needed. ARTmie is capable to do these calculations at a speed that pure Python code can never achieve.

In addition to the standard Mie efficiencies of extinction, scattering and absorption, ARTmie has the ability to calculate the backscatter efficiency by weighted averaging over the total backward scattering angles, as used in lidar-based measurement systems. Nevertheless, ARTmie is designed to be user-friendly with minimal preparation code required.

Furthermore, ARTmie has been used recently to improve the ICON-ART module of the numerical weather and climate prediction model ICON [@Rieger2015; @icon202504] by recalculating all optical properties of aerosols regularly used by ICON-ART and adding additional aerosols with different coatings [@art202504].


# State of the field

There are many Mie scattering libraries available in a variety of programming languages. Python is probably the most accessible programming language for quickly and easily analysing and testing data within the scientific community. ARTmie is not the first Python module for Mie scattering, so it will be compared to three other commonly used Python libraries here:

1. miepython [@swMiepython] :

   miepython is a minimal, pure Python implementation of classical Mie theory that provides efficient scattering, absorption and extinction calculations for homogeneous spheres only. While it is ideal for simple, educational or small-scale use cases, it lacks support for coated particles or size distributions. Its performance is limited due to the absence of a compiled backend accelerator.

2. PyMieScatt [@articlePyMieScatt] :

   PyMieScatt is a feature-rich Python library that supports homogeneous and single-coated spheres. It also offers additional tools for data analysis, angular scattering and inversion methods. However, as with miepython, it is implemented entirely in Python, which makes it significantly slower than ARTmie for coated particles and polydisperse systems, where the computational demand is higher.

3. PyMieSim [@articlePyMieSim,@swPyMieSim] :

   Designed for simulating optical systems, PyMieSim provides a modular interface for defining sources, detectors and particle configurations. It offers support for both uncoated and single-coated spheres. Although it uses a C++ backend to optimise performance, the modular design introduces additional setup overhead. PyMieSim is best suited to integrating Mie scattering into broader optical modelling workflows, but does not yet support polydisperse particles.


# Acknowledgement
This work was done during the funding by ...

# References
