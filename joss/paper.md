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

[//]: # (soopy-doopy summary xD)

[//]: # (What infos would I like to read in a random paper, if it is telling me something about this code?)
[//]: # (- general topic -> Mie scattering calculations)
[//]: # (- core features -> single, coated particles and log-normal distributed particles)
[//]: # (- exceptional feature -> backscattering as seen in athmospheric sciences)
[//]: # (- other things...?)

In atmospheric science, astronomy and similar areas, the optical properties of particles like dust, ash and cloud droplets are necessary when calculating interactions with radiation.
The simplest assumption for these particles is that they are spherical, which leads to the mathematical fully derived theorie of Mie.
The corresponding formulae and first computational algorithms were improved and published by @BH1983.

ARTmie uses these formulae to calculate extinction, scattering, absorption and backscattering efficiencies as well as cross sections by implementing and porting of the relevant functions [@Maetzler2002matlab] from Matlab and the bessel functions [@Amos1986] from Fortran 77 to C++.

ARTmie was developed with focus on speed, especially for particle size distribution and scattering angle weighted backscattering efficiency. ARTmie mainly gains its speed not only by using a C++ backend but also optimizing the amount of necessary calculation when computing combined informations, e.g. for log-normal particle size distributions.


# Statement of need

A lot of research in physics and atmospheric and astronomical science depends on the optical properties of aerosols. The more variability in the constituents and the number of particles there is, the more calculations are needed.  Due to porting of the code to C++, ARTmie is capable to do these calculations at a speed that pure Python code can never achieve.

ARTmie has the ability to calculate the backscatter efficiency by weighted averaging over the total backward scattering angles, as used in lidar-based measurement systems.

Furthermore, ARTmie has been used recently to improve the ICON-ART module of the numerical weather and climate prediction model ICON [@Rieger2015; @icon202504] by recalculating all optical properties of aerosols regularly used by ICON-ART and adding additional aerosols with different coatings [@art202504].


# State of the field

Several other python libraries exist to calculate Mie scattering.

1. Miepython [@swMiepython] :

   ...

2. PyMieScatt [@swPyMieScatt] :

   PyMieScatt provides all the functions (with similar names) as ARTmie does. On top of that, this package has inverse function for retrieving the refractive index from measured Mie efficiencies, which is not planned for ARTmie (at the moment).

   PyMieScatt is written in pure Python too, so ARTmie is way faster, especially for calculation of combined optical properties and those of partical size distributions.

3. PyMieSim [@articlePyMieSim] :

   Like ARTmie, PyMieSim is written in C++. It uses a modular approach with more flexibilities and possibilities. This makes it on the other hand less practical for easy Mie efficiency calculations.

   It is also slower in calculating Mie efficiencies of coated spheres than ARTmie is.


# Acknowledgement
This work was done during the funding by ...

# References
