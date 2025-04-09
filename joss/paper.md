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
authors:
  - name: Enrico P. Metzner
    orcid: 0000-0001-6523-9205
    corresponding: true
    affiliation: '1'
affiliations:
  - index: '1'
    name: Karlsruhe Institute of Technology (KIT), Germany
date: 10 April 2025
bibliography: paper.bib
---

# Summary
The optical properties of particles like dust, ash, droplets and other so called "aerosols" are necessary when calculating interactions with radiation.
The simplest assumption for these particles is that they are spherical, which leads to the mathematical fully derived theorie of Mie.
The corresponding formulae and first computational algorithms were improved and published by @BH1983.

ARTmie uses these formulaes to calculate extinction, scattering, absorption and backscattering efficiencies as well as cross sections by implementing and porting of the relevant functions from matlab, published by @Maetzler2002matlab.

The necessary evaluation of the bessel functions is also ported to C++ based on the well known code by @Amos1986.

ARTmie gains its speed from


# Statement of need
A lot of research in physics and atmospheric and astronomical science depends on the optical properties of aerosols. The more variability in the constituents and the number of particles there is, the more calculations are needed. Due to porting of the code to C++, ARTmie is capable to do these calculations in a fast way which pure python code is never be able to reach.



ARTmie has been used recently to improve the ICON-ART module of the numerical weather and climate prediction model ICON [@Rieger2015; @icon], in the way that all optical properties of aerosols regularly used by ICON-ART were recalculated and additional aerosols with a variaty of coatings were added [@art3].

# Acknowledgement
This work was done during the funding by ...

# References
