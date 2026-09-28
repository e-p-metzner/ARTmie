#ifndef ARTMIE_PHASE_H
#define ARTMIE_PHASE_H

#include <complex>
#include <Python.h>




void scattering_weights(double angres, double *theta, double *dtheta, double *scatwgts);




#define scatfunc_docstring "ScatteringFunction(m, diam, wavelength, theta, /, m_shell=m, fcoat=0.0)\n\n\
Calculates the angle-dependent scattering intensities for parallel, perpendicular polarized and unpolarized light.\n\n\
Parameters\n----------\n\
m : scalar, complex number\n    complex refractive index of the particle (its core, when a coating is given)\n\
diam : scalar, floating point number\n    the diameter of the particle (its core) in nm\n\
wavelength : scalar, floating point number\n    the wavelength of the incident light in nm\n\
theta : array-like, 1dimensional, floating point numbers\n    array of the scattering angles in radians\n\
m_shell : scalar, complex number, optional\n    complex refractive index of the particles coating, defaults to the core's refractive index\n\
fcoat : scalar, floating point number, optional\n    coating fraction as ratio of diameters, default 0.0\n\n\
Results\n-------\n\
sl : array-like, 1dimensional, floating point numbers\n    scattering intensities of perpendicular polarized light\n\
sr : array-like, 1dimensional, floating point numbers\n    scattering intensities of parallel polarized light\n\
su : array-like, 1dimensional, floating point numbers\n    scattering intensities of unpolarized light"
void scattering_function(int anbn_len, std::complex<double> *an, std::complex<double> *bn, int nang, double *theta, int nmax, double *pin, double *taun, double *sl, double *sr, double *su);
PyObject* mie_art_scatfunc(PyObject *self, PyObject *args, PyObject *kwds);




#define cbs_phfunc_docstring "calcVolBackscatteringFromPhFunc(x, theta, phfunc)\n\n\
Calculates the volumetric Mie backscattering efficiency from a provided phase function.\n\n\
Parameters\n----------\n\
x : scalar, floating point number\n    size parameter $x$\n\
theta : array-like, 1dimensional, floating point numbers\n    all angels $\\theta$, for which the phase function was calculated\n\
phfunc : array-like, 1dimensional, floating point numbers\n    phase function values for every angle $\\theta$\n\n\
Returns\n-------\n\
s : scalar, floating point number\n    backscatter coefficient\n    It is not the same as the provided one through MieQ() or ab2mie()!\n\n\
Notes\n-----\n\
If the scattering-angle weighted backscattering coefficient should be calculated for a large number of particles with different  properties,\
the faster way is to use the calcVolBackscattering() and to calculate the backscattering coefficient directly from the series $a_n$, $b_n$, ${\\pi}_n$ and ${\\tau}_n$ with precalculated scattering angles and scatterings weights."
double calc_volbsc_from_phfunc(double x, double *theta, double *phfunc, const int nang);
PyObject* mie_art_calcVolBackScatFromPhFunc(PyObject *self, PyObject *args, PyObject *kwds);




#define sdpf_docstring "Size_Distribution_Phase_Function(m, sizepar1, sizepar2, wavelength, /, nMedium=1.0, fcoat=0.0, mc=mp, density=1.0, resolution=10, effcore=True, normalized=False)\n\n\
Computes the scattering phase function for a log-normal particle size distribution.\n\n\
Parameters\n----------\n\
m : scalar, complex number\n    complex refractive index of the particle (core)\n\
sizepar1 : scalar or 1dimensional array, floating point number(s)\n    mean count diameter (if scalar) or particle sizes (if array) in nanometers\n\
sizepar2 : scalar or 1dimensional array, floating point number(s)\n    geometric std. dev. (if scalar) or dNdlogD in cm$^{-3}$ (if array)\n\
wavelength : scalar, floating point number\n    wavelength of the incident light in nanometers\n\
nMedium : scalar, floating point number, optional\n    refractive index without extinction for the surrounding medium, default 1.0\n\
fcoat : scalar, floating point number, optional\n    coating fraction, ratio of shell thickness to core radius, default 0.0\n\
mc : scalar, complex number, optional\n    complex refractive index of the coating, default m\n\
density : scalar, floating point number, optional\n    particle density in g/cm$^3$, default 1.0\n\
resolution : scalar, floating point number, optional\n    number of bins per power of magnitude within the particle size distribution, default 10\n    ignored when sizepar1 & sizepar2 array-like\n\
effcore : boolean/logical, optional\n    calculates cross-section as nm$^2$/(g of core), default True\n\
normalized : normalized to nm$^2$/g particles, default True\n    setting to False works only with d & dNdlogD (array-like sizeparX)\n\n\
Returns\n-------\n\
theta : array-like, 1dimensional, floating point numbers\n    scattering angles in radians\n\
sl : array-like, 1dimensional, floating point numbers\n    scattering intensities of perpendicular polarized light\n\
sr : array-like, 1dimensional, floating point numbers\n    scattering intensities of parallel polarized light\n\
su : array-like, 1dimensional, floating point numbers\n    scattering intensities of unpolarized light\n\n\
Important Note\n--------------\n\
The size distribution is currently hardcoded to be log-normal. Other distributions may follow in future versions.\n\
1dimensional arguments for sizepar1 and sizepar2 are not implemented yet, they will come in version 0.2.0"
void size_distribution_phase_function(std::complex<double> m_core, double mean_diam, double geom_std, double wavelength, std::complex<double> m_shell, double fcoating, double resolution, double dens, int effcore, int norm2vol, int nang, double *outTheta, double *outSL, double *outSR, double *outSU);
PyObject* mie_art_sdpf(PyObject *self, PyObject *args, PyObject *kwds);




#endif /* ARTMIE_PHASE_H */