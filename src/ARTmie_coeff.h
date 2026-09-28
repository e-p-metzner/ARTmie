#ifndef ARTMIE_COEFF_H
#define ARTMIE_COEFF_H

#include <complex>
#include <Python.h>




#define mieab_docstring "Mie_ab(m, x)\n\n\
Computes external field coefficients $a_n$ and $b_n$ based on inputs of refractive index $m$ and size parameter $x=\\pi\\,d_p/\\lambda$.\n\n\
Parameters\n----------\n\
m : scalar, complex number\n    refractive index of the particle reduced by the refractive index of the surrounding medium\n\
x : scalar, floating point number\n    size parameter of the particle\n\n\
Returns\n-------\n\
an : array-like, 1dimensional, complex numbers\n    external field coefficients $a_n$\n\
bn : array-like, 1dimensional, complex numbers\n    external field coefficients $b_n$"
void mie_ab(std::complex<double> m, double x, std::complex<double> *an, std::complex<double> *bn);
PyObject* mie_art_mieab(PyObject *self, PyObject *args, PyObject *kwds);

#define miecoatedab_docstring "MieCoated_ab(m_core, x_core, m_shell, x_shell)\n\n\
Computes external field coefficients $a_n$ and $b_n$ based on inputs of refractive indices $m_core$ and $m_shell$,\n\
and size parameters $x_core=\\pi\\,d_core/\\lambda$ and $x_shell=\\pi\\,d_shell/\\lambda$.\n\n\
Parameters\n----------\n\
m_core : scalar, complex number\n    refractive index of the particle reduced by the refractive index of the surrounding medium\n\
x_core : scalar, floating point number\n    size parameter of the particle's core without coating\n\
m_shell : scalar, complex number\n    refractive index of the particle's coating reduced by the refractive index of the surrounding medium\n\
x_shell : scalar, floating point number\n    size parameter of the particle including the coating shell\n\n\
Returns\n-------\n\
an : array-like, 1dimensional, complex numbers\n    external field coefficients $a_n$\n\
bn : array-like, 1dimensional, complex numbers\n    external field coefficients $b_n$"
void miecoated_ab(std::complex<double> m_core, double x_core, std::complex<double> m_shell, double x_shell, std::complex<double> *an, std::complex<double> *bn);
PyObject* mie_art_miecoatedab(PyObject *self, PyObject *args, PyObject *kwds);

#define miecd_docstring "Mie_cd(m, x)\n\n\
Computes internal field coefficients $c_n$ and $d_n$ based on inputs of refractive index $m$ and size parameter $x=\\pi\\,d_p/\\lambda$.\n\n\
Parameters\n----------\n\
m : scalar, complex number\n    refractive index of the particle reduced by the refractive index of the surrounding medium\n\
x : scalar, floating point number\n    size parameter of the particle\n\n\
Returns\n-------\n\
cn : array-like, 1dimensional, complex numbers\n    internal field coefficients $c_n$\n\
dn : array-like, 1dimensional, complex numbers\n    internal field coefficients $d_n$"
void mie_cd(std::complex<double> m, double x, std::complex<double> *cn, std::complex<double> *dn);
PyObject* mie_art_miecd(PyObject *self, PyObject *args, PyObject *kwds);

#define miepitau_docstring "Mie_pitau(theta, nmax)\n\n\
Calculates angular functions $\\pi_n$ and $\\tau_n$.\n\n\
Parameters\n----------\n\
theta : scalar or array-like (1dimensional), floating point number\n    the scattering angle(s) $\\theta$\n\
nmax : scalar, integer\n    the maximum number of coefficients to compute.\n    Typically, $nmax = floor\\left(2+x+4x^{1/3}\\right)$, but can be given any integer.\n\n\
Returns\n-------\n\
pin : array-like, floating point numbers\n    coefficient series $\\pi_n$\n    1dimensional if a single value for $\\theta$ was given, else 2dimensional\n\
taun : array-like, floating point numbers\n    coefficient series $\\tau_n$\n    1dimensional if a single value for $\\theta$ was given, else 2dimensional"
void mie_pitau(double theta, int nmax, double *pin, double *taun);
void mie_pitau(int nang, double *theta, int nmax, double *pin, double *taun);
PyObject* mie_art_miepitau(PyObject *self, PyObject *args, PyObject *kwds);




#endif /* ARTMIE_COEFF_H */