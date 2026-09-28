#ifndef ARTMIE_SINGLE_H
#define ARTMIE_SINGLE_H




struct MieResult {
    double qext;
    double qsca;
    double qabs;
    double qback;
    double qpr;
    double qg;
    double qratio;
};




#define abtomie_docstring "ab2mie(an, bn, wavelength, diameter, /, asCrossSection=False, asDict=False)\n\n\
Parameters\n----------\n\
an : array-like, 1dimensional, complex numbers\n    external field coefficients $a_n$\n\
bn : array-like, 1dimensional, complex numbers\n    external field coefficients $b_n$\n\
wavelength : scalar, floating point number\n    the wavelength of incident light, in nm\n\
diameter : scalar, floating point number\n    the diameter of the whole particle, in nm\n\
asCrossSection : scalar, bool, optional\n    if specified and set to True, returns the results as optical cross-sections with units of nm$^2$.\n\
asDict : scalar, bool, optional\n    if specified and set to True, returns the results as a dictionary.\n\n\
Returns\n-------\n\
qext,qsca,qabs,qback,qratio,qpr,g :\n    scalars, floating point numbers\n    Mie efficiencies as described in MieQ()\n\
qext,qsca,qabs,qback,qratio,qpr,g :\n    scalars, floating point numbers\n    Mie efficiencies as optical cross sections, if asCrossSection set to True\n\
q : dict\n   dictionary of the Mie efficiencies, if asDict is set to True\n\
c : dict\n   dictionary of Mie efficiencies as optical cross sections, if asDict and asCrossSection are both set to True"
MieResult ab2mie(int nmax, std::complex<double> an[], std::complex<double> bn[], double wavelength, double diameter, int asCrossSection);
PyObject* mie_art_ab2mie(PyObject *self, PyObject *args, PyObject *kwds);

#define mieq_docstring "MieQ(m, diam, wavelength, /, nMedium=1.0, asCrossSection=False, asDict=False)\n\n\
Computes extinction, scattering, backscattering and absorption efficiencies, radiation pressure and asymmetry parameter\n\n\
Parameters\n----------\n\
m : scalar or array-like, complex number\n    refractive index of the particle reduced by the refractive index of the surrounding medium\n\
diam : scalar, floating point number\n    the diameter of the whole particle, in nm\n\
wavelength : scalar or array-like, floating point number\n    the wavelength of incident light, in nm\n    has to be of the same shape as m\n\
nMedium : scalar, floating point number, optional\n    the refractive index of the surrounding medium without the extinction part\n\
asCrossSection : scalar, bool, optional\n    if specified and set to True, returns the results as optical cross-sections with units of nm$^2$.\n\
asDict : scalar, bool, optional\n    if specified and set to True, returns the results as a dictionary\n\n\
Returns\n-------\n\
qext,qsca,qabs,qback,qratio,qpr,g :\n    scalars or array-like, floating point numbers, same shape as m\n    Mie efficiencies: extinction, scattering, absorption, backscattering, and backscatter-ratio, radiation pressure and asymmetry parameter\n\
cext,csca,cabs,cback,cratio,cpr,g :\n    scalars or array-like, floating point numbers, same shape as m\n    Mie efficiencies as optical cross sections, if asCrossSection set to True\n\
q : dict\n    dictionary of the Mie efficiencies, if asDict is set to True\n    entries have the same shape as m (scalar or array-like)\n\
c : dict\n    dictionary of Mie efficiencies as optical cross sections, if asDict and asCrossSection are both set to True\n    entries have the same shape as m (scalar or array-like)"
PyObject* mie_art_mieq(PyObject *self, PyObject *args, PyObject *kwds);

#define miecoatedq_docstring "MieCoatedQ(m_core, diam_core, m_shell, diam_shell, wavelength, /, nMedium=1.0, asCrossSection=False, asDict=False)\n\n\
As MieQ() but for coated particles\n\n\
Parameters\n----------\n\
m_core : scalar or array-like, complex number\n    refractive index of the particle's core\n\
diam_core : scalar, floating point number\n    the diameter of the particle's core, in nanometers\n\
m_shell : scalar or array-like, complex number\n    refractive index of the coating\n    same shape as m_core\n\
diam_shell : scalar, floating point number\n    the diameter of the particle and its coating, in nanometers\n\
wavelength : scalar or array-like, floating point number\n    the wavelength of incident light, in nanometers\n    same shape as m_core\n\
nMedium : scalar, floating point number, optional\n    the refractive index of the surrounding medium without the extinction part\n\
asCrossSection : scalar, bool, optional\n    if specified and set to True, returns the results as optical cross-sections with units of nm$^2$.\n\
asDict : scalar, bool, optional\n    if specified and set to True, returns the results as a dictionary.\n\n\
Returns\n-------\n\
qext,qsca,qabs,qback,qratio,qpr,g :\n    scalars or array-like, floating point numbers, same shape as m_core\n    Mie efficiencies as described in MieQ()\n\
qext,qsca,qabs,qback,qratio,qpr,g :\n    scalars or array-like, floating point numbers, same shape as m_core\n    Mie efficiencies as optical cross sections, if asCrossSection set to True\n\
q : dict\n    dictionary of the Mie efficiencies, if asDict is set to True\n    entries have the same shape as m_core (scalar or array-like)\n\
c : dict\n    dictionary of Mie efficiencies as optical cross sections, if asDict and asCrossSection are both set to True\n    entries have the same shape as m_core (scalar or array-like)"
PyObject* mie_art_miecoatedq(PyObject *self, PyObject *args, PyObject *kwds);

#define cbs_docstring "calcVolBackscattering(x, an, bn, theta, dtheta, scatwts, pin, taun)\n\n\
Calculates the volumetric Mie backscattering efficiency.\n\n\
Parameters\n----------\n\
x : scalar, floating point value\n    size parameter of the whole particle\n\
an : array-like, 1dimensional, complex numbers\n    external field coefficients $a_n$\n\
bn : array-like, 1dimensional, complex numbers\n    external field coefficients $b_n$\n\
theta : array-like, 1dimensional, floating point numbers\n    all angels $\\theta$, for which the scattering function is calculated\n\
dtheta : array-like, 1dimensional, floating point numbers\n    angle step size for each angle $\\theta$\n\
scatwgts : array-like, 1dimensional, floating point numbers\n    angel dependent scattering weights\n\
pin : array-like, 2dimensional, floating point numbers\n    field coefficients $\\pi_n$ for every angle and every corresponding external field coefficient $a_n$\n\
taun : array-like, 2dimensional, floating point numbers\n    field coefficients $\\tau_n$ for every angle and every corresponding external field coefficient $a_n$\n\n\
Returns\n-------\n\
s : scalar, floating point number\n    backscatter coefficient\n    It is not the same as the provided one through MieQ() or ab2mie()!"
double calcVolBackscattering(double x, int anbn_len, std::complex<double> *an, std::complex<double> *bn, int nang, double *theta, double *dtheta, double *scatwgts, int nmax, double *pin, double *taun);
PyObject* mie_art_calcVolBackScat(PyObject *self, PyObject *args, PyObject *kwds);

#endif /* ARTMIE_SINGLE_H */