#ifndef ARTMIE_SDO_H
#define ARTMIE_SDO_H

#include <Python.h>


struct Mie_tots {
    double bext;
    double bsca;
    double babs;
    double bback;
    double bvolbsc;
    double bssa;
    double bratio;
    double basym;
//    int arr_len;
//    double *ext_arr; //arrays have to be deleted by the caller of size_distribution_optics
//    double *sca_arr;
//    double *abs_arr;
//    double *bck_arr;
//    double *g_arr;
};


#define clnd_docstring "createLogNormalDistribution(mean_diam, stdev_diam, /, fcoat=0.0, res=0.0, norm2core=False, norm2volume=True)\n\n\
Calculates the parameters regarding the log-normal particle size distribution needed internally by Size_Distribution_Optics and Size_Distribution_Phase_Function\n\n\
Parameters\n----------\n\
mean_diam : scalar, floating point number\n    the median count diameter of the particles in nanometers\n\
geom_std : scalar, floating point number\n    the geometric standard deviation of the particle size distribution\n\
fcoat : scalar, floating point number, optional\n    the coating fraction, default 0.0\n\
res : scalar, floating point number, optional\n    bin size of the particle size distribution in log10-space\n    resulting in a t least 16 bins for the range from -4 to +4 geom_std\n    default 1.0\n\
Returns\n-------\n\
x_range : array-like, 1dimensional, floating point numbers\n    diameters of the particle cores\n\
y_range : array-like, 1dimensional, floating point numbers\n    diameters of the particles including the coating\n\
pdf : array-like, 1dimensional, floating point numbers\n    approximated relative fraction of particle count"
void createLogNormalDistribution(double d_gn, double sigma_g, double fcoat, double res, double /*out*/ *x_range, double /*out*/ *y_range, double /*out*/ *pdf);
PyObject* mie_art_createLgNormDist(PyObject *self, PyObject *args, PyObject *kwds);

#define sdo_docstring "Size_Distribution_Optics(m, sizepar1, sizepar2, wavelength, /, nMedium=1.0, fcoat=0.0, mc=mp, density=1.0, resolution=10, effcore=True, normalized=True)\n\n\
Parameters\n----------\n\
m : scalar, complex number\n    complex refractive index of the particle (core)\n\
sizepar1 : scalar or 1dimensional array, floating point number(s)\n    mean count diameter (if scalar) or particle sizes (if array) in nanometers\n\
sizepar2 : scalar or 1dimensional array, floating point number(s)\n    geometric std. dev. (if scalar) or dNdlogD in cm$^{-3}$ (if array)\n\
wavelength : scalar, floating point number\n    wavelength of the incident light in nanometers\n\
nMedium : scalar, floating point number, optional\n    refractive index without extinction for the surrounding medium, default 1.0\n\
fcoat : scalar, floating point number, optional\n    coating fraction, ratio of shell thickness to core radius, default 0.0\n\
mc : scalar, complex number, optional\n    complex refractive index of the coating, default m\n\
density : scalar, floating point number, optional\n    particle density in kg/m$^3$, default 1000.0\n\
resolution : scalar, floating point number, optional\n    number of bins per power of magnitude within the particle size distribution, default 10\n    ignored when sizepar1 & sizepar2 array-like\n\
effcore : boolean/logical, optional\n    gives mass specific cross-section as m$^2$/(g of core), default True\n\
msc : as mass specific cross-section in m$^2$/g , default True\n    setting to False gives results as averaged efficiencies\n\n\
Returns\n-------\n\
mie_tots : dictionary\n    contains the Mie efficiencies of a particle size distribution \"Extinction\", \"Scattering\", \"Absorption\", the \"Asymmetry\" parameter and the \"Backscattering\" efficiency specifically calculated from a weighted average over all scattering angles\n\n\
Important Note\n--------------\n\
The size distribution is currently hardcoded to be log-normal. Other distributions may follow in future versions.\n\
1dimensional arguments for sizepar1 and sizepar2 are not implemented yet, they will come in version 0.2.0"
void size_distribution_optics(std::complex<double> m_core, double mean_diam, double geom_std, double wavelength, std::complex<double> m_shell, double fcoating, double resolution, double dens, int effcore, int msc, int debug, Mie_tots *mie_tots);
PyObject* mie_art_sdo(PyObject *self, PyObject *args, PyObject *kwds);



#endif /* ARTMIE_SDO_H */