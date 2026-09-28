#ifndef ARTMIE_HELPER_H
#define ARTMIE_HELPER_H

#include <exception>
#include <complex>

#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/ndarrayobject.h>
#include <numpy/npy_3kcompat.h>
#include <numpy/npy_math.h>



/* FUNCTIONAL HELPERS */

int parse_arrays(int arg_count, int dtype, PyObject **args, PyObject **arrs);

PyArrayObject* c2py_dblarr(int c_arr_len, double *c_arr);
PyArrayObject* c2py_dblarr(int dim1_len, int dim2_len, double *c_arr);
void py2c_dblarr(PyArrayObject* pyarr, double *carr);

std::complex<double> py2c_cplx(Py_complex np_cplx);
Py_complex c2py_cplx(std::complex<double> c_cplx);
Py_complex nanPyCplx();
PyArrayObject* c2py_cplxarr(int c_arr_len, std::complex<double> *c_arr);
int py2c_cplxarr(PyArrayObject* pyarr, std::complex<double> *carr);

const char* shape2str(int ndim, npy_intp* shape);




/* OPTICS HELPER */

int calc_nmax(double x);
int calc_angles_count(double angres);
int calc_sdo_bin_count(double geom_std, double res);




#endif /* ARTMIE_HELPER_H */