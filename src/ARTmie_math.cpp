#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL ARTmieModule

#include "ARTmie_constants.h"
#include "ARTmie_amos.h"
#include "ARTmie_helper.h"
#include "ARTmie_math.h"




// **** Gamma function

PyObject* mie_art_gamma(PyObject *self, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = { (char*)"x", NULL };

    //use function for integer array and integer value
    double valueX;
    if(PyArg_ParseTupleAndKeywords(args, kwds, "d", kwlist, &valueX)) {
    } else {
        PyErr_SetString(
            PyExc_TypeError,
            "The arguments dtypes do not match, expected (float)"
        );
        return NULL;
    }
    PyObject *res = Py_BuildValue("d",0.0+ std::exp(dgamln(valueX)));
    return res;
}




// **** Bessel functions

PyObject* mie_art_besselj(PyObject *self, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = { (char*)"v", (char*)"z", (char*)"es", (char*)"debug", NULL };

    int valueExpScl = false;
    int valueDebug = false;
    double valueV;
    Py_complex valueNpZ;

    PyObject* arr_cplx_ptr[1] = { NULL };
    PyObject* array_cplx[1]   = { NULL };

    int numArrs = -1;
    int ctype = -1;

    if(PyArg_ParseTupleAndKeywords(args, kwds, "dD|pp", kwlist, &valueV, &valueNpZ, &valueExpScl, &valueDebug)) {
        numArrs = 0;
    } else {
        PyErr_Clear();
    }

    if(numArrs < 0) {
        if(PyArg_ParseTupleAndKeywords(args, kwds, "dO|pp", kwlist, &valueV, &arr_cplx_ptr[0], &valueExpScl, &valueDebug)) {
            if(parse_arrays(1, NPY_COMPLEX64, arr_cplx_ptr, array_cplx))
                ctype = NPY_COMPLEX64;
            if(ctype<0)
            if(parse_arrays(1, NPY_COMPLEX128, arr_cplx_ptr, array_cplx))
                ctype = NPY_COMPLEX128;
            if(ctype<0) {
                PyErr_SetString(
                    PyExc_TypeError,
                    "The z array has to be of type float, double or complex."
                );
                return NULL;
            }
            numArrs = 1;
        } else {
            PyErr_Clear();
        }
    }

    if(numArrs < 0) {
        PyErr_SetString(
            PyExc_TypeError,
            "The arguments dtypes do not match, expected (float,complex)"
        );
        return NULL;
    }
    double bjr[1], bji[1];
    int nz, idum;

    PyObject* res = NULL;
    if(numArrs == 0) {
        std::complex<double> valueZ = py2c_cplx(valueNpZ);
        PySys_WriteStdout("Bessel J: nu=%f  z=%f+i*%f  expscl=%i\n",valueV,valueZ.real(),valueZ.imag(),valueExpScl);
        zbesj(valueZ.real(), valueZ.imag(), valueV, 1+valueExpScl, 1, bjr, bji, &nz, &idum, &valueDebug);
        if(nz>=0 && (idum==0 || idum==3)) {
            valueZ = std::complex<double>(bjr[0],bji[0]);
        } else {
            valueZ = std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
        }
        res = Py_BuildValue("O", PyComplex_FromDoubles(valueZ.real(),valueZ.imag()));
    }

    if(numArrs > 0) {
        // int       ndimZ  = PyArray_NDIM( (PyArrayObject*)array_cplx[0]);
        // npy_intp* shapeZ = PyArray_SHAPE((PyArrayObject*)array_cplx[0]);
        // npy_intp  flatDims[1];
        // flatDims[0]      = PyArray_SIZE( (PyArrayObject*)array_cplx[0]);
        // int       a_len  = (int) flatDims[0];
        // std::complex<double>* valuesZ = new std::complex<double>[a_len];
        // try {
        //     if(ndimZ==1) {
        //         py2c_cplxarr((PyArrayObject*)array_cplx[0], valuesZ);
        //     } else {
        //         PyArray_Dims flatShp = { nullptr, 0 };
        //         flatShp.ptr = flatDims;
        //         flatShp.len = 1;
        //         py2c_cplxarr((PyArrayObject*)PyArray_Newshape((PyArrayObject*)array_cplx[0], &flatShp, NPY_CORDER), valuesZ);
        //     }
        // } catch(const std::exception &e) {
        //     PyErr_SetString(
        //         PyExc_RuntimeError,
        //         e.what()
        //     );
        //     delete[] valuesZ;
        //     return NULL;
        // }
        // for(int i=0; i<a_len; i++) {
        //     PySys_WriteStdout("Bessel J: nu=%f  z=%f+i*%f  expscl=%i\n",valueV,valuesZ[i].real(),valuesZ[i].imag(),valueExpScl);
        //     zbesj(valuesZ[i].real(), valuesZ[i].imag(), valueV, 1+valueExpScl, 1, bjr, bji, &nz, &idum, &valueDebug);
        //     if(nz>=0 && (idum==0 || idum==3)) {
        //         valuesZ[i] = std::complex<double>(bjr[0],bji[0]);
        //     } else {
        //         valuesZ[i] = std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
        //     }
        // }
        // try {
        //     // if(ndimZ==1) {
        //     //     res = Py_BuildValue("O", c2py_cplxarr(a_len, valuesZ));
        //     // } else {
        //     //     PyArray_Dims outShp = { nullptr, 0 };
        //     //     outShp.ptr = shapeZ;
        //     //     outShp.len = ndimZ;
        //     //     res = Py_BuildValue("O", PyArray_Newshape(c2py_cplxarr(a_len, valuesZ), &outShp, NPY_CORDER));
        //     // }
        //     res = Py_BuildValue("O", PyComplex_FromDoubles(valuesZ[0].real(),valuesZ[0].imag()));
        // } catch(const std::exception &e) {
        //     PyErr_SetString(
        //         PyExc_TypeError,
        //         e.what()
        //     );
        //     res = NULL;
        // }
        // delete[] valuesZ;
        res = Py_BuildValue("O", PyComplex_FromDoubles(0.1234, 0.5678));
    }

    return res;
}
PyObject* mie_art_bessely(PyObject *self, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = { (char*)"v", (char*)"z", (char*)"es", (char*)"debug", NULL };

    int valueExpScl = false;
    int valueDebug = false;
    double valueV;
    Py_complex valueNpZ;

    PyObject* arr_cplx_ptr[1] = { NULL };
    PyObject* array_cplx[1]   = { NULL };

    int numArrs = -1;
    int ctype = -1;

    if(PyArg_ParseTupleAndKeywords(args, kwds, "dD|pp", kwlist, &valueV, &valueNpZ, &valueExpScl, &valueDebug)) {
        numArrs = 0;
    } else {
        PyErr_Clear();
    }

    if(numArrs < 0) {
        if(PyArg_ParseTupleAndKeywords(args, kwds, "dO|pp", kwlist, &valueV, &arr_cplx_ptr[0], &valueExpScl, &valueDebug)) {
            if(parse_arrays(1, NPY_COMPLEX64, arr_cplx_ptr, array_cplx))
                ctype = NPY_COMPLEX64;
            if(ctype<0)
            if(parse_arrays(1, NPY_COMPLEX128, arr_cplx_ptr, array_cplx))
                ctype = NPY_COMPLEX128;
            if(ctype<0) {
                PyErr_SetString(
                    PyExc_TypeError,
                    "The z array has to be of type float, double or complex."
                );
                return NULL;
            }
            numArrs = 1;
        } else {
            PyErr_Clear();
        }
    }

    if(numArrs < 0) {
        PyErr_SetString(
            PyExc_TypeError,
            "The arguments dtypes do not match, expected (float,complex)"
        );
        return NULL;
    }

    double byr[1], byi[1], cwr[1], cwi[1];
    int nz, idum;

    PyObject* res = NULL;
    if(numArrs == 0) {
        zbesy(valueNpZ.real, valueNpZ.imag, valueV, 1+valueExpScl, 1, byr, byi, &nz, cwr, cwi, &idum, &valueDebug);
        std::complex<double> valueZ;
        if(nz>=0 && (idum==0 || idum==3)) {
            valueZ = std::complex<double>(byr[0],byi[0]);
        } else {
            valueZ = std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
        }
        res = Py_BuildValue("O", PyComplex_FromDoubles(valueZ.real(),valueZ.imag()));
    }

    if(numArrs > 0) {
        int       ndimZ  = PyArray_NDIM( (PyArrayObject*)array_cplx[0]);
        npy_intp* shapeZ = PyArray_SHAPE((PyArrayObject*)array_cplx[0]);
        npy_intp  flatDims[1];
        flatDims[0]      = PyArray_SIZE( (PyArrayObject*)array_cplx[0]);
        int       a_len  = (int) flatDims[0];
        std::complex<double>* valuesZ = new std::complex<double>[a_len];
        if(ndimZ==1) {
            py2c_cplxarr((PyArrayObject*)array_cplx[0], valuesZ);
        } else {
            PyArray_Dims flatShp = { nullptr, 0 };
            flatShp.ptr = flatDims;
            flatShp.len = 1;
            py2c_cplxarr((PyArrayObject*)PyArray_Newshape((PyArrayObject*)array_cplx[0], &flatShp, NPY_CORDER), valuesZ);
        }
        for(int i=0; i<a_len; i++) {
            zbesy(valuesZ[i].real(), valuesZ[i].imag(), valueV, 1+valueExpScl, 1, byr, byi, &nz, cwr, cwi, &idum, &valueDebug);
            if(nz>=0 && (idum==0 || idum==3)) {
                valuesZ[i] = std::complex<double>(byr[0],byi[0]);
            } else {
                valuesZ[i] = std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
            }
        }
        if(ndimZ==1) {
            res = Py_BuildValue("O", c2py_cplxarr(a_len, valuesZ));
        } else {
            PyArray_Dims outShp = { nullptr, 0 };
            outShp.ptr = shapeZ;
            outShp.len = ndimZ;
            res = Py_BuildValue("O", PyArray_Newshape(c2py_cplxarr(a_len, valuesZ), &outShp, NPY_CORDER));
        }
        delete[] valuesZ;
    }

    return res;
}
PyObject* mie_art_hankel(PyObject *self, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = { (char*)"v", (char*)"z", (char*)"m", (char*)"es", (char*)"debug", NULL };

    int valueExpScl = false;
    int valueDebug = false;
    int valueM;
    double valueV;
    Py_complex valueNpZ;

    PyObject* arr_cplx_ptr[1] = { NULL };
    PyObject* array_cplx[1]   = { NULL };

    int numArrs = -1;
    int ctype = -1;

    if(PyArg_ParseTupleAndKeywords(args, kwds, "dDi|pp", kwlist, &valueV, &valueNpZ, &valueM, &valueExpScl, &valueDebug)) {
        numArrs = 0;
    } else {
        PyErr_Clear();
    }

    if(numArrs < 0) {
        if(PyArg_ParseTupleAndKeywords(args, kwds, "dOi|pp", kwlist, &valueV, &arr_cplx_ptr[0], &valueM, &valueExpScl, &valueDebug)) {
            if(parse_arrays(1, NPY_COMPLEX64, arr_cplx_ptr, array_cplx))
                ctype = NPY_COMPLEX64;
            if(ctype<0)
            if(parse_arrays(1, NPY_COMPLEX128, arr_cplx_ptr, array_cplx))
                ctype = NPY_COMPLEX128;
            if(ctype<0) {
                PyErr_SetString(
                    PyExc_TypeError,
                    "The z array has to be of type float, double or complex."
                );
                return NULL;
            }
            numArrs = 1;
        } else {
            PyErr_Clear();
        }
    }

    double bhr[1], bhi[1];
    int nz, idum;

    PyObject* res = NULL;
    if(numArrs == 0) {
    	zbesh(valueNpZ.real,valueNpZ.imag, valueV, 1+valueExpScl, valueM, 1, bhr,bhi, &nz, &idum, &valueDebug);
        std::complex<double> valueZ;
        if(nz>=0 && (idum==0 || idum==3)) {
            valueZ = std::complex<double>(bhr[0],bhi[0]);
        } else {
            valueZ = std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
        }
        res = Py_BuildValue("O", PyComplex_FromDoubles(valueZ.real(),valueZ.imag()));
    }

    if(numArrs > 0) {
        int       ndimZ  = PyArray_NDIM( (PyArrayObject*)array_cplx[0]);
        npy_intp* shapeZ = PyArray_SHAPE((PyArrayObject*)array_cplx[0]);
        npy_intp  flatDims[1];
        flatDims[0]      = PyArray_SIZE( (PyArrayObject*)array_cplx[0]);
        int       a_len  = (int) flatDims[0];
        std::complex<double>* valuesZ = new std::complex<double>[a_len];
        if(ndimZ==1) {
            py2c_cplxarr((PyArrayObject*)array_cplx[0], valuesZ);
        } else {
            PyArray_Dims flatShp = { nullptr, 0 };
            flatShp.ptr = flatDims;
            flatShp.len = 1;
            py2c_cplxarr((PyArrayObject*)PyArray_Newshape((PyArrayObject*)array_cplx[0], &flatShp, NPY_CORDER), valuesZ);
        }
        for(int i=0; i<a_len; i++) {
            zbesh(valuesZ[i].real(),valuesZ[i].imag(), valueV, 1+valueExpScl, valueM, 1, bhr,bhi, &nz, &idum, &valueDebug);
            if(nz>=0 && (idum==0 || idum==3)) {
                valuesZ[i] = std::complex<double>(bhr[0],bhi[0]);
            } else {
                valuesZ[i] = std::complex<double>(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
            }
        }
        if(ndimZ==1) {
            res = Py_BuildValue("O", c2py_cplxarr(a_len, valuesZ));
        } else {
            PyArray_Dims outShp = { nullptr, 0 };
            outShp.ptr = shapeZ;
            outShp.len = ndimZ;
            res = Py_BuildValue("O", PyArray_Newshape(c2py_cplxarr(a_len, valuesZ), &outShp, NPY_CORDER));
        }
        delete[] valuesZ;
    }

    return res;
}

PyObject* mie_art_besseli(PyObject *self, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = { (char*)"v", (char*)"z", (char*)"es", (char*)"debug", NULL };

    //use function for integer array and integer value
    int valueExpScl = false;
    int valueDebug = false;
    double valueV;
    Py_complex valueNpZ;
    if(PyArg_ParseTupleAndKeywords(args, kwds, "dD|pp", kwlist, &valueV, &valueNpZ, &valueExpScl, &valueDebug)) {
    } else {
        PyErr_SetString(
            PyExc_TypeError,
            "The arguments dtypes do not match, expected (float,complex)"
        );
        return NULL;
    }

    double bir[1], bii[1];
    std::complex<double> bi = std::complex<double>(
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN());;
    int nz, idum;
    zbesi(valueNpZ.real, valueNpZ.imag, valueV, 1+valueExpScl, 1, bir, bii, &nz, &idum, &valueDebug);
    if(nz>=0 && (idum==0 || idum==3)) {
        bi = std::complex<double>(bir[0], bii[0]);
    }
    PyObject *res = Py_BuildValue("O", PyComplex_FromDoubles(bi.real(),bi.imag()));
    return res;
}
PyObject* mie_art_besselk(PyObject *self, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = { (char*)"v", (char*)"z", (char*)"es", (char*)"debug", NULL };

    //use function for integer array and integer value
    int valueExpScl = false;
    int valueDebug = false;
    double valueV;
    Py_complex valueNpZ;
    if(PyArg_ParseTupleAndKeywords(args, kwds, "dD|pp", kwlist, &valueV, &valueNpZ, &valueExpScl, &valueDebug)) {
    } else {
        PyErr_SetString(
            PyExc_TypeError,
            "The arguments dtypes do not match, expected (float,complex)"
        );
        return NULL;
    }

    double bkr[1], bki[1];
    std::complex<double> bk = std::complex<double>(
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN());;
    int nz, idum;
    zbesk(valueNpZ.real, valueNpZ.imag, valueV, 1+valueExpScl, 1, bkr, bki, &nz, &idum, &valueDebug);
    if(nz>=0 && (idum==0 || idum==3)) {
        bk = std::complex<double>(bkr[0], bki[0]);
    }
    PyObject *res = Py_BuildValue("O", PyComplex_FromDoubles(bk.real(),bk.imag()));
    return res;
}




// **** Airy function

PyObject* mie_art_airy(PyObject *self, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = { (char*)"z", (char*)"es", (char*)"debug", NULL };

    //use function for integer array and integer value
    int valueExpScl = false;
    int valueDebug = false;
    Py_complex valueNpZ;
    PyObject* arr_ptr[1] = { NULL };
    PyObject* array[1] =   { NULL };
    int numArrs = -5;

    if(PyArg_ParseTupleAndKeywords(args, kwds, "D|pp", kwlist, &valueNpZ, &valueExpScl, &valueDebug)) {
        numArrs = 0;
    } else {
        PyErr_Clear();
    }
    if(numArrs<0) {
        if(PyArg_ParseTupleAndKeywords(args, kwds, "O|pp", kwlist, &arr_ptr[0], &valueExpScl, &valueDebug)) {
            int dtype = -1;
            if(parse_arrays(1, NPY_COMPLEX64, arr_ptr, array))
                dtype = NPY_COMPLEX64;
            if(dtype<0)
            if(parse_arrays(1, NPY_COMPLEX128, arr_ptr, array))
                dtype = NPY_COMPLEX128;
            if(dtype<0) {
                PyErr_SetString(
                    PyExc_TypeError,
                    "z is expected to be of type complex."
                );
                return NULL;
            }
            numArrs = 1;
        } else {
            PyErr_Clear();
        }
    }

    if(numArrs<0) {
        PyErr_SetString(
            PyExc_TypeError,
            "Arguments do not match function defintion: airy(z, /, es=False)");
        return NULL;
    }

    PyObject* res = NULL;
    if(numArrs==0) {
        double air[1], aii[1];
        std::complex<double> ai = std::complex<double>(
                std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN());
        int nz, idum;
        zairy(valueNpZ.real, valueNpZ.imag, 0, 1+valueExpScl, air, aii, &nz, &idum, &valueDebug);
        if(nz>=0 && (idum==0 || idum==3)) {
            ai = std::complex<double>(air[0], aii[0]);
        }
        res = Py_BuildValue("O", PyComplex_FromDoubles(ai.real(),ai.imag()));
    }
    if(numArrs==1) {
        npy_intp flatDims[1];
        flatDims[0] = PyArray_SIZE((PyArrayObject*)array[0]);
        int a_len = (int) flatDims[0];
        PyArray_Dims flatShp = { nullptr, 0 };
        flatShp.ptr = flatDims;
        flatShp.len = 1;

        std::complex<double>* valuesZ = new std::complex<double>[a_len];
        py2c_cplxarr((PyArrayObject*)PyArray_Newshape((PyArrayObject*)array[0], &flatShp, NPY_CORDER), valuesZ);
        double air[1], aii[1];
        int nz, idum;
        for(int i=0; i<a_len; i++) {
            zairy(valuesZ[i].real(),valuesZ[i].imag(), 0, 1+valueExpScl, air, aii, &nz, &idum, &valueDebug);
            if(nz>=0 && (idum==0 || idum==3)) {
                valuesZ[i] = std::complex<double>(air[0], aii[0]);
            } else {
                valuesZ[i] = std::complex<double>(
                        std::numeric_limits<double>::quiet_NaN(),
                        std::numeric_limits<double>::quiet_NaN());
            }
        }
        PyArray_Dims outShp = { nullptr, 0 };
        outShp.ptr = PyArray_SHAPE((PyArrayObject*)array[0]);
        outShp.len = PyArray_NDIM((PyArrayObject*)array[0]);
        res = Py_BuildValue("O", PyArray_Newshape(c2py_cplxarr(a_len, valuesZ), &outShp, NPY_CORDER));
        delete[] valuesZ;
    }
    return res;
}
