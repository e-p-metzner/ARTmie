#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL ARTmieModule

#include "ARTmie_constants.h"
#include "ARTmie_helper.h"
#include "ARTmie_coeff.h"
#include "ARTmie_single.h"




MieResult ab2mie(int nmax, std::complex<double> an[], std::complex<double> bn[], double wavelength, double diameter, int asCrossSection) {
    MieResult res;
    res.qext = 0.0;
    res.qsca = 0.0;
    res.qg   = 0.0;
    std::complex<double> qbck(0.0,0.0);
    int idx;

    for(idx=1; idx<=nmax; idx++) {
        double i = 1.0 - 2.0*(idx&1);
        double f = 2.0*idx+1.0;
        std::complex<double> a = an[idx-1];
        std::complex<double> b = bn[idx-1];
        std::complex<double> ap(0.0,0.0);
        std::complex<double> bp(0.0,0.0);
        if(idx<nmax) {
            ap = an[idx];
            bp = bn[idx];
        }
        res.qext += f*(a.real() + b.real());
        res.qsca += f*(a.real()*a.real()+a.imag()*a.imag() + b.real()*b.real()+b.imag()*b.imag());
        qbck +=     (a-b)*f*i;
        res.qg +=   (idx*idx+2.0*idx)*(a.real()*ap.real() + a.imag()*ap.imag() + b.real()*bp.real() + b.imag()*bp.imag())/(idx+1.0);
        res.qg +=   f*(a.real()*b.real() + a.imag()*b.imag())/(idx*idx+idx);
    }

    double sizeParam = _PI_ * diameter / wavelength;
    double ix2 = 1.0 / (sizeParam*sizeParam);
    double qbabs2 = qbck.real()*qbck.real()+qbck.imag()*qbck.imag();
    res.qpr    = 2.0*ix2*(res.qext - 2.0*res.qg);
    res.qg    *= 2.0/res.qsca;
    res.qback  = ix2*qbabs2;
    res.qratio = 0.5*qbabs2/res.qsca;
    res.qext  *= 2.0*ix2;
    res.qsca  *= 2.0*ix2;
    res.qabs   = res.qext - res.qsca;

    if(asCrossSection>0) {
        double css = 0.25*_PI_*diameter*diameter;
        res.qext  *= css;
        res.qsca  *= css;
        res.qabs  *= css;
        res.qback *= css;
        res.qpr   *= css;
    }

    return res;
}
PyObject* mie_art_ab2mie(PyObject *self, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = { (char*)"an", (char*)"bn", (char*)"wavelength", (char*)"diameter", (char*)"asCrossSection", (char*)"asDict", NULL };

    PyObject *arr_ptr[2] = { NULL, NULL };
    PyObject *array[2] = { NULL, NULL };
    double valueW;
    double valueD;
    int valueCSS  = false;
    int valueDict = false;

    int dtype = -1;
    if(PyArg_ParseTupleAndKeywords(args, kwds, "OOdd|pp", kwlist,
                &arr_ptr[0], &arr_ptr[1], &valueW, &valueD, &valueCSS, &valueDict)) {
        if (parse_arrays(2, NPY_COMPLEX64, arr_ptr, array)) {
            dtype = NPY_COMPLEX64;
        }
        if (dtype<0)
        if (parse_arrays(2, NPY_COMPLEX128, arr_ptr, array)) {
            dtype = NPY_COMPLEX128;
        }
    } else {
        PyErr_Clear();
    }
    if (dtype<0) {
        PyErr_SetString(
            PyExc_TypeError,
            "The arguments dtypes do not match, expected (complex[:],complex[:],float,float[,bool,bool]) or (float64[:],float64[:],float64,int32[,int32,int32])"
        );
        return NULL;
    }
    int ndim1 = PyArray_NDIM((PyArrayObject*)array[0]);
    int ndim2 = PyArray_NDIM((PyArrayObject*)array[1]);
    if (ndim1!=1 || ndim2!=1) {
        Py_XDECREF(array[0]);
        Py_XDECREF(array[1]);
        PyErr_SetString(PyExc_TypeError, "Arrays must be one-dimensional.");
        return NULL;
    }
    int a_len = (int) PyArray_SIZE((PyArrayObject*)array[0]);
    int b_len = (int) PyArray_SIZE((PyArrayObject*)array[1]);
    if (a_len!=b_len) {
        Py_XDECREF(array[0]);
        Py_XDECREF(array[1]);
        PyErr_SetString(PyExc_TypeError, "Arrays have to be of the same length!");
        return NULL;
    }

    std::complex<double>* an = new std::complex<double>[a_len];
    std::complex<double>* bn = new std::complex<double>[b_len];
    py2c_cplxarr((PyArrayObject*)array[0], an);
    py2c_cplxarr((PyArrayObject*)array[1], bn);
    MieResult mr = ab2mie(a_len, an, bn, valueW, valueD, valueCSS);

    PyObject *res;
    if(valueDict>0) {
        if(valueCSS>0) {
            res = Py_BuildValue("{s:d,s:d,s:d,s:d,s:d,s:d,s:d}",
                                "Cext",mr.qext, "Csca",mr.qsca, "Cabs",mr.qabs, "Cback",mr.qback,
                                "Cratio",mr.qratio, "Cpr",mr.qpr, "g",mr.qg);
        } else {
            res = Py_BuildValue("{s:d,s:d,s:d,s:d,s:d,s:d,s:d}",
                                "Qext",mr.qext, "Qsca",mr.qsca, "Qabs",mr.qabs, "Qback",mr.qback,
                                "Qratio",mr.qratio, "Qpr",mr.qpr, "g",mr.qg);
        }
    } else {
        res = Py_BuildValue("ddddddd",mr.qext,mr.qsca,mr.qabs,mr.qback,mr.qratio,mr.qpr,mr.qg);
    }

    Py_DECREF(array[0]);
    Py_DECREF(array[1]);

    delete[] an;
    delete[] bn;
    return res;
}

PyObject* mie_art_mieq(PyObject *self, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = { (char*)"m", (char*)"diam", (char*)"wavelength", (char*)"nMedium", (char*)"asCrossSection", (char*)"asDict", NULL };

    Py_complex valueNpM;
    double valueW;
    double valueD;
    double valueNmedium = 1.0;
    int valueCSS = false;
    int valueDict = false;

    PyObject* arr_ptr[1] =      { NULL };
    PyObject* array[1] =        { NULL };
    PyObject* arr_cplx_ptr[1] = { NULL };
    PyObject* array_cplx[1] =   { NULL };

    int numArrs = -1;
    int ctype = -1;
    int dtype = -1;

    if(PyArg_ParseTupleAndKeywords(args, kwds, "Ddd|dpp", kwlist, &valueNpM, &valueD, &valueW, &valueNmedium, &valueCSS, &valueDict)) {
        numArrs = 0;
    } else {
        PyErr_Clear();
    }

    if(numArrs<0) {
        if(PyArg_ParseTupleAndKeywords(args, kwds, "OdO|dpp", kwlist, &arr_cplx_ptr[0], &valueD, &arr_ptr[0], &valueNmedium, &valueCSS, &valueDict)) {
            if(parse_arrays(1, NPY_COMPLEX64, arr_cplx_ptr, array_cplx))
                ctype = NPY_COMPLEX64;
            if(ctype<0)
            if(parse_arrays(1, NPY_COMPLEX128, arr_cplx_ptr, array_cplx))
                ctype = NPY_COMPLEX128;
            if(ctype<0) {
                PyErr_SetString(
                    PyExc_TypeError,
                    "The m array has to be of type float, double or complex."
                );
                return NULL;
            }
            if(parse_arrays(1, NPY_FLOAT, arr_ptr, array))
                dtype = NPY_FLOAT;
            if(dtype<0)
            if(parse_arrays(1, NPY_DOUBLE, arr_ptr, array))
                dtype = NPY_DOUBLE;
            if(dtype<0) {
                PyErr_SetString(
                    PyExc_TypeError,
                    "The wavelengths have to be an array of type float or double."
                );
                return NULL;
            }
            numArrs = 2;
        } else {
            PyErr_Clear();
        }
    }

    if(numArrs<0) {
        PyErr_SetString(
            PyExc_TypeError,
            "Arguments do not match function definition: MieQ(m, diameter, wavelength, /, nMedium=1.0, asCrossSection=False, asDict=False)"
        );
    }

    PyObject* res = NULL;
    if(numArrs==0) {
        std::complex<double> valueM = py2c_cplx(valueNpM) / valueNmedium;
        double x = _PI_*valueD/valueW;
        int nmax = calc_nmax(x);
        std::complex<double>* an = new std::complex<double>[nmax];
        std::complex<double>* bn = new std::complex<double>[nmax];
        mie_ab(valueM, x, an, bn);
        MieResult mr = ab2mie(nmax, an, bn, valueW, valueD, valueCSS);

        if(valueDict>0) {
            if(valueCSS>0) {
                res = Py_BuildValue("{s:d,s:d,s:d,s:d,s:d,s:d,s:d}",
                                    "Cext",mr.qext, "Csca",mr.qsca, "Cabs",mr.qabs, "Cback",mr.qback,
                                    "Cratio",mr.qratio, "Cpr",mr.qpr, "g",mr.qg);
            } else {
                res = Py_BuildValue("{s:d,s:d,s:d,s:d,s:d,s:d,s:d}",
                                    "Qext",mr.qext, "Qsca",mr.qsca, "Qabs",mr.qabs, "Qback",mr.qback,
                                    "Qratio",mr.qratio, "Qpr",mr.qpr, "g",mr.qg);
            }
        } else {
            res = Py_BuildValue("ddddddd",mr.qext,mr.qsca,mr.qabs,mr.qback,mr.qratio,mr.qpr,mr.qg);
        }

        delete[] an;
        delete[] bn;
    }
    if(numArrs>1) {
        int ndimM = PyArray_NDIM((PyArrayObject*)array_cplx[0]);
        int ndimW = PyArray_NDIM((PyArrayObject*)array[0]);
        npy_intp* shapeM = PyArray_SHAPE((PyArrayObject*)array_cplx[0]);
        npy_intp* shapeW = PyArray_SHAPE((PyArrayObject*)array[0]);
        int mw_same = true;
        if(ndimM==ndimW) {
            for(int i=0; i<ndimM && mw_same; i++)
                mw_same = (shapeM[i]==shapeW[i]);
        } else {
            mw_same = false;
        }

        if(!mw_same) {
            Py_XDECREF(array_cplx[0]);
            Py_XDECREF(array[0]);
            PyErr_SetString(
                PyExc_IndexError,
                "m and wavelength have to be of the same shape!"
            );
            return res;
        }

        npy_intp flatDims[1];
        flatDims[0] = PyArray_SIZE((PyArrayObject*)array_cplx[0]);
        int a_len = (int) flatDims[0];
        PyArray_Dims flatShp = { nullptr, 0 };
        flatShp.ptr = flatDims;
        flatShp.len = 1;

        std::complex<double>* valuesM = new std::complex<double>[a_len];
        double* valuesW = new double[a_len];
        if(ndimW==1) {
            py2c_cplxarr((PyArrayObject*)array_cplx[0], valuesM);
            py2c_dblarr( (PyArrayObject*)array[0],      valuesW);
        } else {
            py2c_cplxarr((PyArrayObject*)PyArray_Newshape((PyArrayObject*)array_cplx[0], &flatShp, NPY_CORDER), valuesM);
            py2c_dblarr( (PyArrayObject*)PyArray_Newshape((PyArrayObject*)array[0],      &flatShp, NPY_CORDER), valuesW);
        }

        double* bext = new double[a_len];
        double* bsca = new double[a_len];
        double* babs = new double[a_len];
        double* bbck = new double[a_len];
        double* brat = new double[a_len];
        double* bpr  = new double[a_len];
        double* bg   = new double[a_len];
        for(int i=0; i<a_len; i++) {
            double x = _PI_*valueD/valuesW[i];
            int nmax = calc_nmax(x);
            std::complex<double>* an = new std::complex<double>[nmax];
            std::complex<double>* bn = new std::complex<double>[nmax];
            mie_ab(valuesM[i], x, an, bn);
            MieResult mr = ab2mie(nmax, an, bn, valuesW[i], valueD, valueCSS);
            bext[i] = mr.qext;
            bsca[i] = mr.qsca;
            babs[i] = mr.qabs;
            bbck[i] = mr.qback;
            brat[i] = mr.qratio;
            bpr[i]  = mr.qpr;
            bg[i]   = mr.qg;
            delete[] an;
            delete[] bn;
        }

        PyArray_Dims outShp = { nullptr, 0 };
        outShp.ptr = shapeM;
        outShp.len = ndimM;
        PyObject* pyext_arr = NULL;
        PyObject* pysca_arr = NULL;
        PyObject* pyabs_arr = NULL;
        PyObject* pybck_arr = NULL;
        PyObject* pyrat_arr = NULL;
        PyObject* pypr_arr  = NULL;
        PyObject* pyasy_arr = NULL;
        if(ndimW==1) {
            pyext_arr = (PyObject*)c2py_dblarr(a_len, bext);
            pysca_arr = (PyObject*)c2py_dblarr(a_len, bsca);
            pyabs_arr = (PyObject*)c2py_dblarr(a_len, babs);
            pybck_arr = (PyObject*)c2py_dblarr(a_len, bbck);
            pyrat_arr = (PyObject*)c2py_dblarr(a_len, brat);
            pypr_arr  = (PyObject*)c2py_dblarr(a_len, bpr);
            pyasy_arr = (PyObject*)c2py_dblarr(a_len, bg);
        } else {
            pyext_arr = PyArray_Newshape(c2py_dblarr(a_len, bext), &outShp, NPY_CORDER);
            pysca_arr = PyArray_Newshape(c2py_dblarr(a_len, bsca), &outShp, NPY_CORDER);
            pyabs_arr = PyArray_Newshape(c2py_dblarr(a_len, babs), &outShp, NPY_CORDER);
            pybck_arr = PyArray_Newshape(c2py_dblarr(a_len, bbck), &outShp, NPY_CORDER);
            pyrat_arr = PyArray_Newshape(c2py_dblarr(a_len, brat), &outShp, NPY_CORDER);
            pypr_arr  = PyArray_Newshape(c2py_dblarr(a_len, bpr),  &outShp, NPY_CORDER);
            pyasy_arr = PyArray_Newshape(c2py_dblarr(a_len, bg),   &outShp, NPY_CORDER);
        }
        if(valueDict>0) {
            if(valueCSS>0) {
                res = Py_BuildValue("{s:O,s:O,s:O,s:O,s:O,s:O,s:O}",
                                    "Cext",pyext_arr, "Csca",pysca_arr, "Cabs",pyabs_arr, "Cback",pybck_arr,
                                    "Cratio",pyrat_arr, "Cpr",pypr_arr, "g",pyasy_arr);
            } else {
                res = Py_BuildValue("{s:O,s:O,s:O,s:O,s:O,s:O,s:O}",
                                    "Qext",pyext_arr, "Qsca",pysca_arr, "Qabs",pyabs_arr, "Qback",pybck_arr,
                                    "Qratio",pyrat_arr, "Qpr",pypr_arr, "g",pyasy_arr);
            }
        } else {
            res = Py_BuildValue("OOOOOOO",
                                pyext_arr, pysca_arr, pyabs_arr,pybck_arr,
                                pyrat_arr, pypr_arr, pyasy_arr);
        }

        delete[] valuesM;
        delete[] valuesW;
        delete[] bext;
        delete[] bsca;
        delete[] babs;
        delete[] bbck;
        delete[] brat;
        delete[] bpr;
        delete[] bg;
    }

    return res;
}
PyObject* mie_art_miecoatedq(PyObject *self, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = { (char*)"m_core", (char*)"diam_core", (char*)"m_shell", (char*)"diam_shell", (char*)"wavelength", (char*)"nMedium", (char*)"asCrossSection", (char*)"asDict", NULL };

    Py_complex valueNpMcore;
    Py_complex valueNpMshell;
    double valueW;
    double valueDcore;
    double valueDshell;
    double valueNmedium = 1.0;
    int valueCSS = false;
    int valueDict = false;

    PyObject* arr_ptr[1] =      { NULL };
    PyObject* array[1] =        { NULL };
    PyObject* arr_cplx_ptr[2] = { NULL, NULL };
    PyObject* array_cplx[2] =   { NULL, NULL };

    int numArrs = -1;
    int ctype = -1;
    int dtype = -1;

    if(PyArg_ParseTupleAndKeywords(args, kwds, "DdDdd|dpp", kwlist, &valueNpMcore, &valueDcore, &valueNpMshell, &valueDshell, &valueW, &valueNmedium, &valueCSS, &valueDict)) {
        numArrs = 0;
    } else {
        PyErr_Clear();
    }

    if(numArrs<0) {
        if(PyArg_ParseTupleAndKeywords(args, kwds, "OdOdO|dpp", kwlist, &arr_cplx_ptr[0], &valueDcore, &arr_cplx_ptr[1], &valueDshell, &arr_ptr[0], &valueNmedium, &valueCSS, &valueDict)) {
            if(parse_arrays(2, NPY_COMPLEX64, arr_cplx_ptr, array_cplx))
                ctype = NPY_COMPLEX64;
            if(ctype<0)
            if(parse_arrays(2, NPY_COMPLEX128, arr_cplx_ptr, array_cplx))
                ctype = NPY_COMPLEX128;
            if(ctype<0) {
                PyErr_SetString(
                    PyExc_TypeError,
                    "The arrays m_core and m_shell have to be of type float, double or complex."
                );
                return NULL;
            }
            if(parse_arrays(1, NPY_FLOAT, arr_ptr, array))
                dtype = NPY_FLOAT;
            if(dtype<0)
            if(parse_arrays(1, NPY_DOUBLE, arr_ptr, array))
                dtype = NPY_DOUBLE;
            if(dtype<0) {
                PyErr_SetString(
                    PyExc_TypeError,
                    "The wavelengths have to be an array of type float or double."
                );
                return NULL;
            }
            numArrs = 3;
        } else {
            PyErr_Clear();
        }
    }

    if(numArrs<0) {
        PyErr_SetString(
            PyExc_TypeError,
            "Arguments do not match function definition: MieCoatedQ(m_core, diam_core, m_shell, diam_shell, wavelength, /, nMedium=1.0, asCrossSection=False, asDict=False)"
        );
    }

    PyObject* res = NULL;
    if(numArrs==0) {
        std::complex<double> valueMcore = py2c_cplx(valueNpMcore) / valueNmedium;
        std::complex<double> valueMshell = py2c_cplx(valueNpMshell) / valueNmedium;
        double x = _PI_*valueDcore/valueW;
        double y = _PI_*valueDshell/valueW;
        int nmax = calc_nmax(y);
        std::complex<double>* an = new std::complex<double>[nmax];
        std::complex<double>* bn = new std::complex<double>[nmax];
        miecoated_ab(valueMcore, x, valueMshell, y, an, bn);
        MieResult mr = ab2mie(nmax, an, bn, valueW, valueDshell, valueCSS);

        if(valueDict>0) {
            if(valueCSS>0) {
                res = Py_BuildValue("{s:d,s:d,s:d,s:d,s:d,s:d,s:d}",
                                    "Cext",mr.qext, "Csca",mr.qsca, "Cabs",mr.qabs, "Cback",mr.qback,
                                    "Cratio",mr.qratio, "Cpr",mr.qpr, "g",mr.qg);
            } else {
                res = Py_BuildValue("{s:d,s:d,s:d,s:d,s:d,s:d,s:d}",
                                    "Qext",mr.qext, "Qsca",mr.qsca, "Qabs",mr.qabs, "Qback",mr.qback,
                                    "Qratio",mr.qratio, "Qpr",mr.qpr, "g",mr.qg);
            }
        } else {
            res = Py_BuildValue("ddddddd",mr.qext,mr.qsca,mr.qabs,mr.qback,mr.qratio,mr.qpr,mr.qg);
        }

        delete[] an;
        delete[] bn;
    }
    if(numArrs>1) {
        int ndimMc = PyArray_NDIM((PyArrayObject*)array_cplx[0]);
        int ndimMs = PyArray_NDIM((PyArrayObject*)array_cplx[1]);
        int ndimW = PyArray_NDIM((PyArrayObject*)array[0]);
        npy_intp* shapeMc = PyArray_SHAPE((PyArrayObject*)array_cplx[0]);
        npy_intp* shapeMs = PyArray_SHAPE((PyArrayObject*)array_cplx[1]);
        npy_intp* shapeW = PyArray_SHAPE((PyArrayObject*)array[0]);
        int mw_same = (ndimMc==ndimMs && ndimMc==ndimW);
        if(mw_same) {
            for(int i=0; i<ndimMc && mw_same; i++)
                mw_same = (shapeMc[i]==shapeMs[i] && shapeMc[i]==shapeW[i]);
        }

        if(!mw_same) {
            Py_XDECREF(array_cplx[0]);
            Py_XDECREF(array_cplx[1]);
            Py_XDECREF(array[0]);
            PyErr_Format(
                PyExc_ValueError,
                "m_core, m_shell and wavelength cannot be broadcast together, found shapes %s, %s and %s.",
                shape2str(ndimMc,shapeMc), shape2str(ndimMs,shapeMs), shape2str(ndimW,shapeW)
            );
            return res;
        }

        npy_intp flatDims[1];
        flatDims[0] = PyArray_SIZE((PyArrayObject*)array_cplx[0]);
        int a_len = (int) flatDims[0];
        PyArray_Dims flatShp = { nullptr, 0 };
        flatShp.ptr = flatDims;
        flatShp.len = 1;

        std::complex<double>* valuesMcore  = new std::complex<double>[a_len];
        std::complex<double>* valuesMshell = new std::complex<double>[a_len];
        double*               valuesW      = new double[a_len];
        if(ndimW==1) {
            py2c_cplxarr((PyArrayObject*)array_cplx[0], valuesMcore);
            py2c_cplxarr((PyArrayObject*)array_cplx[1], valuesMshell);
            py2c_dblarr( (PyArrayObject*)array[0],      valuesW);
        } else {
            py2c_cplxarr((PyArrayObject*)PyArray_Newshape((PyArrayObject*)array_cplx[0], &flatShp, NPY_CORDER), valuesMcore);
            py2c_cplxarr((PyArrayObject*)PyArray_Newshape((PyArrayObject*)array_cplx[1], &flatShp, NPY_CORDER), valuesMshell);
            py2c_dblarr( (PyArrayObject*)PyArray_Newshape((PyArrayObject*)array[0],      &flatShp, NPY_CORDER), valuesW);
        }

        double* bext = new double[a_len];
        double* bsca = new double[a_len];
        double* babs = new double[a_len];
        double* bbck = new double[a_len];
        double* brat = new double[a_len];
        double* bpr  = new double[a_len];
        double* bg   = new double[a_len];
        for(int i=0; i<a_len; i++) {
            valuesMcore[i]  /= valueNmedium;
            valuesMshell[i] /= valueNmedium;
            double x = _PI_*valueDcore/valuesW[i];
            double y = _PI_*valueDshell/valuesW[i];
            int nmax = calc_nmax(y);
            std::complex<double>* an = new std::complex<double>[nmax];
            std::complex<double>* bn = new std::complex<double>[nmax];
            miecoated_ab(valuesMcore[i], x, valuesMshell[i], y, an, bn);
            MieResult mr = ab2mie(nmax, an, bn, valuesW[i], valueDshell, valueCSS);
            bext[i] = mr.qext;
            bsca[i] = mr.qsca;
            babs[i] = mr.qabs;
            bbck[i] = mr.qback;
            brat[i] = mr.qratio;
            bpr[i]  = mr.qpr;
            bg[i]   = mr.qg;

            delete[] an;
            delete[] bn;
        }

        PyArray_Dims outShp = { nullptr, 0 };
        outShp.ptr = shapeMc;
        outShp.len = ndimMc;
        PyObject* pyext_arr = NULL;
        PyObject* pysca_arr = NULL;
        PyObject* pyabs_arr = NULL;
        PyObject* pybck_arr = NULL;
        PyObject* pyrat_arr = NULL;
        PyObject* pypr_arr  = NULL;
        PyObject* pyasy_arr = NULL;
        if(ndimW==1) {
            pyext_arr = (PyObject*)c2py_dblarr(a_len, bext);
            pysca_arr = (PyObject*)c2py_dblarr(a_len, bsca);
            pyabs_arr = (PyObject*)c2py_dblarr(a_len, babs);
            pybck_arr = (PyObject*)c2py_dblarr(a_len, bbck);
            pyrat_arr = (PyObject*)c2py_dblarr(a_len, brat);
            pypr_arr  = (PyObject*)c2py_dblarr(a_len, bpr);
            pyasy_arr = (PyObject*)c2py_dblarr(a_len, bg);
        } else {
            pyext_arr = PyArray_Newshape(c2py_dblarr(a_len, bext), &outShp, NPY_CORDER);
            pysca_arr = PyArray_Newshape(c2py_dblarr(a_len, bsca), &outShp, NPY_CORDER);
            pyabs_arr = PyArray_Newshape(c2py_dblarr(a_len, babs), &outShp, NPY_CORDER);
            pybck_arr = PyArray_Newshape(c2py_dblarr(a_len, bbck), &outShp, NPY_CORDER);
            pyrat_arr = PyArray_Newshape(c2py_dblarr(a_len, brat), &outShp, NPY_CORDER);
            pypr_arr  = PyArray_Newshape(c2py_dblarr(a_len, bpr),  &outShp, NPY_CORDER);
            pyasy_arr = PyArray_Newshape(c2py_dblarr(a_len, bg),   &outShp, NPY_CORDER);
        }
        if(valueDict>0) {
            if(valueCSS>0) {
                res = Py_BuildValue("{s:O,s:O,s:O,s:O,s:O,s:O,s:O}",
                                    "Cext",pyext_arr, "Csca",pysca_arr, "Cabs",pyabs_arr, "Cback",pybck_arr,
                                    "Cratio",pyrat_arr, "Cpr",pypr_arr, "g",pyasy_arr);
            } else {
                res = Py_BuildValue("{s:O,s:O,s:O,s:O,s:O,s:O,s:O}",
                                    "Qext",pyext_arr, "Qsca",pysca_arr, "Qabs",pyabs_arr, "Qback",pybck_arr,
                                    "Qratio",pyrat_arr, "Qpr",pypr_arr, "g",pyasy_arr);
            }
        } else {
            res = Py_BuildValue("OOOOOOO",
                                pyext_arr, pysca_arr, pyabs_arr, pybck_arr,
                                pyrat_arr, pypr_arr, pyasy_arr);
        }

        delete[] valuesMcore;
        delete[] valuesMshell;
        delete[] valuesW;
        delete[] bext;
        delete[] bsca;
        delete[] babs;
        delete[] bbck;
        delete[] brat;
        delete[] bpr;
        delete[] bg;
    }

    return res;
}

double calcVolBackscattering(double x, int anbn_len, std::complex<double> *an, std::complex<double> *bn, int nang, double *theta, double *dtheta, double *scatwgts, int nmax, double *pin, double *taun) {
    //an,bn are complex[anbn_len]
    //theta,dtheta,scatwgts are double[nang]
    //pin,taun are double[nang][nmax] flattened into double[nang*nmax]
    double ssq = 0.0;
    int a,n,t;
    for(a=0; a<nang; a++) {
        std::complex<double> S1(0.0,0.0);
        std::complex<double> S2(0.0,0.0);
        t = a*nmax;
        for(n=anbn_len-1; n>=0; n--) {
            double n2 = (2.0+1.0/(n+1))/(n+2.0);
            std::complex<double> s1 = an[n]*pin[t+n] + bn[n]*taun[t+n];
            std::complex<double> s2 = an[n]*taun[t+n] + bn[n]*pin[t+n];
            S1 = S1 + n2*s1;
            S2 = S2 + n2*s2;
        }
        double S1mag2 = S1.real()*S1.real() + S1.imag()*S1.imag();
        double S2mag2 = S2.real()*S2.real() + S2.imag()*S2.imag();
        ssq += (S1mag2+S2mag2) * dtheta[a] * scatwgts[a];
    }
    return ssq/(x*x);
}
PyObject* mie_art_calcVolBackScat(PyObject *self, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = { (char*)"x", (char*)"an", (char*)"bn", (char*)"theta", (char*)"dtheta", (char*)"scatwgts", (char*)"pin", (char*)"taun", NULL };

    double valueX;
    PyObject *ab_ptr[] = { NULL, NULL };
    PyObject *th_ptr[] = { NULL, NULL, NULL };
    PyObject *pt_ptr[] = { NULL, NULL };
    PyObject *ab_array[] = { NULL, NULL };
    PyObject *th_array[] = { NULL, NULL, NULL };
    PyObject *pt_array[] = { NULL, NULL };
    int dtype_ab = -1;
    int dtype_th = -1;
    int dtype_pt = -1;
    if(PyArg_ParseTupleAndKeywords(args, kwds, "dOOOOOOO", kwlist, &valueX, &ab_ptr[0], &ab_ptr[1], &th_ptr[0], &th_ptr[1], &th_ptr[2], &pt_ptr[0], &pt_ptr[1])) {
        {
            if(parse_arrays(2, NPY_COMPLEX64, ab_ptr, ab_array)) dtype_ab = NPY_COMPLEX64;
        }
        if(dtype_ab<0) {
            if(parse_arrays(2, NPY_COMPLEX128, ab_ptr, ab_array)) dtype_ab = NPY_COMPLEX128;
        }
        {
            if(parse_arrays(3, NPY_FLOAT, th_ptr, th_array)) dtype_th = NPY_FLOAT;
        }
        if(dtype_th<0) {
            if(parse_arrays(3, NPY_DOUBLE, th_ptr, th_array)) dtype_th = NPY_DOUBLE;
        }
        {
            if(parse_arrays(2, NPY_FLOAT, pt_ptr, pt_array)) dtype_pt = NPY_FLOAT;
        }
        if(dtype_pt<0) {
            if(parse_arrays(2, NPY_DOUBLE, pt_ptr, pt_array)) dtype_pt = NPY_DOUBLE;
        }
    } else {
        PyErr_Clear();
    }
    if(dtype_ab<0 || dtype_th<0 || dtype_pt<0) {
        PyErr_SetString(
            PyExc_TypeError,
            "The arguments dtypes do not match, expected x to be float or double, an and bn of same type complex[], theta and dtheta and scatwgts of same type float[] or double[], pin and taun of same type float[][] or double[][]"
        );
        Py_XDECREF(ab_array[0]); Py_XDECREF(ab_array[1]);
        Py_XDECREF(th_array[0]); Py_XDECREF(th_array[1]); Py_XDECREF(th_array[2]);
        Py_XDECREF(pt_array[0]); Py_XDECREF(pt_array[1]);
        return NULL;
    }

    int ndimA = PyArray_NDIM((PyArrayObject *)ab_array[0]);
    int ndimB = PyArray_NDIM((PyArrayObject *)ab_array[1]);
    int sizeA = (int) PyArray_SIZE((PyArrayObject *)ab_array[0]);
    int sizeB = (int) PyArray_SIZE((PyArrayObject *)ab_array[1]);
    int nmax = sizeA;
    if(ndimA!=1 || ndimB!=1 || sizeA!=sizeB) {
        PyErr_SetString(
            PyExc_IndexError,
            "The arrays an and bn have to be both 1dimensional and they have be the same length"
        );
        Py_XDECREF(ab_array[0]); Py_XDECREF(ab_array[1]);
        Py_XDECREF(th_array[0]); Py_XDECREF(th_array[1]); Py_XDECREF(th_array[2]);
        Py_XDECREF(pt_array[0]); Py_XDECREF(pt_array[1]);
        return NULL;
    }
    int ndimH = PyArray_NDIM((PyArrayObject *)th_array[0]);
    int ndimD = PyArray_NDIM((PyArrayObject *)th_array[1]);
    int ndimS = PyArray_NDIM((PyArrayObject *)th_array[2]);
    int sizeH = (int) PyArray_SIZE((PyArrayObject *)th_array[0]);
    int sizeD = (int) PyArray_SIZE((PyArrayObject *)th_array[1]);
    int sizeS = (int) PyArray_SIZE((PyArrayObject *)th_array[2]);
    int nang = sizeH;
    if(ndimH!=1 || ndimD!=1 || ndimS!=1 || sizeH!=sizeD || sizeH!=sizeS) {
        PyErr_SetString(
            PyExc_IndexError,
            "The arrays theta, dtheta and scatwgts have to be all 1dimensional and they have to be the same length"
        );
        Py_XDECREF(ab_array[0]); Py_XDECREF(ab_array[1]);
        Py_XDECREF(th_array[0]); Py_XDECREF(th_array[1]); Py_XDECREF(th_array[2]);
        Py_XDECREF(pt_array[0]); Py_XDECREF(pt_array[1]);
        return NULL;
    }
    int ndimP = PyArray_NDIM((PyArrayObject *)pt_array[0]);
    int ndimT = PyArray_NDIM((PyArrayObject *)pt_array[1]);
    npy_intp *shapeP = PyArray_DIMS((PyArrayObject *)pt_array[0]);
    npy_intp *shapeT = PyArray_DIMS((PyArrayObject *)pt_array[1]);
    if (ndimP!=2 || ndimT!=2) {
        PyErr_SetString(
            PyExc_IndexError,
            "Both arrays pin and taun have to be two-dimensional"
        );
        Py_XDECREF(ab_array[0]); Py_XDECREF(ab_array[1]);
        Py_XDECREF(th_array[0]); Py_XDECREF(th_array[1]); Py_XDECREF(th_array[2]);
        Py_XDECREF(pt_array[0]); Py_XDECREF(pt_array[1]);
        return NULL;
    }
    if (nang!=(int)shapeP[0] || nang!=(int)shapeT[0] || nmax>shapeP[1] || shapeP[1]!=shapeT[1]) {
        PyErr_SetString(
            PyExc_IndexError,
            "Both arrays pin and taun have to be of the same shape (len(theta),>=len(an))"
        );
        Py_XDECREF(ab_array[0]); Py_XDECREF(ab_array[1]);
        Py_XDECREF(th_array[0]); Py_XDECREF(th_array[1]); Py_XDECREF(th_array[2]);
        Py_XDECREF(pt_array[0]); Py_XDECREF(pt_array[1]);
        return NULL;
    }

    std::complex<double>* an = new std::complex<double>[nmax];
    std::complex<double>* bn = new std::complex<double>[nmax];
    py2c_cplxarr((PyArrayObject *)ab_array[0], an);
    py2c_cplxarr((PyArrayObject *)ab_array[1], bn);
    double* theta = new double[nang];
    double* dtheta = new double[nang];
    double* scatwgts = new double[nang];
    py2c_dblarr((PyArrayObject *)th_array[0], theta);
    py2c_dblarr((PyArrayObject *)th_array[1], dtheta);
    py2c_dblarr((PyArrayObject *)th_array[2], scatwgts);
    int pitau_len = nang*nmax;
    double *pin = new double[pitau_len];
    double *taun = new double[pitau_len];
    py2c_dblarr((PyArrayObject *)pt_array[0], pin);
    py2c_dblarr((PyArrayObject *)pt_array[1], taun);

    double backscat = calcVolBackscattering(valueX, nmax, an, bn, nang, theta, dtheta, scatwgts, (int)shapeP[1], pin, taun);

    PyObject *res = Py_BuildValue("d", backscat);
    Py_DECREF(ab_array[0]);
    Py_DECREF(ab_array[1]);
    Py_DECREF(th_array[0]);
    Py_DECREF(th_array[1]);
    Py_DECREF(th_array[2]);
    Py_DECREF(pt_array[0]);
    Py_DECREF(pt_array[1]);
    delete[] an;
    delete[] bn;
    delete[] theta;
    delete[] dtheta;
    delete[] scatwgts;
    delete[] pin;
    delete[] taun;
    return res;
}
