#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL ARTmieModule

#include "ARTmie_constants.h"
#include "ARTmie_helper.h"
#include "ARTmie_coeff.h"
#include "ARTmie_phase.h"




void scattering_weights(double angres, double *theta, double *dtheta, double *scatwgts) {
    int nang = calc_angles_count(angres);
    double ares = angres*_PI_/180.0; // convert from degrees to radians
    int idx;
    for(idx=0; idx<nang; idx++) {
        theta[idx]  = idx*ares;
        dtheta[idx] = (idx==0 || idx==nang-1) ? 0.5*ares : ares;
        double bsflag = 2*idx<nang ? 0.0 : 2*idx>nang ? 1.0 : 0.5;
        double wgt = std::sin(theta[idx]);
        scatwgts[idx] = wgt*bsflag;
    }
}

void scattering_function(int anbn_len, std::complex<double> *an, std::complex<double> *bn, int nang, double *theta, int nmax, double *pin, double *taun, double *sl, double *sr, double *su) {
    int a,n,t;
    for(a=0; a<nang; a++) {
        std::complex<double> S1(0.0,0.0);
        std::complex<double> S2(0.0,0.0);
        t = a*nmax-1;
        for(n=anbn_len; n>0; n--) {
            double n2 = (2.0+1.0/n)/(n+1.0);
            std::complex<double> s1 = an[n-1]*pin[t+n] + bn[n-1]*taun[t+n];
            std::complex<double> s2 = an[n-1]*taun[t+n] + bn[n-1]*pin[t+n];
            S1 = S1 + n2*s1;
            S2 = S2 + n2*s2;
        }
        sl[a] = S1.real()*S1.real() + S1.imag()*S1.imag();
        sr[a] = S2.real()*S2.real() + S2.imag()*S2.imag();
        su[a] = 0.5*(sl[a]+sr[a]);
    }
}
PyObject* mie_art_scatfunc(PyObject *self, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = { (char*)"m", (char*)"diam", (char*)"wavelength", (char*)"theta", (char*)"m_shell", (char*)"fcoat", NULL };

    Py_complex valueNpMcore;
    Py_complex valueNpMshell = nanPyCplx();
    double valueD;
    double valueW;
    double valueFcoat = 0.0;
    PyObject *th_ptr[] = { NULL };
    PyObject *th_arr[] = { NULL };
    int dtype_th = -1;
    if(PyArg_ParseTupleAndKeywords(args, kwds, "DddO|Dd", kwlist, &valueNpMcore, &valueD, &valueW, &th_ptr, &valueNpMshell, &valueFcoat)) {
        {
            if(parse_arrays(1, NPY_FLOAT, th_ptr, th_arr)) dtype_th = NPY_FLOAT;
        }
        if(dtype_th<0) {
            if(parse_arrays(1, NPY_DOUBLE, th_ptr, th_arr)) dtype_th = NPY_DOUBLE;
        }
    } else {
        PyErr_Clear();
    }

    if(dtype_th<0) {
        PyErr_SetString(
            PyExc_TypeError,
            "The arguments dtypes do not match, expected m and m_shell to be complex, dia, wavelength and fcoat to be float or double, theta to be of the type float[] or double[]"
        );
        Py_XDECREF(th_arr[0]);
        return NULL;
    }

    int ndimH = PyArray_NDIM((PyArrayObject *)th_arr[0]);
    if(ndimH!=1) {
        PyErr_SetString(
            PyExc_IndexError,
            "Theta has to be 1dimensional"
        );
        Py_XDECREF(th_arr[0]);
        return NULL;
    }
    int nang = (int) PyArray_SIZE((PyArrayObject *)th_arr[0]);

    std::complex<double> m_core  = py2c_cplx(valueNpMcore);
    std::complex<double> m_shell = py2c_cplx(valueNpMshell);
    if(std::isnan(valueNpMshell.real) || std::isnan(valueNpMshell.imag)) {
        m_shell = py2c_cplx(valueNpMcore);
    }
    double xval = valueD*_PI_/valueW;
    double yval = valueD*(1.0+valueFcoat)*_PI_/valueW;
    int nmax = calc_nmax(yval);
    std::complex<double>* an = new std::complex<double>[nmax];
    std::complex<double>* bn = new std::complex<double>[nmax];

    int is_coated = (valueFcoat>EPS && (m_core!=m_shell));

    if (is_coated) {
        miecoated_ab(m_core, xval, m_shell, yval, an, bn);
    } else {
        mie_ab(m_core, yval, an, bn);
    }

    double* theta = new double[nang];
    py2c_dblarr((PyArrayObject *)th_arr[0], theta);

    int pitau_len = nang*nmax;
    double *pin = new double[pitau_len];
    double *taun = new double[pitau_len];
    mie_pitau(nang, theta, nmax, pin, taun);

    double* outSL = new double[nang];
    double* outSR = new double[nang];
    double* outSU = new double[nang];
    scattering_function(nmax, an, bn, nang, theta, nmax, pin, taun, outSL, outSR, outSU);

    PyObject *res = Py_BuildValue("OOO",
        c2py_dblarr(nang, outSL),
        c2py_dblarr(nang, outSR),
        c2py_dblarr(nang, outSU)
        );
    Py_DECREF(th_arr[0]);

    delete[] an;
    delete[] bn;
    delete[] theta;
    delete[] pin;
    delete[] taun;
    delete[] outSL;
    delete[] outSR;
    delete[] outSU;

    return res;
}

double calc_volbsc_from_phfunc(double x, double *theta, double *phfunc, const int nang) {
    double volbsc = 0.0;
    for(int idx=0; idx<nang; idx++) {
        if ( theta[idx] <= 0.5*_PI_-EPS ) {
            continue;
        }
        double dtheta = 0.0;
        if(idx==0) {
            dtheta = theta[idx+1] - theta[idx];
        } else if(idx==nang-1) {
            dtheta = theta[idx] - theta[idx-1];
        } else {
            dtheta = 0.5 * (theta[idx+1] - theta[idx-1]);
        }
        double wgt = std::sin(theta[idx]);
        if (std::abs(theta[idx]-0.5*_PI_)<EPS) {
            wgt *= 0.5;
        }
        volbsc += phfunc[idx] * dtheta * wgt;
    }
    return volbsc * 2.0 / (x*x);
}
PyObject* mie_art_calcVolBackScatFromPhFunc(PyObject *self, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = { (char*)"x", (char*)"theta", (char*)"phfunc", NULL };

    double valueX;
    PyObject *data_ptr[] = { NULL, NULL };
    PyObject *data_arr[] = { NULL, NULL };
    int dtype = -1;
    if(PyArg_ParseTupleAndKeywords(args, kwds, "dOO", kwlist, &valueX, &data_ptr[0], &data_ptr[1])) {
        {
            if(parse_arrays(2, NPY_FLOAT, data_ptr, data_arr)) dtype = NPY_FLOAT;
        }
        if(dtype<0) {
            if(parse_arrays(2, NPY_DOUBLE, data_ptr, data_arr)) dtype = NPY_DOUBLE;
        }
    } else {
        return NULL;
    }
    if(dtype < 0) {
        PyErr_SetString(
            PyExc_TypeError,
            "The arguments dtypes do not match, expected theta and phfunc to be of the type float[] or double[]!"
        );
        Py_XDECREF(data_arr[0]);
        Py_XDECREF(data_arr[1]);
        return NULL;
    }

    int sizeT = (int) PyArray_SIZE((PyArrayObject *)data_arr[0]);
    int sizeP = (int) PyArray_SIZE((PyArrayObject *)data_arr[1]);
    if(sizeT!=sizeP) {
        PyErr_SetString(
            PyExc_IndexError,
            "The length of the arrays theta and phfunc have to be the same!"
        );
        Py_XDECREF(data_arr[0]);
        Py_XDECREF(data_arr[1]);
        return NULL;
    }

    double* theta = new double[sizeT];
    double* phfunc = new double[sizeP];
    py2c_dblarr((PyArrayObject *)data_arr[0], theta);
    py2c_dblarr((PyArrayObject *)data_arr[1], phfunc);
    double volbsc = calc_volbsc_from_phfunc(valueX, theta, phfunc, sizeT);

    PyObject *res = Py_BuildValue("d", volbsc);
    Py_DECREF(data_arr[0]);
    Py_DECREF(data_arr[1]);
    delete[] theta;
    delete[] phfunc;
    return res;
}
