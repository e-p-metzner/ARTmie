#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL ARTmieModule

#include "ARTmie_constants.h"
#include "ARTmie_helper.h"
#include "ARTmie_amos.h"
#include "ARTmie_coeff.h"

void mie_ab(std::complex<double> m, double x, std::complex<double> *an, std::complex<double> *bn) {
    int nmax = calc_nmax(x);
    std::complex<double> mx = m*x;
    int nmx  = 16 + std::max(nmax, (int)(std::abs(mx)+0.5));
    double sx = std::sqrt(HPI*x);
    int idx;

    std::complex<double>* Dn = new std::complex<double>[nmx];
    Dn[nmx-1] = std::complex<double>(0.0,0.0);
    for(idx=nmx-2; idx>=0; idx--) {
        std::complex<double> invM = (2.0+idx) / mx;
        Dn[idx] = invM - 1.0/(Dn[idx+1] + invM);
    }

    std::complex<double>* jv = new std::complex<double>[nmax];
    std::complex<double>* yv = new std::complex<double>[nmax];

    double jvr[1], jvi[1], yvr[1], yvi[1], cwr[1], cwi[1];
    double v_start = 1.5;
    int success, ierr;
    int debug = false;
    for(idx=0; idx<nmax; idx++) {
        double nu = v_start+idx;
        zbesj(x, 0.0, nu, 1, 1, jvr, jvi, &success,           &ierr, &debug);
        zbesy(x, 0.0, nu, 1, 1, yvr, yvi, &success, cwr, cwi, &ierr, &debug);
        jv[idx] = std::complex<double>(jvr[0],jvi[0]);
        yv[idx] = std::complex<double>(yvr[0],yvi[0]);
    }

    std::complex<double> cx(x,0.0);
    std::complex<double> p1x(std::sin(x),0.0);
    std::complex<double> ch1x(std::cos(x),0.0);

    for(idx=0; idx<nmax; idx++) {
        std::complex<double> px   =  sx*jv[idx];
        std::complex<double> chx  = -sx*yv[idx];
        std::complex<double> gsx  = px - chx*cplxJ;
        std::complex<double> gs1x = p1x - ch1x*cplxJ;
        std::complex<double> af   = Dn[idx]/m + (1.0+idx)/x;
        std::complex<double> bf   = m*Dn[idx] + (1.0+idx)/x;

        an[idx] = (px*af-p1x) / (gsx*af-gs1x);
        bn[idx] = (px*bf-p1x) / (gsx*bf-gs1x);

        p1x  = px;
        ch1x = chx;
    }

    delete[] Dn;
    delete[] jv;
    delete[] yv;
}
PyObject* mie_art_mieab(PyObject *self, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = { (char*)"m", (char*)"x", NULL };

    //use function for integer array and integer value
    Py_complex valueNpM;
    double valueX;
    PyObject *res = NULL;
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "Dd", kwlist, &valueNpM, &valueX)) {
        PyErr_SetString(
            PyExc_TypeError,
            "The arguments dtypes do not match, expected (complex,float)"
        );
        return res;
    }
    int nmax = calc_nmax(valueX);
    std::complex<double> valueM = py2c_cplx(valueNpM);
    std::complex<double>* an = new std::complex<double>[nmax];
    std::complex<double>* bn = new std::complex<double>[nmax];
    try {
        mie_ab(valueM, valueX, an, bn);
    } catch(const std::exception &e) {
        delete[] an;
        delete[] bn;
        PyErr_SetString(
            PyExc_ArithmeticError,
            e.what()
        );
        return res;
    }
    try {
        PyArrayObject *pyan = c2py_cplxarr(nmax, an);
        PyArrayObject *pybn = c2py_cplxarr(nmax, bn);
        res = Py_BuildValue("OO", pyan, pybn);
    } catch(const std::exception& e) {
        PyErr_SetString(PyExc_Exception, e.what());
        res = NULL;
    }
    delete[] an;
    delete[] bn;
    return res;
}

void miecoated_ab(std::complex<double> m_core, double x_core, std::complex<double> m_shell, double x_shell, std::complex<double> *an, std::complex<double> *bn) {
    int nmax = calc_nmax(x_shell);
    std::complex<double> cx_shell(x_shell,0.0);
    std::complex<double> m = m_shell / m_core;
    std::complex<double> u = m_core * x_core; // #arguments for the Bessel functions
    std::complex<double> v = m_shell * x_core;
    std::complex<double> w = m_shell * x_shell;
    double mx = std::max(std::abs(m_core)*x_shell, std::abs(m_shell)*x_shell);
    int nmx  = 16 + std::max(nmax, (int)(mx+0.5));
    int idx;

    std::complex<double>* dnu = new std::complex<double>[nmx];
    std::complex<double>* dnv = new std::complex<double>[nmx];
    std::complex<double>* dnw = new std::complex<double>[nmx];
    dnu[nmx-1] = std::complex<double>(0.0,0.0);
    dnv[nmx-1] = std::complex<double>(0.0,0.0);
    dnw[nmx-1] = std::complex<double>(0.0,0.0);
    dnu[nmx-2] = std::complex<double>(0.0,0.0);
    dnv[nmx-2] = std::complex<double>(0.0,0.0);
    dnw[nmx-2] = std::complex<double>(0.0,0.0);
    for(idx=nmx-3; idx>=0; idx--) {
        std::complex<double> invU = (2.0+idx)/u;
        std::complex<double> invV = (2.0+idx)/v;
        std::complex<double> invW = (2.0+idx)/w;
        dnu[idx] = invU - 1.0/(dnu[idx+1]+invU);
        dnv[idx] = invV - 1.0/(dnv[idx+1]+invV);
        dnw[idx] = invW - 1.0/(dnw[idx+1]+invW);
    }

    std::complex<double> sv = std::sqrt(HPI*v);
    std::complex<double> sw = std::sqrt(HPI*w);
    std::complex<double> sy(std::sqrt(HPI*x_shell),0.0);
    std::complex<double> p1y(std::sin(x_shell),0.0);
    std::complex<double> ch1y(std::cos(x_shell),0.0);
    std::complex<double> gs1y(0.0,0.0);

    std::complex<double>* jv = new std::complex<double>[nmax];
    std::complex<double>* yv = new std::complex<double>[nmax];
    std::complex<double>* jw = new std::complex<double>[nmax];
    std::complex<double>* yw = new std::complex<double>[nmax];
    std::complex<double>* jy = new std::complex<double>[nmax];
    std::complex<double>* yy = new std::complex<double>[nmax];

    double jvr[1], jvi[1], yvr[1], yvi[1], cwr[1], cwi[1];
    double v_start = 1.5;
    int success, ierrj, ierry;
    int debug = false;
    for(idx=0; idx<nmax; idx++) {
        double nu = v_start+idx;
        zbesj(v.real(),v.imag(), nu, 2, 1, jvr, jvi, &success,           &ierrj, &debug);
        zbesy(v.real(),v.imag(), nu, 2, 1, yvr, yvi, &success, cwr, cwi, &ierry, &debug);
        jv[idx] = std::complex<double>(jvr[0],jvi[0]);
        yv[idx] = std::complex<double>(yvr[0],yvi[0]);
        zbesj(w.real(),w.imag(), nu, 2, 1, jvr, jvi, &success,           &ierrj, &debug);
        zbesy(w.real(),w.imag(), nu, 2, 1, yvr, yvi, &success, cwr, cwi, &ierry, &debug);
        jw[idx] = std::complex<double>(jvr[0],jvi[0]);
        yw[idx] = std::complex<double>(yvr[0],yvi[0]);
        zbesj(x_shell,0.0, nu, 1, 1, jvr, jvi, &success,           &ierrj, &debug);
        zbesy(x_shell,0.0, nu, 1, 1, yvr, yvi, &success, cwr, cwi, &ierry, &debug);
        jy[idx] = std::complex<double>(jvr[0],jvi[0]);
        yy[idx] = std::complex<double>(yvr[0],yvi[0]);
    }
    double aiv = std::abs(v.imag());
    double aiw = std::abs(w.imag());
    double ewvv = std::exp(aiw-aiv-aiv);
    double ew   = std::exp(aiw);

    for(idx=0; idx<nmax; idx++) {
        std::complex<double> pv   = jv[idx]*sv;
        std::complex<double> pw   = jw[idx]*sw;
        std::complex<double> py   = jy[idx]*sy;
        std::complex<double> chv  = -yv[idx]*sv;
        std::complex<double> chw  = -yw[idx]*sw;
        std::complex<double> chy  = -yy[idx]*sy;
        std::complex<double> gsy  = py - chy*cplxJ;
        std::complex<double> gs1y = p1y - ch1y*cplxJ;

        std::complex<double> uu   = dnu[idx]*m - dnv[idx];
        std::complex<double> vv   = dnu[idx]/m - dnv[idx];
        std::complex<double> fv   = pv/chv;
        std::complex<double> pw_chw_fv = (pw-chw*fv)*ew;
        std::complex<double> pw_pv_chv = pw/(pv*chv)*ewvv;
        std::complex<double> ku1  = (uu*fv/pw)/ew;
        std::complex<double> kv1  = (vv*fv/pw)/ew;
        std::complex<double> ku2  = uu*pw_chw_fv + pw_pv_chv;
        std::complex<double> kv2  = vv*pw_chw_fv + pw_pv_chv;
        std::complex<double> dns  = ku1 / ku2 + dnw[idx];
        std::complex<double> gns  = kv1 / kv2 + dnw[idx];

        std::complex<double> af   = dns/m_shell + (1.0+idx)/x_shell;
        std::complex<double> bf   = gns*m_shell + (1.0+idx)/x_shell;
        an[idx] = (py*af-p1y) / (gsy*af-gs1y);
        bn[idx] = (py*bf-p1y) / (gsy*bf-gs1y);

        p1y = py;
        ch1y = chy;
    }

    delete[] dnu;
    delete[] dnv;
    delete[] dnw;
    delete[] jv;
    delete[] yv;
    delete[] jw;
    delete[] yw;
    delete[] jy;
    delete[] yy;
}
PyObject* mie_art_miecoatedab(PyObject *self, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = { (char*)"m_core", (char*)"x_core", (char*)"m_shell", (char*)"x_shell", NULL };

    //use function for integer array and integer value
    Py_complex valueNpMcore;
    Py_complex valueNpMshell;
    double valueXcore;
    double valueXshell;
    if(PyArg_ParseTupleAndKeywords(args, kwds, "DdDd", kwlist, &valueNpMcore, &valueXcore, &valueNpMshell, &valueXshell)) {
    } else {
        PyErr_SetString(
            PyExc_TypeError,
            "The arguments dtypes do not match, expected (complex,float,complex,float)"
        );
        return NULL;
    }
    int nmax = calc_nmax(valueXshell);
    std::complex<double> valueMcore = py2c_cplx(valueNpMcore);
    std::complex<double> valueMshell = py2c_cplx(valueNpMshell);
    std::complex<double>* an = new std::complex<double>[nmax];
    std::complex<double>* bn = new std::complex<double>[nmax];
    miecoated_ab(valueMcore, valueXcore, valueMshell, valueXshell, an, bn);
    PyArrayObject *pyan = c2py_cplxarr(nmax, an);
    PyArrayObject *pybn = c2py_cplxarr(nmax, bn);
    PyObject *res = Py_BuildValue("OO", pyan, pybn);
    delete[] an;
    delete[] bn;
    return res;
}

void mie_cd(std::complex<double> m, double x, std::complex<double> *cn, std::complex<double> *dn) {
    std::complex<double> cx(x,0.0);
    std::complex<double> mx = m*x;
    std::complex<double> m2 = m*m;
    std::complex<double> mx2 = mx*mx;
    int nmax = calc_nmax(x);
    int nmx  = 16+std::max(nmax, (int)(std::abs(mx)+0.5));
    int idx;

    std::complex<double>* cnx = new std::complex<double>[nmx];
    cnx[nmx-1] = std::complex<double>(0.0,0.0);
    for(idx=nmx; idx>1; idx--) {
        cnx[idx-1] = (double)idx - mx2/(cnx[idx]+(double)idx);
    }

    double rx1 = std::sqrt(HPI/x);
    std::complex<double> rx2 = std::sqrt(mx/HPI);
    std::complex<double> b1x(std::sin(x)/x,0.0);
    std::complex<double> y1x(std::cos(x)/x,0.0);

    std::complex<double>* jv = new std::complex<double>[nmax];
    std::complex<double>* yv = new std::complex<double>[nmax];
    double jvr[1], jvi[1], yvr[1], yvi[1], cwr[1], cwi[1];
    double v_start = 1.5;
    int success, ierr;
    int  debug = false;
    for(idx=0; idx<nmax; idx++) {
        double nu = v_start + idx;
        zbesj(x,0.0, nu, 1, 1, jvr, jvi, &success,           &ierr, &debug);
        zbesy(x,0.0, nu, 1, 1, yvr, yvi, &success, cwr, cwi, &ierr, &debug);
        jv[idx] = std::complex<double>(jvr[0],jvi[0]);
        yv[idx] = std::complex<double>(yvr[0],yvi[0]);
    }
    for(idx=0; idx<nmax; idx++) {
        std::complex<double> jnx  = rx1*jv[idx];
        std::complex<double> jnmx = rx2/jv[idx];
        std::complex<double> yx   = rx1*yv[idx];
        std::complex<double> hx   = jnx + yx*cplxJ;
        std::complex<double> hn1x = b1x + y1x*cplxJ;
        std::complex<double> ax   = b1x*x - jnx*(1.0+idx);
        std::complex<double> ahx  = hn1x*x - hx*(1.0*idx);
        std::complex<double> numerator = jnx*ahx - hx*ax;
        std::complex<double> c_denom =   ahx - hx*cnx[idx];
        std::complex<double> d_denom =   m2*ahx - hx*cnx[idx];

        cn[idx] = jnmx * numerator/c_denom;
        dn[idx] = jnmx*m *numerator/d_denom;

        b1x = jnx;
        y1x = yx;
    }

    delete[] cnx;
    delete[] jv;
    delete[] yv;
}
PyObject* mie_art_miecd(PyObject *self, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = { (char*)"m", (char*)"x", NULL };

    //use function for integer array and integer value
    Py_complex valueNpM;
    double valueX;
    if(PyArg_ParseTupleAndKeywords(args, kwds, "Dd", kwlist, &valueNpM, &valueX)) {
    } else {
        PyErr_SetString(
            PyExc_TypeError,
            "The arguments dtypes do not match, expected (complex,float)"
        );
        return NULL;
    }
    int nmax = calc_nmax(valueX);
    std::complex<double> valueM = py2c_cplx(valueNpM);
    std::complex<double>* cn = new std::complex<double>[nmax];
    std::complex<double>* dn = new std::complex<double>[nmax];
    mie_cd(valueM, valueX, cn, dn);
    PyArrayObject *pycn = c2py_cplxarr(nmax, cn);
    PyArrayObject *pydn = c2py_cplxarr(nmax, dn);
    PyObject *res = Py_BuildValue("OO", pycn, pydn);
    delete[] cn;
    delete[] dn;
    return res;
}

void mie_pitau(double theta, int nmax, double *pin, double *taun) {
    double mu = std::cos(theta);
    pin[0] = 1.0;
    pin[1] = 3*mu;
    taun[0] = mu;
    taun[1] = 3*std::cos(2*theta);
    int idx;
    for(idx=2; idx<nmax; idx++) {
        pin[idx]  = ((2.0*idx+1.0)*mu*pin[idx-1] - (1.0+idx)*pin[idx-2]) / idx;
        taun[idx] = (1.0+idx)*mu*pin[idx] - (2.0+idx)*pin[idx-1];
    }
}
void mie_pitau(int nang, double *theta, int nmax, double *pin, double *taun) {
    int a,n;
    double* sub_pi = new double[nmax];
    double* sub_tau = new double[nmax];
    for(a=0; a<nang; a++) {
        mie_pitau(theta[a], nmax, sub_pi, sub_tau);
        int off = a*nmax;
        for(n=0; n<nmax; n++) {
            pin[off+n]  = sub_pi[n];
            taun[off+n] = sub_tau[n];
        }
    }
    delete[] sub_pi;
    delete[] sub_tau;
}
PyObject* mie_art_miepitau(PyObject *self, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = { (char*)"theta", (char*)"nmax", NULL };

    //use function for integer array and integer value
    double valueTheta;
    int valueNmax;
    if(PyArg_ParseTupleAndKeywords(args, kwds, "di", kwlist, &valueTheta, &valueNmax)) {
        double* pin = new double[valueNmax];
        double* taun = new double[valueNmax];
        mie_pitau(valueTheta, valueNmax, pin, taun);
        PyArrayObject *pypin = c2py_dblarr(valueNmax, pin);
        PyArrayObject *pytaun = c2py_dblarr(valueNmax, taun);
        PyObject *res = Py_BuildValue("OO", pypin, pytaun);
        delete[] pin;
        delete[] taun;
        return res;
    } else {
        PyErr_Clear();
    }

    PyObject* arr_ptr[] = { NULL };
    PyObject* array[]   = { NULL };
    int dtype = -1;
    if(PyArg_ParseTupleAndKeywords(args, kwds, "Oi", kwlist, &arr_ptr[0], &valueNmax)) {
        {
            if(parse_arrays(1, NPY_FLOAT, arr_ptr, array)) {
                dtype = NPY_FLOAT;
            }
        }
        if(dtype==-1) {
            if (parse_arrays(1, NPY_DOUBLE, arr_ptr, array)) {
                dtype = NPY_DOUBLE;
            }
        }
    } else {
        PyErr_Clear();
    }

    if(dtype<0) {
        PyErr_SetString(
            PyExc_TypeError,
            "The arguments dtypes do not match, expected (float,int), (double,int), (float[],int) or (double[],int)"
        );
        return NULL;
    }
    int ndimT = PyArray_NDIM((PyArrayObject *)array[0]);
    if(ndimT!=1) {
        PyErr_SetString( PyExc_TypeError, "First argument must be scalar or one-dimensional!");
        Py_XDECREF(array[0]);
        return NULL;
    }
    int nang = (int) PyArray_SIZE((PyArrayObject *)array[0]);
    double *theta = new double[nang];
    py2c_dblarr((PyArrayObject *)array[0], theta);
    int arr_len = nang*valueNmax;
    double *pin = new double[arr_len];
    double *taun = new double[arr_len];
    mie_pitau(nang, theta, valueNmax, pin, taun);
    PyArrayObject *pypin = c2py_dblarr(nang, valueNmax, pin);
    PyArrayObject *pytaun = c2py_dblarr(nang, valueNmax, taun);
    PyObject *res = Py_BuildValue("OO", pypin, pytaun);
    Py_DECREF(array[0]);
    delete[] theta;
    delete[] pin;
    delete[] taun;
    return res;
}
