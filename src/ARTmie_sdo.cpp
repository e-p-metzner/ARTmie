#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL ARTmieModule

#include "ARTmie_constants.h"
#include "ARTmie_helper.h"
#include "ARTmie_coeff.h"
#include "ARTmie_single.h"
#include "ARTmie_phase.h"
#include "ARTmie_sdo.h"

void createLogNormalDistribution(double d_gn, double sigma_g, double fcoat, double res, double /*out*/ *x_range, double /*out*/ *y_range, double /*out*/ *pdf) {
    //Helper function: createLogNormalDistribution
    //Input
    //  d_gn, sigma_g = count median diameter, geometric standart deviation
    //  fcoat = coating fraction for coated particles
    //  res = resolution of the particle size distribution
    //  dens = particle density in kg/m3 (1g/cm3 = (0.001kg / 0.000001m3) = 1000kg/m3 ?)
    int nbin = calc_sdo_bin_count(sigma_g, res);
    double loggsd = std::log(sigma_g);
    double lgg2 = 2.0*loggsd*loggsd;

    double dlogd = std::min(res*LN10, std::log(sigma_g)*0.25);
    double cfp1 = 1.0+fcoat;
    double gaussnorm = dlogd / (loggsd * std::sqrt(2.0*_PI_));

    int idx;
    int hcnt = nbin/2;
    for(idx=0; idx<nbin; idx++) {
        double lnx = (idx-hcnt)*dlogd;
        x_range[idx] = d_gn * std::exp(lnx);
        y_range[idx] = cfp1*x_range[idx];
        pdf[idx] =  std::exp(-lnx*lnx/lgg2) * gaussnorm;
    }
}
PyObject* mie_art_createLgNormDist(PyObject *self, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = { (char*)"mean_diam", (char*)"geom_std", (char*)"fcoat", (char*)"res", NULL };

    double valueMu;
    double valueStd;
    double valueFcoat = 0.0;
    double valueRes =   1.0;
    int valueN2core =   false;
    int valueN2vol =    true;
    if(!PyArg_ParseTupleAndKeywords(args, kwds, "dd|dd", kwlist, &valueMu, &valueStd, &valueFcoat, &valueRes)) {
        PyErr_SetString(
            PyExc_TypeError,
            "The arguments dtypes do not match, expected (float,float[,float,float])"
        );
        return NULL;
    }

    int nbin = calc_sdo_bin_count(valueStd, valueRes);
    double* core_diams  = new double[nbin];
    double* shell_diams = new double[nbin];
    double* pdf         = new double[nbin];
    createLogNormalDistribution(valueMu, valueStd, valueFcoat, valueRes, core_diams, shell_diams, pdf);

    // PyObject *res = Py_BuildValue("OOO",
    //                               c2py_dblarr(nbin, core_diams),
    //                               c2py_dblarr(nbin, shell_diams),
    //                               c2py_dblarr(nbin, pdf));

    PyArrayObject* arr = c2py_dblarr(nbin, pdf);

    PyObject *res = Py_BuildValue("ddd", 1.0, 2.3, 6.7);

    delete[] core_diams;
    delete[] shell_diams;
    delete[] pdf;

    return res;
}

void size_distribution_optics(std::complex<double> m_core, double mean_diam, double geom_std, double wavelength, std::complex<double> m_shell, double fcoating, double resolution, double dens, int effcore, int msc, int debug, Mie_tots *mie_tots) {
    if(debug) {
        PySys_WriteStdout("[DEBUG] SDO-Input: m_core=%.7f+i*%14.7e  m_shell=%.7f+i*%14.7e  psd=%.1f+/-%4f  fcoat=%.4f  wl=%14.7e\n",
            m_core.real(),m_core.imag(),m_shell.real(),m_shell.imag(),mean_diam,geom_std,fcoating,wavelength);
        PySys_WriteStdout("[DEBUG]            resolution=%.7f  density=%.7f  effcore=%s  msc=%s\n",
            resolution,dens,effcore?"true":"false",msc?"true":"false");
    }
    double res = 1.0/resolution;
    int nbin = calc_sdo_bin_count(geom_std, res);
    double* core_diams  = new double[nbin];
    double* shell_diams = new double[nbin];
    double* pdf         = new double[nbin];
    createLogNormalDistribution(mean_diam, geom_std, fcoating, res, core_diams, shell_diams, pdf);
    double area_factor =   1.0e-18 /*nm2->m2*/ * _PI_ / 4.0;
    double volume_factor = 1.0e-27 /*nm3->m3*/ * _PI_ / 6.0;
    double mean_area = 0.0;
    double mean_volume = 0.0;
    double max_shell_diam = 0.0;
    int idx;
    for(idx=0; idx<nbin; idx++) {
        if(max_shell_diam < shell_diams[idx])
            max_shell_diam = shell_diams[idx];
        double area = area_factor * shell_diams[idx] * shell_diams[idx] * pdf[idx];
        mean_area += area;
        double diam = effcore ? core_diams[idx] : shell_diams[idx];
        double volume = volume_factor * diam * diam * diam * pdf[idx];
        mean_volume += volume;
        //optical properties of size distributions are weighted by cross sections
        pdf[idx] = area;
    }
    double normWeight = 1.0 / mean_area;
    if (msc) {
        double mean_mass = mean_volume * 1000.0 /*kg -> g*/ * dens;
        normWeight = 1.0 / mean_mass;
    }

    double angres = 0.25; //degrees
    int nang = calc_angles_count(angres);
    double* theta    = new double[nang];
    double* dtheta   = new double[nang];
    double* scatwgts = new double[nang];
    scattering_weights(angres, theta, dtheta, scatwgts);

    double maxy = _PI_*max_shell_diam/wavelength;
    int nmax    = calc_nmax(maxy);
    int arr_len = nang*nmax;
    double *pin  = new double[arr_len];
    double *taun = new double[arr_len];
    mie_pitau(nang, theta, nmax, pin, taun);
    std::complex<double> *an = new std::complex<double>[nmax];
    std::complex<double> *bn = new std::complex<double>[nmax];

    int is_coated = (fcoating>EPS) && (m_core!=m_shell);

    mie_tots->bext = 0.0;
    mie_tots->bsca = 0.0;
    mie_tots->babs = 0.0;
    mie_tots->bback = 0.0;
    mie_tots->bssa = 0.0;
    mie_tots->basym = 0.0;
    mie_tots->bratio = 0.0;

    for(idx=0; idx<nbin; idx++) {
        try {
            double xval = _PI_*core_diams[idx]/wavelength;
            double yval = _PI_*shell_diams[idx]/wavelength;
            int nmax2 = calc_nmax(yval);
            if (is_coated) {
                miecoated_ab(m_core, xval, m_shell, yval, an, bn);
            } else {
                mie_ab(m_core, yval, an, bn);
            }

            MieResult one_result = ab2mie(nmax2, an, bn, wavelength, shell_diams[idx], 0);
            double backscat = calcVolBackscattering(yval, nmax2, an, bn, nang, theta, dtheta, scatwgts, nmax, pin, taun);
            if(debug) {
                PySys_WriteStdout("[DEBUG]    (cdia=%.2f) -> one_res{ qext=%.6f, qsca=%.6f, qabs=%.6f, qratio=%.6f, qg=%.6f, ... } + volbsc=%.6e\n",
                                  core_diams[idx], one_result.qext,one_result.qsca,one_result.qabs,one_result.qratio,one_result.qg,backscat);
            }

            mie_tots->bext   += one_result.qext * pdf[idx];
            mie_tots->bsca   += one_result.qsca * pdf[idx];
            mie_tots->babs   += one_result.qabs * pdf[idx];
            mie_tots->bback  += one_result.qback * pdf[idx];
            mie_tots->bvolbsc += backscat * pdf[idx];
            mie_tots->basym  += one_result.qg * one_result.qsca * pdf[idx];
        } catch(const std::exception& e) {
            PySys_WriteStdout("[exception] %s\n", e.what());
        }
    }
    if(debug) {
        PySys_WriteStdout("[DEBUG]    ext=%.6e, sca=%.6e, abs=%.6e, asy=%.6e, normWeight=%.6e\n",
                          mie_tots->bext, mie_tots->bsca, mie_tots->babs, mie_tots->basym, normWeight);
    }
//    mie_tots->arr_len = dcount;
    mie_tots->basym   /= mie_tots->bsca;
    mie_tots->bratio   = mie_tots->bback / mie_tots->bsca;
    mie_tots->bssa     = mie_tots->bsca / mie_tots->bext;
    mie_tots->bext    *= normWeight;
    mie_tots->bsca    *= normWeight;
    mie_tots->babs    *= normWeight;
    mie_tots->bback   *= normWeight;
    mie_tots->bvolbsc *= normWeight;

    delete[] core_diams;
    delete[] shell_diams;
    delete[] pdf;
    delete[] theta;
    delete[] dtheta;
    delete[] scatwgts;
    delete[] pin;
    delete[] taun;
    delete[] an;
    delete[] bn;
}
PyObject* mie_art_sdo(PyObject *self, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = { (char*)"m", (char*)"sizepar1", (char*)"sizepar2", (char*)"wavelength", (char*)"nMedium", (char*)"fcoat", (char*)"mc", (char*)"density", (char*)"resolution", (char*)"effcore", (char*)"msc", (char*)"debug", NULL };

    Py_complex valueNpMcore;
    Py_complex valueNpMshell = nanPyCplx();
    double valueDmu;
    double valueDstd;
    double valueW;
    double valueNmedium;
    double valueFcoat;
    double valueDens = 1.0;
    double valueRes = 10.0;
    int valueN2core = true;
    int valueAsMSC = true;
    int valueDebug = false;
    Mie_tots *mie_tots = new Mie_tots();
//    int array_sizepar = 0;
    if(PyArg_ParseTupleAndKeywords(args, kwds, "Dddd|ddDddppp", kwlist, &valueNpMcore, &valueDmu, &valueDstd, &valueW, &valueNmedium, &valueFcoat, &valueNpMshell, &valueDens, &valueRes, &valueN2core, &valueAsMSC, &valueDebug)) {
        std::complex<double> valueMcore  = py2c_cplx(valueNpMcore);
        std::complex<double> valueMshell = py2c_cplx(valueNpMshell);
        if(std::isnan(valueNpMshell.real) || std::isnan(valueNpMshell.imag)) {
            valueMshell = py2c_cplx(valueNpMcore);
        }
        size_distribution_optics(valueMcore, valueDmu, valueDstd, valueW, valueMshell, valueFcoat, valueRes, valueDens, valueN2core, valueAsMSC, valueDebug, mie_tots);
    } else {
//        array_sizepar = 1;
//        PyErr_Clear();
        return NULL;
    }

//    if(array_sizepar) {
//        PyObject *arr_ptr[] = { NULL, NULL };
//        PyObject *array[] =   { NULL, NULL };
//        int dtype = -1;
//        if(PyArg_ParseTupleAndKeywords(args, kwds, "DOOd|ddDddpp", kwlist, &valueNpMcore, &arr_ptr[0], &arr_ptr[1], &valueW, &valueNmedium, &valueFcoat, &valueNpMshell, &valueDens, &valueRes, &valueN2core, &valueAsCrossSec)) {
//            {
//                if(parse_arrays(2, NPY_FLOAT, arr_ptr, array)) dtype = NPY_FLOAT;
//            }
//            if(dtype<0) {
//                if(parse_arrays(2, NPY_DOUBLE, arr_ptr, array)) dtype = NPY_DOUBLE;
//            }
//        } else {
//            PyErr_Clear();
//        }
//
//        int ndim1 = PyArray_NDIM((PyArrayObject *)array[0]);
//        int ndim2 = PyArray_NDIM((PyArrayObject *)array[1]);
//        if(ndim1!=1 || ndim2!=1) {
//            PyErr_SetString(
//                PyExc_TypeError,
//                "sizepar1 and sizepar2 have to be both 1dimensional arrays of type float[] or double[]"
//            );
//            Py_XDECREF(array[0]);
//            Py_XDECREF(array[1]);
//            return NULL;
//        }
//        int size1 = PyArray_SIZE((PyArrayObject *)array[0]);
//        int size2 = PyArray_SIZE((PyArrayObject *)array[1]);
//        if(size1!=size2) {
//            PyErr_SetString(
//                PyExc_IndexError,
//                "sizepar1 and sizepar2 have to be arrays of the same length"
//            );
//            Py_XDECREF(array[0]);
//            Py_XDECREF(array[1]);
//            return NULL;
//        }
//
////        std::complex<double> valueMcore  = py2c_cplx(valueNpMcore);
////        std::complex<double> valueMshell = py2c_cplx(valueNpMshell);
////        if(std::isnan(valueNpMshell.real) || std::isnan(valueNpMshell.imag)) {
////            valueMshell = py2c_cplx(valueNpMcore);
////        }
//
//        PyErr_SetString(
//            PyExc_NotImplementedError,
//            "array-like input for sizepar1 and sizepar2 not implemented yet."
//        );
//        Py_XDECREF(array[0]);
//        Py_XDECREF(array[1]);
//        return NULL;
//    }

//    PyObject *res = Py_BuildValue("{s:d,s:d,s:d,s:d,s:d,s:d,s:O,s:O,s:O,s:O,s:O}",
//        "Extinction", mie_tots.bext,
//        "Scattering", mie_tots.bsca,
//        "Absorption", mie_tots.babs,
//        "Backscattering", mie_tots.bback,
//        "SSA", mie_tots.bssa,
//        "Asymmetry", mie_tots.basym,
//        "Extinction Coefficients", c2py_dblarr(mie_tots.arr_len, mie_tots.ext_arr),
//        "Scattering Coefficients", c2py_dblarr(mie_tots.arr_len, mie_tots.sca_arr),
//        "Absorption Coefficients", c2py_dblarr(mie_tots.arr_len, mie_tots.abs_arr),
//        "Backscattering Coefficients", c2py_dblarr(mie_tots.arr_len, mie_tots.bck_arr),
//        "Asymmetry Coefficients", c2py_dblarr(mie_tots.arr_len, mie_tots.g_arr)
//    );
    PyObject *res = Py_BuildValue("{s:d,s:d,s:d,s:d,s:d,s:d,s:d,s:d}",
        "Extinction", 0.0+mie_tots->bext,
        "Scattering", 0.0+mie_tots->bsca,
        "Absorption", 0.0+mie_tots->babs,
        "Backscattering", 0.0+mie_tots->bback,
        "AngWgtBackscat", 0.0+mie_tots->bvolbsc,
        "SSA", 0.0+mie_tots->bssa,
        "BackscatterRatio", 0.0+mie_tots->bratio,
        "Asymmetry", 0.0+mie_tots->basym
    );

//    delete[] mie_tots.ext_arr;
//    delete[] mie_tots.sca_arr;
//    delete[] mie_tots.abs_arr;
//    delete[] mie_tots.bck_arr;
//    delete[] mie_tots.g_arr;
    delete mie_tots;
    return res;
}

void size_distribution_phase_function(std::complex<double> m_core, double mean_diam, double geom_std, double wavelength, std::complex<double> m_shell, double fcoating, double resolution, double dens, int effcore, int norm2vol, int nang, double *outTheta, double *outSL, double *outSR, double *outSU) {
    double res = 1.0/resolution;
    int nbin = calc_sdo_bin_count(geom_std, res);
    double* core_diams  = new double[nbin];
    double* shell_diams = new double[nbin];
    double* pdf         = new double[nbin];
    double* crossArea   = new double[nbin];
    createLogNormalDistribution(mean_diam, geom_std, fcoating, res, core_diams, shell_diams, pdf);
    double max_shell_diam = 0.0;
    double total_area = 0.0;
    int idx, a;
    for(idx=0; idx<nbin; idx++) {
        if(max_shell_diam < shell_diams[idx])
            max_shell_diam = shell_diams[idx];
        // phase function has to be weighted by cross-section analog
        pdf[idx] *= 0.25 * _PI_ * shell_diams[idx]*shell_diams[idx];
        total_area += pdf[idx];
    }
    double normWeight = 1.0/total_area;

    double* sl = new double[nang];
    double* sr = new double[nang];
    double* su = new double[nang];

    double maxy = _PI_*max_shell_diam/wavelength;
    int nmax = calc_nmax(maxy);
    //std::complex<double> maxmy = m_core*maxy;
    //PySys_WriteStdout("%12.6f+i*%12.6f\n",maxmy.real(),maxmy.imag());
    int arr_len = nang*nmax;
    double *pin = new double[arr_len];
    double *taun = new double[arr_len];
    mie_pitau(nang, outTheta, nmax, pin, taun);
    std::complex<double> *an = new std::complex<double>[nmax];
    std::complex<double> *bn = new std::complex<double>[nmax];

    int is_coated = (fcoating>EPS) && (m_core!=m_shell);

    for(a=0; a<nang; a++) {
        outSL[a] = 0.0;
        outSR[a] = 0.0;
        outSU[a] = 0.0;
    }

    for(idx=0; idx<nbin; idx++) {
        double xval = _PI_*core_diams[idx]/wavelength;
        double yval = _PI_*shell_diams[idx]/wavelength;
        int anbn_used_len = calc_nmax(yval);
        if (is_coated) {
            miecoated_ab(m_core, xval, m_shell, yval, an, bn);
        } else {
            mie_ab(m_core, yval, an, bn);
        }

        scattering_function(anbn_used_len, an, bn, nang, outTheta, nmax, pin, taun, sl, sr, su);

        for(a=0; a<nang; a++) {
            outSL[a] += sl[a]*pdf[idx];
            outSR[a] += sr[a]*pdf[idx];
            outSU[a] += su[a]*pdf[idx];
        }
    }

    for(a=0; a<nang; a++) {
        outSL[a] *= normWeight;
        outSR[a] *= normWeight;
        outSU[a] *= normWeight;
    }

    delete[] core_diams;
    delete[] shell_diams;
    delete[] pdf;
    delete[] sl;
    delete[] sr;
    delete[] su;
    delete[] pin;
    delete[] taun;
    delete[] an;
    delete[] bn;
}
PyObject* mie_art_sdpf(PyObject *self, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = { (char*)"m", (char*)"sizepar1", (char*)"sizepar2", (char*)"wavelength", (char*)"nMedium", (char*)"fcoat", (char*)"mc", (char*)"density", (char*)"resolution", (char*)"effcore", (char*)"normalized", NULL };

    Py_complex valueNpMcore;
    Py_complex valueNpMshell = nanPyCplx();
    double valueDmu;
    double valueDstd;
    double valueW;
    double valueNmedium;
    double valueFcoat;
    double valueDens = 1.0;
    double valueRes = 10.0;
    int valueN2core = true;
    int valueAsCrossSec = true;
    PyObject *res = NULL;
    //    int array_sizepar = 0;

    if(PyArg_ParseTupleAndKeywords(args, kwds, "Dddd|ddDddpp", kwlist, &valueNpMcore, &valueDmu, &valueDstd, &valueW, &valueNmedium, &valueFcoat, &valueNpMshell, &valueDens, &valueRes, &valueN2core, &valueAsCrossSec)) {
        std::complex<double> valueMcore  = py2c_cplx(valueNpMcore);
        std::complex<double> valueMshell = py2c_cplx(valueNpMshell);
        if(std::isnan(valueNpMshell.real) || std::isnan(valueNpMshell.imag)) {
            valueMshell = py2c_cplx(valueNpMcore);
        }
        double angres = 0.25; //degrees
        int nang = calc_angles_count(angres);
        double* theta = new double[nang];
        double* pf_sl = new double[nang];
        double* pf_sr = new double[nang];
        double* pf_su = new double[nang];
        double ares = angres*_PI_/180.0;
        for(int a=0; a<nang; a++) {
            theta[a] = a*ares;
        }
        size_distribution_phase_function(valueMcore, valueDmu, valueDstd, valueW, valueMshell, valueFcoat, valueRes, valueDens, valueN2core, valueAsCrossSec, nang, theta, pf_sl, pf_sr, pf_su);
        res = Py_BuildValue("OOOO",
                c2py_dblarr(nang, theta),
                c2py_dblarr(nang, pf_sl),
                c2py_dblarr(nang, pf_sr),
                c2py_dblarr(nang, pf_su)
            );
        delete[] theta;
        delete[] pf_sl;
        delete[] pf_sr;
        delete[] pf_su;
    } else {
//        array_sizepar = 1;
//        PyErr_Clear();
    }

//    if(array_sizepar) {
//        PyObject *arr_ptr[] = { NULL, NULL };
//        PyObject *array[] =   { NULL, NULL };
//        int dtype = -1;
//        if(PyArg_ParseTupleAndKeywords(args, kwds, "DOOd|ddDddpp", kwlist, &valueNpMcore, &arr_ptr[0], &arr_ptr[1], &valueW, &valueNmedium, &valueFcoat, &valueNpMshell, &valueDens, &valueRes, &valueN2core, &valueAsCrossSec)) {
//            {
//                if(parse_arrays(2, NPY_FLOAT, arr_ptr, array)) dtype = NPY_FLOAT;
//            }
//            if(dtype<0) {
//                if(parse_arrays(2, NPY_DOUBLE, arr_ptr, array)) dtype = NPY_DOUBLE;
//            }
//        } else {
//            PyErr_Clear();
//        }
//
//        int ndim1 = PyArray_NDIM((PyArrayObject *)array[0]);
//        int ndim2 = PyArray_NDIM((PyArrayObject *)array[1]);
//        if(ndim1!=1 || ndim2!=1) {
//            PyErr_SetString(
//                PyExc_TypeError,
//                "sizepar1 and sizepar2 have to be both 1dimensional arrays of type float[] or double[]"
//            );
//            Py_XDECREF(array[0]);
//            Py_XDECREF(array[1]);
//            return NULL;
//        }
//        int size1 = PyArray_SIZE((PyArrayObject *)array[0]);
//        int size2 = PyArray_SIZE((PyArrayObject *)array[1]);
//        if(size1!=size2) {
//            PyErr_SetString(
//                PyExc_IndexError,
//                "sizepar1 and sizepar2 have to be arrays of the same length"
//            );
//            Py_XDECREF(array[0]);
//            Py_XDECREF(array[1]);
//            return NULL;
//        }
//
////        std::complex<double> valueMcore  = py2c_cplx(valueNpMcore);
////        std::complex<double> valueMshell = py2c_cplx(valueNpMshell);
////        if(std::isnan(valueNpMshell.real) || std::isnan(valueNpMshell.imag)) {
////            valueMshell = py2c_cplx(valueNpMcore);
////        }
//
//        PyErr_SetString(
//            PyExc_NotImplementedError,
//            "array-like input for sizepar1 and sizepar2 not implemented yet."
//        );
//        Py_XDECREF(array[0]);
//        Py_XDECREF(array[1]);
//        return NULL;
//    }

    return res;
}
