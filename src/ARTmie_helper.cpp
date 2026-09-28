#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL ARTmieModule

#include "ARTmie_constants.h"
#include "ARTmie_helper.h"

int parse_arrays(int arg_count, int dtype, PyObject **args, PyObject **arrs) {
    int i, e=0;

    for (i=0; i<arg_count; i++) {
        arrs[i] = PyArray_FROM_OTF(args[i], dtype, NPY_ARRAY_IN_ARRAY);
        if(arrs[i] == NULL) e++;
    }

    if( e>0 ) {
        for (i=0; i<arg_count; i++) {
            Py_XDECREF(arrs[i]);
        }
        PyErr_Clear();
        return 0; //some error occured
    }

    return 1; //no error
}

PyArrayObject* c2py_dblarr(int c_arr_len, double *c_arr) {
    const npy_intp dims[1] = { (npy_intp)c_arr_len };
    PyArrayObject *pyarr = (PyArrayObject *) PyArray_Zeros(1, dims, PyArray_DescrFromType(NPY_DOUBLE), 0);
    double* pydata = (double *) PyArray_DATA(pyarr);
    for(int idx=0; idx<c_arr_len; idx++) {
        pydata[idx] = c_arr[idx];
    }
    Py_INCREF(pyarr);
    return pyarr;
}
PyArrayObject* c2py_dblarr(int dim1_len, int dim2_len, double *c_arr) {
    npy_intp dims[2];
    dims[0] = (npy_intp)dim1_len;
    dims[1] = (npy_intp)dim2_len;
    PyArrayObject *pyarr = (PyArrayObject *) PyArray_Zeros(2, dims, PyArray_DescrFromType(NPY_DOUBLE), 0);

    PyArrayObject *op[1] = { pyarr };
    npy_uint32 op_flags[1] = { NPY_ITER_WRITEONLY };
    npy_uint32 flags = NPY_ITER_EXTERNAL_LOOP | NPY_ITER_BUFFERED | NPY_ITER_GROWINNER;
    NpyIter *iter = NpyIter_MultiNew(1, op, flags, NPY_KEEPORDER, NPY_NO_CASTING, op_flags, NULL);
    if( iter==NULL ) {
        Py_XDECREF(pyarr);
        return NULL;
    }
    NpyIter_IterNextFunc * iternext = NpyIter_GetIterNext(iter, NULL);
    if (iternext == NULL) {
        NpyIter_Deallocate(iter);
        Py_XDECREF(pyarr);
        return NULL;
    }

    // -- iterate ------------------
    npy_intp count;
    int idx=0;
    char ** dataptr = NpyIter_GetDataPtrArray(iter);
    npy_intp * strideptr = NpyIter_GetInnerStrideArray(iter);
    npy_intp * innersizeptr = NpyIter_GetInnerLoopSizePtr(iter);
    do {
        count = *innersizeptr;

        while(count--) {
            *(double *)dataptr[0] = c_arr[idx];

            dataptr[0] += strideptr[0];
            idx++;
        }
    } while (iternext(iter));

    if(NpyIter_Deallocate(iter) != NPY_SUCCEED) {
        Py_XDECREF(pyarr);
        return NULL;
    }

//    pydata = (double *) PyArray_DATA(pyarr);
//    for(int idx=0; idx<dim1_len*dim2_len; idx++) {
////        int j = idx/dim2_len;
////        int i = idx%dim2_len;
////        pydata[idx] = c_arr[j][i];
//        pydata[idx] = c_arr[idx];
//    }
//    Py_INCREF(pyarr);
    return pyarr;
}
void py2c_dblarr(PyArrayObject* pyarr, double *carr) {
    PyArrayObject* inarr = nullptr;
    if(PyArray_IS_C_CONTIGUOUS(pyarr)) {
        inarr = pyarr;
    } else {
        inarr = (PyArrayObject*)PyArray_NewCopy(pyarr, NPY_CORDER);
    }
    //inarr is C-contiguous, so we use the flat representation of the data
    int arr_len = (int) PyArray_SIZE(inarr);
    int dtype = PyArray_DTYPE(inarr)->type_num;
    int idx;
    if (dtype==NPY_FLOAT) {
        float *pydata = (float *) PyArray_DATA(inarr);
        for(idx=0; idx<arr_len; idx++)
            carr[idx] = (double)pydata[idx];
        return;
    }
    if (dtype==NPY_DOUBLE) {
        double *pydata = (double *) PyArray_DATA(inarr);
        for(idx=0; idx<arr_len; idx++)
            carr[idx] = pydata[idx];
        return;
    }
    for(idx=0; idx<arr_len; idx++)
        carr[idx] = std::numeric_limits<double>::quiet_NaN();
}

std::complex<double> py2c_cplx(Py_complex np_cplx) {
    std::complex<double> c_cplx(np_cplx.real, np_cplx.imag);
    return c_cplx;
}
Py_complex c2py_cplx(std::complex<double> c_cplx) {
    Py_complex np_cplx;
    np_cplx.real = c_cplx.real();
    np_cplx.imag = c_cplx.imag();
    return np_cplx;
}
Py_complex nanPyCplx() {
    Py_complex py_cplx;
    py_cplx.real = std::numeric_limits<double>::quiet_NaN();
    py_cplx.imag = std::numeric_limits<double>::quiet_NaN();
    return py_cplx;
}
PyArrayObject* c2py_cplxarr(int c_arr_len, std::complex<double> *c_arr) {
    npy_intp dims[1];
    dims[0] =  (npy_intp)c_arr_len;
    PyArrayObject *pyarr = (PyArrayObject *) PyArray_Zeros(1, dims, PyArray_DescrFromType(NPY_COMPLEX128), 0);
    Py_complex *pydata = (Py_complex *) PyArray_DATA(pyarr);
    for(int idx=0; idx<c_arr_len; idx++) {
        pydata[idx].real = c_arr[idx].real();
        pydata[idx].imag = c_arr[idx].imag();
    }
    Py_INCREF(pyarr);
    return pyarr;
}
int py2c_cplxarr(PyArrayObject* pyarr, std::complex<double> *carr) {
    PyArrayObject* inarr = nullptr;
    if(PyArray_IS_C_CONTIGUOUS(pyarr)) {
        inarr = pyarr;
    } else {
        inarr = (PyArrayObject*)PyArray_NewCopy(pyarr, NPY_CORDER);
    }
//    int dtype = PyArray_DTYPE(inarr)->type_num;
//    //if(!(dtype==NPY_FLOAT || dtype==NPY_DOUBLE || dtype==NPY_COMPLEX64 || dtype==NPY_COMPLEX128)) {
//    //    inarr = (PyArrayObject*)PyArray_Cast(inarr, NPY_COMPLEX128);
//    //    return 0;
//    //}
    int arr_len = (int) PyArray_SIZE(inarr);
    int idx;
//    if(dtype==NPY_FLOAT) {
//        float* pydata = (float *) PyArray_DATA(inarr);
//        for(idx=0; idx<arr_len; idx++)
//            carr[idx] = std::complex<double>((double)pydata[idx], 0.0);
//        return 1;
//    }
//    if(dtype==NPY_DOUBLE) {
//        double* pydata = (double *) PyArray_DATA(inarr);
//        for(idx=0; idx<arr_len; idx++)
//            carr[idx] = std::complex<double>(pydata[idx], 0.0);
//        return 1;
//    }
//    if(dtype==NPY_COMPLEX64 || dtype==NPY_COMPLEX128) {
        Py_complex *pydata = (Py_complex *) PyArray_DATA(inarr);
        for(idx=0; idx<arr_len; idx++)
            carr[idx] = std::complex<double>(pydata[idx].real, pydata[idx].imag);
        return 1;
//    }
//    for(idx=0; idx<arr_len; idx++)
//        carr[idx] = std::complex<double>(
//                std::numeric_limits<double>::quiet_NaN(),
//                std::numeric_limits<double>::quiet_NaN());
//    return 1;
}

const char* shape2str(int ndim, npy_intp* shape) {
    std::string s = "(";
    for(int i=0; i<ndim; i++) {
        if(i>0)
            s += ",";
        s += std::to_string((long)shape[i]);
    }
    return (s+")").c_str();
}




int calc_nmax(double x) {
    return (int)(x + 4.0*std::cbrt(x) + 2.5);
}
int calc_angles_count(double angres) {
    return 1 + (int) (180.0 / angres + EPS);
}
int calc_sdo_bin_count(double geom_std, double res) {
    double dlogd = std::min(res*LN10, std::log(geom_std)*0.25);
    int half_bin_count = (int)(4.0*std::log(geom_std)/dlogd);
    return 2*half_bin_count + 1; // +1 for the bin containing the median
}
