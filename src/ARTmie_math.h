#ifndef ARTMIE_MATH_H
#define ARTMIE_MATH_H

#include <Python.h>




// **** Gamma function

#define gm_docstring "gamma(x)\n\n\
Calculates the gamma function\n\n\
Parameters\n----------\n\
x : scalar, floating point number\n    argument\n\n\
Returns\n-------\n\
g : scalar, floating point number\n    result of the gamma function"
PyObject* mie_art_gamma(PyObject *self, PyObject *args, PyObject *kwds);




// **** Bessel functions

#define bj_docstring "besselj(v, z, /, es=False)\n\n\
Calculates the Bessel function of the first kind\n\n\
Parameters\n----------\n\
v : scalar, float\n    order of the Bessel function\n\
z : scalar or array-like, complex\n    the argument/location, where the Bessel function has to be evaluated\n\
es : scalar, boolean, optional\n    exponentially scales the result by exp(2/3*z**1.5) if set to True, default: False\n\n\
Returns\n-------\n\
r : scalar or array-like, complex\n    result of the Bessel function of the first kind and order v at complex value z, same shape as z"
PyObject* mie_art_besselj(PyObject *self, PyObject *args, PyObject *kwds);
#define by_docstring "bessely(v, z, /, es=False)\n\n\
Calculates the Bessel function of the second kind\n\n\
Parameters\n----------\n\
v : scalar, float\n    order of the Bessel function\n\
z : scalar or array-like, complex\n    the argument/location, where the Bessel function has to be evaluated\n\
es : scalar, boolean, optional\n    exponentially scales the result by exp(2/3*z**1.5) if set to True, default: False\n\n\
Returns\n-------\n\
r : scalar or array-like, complex number\n    result of the Bessel function of the second kind and order v at complex value z, same shape as z"
PyObject* mie_art_bessely(PyObject *self, PyObject *args, PyObject *kwds);
#define hv_docstring "hankel(v, z, m, /, es=False)\n\n\
Calculates the Bessel function of the third kind, also known as Hankel function\n\n\
Parameters\n----------\n\
v : scalar, float\n    order of the Bessel function\n\
z : scalar or array-like, complex\n    the argument/location, where the Bessel function has to be evaluated\n\
m : scalar, integer\n    kind of the Hankel function, possible values: 1, 2\n\
es : scalar, boolean, optional\n    exponentially scales the result by exp(2/3*z**1.5) if set to True, default: False\n\n\
Returns\n-------\n\
r : scalar or array-like, complex\n    result of the bessel function of the second kind and order v at complex value z, same shape as z"
PyObject* mie_art_hankel(PyObject *self, PyObject *args, PyObject *kwds);

#define bi_docstring "besseli(v, z, m, /, es=False)\n\n\
Calculates the modified Bessel function of the first kind\n\n\
Parameters\n----------\n\
v : scalar, floating point number\n    order of the Bessel function\n\
z : scalar, complex number\n    the argument/location, where the Bessel function has to be evaluated\n\
es : scalar, boolean, optional\n    exponentially scales the result by exp(2/3*z**1.5) if set to True, default: False\n\n\
Returns\n-------\n\
r : scalar, complex number\n    result of the modified Bessel function of the first kind and order v at complex value z"
PyObject* mie_art_besseli(PyObject *self, PyObject *args, PyObject *kwds);
#define bk_docstring "besselk(v, z, m, /, es=False)\n\n\
Calculates the modified Bessel function of the second kind\n\n\
Parameters\n----------\n\
v : scalar, floating point number\n    order of the Bessel function\n\
z : scalar, complex number\n    the argument/location, where the Bessel function has to be evaluated\n\
es : scalar, boolean, optional\n    exponentially scales the result by exp(2/3*z**1.5) if set to True, default: False\n\n\
Returns\n-------\n\
r : scalar, complex number\n    result of the modified Bessel function of the second kind and order v at complex value z"
PyObject* mie_art_besselk(PyObject *self, PyObject *args, PyObject *kwds);



// **** Airy function

#define ai_docstring "airy(z, /, es=False)\n\n\
Calculates the Airy function\n\n\
Parameters\n----------\n\
z : scalar, complex number\n    the argument/location, where the airy function has to be evaluated\n\
es : scalar, boolean, optional\n    exponentially scales the result by exp(2/3*z**1.5) if set to True, default: False\n\n\
Returns\n-------\n\
r : scalar, complex number\n    result of the Airy function at complex value z"
PyObject* mie_art_airy(PyObject *self, PyObject *args, PyObject *kwds);




#endif /* ARTMIE_MATH_H */