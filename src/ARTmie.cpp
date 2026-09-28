#include <vector>
#include <string>

#define PY_ARRAY_UNIQUE_SYMBOL ARTmieModule
#include "ARTmie_helper.h"
#include "ARTmie_math.h"
#include "ARTmie_coeff.h"
#include "ARTmie_single.h"
#include "ARTmie_phase.h"
#include "ARTmie_sdo.h"

// **** Module definition ****

PyMethodDef mie_methods[] = {
    {"gamma",   (PyCFunction)(void(*)(void))mie_art_gamma,   METH_VARARGS|METH_KEYWORDS, gm_docstring},
    {"besselj", (PyCFunction)(void(*)(void))mie_art_besselj, METH_VARARGS|METH_KEYWORDS, bj_docstring},
    {"bessely", (PyCFunction)(void(*)(void))mie_art_bessely, METH_VARARGS|METH_KEYWORDS, by_docstring},
    {"hankel",  (PyCFunction)(void(*)(void))mie_art_hankel,  METH_VARARGS|METH_KEYWORDS, hv_docstring},
    {"besseli", (PyCFunction)(void(*)(void))mie_art_besseli, METH_VARARGS|METH_KEYWORDS, bi_docstring},
    {"besselk", (PyCFunction)(void(*)(void))mie_art_besselk, METH_VARARGS|METH_KEYWORDS, bk_docstring},
    {"airy",    (PyCFunction)(void(*)(void))mie_art_airy,    METH_VARARGS|METH_KEYWORDS, ai_docstring},

    {"MieQ",               (PyCFunction)(void(*)(void))mie_art_mieq,        METH_VARARGS|METH_KEYWORDS, mieq_docstring},
    {"MieCoatedQ",         (PyCFunction)(void(*)(void))mie_art_miecoatedq,  METH_VARARGS|METH_KEYWORDS, miecoatedq_docstring},
    {"Mie_ab",             (PyCFunction)(void(*)(void))mie_art_mieab,       METH_VARARGS|METH_KEYWORDS, mieab_docstring},
    {"MieCoated_ab",       (PyCFunction)(void(*)(void))mie_art_miecoatedab, METH_VARARGS|METH_KEYWORDS, miecoatedab_docstring},
    {"Mie_cd",             (PyCFunction)(void(*)(void))mie_art_miecd,       METH_VARARGS|METH_KEYWORDS, miecd_docstring},
    {"Mie_pitau",          (PyCFunction)(void(*)(void))mie_art_miepitau,    METH_VARARGS|METH_KEYWORDS, miepitau_docstring},
    {"ab2mie",             (PyCFunction)(void(*)(void))mie_art_ab2mie,      METH_VARARGS|METH_KEYWORDS, abtomie_docstring},
    {"ScatteringFunction", (PyCFunction)(void(*)(void))mie_art_scatfunc,    METH_VARARGS|METH_KEYWORDS, scatfunc_docstring},

    {"createLogNormalDistribution",      (PyCFunction)(void(*)(void))mie_art_createLgNormDist,          METH_VARARGS|METH_KEYWORDS, clnd_docstring},
    {"calcVolBackscattering",            (PyCFunction)(void(*)(void))mie_art_calcVolBackScat,           METH_VARARGS|METH_KEYWORDS, cbs_docstring},
    {"calcVolBackscatteringFromPhFunc",  (PyCFunction)(void(*)(void))mie_art_calcVolBackScatFromPhFunc, METH_VARARGS|METH_KEYWORDS, cbs_phfunc_docstring},
    {"Size_Distribution_Optics",         (PyCFunction)(void(*)(void))mie_art_sdo,                       METH_VARARGS|METH_KEYWORDS, sdo_docstring},
    {"Size_Distribution_Phase_Function", (PyCFunction)(void(*)(void))mie_art_sdpf,                      METH_VARARGS|METH_KEYWORDS, sdpf_docstring},

    {NULL, NULL, 0, NULL} /* sentinel */
};

PyModuleDef artmiemodule = { // @suppress("Invalid arguments")
    PyModuleDef_HEAD_INIT,
    "ARTmie",
    NULL,
    -1,
    mie_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit_ARTmie(void) {
    
    PyObject *m;
    
    import_array();
    
    m = PyModule_Create(&artmiemodule);
    if (!m) {
        return NULL;
    }
    
    return m;
}
