#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include "tsne.h"

/* Must define Py_TYPE for Python 2.5 or older */
#ifndef Py_TYPE
#  define Py_TYPE(o) ((o)->ob_type)
#endif

/* Must define PyVarObject_HEAD_INIT for Python 2.5 or older */
#ifndef PyVarObject_HEAD_INIT
#define PyVarObject_HEAD_INIT(type, size)       \
        PyObject_HEAD_INIT(type) size,
#endif


