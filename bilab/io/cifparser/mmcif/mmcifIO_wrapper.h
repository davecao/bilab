//
//  mmcifIO_wrapper.h
// 
//
//  Created by 曹巍 on 2018/06/19.
//  Copyright © 2018年 巍 曹. All rights reserved.
//
#ifndef MMCIFIO_WRAPPER_H
#define MMCIFIO_WRAPPER_H

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <stdio.h>
#include <iostream>

/* Must define Py_TYPE for Python 2.5 or older */
#ifndef Py_TYPE
#  define Py_TYPE(o) ((o)->ob_type)
#endif

/* Must define PyVarObject_HEAD_INIT for Python 2.5 or older */
#ifndef PyVarObject_HEAD_INIT
#define PyVarObject_HEAD_INIT(type, size)       \
        PyObject_HEAD_INIT(type) size,
#endif

#include "mmcif.h"


namespace mmcif = bilab::mmcif;

bool reader_mmcif(std::string& filename, mmcif::CIFDocument& doc);

#endif  //! MMCIFIO_WRAPPER_H