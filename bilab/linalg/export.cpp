#ifndef EXPORT_H
#define EXPORT_H

#include "common.h"
#include "matrix.h"

#include <boost/python.hpp>

namespace bpy = boost::python;

template<class T>
void export_matrix() {
  bpy::class_< CWmatrix<T> >("CWmatrix", bpy::init<T>())
//      .def(init<float, float, float, std::string>())
    .def("show", &CWmatrix<T>::show);
}

BOOST_PYTHON_MODULE(matrix)
{
  export_matrix<int>();
}

#endif
