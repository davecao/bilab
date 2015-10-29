#ifndef EXPORT_H
#define EXPORT_H

#include "CIsoSurface.h"
#include <boost/python.hpp>
#include <boost/python/class.hpp>
#include <boost/python/raw_function.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/module.hpp>

namespace bpy = boost::python;

template<class T>
void export_cisosurface() {
  bpy::class_< CIsoSurface<T> >("CIsoSurface", bpy::init<T>())
    .def("GenerateSurface", &CIsoSurface<T>::GenerateSurface)

    //(const T* ptScalarField, T tIsoLevel, unsigned int nCellsX, unsigned int nCellsY,  unsigned int nCellsZ, float fCellLengthX, float fCellLengthY, float fCellLengthZ);

    // Returns true if a valid surface has been generated.
    .def("IsSurfaceValid", &CIsoSurface<T>::IsSurfaceValid)

    // Deletes the isosurface.
    .def("DeleteSurface", &CIsoSurface<T>::DeleteSurface)

    // Returns the length, width, and height of the volume in which the
    // isosurface in enclosed in.  Returns -1 if the surface is not
    // valid.
    .def("GetVolumeLengths", &CIsoSurface<T>::GetVolumeLengths);
    //(float& fVolLengthX, float& fVolLengthY, float& fVolLengthZ);
    //  .def(init<float, float, float, std::string>())
    //.def("show", &CWmatrix<T>().show);
}

BOOST_PYTHON_MODULE(matrix)
{
  export_cisosurface<int>();
}

#endif
