#ifndef BHTSNE_EXPORT_
#define BHTSNE_EXPORT_

#include <Python.h>
#include <numpy/ndarrayobject.h>

#include <boost/python.hpp>
#include <boost/python/class.hpp>
#include <boost/python/raw_function.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/module.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/tuple.hpp>

#include "bhtsne.h"

namespace bpy = boost::python;

namespace num_util {
  PyArray_TYPES type(bpy::numeric::array& arr){
    return PyArray_TYPES(PyArray_TYPE(arr.ptr()));
  }
  void check_type(bpy::numeric::array& arr, PyArray_TYPES tp) {
    PyArray_TYPES pyType = type(arr);
    if (pyType != tp){
      PyErr_SetString(PyExc_TypeError, "Input type is wrong");
      bpy::throw_error_already_set();
    }
    return;
  }
  void* data(bpy::numeric::array& arr) {
    if (!PyArray_Check(arr.ptr())){
      PyErr_SetString(PyExc_ValueError, "Error: input argument is not a PyArrayObject.");
      bpy::throw_error_already_set();
    }
    return PyArray_DATA(arr.ptr());
  }


}

void BHTSNE_wrapper(bpy::numeric::array& data, bpy::numeric::array& Y, int numOfsamples, int dims, int mapped_D, double perplex, double th, int rseed, bool verbose) {
  // convert data of numpy array to double pointer
  double* samples = (double*)num_util::data(data);

  //create an object
  BHTSNE bhtsne(samples, numOfsamples, dims, mapped_D, perplex, th, rseed, verbose);
  bhtsne.run();
  // convert double** to numeric array
  for (int i = 0; i < numOfsamples; ++i){
    int ind1 = i*mapped_D;
    for (int j = 0; j < mapped_D; ++j){
      Y[bpy::make_tuple(i,j)] = bhtsne.Y[ind1+j];
    }
  }
}

void export_bhtsne() {
//  bpy::class_< BHTSNE_wrapper > ("BHTSNE", 
//    bpy::init<double*, int, int, double, double, int, bool>())
//    .def("run", &BHTSNE::run);
  def("bhtsne", BHTSNE_wrapper);
}

BOOST_PYTHON_MODULE(_bhtsne_wrap)
{
  import_array();
  bpy::numeric::array::set_module_and_type("numpy", "ndarray");
  export_bhtsne();
}

#endif
