#ifndef BHTSNE_EXPORT_
#define BHTSNE_EXPORT_

#include <boost/python.hpp>
#include <boost/python/class.hpp>
#include <boost/python/raw_function.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/module.hpp>

#include "bhtsne.h"

namespace bpy = boost::python;

void export_bhtsne() {
  bpy::class_< BHTSNE > ("BHTSNE", 
    bpy::init<double*, int, int, double, double, int, bool>())
    .def("run", &BHTSNE::run);
}

BOOST_PYTHON_MODULE(_bhtsne_wrap)
{
  export_bhtsne();
}

#endif
