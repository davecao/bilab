//
//  NSCSurface.c
//
//
//  Created by 曹巍 on 2018/12/25.
//  Copyright © 2018年 巍 曹. All rights reserved.
//
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include "nsc.h"

/* Must define Py_TYPE for Python 2.5 or older */
#ifndef Py_TYPE
#  define Py_TYPE(o) ((o)->ob_type)
#endif

/* Must define PyVarObject_HEAD_INIT for Python 2.5 or older */
#ifndef PyVarObject_HEAD_INIT
#define PyVarObject_HEAD_INIT(type, size)       \
        PyObject_HEAD_INIT(type) size,
#endif

#define VARNAME(var) #var
// void printPyObject(PyObject *myPyObject) {
//   // @todo Get this to actually work for complex data types!
//   static size_t depth = 0;    //track how far we descend into the object
//   Py_ssize_t pos;
//   PyObject* myData;
//   const char* s;
//
//   //cout << string(depth, ' ') + "myPyObject address: " << (int)myPyObject << endl;
//   //process based on the object type
//   //Note: really bad case/switch statement :-(
//
//   //null
//   if (myPyObject == NULL || myPyObject == Py_None) {
// 	//do nothing, empty object
//   }
//   //string
//   else if (PyUnicode_Check(myPyObject)) {
// 	string myString;
// 	myString = PyUnicode_AsUTF8(myPyObject);
// 	cout << string(depth, ' ') + "string: " << myString << endl;
//   }
//   //bool
//   else if (PyBool_Check(myPyObject)) {
//     bool myBool;
//     myBool = (myPyObject == Py_True);
//     cout << string(depth, ' ') + "bool: " << myBool << endl;
//   }
//   //float
//   else if (PyFloat_Check(myPyObject)) {
//     double myFloat;
//     myFloat = PyFloat_AsDouble(myPyObject);
//     cout << string(depth, ' ') + "float: " << myFloat << endl;
//   }
//   //long (integer)
//   else if (PyLong_Check(myPyObject)) {
//     long myInt;
//     myInt = PyLong_AsLong(myPyObject);
//     cout << string(depth, ' ') + "long: " << myInt << endl;
//   }
//   //tuple
//   else if (PyTuple_Check(myPyObject)) {
//     // @todo Get tuples to work
//     depth++;
//     cout << string(depth, ' ') + "tuple: " << endl;
//     //iterate through all items in the tuple
//     int tsize = PyTuple_Size(myPyObject);
//     for (int pos = 0; pos < tsize; pos++) {
//       myData = PyTuple_GetItem(myPyObject, pos);
//       //int thisInt = PyLong_AsLong(myData);
//       //cout << "item " << pos << ": " << thisInt << endl;
//       Py_INCREF(myData);
//       printPyObject(myData);
//       Py_DECREF(myData);
//     }
// 	depth--;
//   }
//   //list
//   else if (PyList_Check(myPyObject)) {
//     // @todo get lists to work
//     depth++;
//     cout << string(depth, ' ') + "list: " << endl;
//     //iterate through all items in the list
//     int tsize = PyList_Size(myPyObject);
//     for (int pos = 0; pos < tsize; pos++) {
//       myData = PyList_GetItem(myPyObject, pos);
//       //int thisInt = PyLong_AsLong(myData);
//       //cout << "item " << pos << ": " << thisInt << endl;
//       Py_INCREF(myData);
//       printPyObject(myData);
//       Py_DECREF(myData);
//     }
//     depth--;
//   }
//   //dict
//   else if (PyDict_Check(myPyObject)) {
//     // @todo get dicts to work
// 	depth++;
//     cout << string(depth, ' ') + "dict: " << endl;
//     //iterate through all items in the dict
//     PyObject *key, *value;
//     Py_ssize_t dpos = 0;
//
//     while (PyDict_Next(myPyObject, &dpos, &key, &value)) {
//       Py_INCREF(key);
//       Py_INCREF(value);
//
//       printPyObject(key);
//       printPyObject(value);
//
//       Py_DECREF(key);
//       Py_DECREF(value);
// 	}
//     depth--;
//   }
// }

int not_double_matrixs(PyArrayObject *vectors) {
  int dims = PyArray_NDIM(vectors);
  PyArray_Descr *descr = PyArray_DESCR(vectors);
  if (descr->type_num != NPY_DOUBLE || dims != 2)  {
    PyErr_SetString(PyExc_ValueError,
      "In not2Ddouble_vectors: array must be of type Float and 2 dimensional (n).");
      return 1;
  }
  return 0;
}

int not_doublevector(PyArrayObject *vec)  {
  int dims = PyArray_NDIM(vec);
  PyArray_Descr *descr = PyArray_DESCR(vec);
  if (descr->type_num != NPY_DOUBLE || dims != 1)  {
    PyErr_SetString(PyExc_ValueError,
      "In not_doublevector: array must be of type Float and 1 dimensional (n).");
      return 1;
  }
  return 0;
}

double **ptrvector(long n){
  double **v;
  v = (double **) malloc((size_t)(n*sizeof(double)));
  if(!v){
    printf("ptrvector: Allocation of memory for double array failed.\n");
	exit(0);
  }
  return v;
}

double **pymatrix_to_Carrayptrs(PyArrayObject *arrayin){
	double **c, *a;
	int i, n, m;
	npy_intp *dims = PyArray_DIMS(arrayin);
	n = dims[0]; //arrayin->dimensions[0];
	m = dims[1]; //arrayin->dimensions[1];
	c = ptrvector(n);
	//a = (double *) arrayin->data;
	a = PyArray_DATA(arrayin);
	for(i=0; i<n; i++){
		c[i] = a + i*m;
	}
	return c;
}

// Create 1D Caray from PyArray
double* pyvector_to_Carrayptrs(PyArrayObject *arrayin){
	int n;
	npy_intp *dims = PyArray_DIMS(arrayin);
	n = dims[0];//arrayin->dimensions[0];
	//return (double *) arrayin->data;
	__attribute__((unused)) double* data = PyArray_DATA(arrayin);
	return data;
}

PyArrayObject *pyvector(PyObject *objin) {
	return (PyArrayObject *) PyArray_ContiguousFromObject(objin,
		NPY_DOUBLE, 1, 1);
}

PyArrayObject *pymat(PyObject *objin) {
	return (PyArrayObject *) PyArray_ContiguousFromObject(objin,
		NPY_DOUBLE, 2, 2);
}

void free_Carrayptrs(double **v) {
	free((char*) v);
}

static
PyObject* NSC_wrapper(PyObject *self, PyObject *args){

	// Input:
	// 	Coord (numpy array),
	//  Radii (numpy array),
	//	TotalAt (int),
	//  DotsSphere (int),
	//	OverallArea (int)
	PyArrayObject *Coord=NULL, *Radii=NULL;
	PyTupleObject *Result_tuple=NULL;
	// Output:
	//	 OverallArea
	PyArrayObject *areaPerAtom=NULL;
	PyObject *Coord_Flatten=NULL;
	// C parameters
	double *Coord_c, *Radii_c, overallArea_c;
	double *areaPerAtom_c;

	int Coord_length, Radii_length;
	int totalAt_c, dotsPerSphere_c;
	int dims[2];

	totalAt_c = 0;
	dotsPerSphere_c = 600;
	overallArea_c = 0;

    if (!PyArg_ParseTuple(
		args, "O!O!ii",
		&PyArray_Type, &Coord,
		&PyArray_Type, &Radii,
		&totalAt_c, &dotsPerSphere_c)) {
		return NULL;
    }
	if (NULL == Coord)  return NULL;
	if (NULL == Radii)  return NULL;
	/* Not needed if python wrapper function checks
	   before call to this routine */
	if (not_double_matrixs(Coord)) return NULL;
	if (not_doublevector(Radii)) return NULL;

	// Get 2D vectors dimensions: PyArray_Flatten
	Coord_length = (int) PyArray_DIM(Coord, 0);
	Radii_length = (int) PyArray_DIM(Radii, 0);

	if (Coord_length != Radii_length || Coord_length != totalAt_c ||
		Radii_length != totalAt_c) {
		printf("The length of Coord or Radii are not equal.\n");
		printf("Or Coord and Radii are not same size.\n");
		exit(0);
	}

	/* Change contiguous arrays into C * arrays   */
	Coord_Flatten = PyArray_Flatten(Coord, NPY_CORDER);
    Coord_c = pyvector_to_Carrayptrs(Coord_Flatten);
    Radii_c = pyvector_to_Carrayptrs(Radii);

	// Call NSC accessible surface area
	NSC(Coord_c, Radii_c, totalAt_c, dotsPerSphere_c,
		FLAG_ATOM_AREA, &overallArea_c, &areaPerAtom_c,
		NULL, NULL, NULL);

	/* Return: Make a new double array of same dims */
	dims[0] = Coord_length;
	dims[1] = 0;
	//areaPerAtom = (PyArrayObject *) PyArray_FromDims(1, dims, NPY_DOUBLE);
	areaPerAtom = (PyArrayObject *) PyArray_SimpleNewFromData(
		1, dims, NPY_DOUBLE, areaPerAtom_c);
	if (areaPerAtom == NULL) {
	  printf("Failed to create an array, areaPerAtom, dims=%dx%d\n",
			  dims[0],dims[1]);
	  return NULL;
	}
	// Return a tuple, (overall, areaPerAtom)
	// OverallArea = PyFloat_FromDouble(overallArea_c);
	Result_tuple = PyTuple_New(2);
	if (NULL == Result_tuple) {
	  printf("Failed to create a tuple for returning.\n");
	  exit(0);
	}
	PyTuple_SetItem(Result_tuple, 0, PyFloat_FromDouble(overallArea_c));
	PyTuple_SetItem(Result_tuple, 1, areaPerAtom);
	// free memory
	//free_Carrayptrs(Coord_c);
	//free_Carrayptrs(Radii_c);
	//free_Carrayptrs(Coord_Flatten);
	//free_Carrayptrs(areaPerAtom);

	return Result_tuple;
	//return PyArray_Return(areaPerAtom);
}

static PyMethodDef Py_methods[] = {
    {"_nsc", NSC_wrapper, METH_VARARGS, "C extension for calculate accessible surface area"},
    {NULL, NULL, 0, NULL}  /* Sentinel */
};
/* ========================================================================== */
/* -- Initialization -------------------------------------------------------- */
/* ========================================================================== */


#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_nsc",                     // name of this module
        " C extension for NSC ASA", // Doc string
        -1,                         // m_size
        Py_methods,                 // m_methods
        NULL,                       // m_reload
        NULL,                       // m_traverse
        NULL,                       // m_clear
        NULL                        // m_free
};
#endif

PyMODINIT_FUNC
#if PY_MAJOR_VERSION >= 3
PyInit__nsc(void)
#else
init_nsc(void)
#endif
{
  PyObject *module;
  import_array();

#if PY_MAJOR_VERSION >= 3
  module = PyModule_Create(&moduledef);
  if (module==NULL) return NULL;
#else
  module = Py_InitModule("_nsc", Py_methods);
  if (module==NULL) return;
#endif

  if (PyErr_Occurred()) Py_FatalError("can't initialize module _nsc");
#if PY_MAJOR_VERSION >= 3
  return module;
#endif
}
