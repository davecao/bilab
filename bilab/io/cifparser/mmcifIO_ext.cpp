//
//  mmcifIO_ext.cpp
// 
//
//  Created by 曹巍 on 2018/06/19.
//  Copyright © 2018年 巍 曹. All rights reserved.
//
#include "mmcifIO_wrapper.h"

static 
PyObject* reader_mmcif_wrapper(PyObject *self, PyObject *args){
  
  PyObject *pDict = nullptr;

  const char *file;
  bool state;

  std::string filename;
  std::string source_name;
  std::string category_name;
  std::string item_name;
  

  // parse arguments
  if (!PyArg_ParseTuple(args, "s", &file)) {
    return NULL;
  }
  filename = file;
  mmcif::CIFDocument* doc = new mmcif::CIFDocument();
  // allocate for pDict
  pDict = PyDict_New(); // new reference
  assert(PyDict_Check(pDict));
  
  // Parse the document
  // call
  state = reader_mmcif(filename, *doc);
  if(!state){
    std::cout<<"Failed to parse the file: "<< filename << std::endl; 
  }
  // Source: PDBID
  source_name = doc->source;
  //doc -> show_block_names();
  // Convert to python's dictionary
  for (unsigned int i = 0; i < doc->getNumBlocks(); i++){
    auto block = doc->blocks[i];
    // Example: _atom_site.id:xxxx
    // category name: atom_site
    // item name: id
    category_name = block->name;
    PyObject* pCategoryDict = PyDict_New(); // new reference
    assert(PyDict_Check(pCategoryDict));
    
    // get the table
    auto tbl = block->table;
    int rows = tbl->GetNumTuples();
    int cols = tbl->GetNumColumns();
    
    for (unsigned int j = 0; j < rows; j++) {
      // item name : (j, 0) i.e., k=0
      // e.g., group_PDB 
      std::string pItemKey = (*tbl)(j, 0);
      // item values
      
      PyObject *pList = PyList_New(0);
      if(!pList) {
        Py_RETURN_NONE;
      }
      // values : ["ATOM", "ATOM", ..., "ATOM"]
      for (unsigned int k = 1; k < cols; k++) {
        // append
        PyList_Append(pList, Py_BuildValue("s", (*tbl)(j, k).c_str()));
        //std::cout<< (*tbl)(j, k) <<"@@";
      }
      //std::cout<<"\n";
      //std::cout << category_name <<"."<<pItemKey <<"("<<cols-1<<")-" <<PyList_Size(pList) <<"\n";
      if (category_name == "atom_site") {
        for(unsigned int l=0; l<PyList_Size(pList); ++l){
          PyObject *pValue = PyList_GetItem(pList, l);
          std::cout << PyString_AsString(pValue) << "\n";
        }
      }
      // save to pCategoryDict
      int state = PyDict_SetItem(pCategoryDict, PyString_FromString(pItemKey.c_str()), pList);
      if (state == 0) {
        // Release or clear
        Py_DECREF(pList);
      }else {
        std::cout<<"Failed to append."<<"\n";
      }
    }// End of category
    PyDict_SetItem(pDict, PyString_FromString(category_name.c_str()), pCategoryDict);
  } // End of blocks
  PyDict_SetItemString(pDict, "source", Py_BuildValue("s", source_name.c_str()));
  
  
  // return dict
  // Release
  if (doc != nullptr) {
    delete doc;
  }
  if (pDict == nullptr){
    Py_RETURN_NONE;
  }
  return pDict;
} 


static PyMethodDef Py_methods[] = {
    {"reader_mmcif", reader_mmcif_wrapper, METH_VARARGS, "C++ extension for reading mmcif format"},
    {NULL, NULL, 0, NULL}  /* Sentinel */
};

/* ========================================================================== */
/* -- Initialization -------------------------------------------------------- */
/* ========================================================================== */

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_mmcifo",                          // name of this module
        " C++ extension for mmcif reading ", // Doc string       
        -1,
        Py_methods,
        NULL,
        NULL,
        NULL,
        NULL
};

PyObject *
PyInit__mmcifio(void)

#else
PyMODINIT_FUNC
init_mmcifio(void)
#endif
{
  PyObject *module;

  import_array();

#if PY_MAJOR_VERSION >= 3
  module = PyModule_Create(&moduledef);
  if (module==NULL) return NULL;
#else
  module = Py_InitModule("_mmcifio", Py_methods);
  if (module==NULL) return;
  
#endif

  if (PyErr_Occurred()) Py_FatalError("can't initialize module _mmcifio");
#if PY_MAJOR_VERSION >= 3
  return module;
#endif
}