/* 
pycapi_utils.h

Containers of C utilities routines for conversion of data objects between
Python and C (usually through weave)

NOTE: requires -std=c++11 compiler extension to be compiled

Maurizio De Pitta', The University of Chicago, October 2nd, 2014.
*/

#if ! defined _PYCAPI_UTILS_H
#define _PYCAPI_UTILS_H 1

#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>    // intptr_t type
#include <Python.h>
#include <numpy/arrayobject.h>
#include <tuple>
#include <list>

typedef std::tuple<const char*,PyObject*> pykeyval; //Tuple type (string,Pyobj*) as dictionary entry (key,val)
typedef std::list<pykeyval> kvlist;                 //Type List of tuples for all dictionary entries

void init_numpy(); // Initialize Python/Numpy interface (must be called within each module using Python data)
PyObject* array_double_to_pyobj(double* v_c, intptr_t NUMEL);        // Convert from array to Python list (double)
PyObject* array_int_to_pyobj(int* v_c, intptr_t NUMEL);              // Convert from array to Python list (int)
PyObject* build_PyDict(kvlist & items);                              // Convert from kvlist of PyObjects to PyDict

// Call Python from C
PyObject* callPyFunc(char *module,char *function,double *args, intptr_t NARGS);              // Call function in module with arguments args

#endif