//pycapi_utils.cpp

#include "pycapi_utils.h"

void init_numpy(){
    Py_Initialize();
    if(PyArray_API==NULL){
        import_array();
    }
}

PyObject* array_double_to_pyobj(double* v_c, intptr_t NUMEL){
    init_numpy();
    PyObject* out_array = PyArray_SimpleNew(1, &NUMEL, NPY_DOUBLE);
    double* v_b = (double*) ((PyArrayObject*) out_array)->data;
    for (int i=0;i<NUMEL;i++) v_b[i] = v_c[i];
    free(v_c);
    return out_array;
    //Py_XDECREF(out_array);
}

PyObject* array_int_to_pyobj(int* v_c, intptr_t NUMEL){
    init_numpy();
    PyObject* out_array = PyArray_SimpleNew(1, &NUMEL, NPY_INT);
    int* v_b = (int*) ((PyArrayObject*) out_array)->data;
    for (int i=0;i<NUMEL;i++) v_b[i] = v_c[i];
    free(v_c);
    return out_array;
    //Py_XDECREF(out_array);
}

PyObject* build_PyDict(kvlist & items){
/*
Convert a kvlist into a Python dictionary object
*/
    init_numpy();
    PyObject* pyd = PyDict_New();
    for(std::list<pykeyval>::iterator it = items.begin(); it != items.end(); it++){
            PyDict_SetItem(pyd,PyString_FromString(std::get<0>(*it)),std::get<1>(*it));
    }
    return pyd;
}

PyObject* callPyFunc(char *module, char *function, double *args, intptr_t NARGS){
/*
  Call Python function in module with NARGS arguments.

  Adapted from https://docs.python.org/release/2.6.5/extending/embedding.html#pure-embedding
*/

    PyObject *pName, *pModule, *pDict, *pFunc;
    PyObject *pArgs, *pValue;
    int i;

    Py_Initialize();
    pName = PyString_FromString(module);
    /* Error checking of pName left out */

    pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    if (pModule != NULL) {
        pFunc = PyObject_GetAttrString(pModule, function);
        /* pFunc is a new reference */

        if (pFunc && PyCallable_Check(pFunc)) {
            pArgs = PyTuple_New(NARGS);
            for (i=0;i<NARGS;++i) {
                pValue = PyFloat_FromDouble(args[i]);
                if (!pValue) {
                    Py_DECREF(pArgs);
                    Py_DECREF(pModule);
                    fprintf(stderr, "Cannot convert argument\n");
                    return NULL;
                }
                /* pValue reference stolen here: */
                PyTuple_SetItem(pArgs, i, pValue);
            }
            pValue = PyObject_CallObject(pFunc, pArgs);
            Py_DECREF(pArgs);
            // TODO complete
            if (pValue != NULL) {
//                Py_Finalize();
                return pValue;
//                Py_DECREF(pValue);
            }
            else {
                Py_DECREF(pFunc);
                Py_DECREF(pModule);
                PyErr_Print();
                fprintf(stderr,"Call failed\n");
//                Py_Finalize();
                return NULL;
            }
        }
        else {
            if (PyErr_Occurred())
                PyErr_Print();
            fprintf(stderr, "Cannot find function \"%s\"\n", function);
        }
        Py_XDECREF(pFunc);
        Py_DECREF(pModule);
//        Py_Finalize();
        return NULL;
    }
    else {
        PyErr_Print();
        fprintf(stderr, "Failed to load \"%s\"\n", module);
//        Py_Finalize();
        return NULL;
    }
}