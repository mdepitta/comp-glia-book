// solver_options.cpp
// Source file for solver_options.h

#include "solver_options.h"

//-----------------------------------------------------------------------------
// Defintion of class set_solver
//-----------------------------------------------------------------------------
set_solver::set_solver(){
    solver_id = "euler";
    t0        = 0.0;
    tfin      = 1.0;
    h         = 1e-3;
    transient = 0.0;
    tbin      = 1e-2;
    nout      = (int)(tbin/h);
}

set_solver::set_solver(PyObject* opts_dict){
    solver_id = PyString_AsString(PyDict_GetItemString(opts_dict,"solver"));
    t0        = PyFloat_AS_DOUBLE(PyDict_GetItemString(opts_dict,"t0"));
    tfin      = PyFloat_AS_DOUBLE(PyDict_GetItemString(opts_dict,"tfin"));
    h         = PyFloat_AS_DOUBLE(PyDict_GetItemString(opts_dict,"dt"));
    transient = PyFloat_AS_DOUBLE(PyDict_GetItemString(opts_dict,"transient"));
    tbin      = PyFloat_AS_DOUBLE(PyDict_GetItemString(opts_dict,"tbin"));
    nout      = (int)(tbin/h);
}

set_solver::~set_solver(){}

set_solver_gsl::set_solver_gsl(){
    solver_id = "gsl_rk8pd";
    t0        = 0.0;
    tfin      = 1.0;
    dt        = 1e-3;
    atol      = 1e-8;
    rtol      = 1e-6;
    nstep     = (long int)((tfin-t0)/dt + 1);
}

set_solver_gsl::set_solver_gsl(PyObject* opts_dict){
    solver_id = PyString_AsString(PyDict_GetItemString(opts_dict,"solver"));
    t0        = PyFloat_AS_DOUBLE(PyDict_GetItemString(opts_dict,"t0"));
    tfin      = PyFloat_AS_DOUBLE(PyDict_GetItemString(opts_dict,"tfin"));
    dt        = PyFloat_AS_DOUBLE(PyDict_GetItemString(opts_dict,"dt"));
    atol      = PyFloat_AS_DOUBLE(PyDict_GetItemString(opts_dict,"atol"));
    rtol      = PyFloat_AS_DOUBLE(PyDict_GetItemString(opts_dict,"rtol"));
    nstep     = (long int)((tfin-t0)/dt + 1);
}

set_solver_gsl::~set_solver_gsl(){}