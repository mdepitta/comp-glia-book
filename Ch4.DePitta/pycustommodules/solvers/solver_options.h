// solver_options.h
// Header file containing class and structures for solver options (custom C-based solvers)
// and GSL-based solvers). This header is included in ode_solvers.h

#include <Python.h>

#if ! defined _SOLVER_OPTIONS_H
#define _SOLVER_OPTIONS_H 1

// Solver Options class (for custom C-based solvers)
class set_solver{
    public:
    const char*  solver_id;
    double t0,tfin,transient;
    double h;   //dt
    double tbin;
    int nout;
    set_solver();
    set_solver(PyObject* opts_dict);
    ~set_solver();
};

// Options for GSL solvers
class set_solver_gsl{
    public:
    const char* solver_id;
    double t0,tfin;
    double dt;
    double atol,rtol;
    long int nstep;
    set_solver_gsl();
    set_solver_gsl(PyObject* opts_dict);
    ~set_solver_gsl();
};

// Output struct for the integrator
struct out_solver{
    long int NUMEL;
    double* ysol;
};

#endif