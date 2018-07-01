#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <functional> //std::bind
#include <cstring>
//#include <memory>
#include <iterator>   //noise generator
#include <random>     //noise generator
#include <gsl/gsl_errno.h>  // GSL error library
#include <gsl/gsl_matrix.h> // GSL matrix library for Jacobian
#include <gsl/gsl_odeiv2.h> // GSL ODE solver library

#include "solver_options.h"
//#define DBL_MEMCPY(dest,src,n) memcpy((dest),(src),(n)*sizeof(double))

#if ! defined _STOCHASTIC_SOLVERS_H
#define _STOCHASTIC_SOLVERS_H 1

// Typedef Pointer to Jacobian function;
typedef int(*jacPtr)(double t, const double y[], double * dfdy, double dfdt[], void * params);

//----------------------------------------------------------------------------------------------------------------------
// Stochastic ODE system definition
//----------------------------------------------------------------------------------------------------------------------
class sde{
    public:
        int noise = 0;  // noise flag
        size_t N_d,N_s,NEQ;     // # of deterministic variables, # of stochastic variables, # equations
        double *y;      // variables
        double *ddy;    // deterministic dy
        double *sdy;    // stochastic dy
        double *sdy_d;  // derivative of stochastic dy
//        void   *params; // Model parameters
        sde(size_t Nd_,size_t Ns_,double *y0,
            int (*dfp)(double, const double*, double*, void*),
            int (*jcb)(double, const double*, double*, double*, void*),
            int (*sfp)(double, double, const double*, double*, void*),
            int (*sdfp)(double, double, const double*, double*, void*));
        ~sde();
        int (*drift)(double t, const double y_[], double dy_[], void * params);
        int (*jac)(double t, const double y_[], double *dfdy_, double dfdt_[], void *params);
        int (*diffusion)(double t, double dt, const double y[], double dydt[], void *params);
        int (*diffusion_derivative)(double t, double dt, const double y[], double dydt[], void *params);
        static int drift_wrapper(double t, const double y_[], double dy_[], void * params);
        static int jac_wrapper(double t, const double y_[], double *dfdy_, double dfdt_[], void *params);
        static int diffusion_wrapper(double t, double dt, const double y[], double dydt[], void *params);
        static int diffusion_derivative_wrapper(double t, double dt, const double y[], double dydt[], void *params);
};

//----------------------------------------------------------------------------------------------------------------------
// Generic stochastic solver
//----------------------------------------------------------------------------------------------------------------------
class de_solver{
    size_t N_d,N_s,NEQ; // System dimensions
    int noise;          // noise flag
    public:
        sde *sys;        // SDE/ODE system to integrate
        double *y_old;   // value of y before an integration step
        double t;        // time variable
        int event = 0;   // Flag for event-related solutions
        de_solver(sde *sys_);
        ~de_solver();
        void milstein(double t, double dt, void *params);
        void integrate(set_solver_gsl opts, void *params, void (*fpa)(double*, double p[]),
                       double *f_var, double *f_par,double *ts, double **sol);
        int integrate_event(set_solver_gsl opts, void *params, void (*fpa)(void**, void **p),
                            void **f_var, void **f_par,double *ts, double **sol);
};
#endif