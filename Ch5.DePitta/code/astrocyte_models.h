/*
astrocyte_models.h
HEADER file for astrocyte_models.cpp
Contains class, structure and function declaration to simulate several astrocyte
models.

v2.0
Implemented integration by GSL routines, and added LR model + PKC model.
Maurizio De Pitta', INRIA Rhone-Alpes, November 3rd, 2017.

v1.0
Maurizio De Pitta', The University of Chicago, February 28th, 2015.
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Python.h>
#include <numpy/arrayobject.h>  /* NumPy as seen from C */
#include <type_traits>
#include <gsl/gsl_errno.h>  // GSL error library
#include <gsl/gsl_matrix.h> // GSL matrix library for Jacobian
#include <gsl/gsl_odeiv2.h> // GSL ODE solver library
#include <cstddef>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <memory>

#include "pycapi_utils.h"
#include "solver_options.h"

//Define Hill function
#define HILL(x,k,n) (pow(x,n)/(pow(x,n)+pow(k,n)))
// Derivative Hill wrt x
#define Hx(x,k,n) (n/x*HILL(x,k,n)*(1.0-HILL(x,k,n)))
// Q2
#define Q2(I,d1,d2,d3) (d2*(I+d1)/(I+d3))
// Derivative of Q2
#define Q2x(I,d1,d2,d3) (d2*(d1-d3)/d3*HILL(I,d3,1)*(1.0-HILL(I,d3,1)))
//// Heaviside function
//#define HEAV(x) (x>0.0)?1.0:0.0
//// Square pulse
//#define SQRP(x,d) (HEAV(x)-HEAV(x-d))

#if ! defined _ASTROCYTE_MODELS_H
#define _ASTROCYTE_MODELS_H 1

/*-----------------------------------------------------------------------------
DATA TYPES
-----------------------------------------------------------------------------*/
//Li-Rinzel par data struct
struct lra_pars{
    double d1,d2,d3,d5,a2;
    double c1,c0,rc,rl;
    double ver,ker;
    double ip3;
};

//CHI Par data struct
struct chi_pars{
    double vbias,vbeta;
    double vdelta,kappad,kdelta;
    double v3k,kd,k3,r5p;
    double d1,d2,d3,d5,a2;
    double c1,c0,rc,rl;
    double ver,ker;
};

//GCHI Par data struct
struct gchi_pars{
    double yrel,vp,omp,zeta,kkc;
    double *stimulus;
    double vbias,vbeta;
    double vdelta,kappad,kdelta;
    double v3k,kd,k3,r5p;
    double d1,d2,d3,d5,a2;
    double c1,c0,rc,rl;
    double ver,ker;
};

//eGCHI Par data struct
struct egchi_pars{
    double yrel,vp,omp;
    double vkd,kkc,vk,omkd;
    double vd,kdc,kdd,omd;
//    double vkp,rtot;
//    double vk,kkc,kkd;
//    double vd,kdc,kdd,omd;
//    double alpha,beta;
    double vbias,vbeta;
    double vdelta,kappad,kdelta;
    double v3k,kd,k3,r5p;
    double d1,d2,d3,d5,a2;
    double c1,c0,rc,rl;
    double ver,ker;
};

/*-----------------------------------------------------------------------------
CLASS DECLARATION
-----------------------------------------------------------------------------*/
// LR model class
class lra{
    //Parameters
    lra_pars pars;
    public:
        size_t NEQ = 2;
        // Variables
        double y[2],dy[2];  //C,h,I state variables; their derivatives; and the old values
        lra();  //constructor
        ~lra(){}; //destructor
        void set_pars(lra_pars pars);      //Initializer of custom parameters
        void set_pars(PyObject* pardict);  //Initializer of custom parameters (through Pydictionary)
        void set_ics(double* ics);         //Initializer for ICs
        void set_ics(PyObject* pardict);   //Initializer for ICs (through Pydictionary)
        double bias(double t,double dt);   //time-dependent input signal in the ode
        void ode(double t,double dt);      //ODE (used for integration) (standard form)
        int ode_gsl(double t, const double y_[], double dy_[]); // GSL wrapper
        static int ode_gsl_wrapper(double t, const double y_[], double dy_[], void * params);
};

// CHI model class
class chi{
    //Parameters
    chi_pars pars;
    public:
        size_t NEQ = 3;
        // Variables
        double y[3],dy[3];  //C,h,I state variables; their derivatives; and the old values
        chi();  //constructor
        ~chi(){}; //destructor
        void set_pars(chi_pars pars);      //Initializer of custom parameters
        void set_pars(PyObject* pardict);  //Initializer of custom parameters (through Pydictionary)
        void set_ics(double* ics);         //Initializer for ICs
        void set_ics(PyObject* pardict);   //Initializer for ICs (through Pydictionary)
        double bias(double t,double dt);   //time-dependent input signal in the ode         
        void ode(double t,double dt);      //ODE (used for integration) (standard form)
        int ode_gsl(double t, const double y_[], double dy_[]); // GSL wrapper
        static int ode_gsl_wrapper(double t, const double y_[], double dy_[], void * params);
};

// G-CHI model class
class gchi{
    //Parameters
    gchi_pars pars;
    public:
        size_t NEQ = 4;
        // Variables
        double y[4],dy[4];  // C,h,I state variables; their derivatives; and the old values
        gchi();  // constructor
        ~gchi(){free(pars.stimulus);}; // destructor
        void set_pars(gchi_pars pars);     // Initializer of custom parameters
        void set_pars(PyObject* pardict);  // Initializer of custom parameters (through Pydictionary)
        void set_ics(double* ics);         // Initializer for ICs
        void set_ics(PyObject* pardict);   // Initializer for ICs (through Pydictionary)
        double bias(double t,double dt);   // time-dependent input signal in the ode
        void ode(double t,double dt);      // ODE (used for integration) (standard form)
        int ode_gsl(double t, const double y_[], double dy_[]); // GSL wrapper
        static int ode_gsl_wrapper(double t, const double y_[], double dy_[], void * params);
};

// eG-CHI model class
class egchi{
    //Parameters
    egchi_pars pars;
    public:
        size_t NEQ = 6;
        // Variables
        double y[6],dy[6];  // C,h,I state variables; their derivatives; and the old values
        double dfdy[36];    // Jacobian, row-ordered as required by GSL
        egchi();    // constructor
        ~egchi(){}; // destructor
        void set_pars(egchi_pars pars);     // Initializer of custom parameters
        void set_pars(PyObject* pardict);  // Initializer of custom parameters (through Pydictionary)
        void set_ics(double* ics);         // Initializer for ICs
        void set_ics(PyObject* pardict);   // Initializer for ICs (through Pydictionary)
        double bias(double t,double dt);   // time-dependent input signal in the ode
        void ode(double t,double dt);      // ODE (used for integration) (standard form)
        int ode_gsl(double t, const double y_[], double dy_[]); // GSL wrapper
        static int ode_gsl_wrapper(double t, const double y_[], double dy_[], void * params);
        void jacobian(double t, double dt);     // Jacobian
        int jacobian_gsl(double t, const double y_[], double * dfdy_, double dfdt_[]);
        static int jacobian_gsl_wrapper(double t, const double y_[], double * dfdy_, double dfdt_[], void * params);
};

/*-----------------------------------------------------------------------------
//Output data struct
-----------------------------------------------------------------------------*/
class out_lra{
    long int N = 1;
    public:
    double* ts; //time of solution
    double* ca;
    double* h;
    out_lra();
    out_lra(long int N);
    ~out_lra(){}; // Empty destructor to avoid deallocated variables useful data to pass to python
    PyObject* make_PyDict();   //Convert output to PyDictionary
};

class out_chi{
    long int N = 1;
    public:
    double* ts; //time of solution    
    double* ip3;
    double* ca; 
    double* h;
    out_chi();    
    out_chi(long int N);
    ~out_chi(){}; // Empty destructor to avoid deallocated variables useful data to pass to python
    PyObject* make_PyDict();   //Convert output to PyDictionary
};

class out_gchi{
    long int N = 1;
    public:
    double* ts; //time of solution
    double* rec;    
    double* ip3;
    double* ca; 
    double* h;
    out_gchi();    
    out_gchi(long int N);
    ~out_gchi(){}; // Empty destructor to avoid deallocated variables useful data to pass to python
    PyObject* make_PyDict();   //Convert output to PyDictionary
};

class out_egchi{
    long int N = 1;
    public:
    double* ts; //time of solution
    double* rec;
    double* ip3;
    double* ca;
    double* h;
    double* dag;
    double* pkc;
    out_egchi();
    out_egchi(long int N);
    ~out_egchi(){}; // Empty destructor to avoid deallocated variables useful data to pass to python
    PyObject* make_PyDict();   //Convert output to PyDictionary
};

/*-----------------------------------------------------------------------------
Functions
-----------------------------------------------------------------------------*/
// GSL wrappers
out_lra lra_simulator_gsl(PyObject* pars,PyObject* solver_options);    //LRA model simulator
out_chi chi_simulator_gsl(PyObject* pars, PyObject* solver_opts);      //CHI model simulator
out_gchi gchi_simulator_gsl(PyObject* pars, PyObject* solver_options); //G-CHI model simulator
PyObject* gchi_simulator_py_gsl(PyObject* pars, PyObject* solver_options);  // NOT implemented
out_egchi egchi_simulator_gsl(PyObject* pars, PyObject* solver_options);    // ExG-CHI model simulator
PyObject* egchi_simulator_py_gsl(PyObject* pars, PyObject* solver_options); // NOT implemented
#endif