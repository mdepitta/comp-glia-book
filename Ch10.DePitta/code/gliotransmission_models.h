/*
gliotransmission_models.h
HEADER file for gliotransmission_models.cpp
Contains class, structure and function declaration to simulate several models of gliotransmission.

v1.0
Maurizio De Pitta', Basque Center of Applied Mathematics, Bilbao, Bizkaia, Spain
June  2018.
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
#include "stochastic_solvers.h"

// Define Hill function
#define HILL(x,k,n) (pow(x,n)/(pow(x,n)+pow(k,n)))
// Derivative Hill wrt x
#define Hx(x,k,n) (n*pow(x,n-1)/(pow(x,n)+pow(k,n))*HILL(k,x,n))
// U0
#define U0(u0_star,alpha_,gs_) ((1.0-gs_)*u0_star + alpha_*gs_)

#if ! defined _ASTROCYTE_MODELS_H
#define _ASTROCYTE_MODELS_H 1

/*-----------------------------------------------------------------------------
GENERAL UTILITY FUNCTION
-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------
Auxiliary functions (event handlers)
-----------------------------------------------------------------------------*/
void clip_h(double *x, double lims[]);
void c_thr(void **x, void **vals);
void clip_thr(void **x, void **vals);
void u0_update(double ts,void **x, void **gs_);

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

struct release_pars{
    double u0,taud,tauf;
    double taui,psc0;
};

struct gtrelease_pars{
    lra_pars calcium;
    // Gt. exocytosis parameters
    double cthr;        // Ca2+ threshold for release
    double ua, taua;
};

struct asn_pars{
    gtrelease_pars gliot;
    release_pars synapse;
    // Additional parameters;
    double rhoe,Gtot,taue,Ou,Ku,js; // Gt. clearance parameters
    double Og,taug,alpha; // presynaptic receptors
};
/*-----------------------------------------------------------------------------
CLASS DECLARATION
-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------
//Output data struct
-----------------------------------------------------------------------------*/
class out_lra{
    long int N = 1;
    public:
    double* ts; //time of solution
    double* ca;
    double* h;
    out_lra(){};
    out_lra(long int N);
    ~out_lra(){}; // Empty destructor to avoid deallocated variables useful data to pass to python
    PyObject* make_PyDict();   //Convert output to PyDictionary
};

// LR model class
class lra{
    //Parameters
    lra_pars pars;
    int noise = 0;  // noise flag (default no noise)
    public:
        size_t NEQ = 2;
        size_t N_d, N_s;   //Specification of N_deterministic vs. N_stochastic with N_d + N_s = NEQ
        // Variables
        double *y, *dy; //C,h,I state variables; their derivatives; and the old values
        double *dfdy;    //Jacobian, row-ordered as required by GSL
        lra();  //constructor
        ~lra(); //destructor
        void set_pars(lra_pars parameters);      //Initializer of custom parameters
        void set_pars(PyObject* pardict);  //Initializer of custom parameters (through Pydictionary)
        lra_pars get_pars();               //Retrieve model parameters
        void set_ics(double* ics);         //Initializer for ICs
        void set_ics(PyObject* pardict);   //Initializer for ICs (through Pydictionary)
        double tauh(double t,double dt);   //Provide tau_h (needed in stochastic part)
        double bias(double t,double dt);   //time-dependent input signal in the ode
        // RHS (deterministic version)
        void ode(double t,double dt);      //ODE (used for integration) (standard form)
        int ode_gsl(double t, const double y_[], double dy_[]); // GSL wrapper
        static int ode_gsl_wrapper(double t, const double y_[], double dy_[], void * params);
        // Jacobian
        void jacobian(double t, double dt);     // Jacobian
        int jacobian_gsl(double t, const double y_[], double * dfdy_, double dfdt_[]);
        static int jacobian_gsl_wrapper(double t, const double y_[], double * dfdy_, double dfdt_[], void * params);
        // RHS (Stochastic version of the model)
        // member functions must be in the form of sde class member functions)
        int ode_s(double t,double dt,const double y_[], double dy_[]);    //Determininistic part of ODEs (used for integration) (standard form)
        int ode_sd(double t,double dt,const double y_[], double dy_[]);   //Stochastic part of ODE (used for integration) (standard form)
        static int ode_s_wrapper(double t,double dt,const double y_[], double dy_[], void *params);
        static int ode_sd_wrapper(double t,double dt,const double y_[], double dy_[], void *params);
        // Simulator
        out_lra simulate();
        out_lra simulate(PyObject* pars,PyObject* solver_options);
};

/*-----------------------------------------------------------------------------
//Exocytosis models
-----------------------------------------------------------------------------*/
class out_release{
    long int N_spk;
    size_t NVAR, N_out;  // N_out is used only in mean-field output and specify which variables should be provided in the output
    int mean_field = 0;  // Will inform ::make_PyDict() on how to handle data.
    public:
        double* ts;
        double* u;
        double* x;
        double* y; // PSC (may not be needed)
        out_release(){}; // Empty constructor (for initialization of class inside other classes)
        out_release(long int Ns_, int Nv_); // Spiking model only
        out_release(long int NUMEL, int Nv_, int Nout_);        // Mean-field model only
        ~out_release(){}; // Empty destructor to avoid deallocated variables useful data to pass to python
        void get_dims();  // Print system dimenstions
        PyObject* make_PyDict();   //Convert output to PyDictionary
};

// Spiking model (Tsodyks and Markram)
class release{
    size_t NVAR;
    size_t N_syn;
    long int N_spk;
    long int counter;
    int create_ics = 0;     // If ics are created internally (-->1) (or passed -->0)
    int allocate_memory = 0;// Regulate freeing in case empty class is initialized
    public:
        double *ics;        // x(0),y(0),u(0)
        release_pars pars;  // model parameters must be made public in this case if we want to change them by gliotransmission
        double *u, *x, *y;
        double *last_spike;
        long int *spike_index;
        release(){};
        release(int Nv_,int Nsyn_,long int Nspk_);
        release(int Nv_,int Nsyn_,long int Nspk_,PyObject* pardict);
        ~release();
        void set_pars(release_pars parameters);// Set parameters
        void set_pars(PyObject* pardict);// Set parameters
        void set_ics(double *ics_);   // Initializer of ICs
        void set_ics(PyObject* ics_); // Set ICs from python dictionary
        void get_dims(); // Print private variables for size
        // Spiking (TM) model
        void x_sol(double x0_,double *x_,double dt_);  // Analytical solution of x starting from x(0) for [0,dt]
        void uy_sol(double u0_,double *u_,double dt_);  // Analytical solution of u (or y) starting from u(0) (y(0)) for [0,dt]
        void update(double* state, double ts, int isyn);
        // Simulator
        out_release simulate(double* spikes); // Plain synaptic model
        out_release simulate(double *spikes, void (*fpa)(double,void**,void**),void **f_var, void **f_par); // Allows operations by (fp) before update
};

// Mean-field model
class release_ave{
    size_t NVAR, NEQ;
    public:
        release_pars pars;  // model parameters must be made public in this case if we want to change them by gliotransmission
        double nu = 0;      // Rate of stimulation (currently passed as separate parameter)
        double *y,*dy,*dfdy;// state variables, derivatives and Jacobian are of size specified by NVAR
        release_ave(int Nv_);
        ~release_ave();
        void set_pars(PyObject* pardict);
        void set_ics(double* ics);
        void set_ics(PyObject* pardict);
        // RHS (deterministic version)
        double bias(double t,double dt);
        void ode(double t,double dt);      //ODE (used for integration) (standard form)
        int ode_gsl(double t, const double y_[], double dy_[]); // GSL wrapper
        static int ode_gsl_wrapper(double t, const double y_[], double dy_[], void * params);
        // Jacobian
        void jacobian(double t, double dt);     // Jacobian
        int jacobian_gsl(double t, const double y_[], double * dfdy_, double dfdt_[]);
        static int jacobian_gsl_wrapper(double t, const double y_[], double * dfdy_, double dfdt_[], void * params);
        // Simulator
        out_release simulate(double rate_,PyObject *pardict,PyObject* solver_options);
};

/*-----------------------------------------------------------------------------
// Composed models
-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------
// Calcium-dependent gliotransmitter exocytosis
-----------------------------------------------------------------------------*/
class out_gtrelease{
    long int N;
    public:
        out_lra calcium;
        double* x;
        out_gtrelease(){};
        out_gtrelease(long int NUMEL) : N(NUMEL),
                                        calcium(NUMEL){
            x = (double*)calloc(N,sizeof(double));
        };                  // Will use values in the bins to also deduce time values (in the python code)
        ~out_gtrelease(){}; // Empty destructor to avoid deallocated variables useful data to pass to python
        PyObject* make_PyDict();   //Convert output to PyDictionary
};

// Spiking gliotransmitter release model where release instant are computed in real time based on calcium dynamics
class gtrelease{
    public:
        lra calcium;               // Calcium signal generator
        gtrelease_pars pars;       // model parameters must be made public in this case if we want to change them by gliotransmission
        double x, r = 0;           // x and release r
        double last_spike = 0;     // Last spike value (set to 0) for proper handling of ICs
        double next_spike = 0;     // Next spike value (used in asn class)
        gtrelease(){};
        ~gtrelease(){};
        void set_pars(PyObject* pardict);
        void set_ics(double *ics_);   // Initializer of ICs
        void set_ics(PyObject* ics_); // Set ICs from python dictionary
        // Spiking (TM) model
        void update(double ts);
        // Simulator
        out_gtrelease simulate(PyObject* pardict,PyObject* solver_options);
};

/*-----------------------------------------------------------------------------
// Tripartite synapse class
-----------------------------------------------------------------------------*/
class out_asn{
    long int N, N_spk, N_gre;
    size_t NEQ,NVAR;
    int mean_field = 0;  // Mean-field flag
    public:
        out_gtrelease gliot;
        out_release synapse;
        double *ts;
        double *xa, *ga, *gammas;
        out_asn(){};
        out_asn(size_t Neq_,size_t Nv_,long int NUMEL,long int Ns_,long int Ng_) : N(NUMEL),
                                                                                   N_spk(Ns_),
                                                                                   N_gre(Ng_),
                                                                                   NEQ(Neq_),
                                                                                   NVAR(Nv_),
                                                                                   synapse(Ns_,Nv_),
                                                                                   gliot(NEQ>2 ? NUMEL : Ng_){
            // Additional variables
            ts = (double*)calloc(N,sizeof(double));
            ga = (double*)calloc(N,sizeof(double));
            gammas = (double*)calloc(N,sizeof(double));
        };
        out_asn(size_t Neq_,size_t Nv_,long int NUMEL) : N(NUMEL),
                                                         NEQ(Neq_),
                                                         NVAR(Nv_),
                                                         synapse(NUMEL,Nv_,Nv_+1){
            // Additional variables
            ts = (double*)calloc(N,sizeof(double));
            xa = (double*)calloc(N,sizeof(double));
            ga = (double*)calloc(N,sizeof(double));
            gammas = (double*)calloc(N,sizeof(double));
            mean_field = 1;
        };
        ~out_asn(){};
        PyObject* make_PyDict();
};

// Spiking tripartite synapse model
class asn{
    asn_pars pars;
    size_t NEQ,N_d,N_s;
    size_t NVAR;
    size_t N_syn;
    long int N, N_spk, N_gre;
    public:
        // Continuous variables
        gtrelease gliot;       // Gliotransmitter release
        release synapse;       // Synaptic exocytosis
        release astro;         // Gliotransmitter release (if spikes are provided)
        double y[2], dy[4], dfdy[16];  // Continuous variables // There is a minor overuse of memory when NEQ<3 since
                                       // it should just be dy[2] and dfdy[4]. But allocating on stack minimizes leaks
        asn(int Neq_,int Nv_,int Nsyn_, long int Ns_,long int Ng_) : NEQ((size_t)Neq_),
                                                                     NVAR((size_t)Nv_),
                                                                     N_syn((size_t)Nsyn_),
                                                                     N_spk(Ns_),
                                                                     N_gre(Ng_),
                                                                     synapse(Nv_,Nsyn_,Ns_),
                                                                     astro(1,1,Ng_){};
        ~asn(){};
        void set_pars(PyObject* pardict);
        void set_ics(PyObject* pardict);
        // RHS (deterministic version)
        void ode(double t,double dt);      //ODE (used for integration) (standard form)
        int ode_gsl(double t, const double y_[], double dy_[]); // GSL wrapper
        static int ode_gsl_wrapper(double t, const double y_[], double dy_[], void * params);
        // Jacobian
        void jacobian(double t, double dt);     // Jacobian
        int jacobian_gsl(double t, const double y_[], double * dfdy_, double dfdt_[]);
        static int jacobian_gsl_wrapper(double t, const double y_[], double * dfdy_, double dfdt_[], void * params);
        // Because we are using gtrelease class to simulate calcium-dependent Gt. release, we need a wrapper around this
        // class to access the static routines of gtrelease::calcium responsible of the noisy part
        static int ode_s_wrapper(double t,double dt,const double y_[], double dy_[], void *params);
        static int ode_sd_wrapper(double t,double dt,const double y_[], double dy_[], void *params);
        // Simulation routines
        out_asn simulate(PyObject* pardict,PyObject* solver_options,double *spikes, double *gres);
};

// Mean field tripartite synapse model
class asn_ave{
    size_t NVAR, NEQ;
    public:
        asn_pars pars;                 // Model parameters
        release_ave synapse;           // Synaptic dynamics
        release_ave gliot;          // Gt. dynamics
        double y[6], dy[6], dfdy[36];  // State variables, derivatives and Jacobian are of size specified by NVAR
        asn_ave(int Nv_) : NVAR((size_t)Nv_),
                           synapse(Nv_),
                           gliot(1){NEQ = (NVAR+1)+3;};
        ~asn_ave(){};
        void set_pars(PyObject* pardict);
        void set_ics(PyObject* pardict);
        // RHS (deterministic version)
        void ode(double t,double dt);  // ODE (used for integration) (standard form)
        int ode_gsl(double t, const double y_[], double dy_[]); // GSL wrapper
        static int ode_gsl_wrapper(double t, const double y_[], double dy_[], void * params);
        // Jacobian
        void jacobian(double t, double dt);     // Jacobian
        int jacobian_gsl(double t, const double y_[], double * dfdy_, double dfdt_[]);
        static int jacobian_gsl_wrapper(double t, const double y_[], double * dfdy_, double dfdt_[], void * params);
        // Simulator
        out_asn simulate(double rate_s,double rate_a,PyObject *pardict,PyObject* solver_options);
};

#endif