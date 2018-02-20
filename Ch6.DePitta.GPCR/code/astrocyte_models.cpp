// Compile as
// gcc -std=c++11 -lpython2.7 -lstdc++ -lm -lgsl -lgslcblas -I/usr/include/python2.7 -I/usr/include/numpy -I/home/maurizio/Dropbox/Ongoing.Projects/pycustommodules -I/home/maurizio/Dropbox/Ongoing.Projects/pycustommodules/solvers -I/home/maurizio/Editorial.Projects/Chapter.DePitta.GPCR/code /home/maurizio/Dropbox/Ongoing.Projects/pycustommodules/pycapi_utils.cpp /home/maurizio/Editorial.Projects/Chapter.DePitta.GPCR/code/astro_models.cpp -o test.o

#include "astrocyte_models.h"

/*-----------------------------------------------------------------------------
Definitions
-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------
LRA model
-----------------------------------------------------------------------------*/
// LRA Model Constructor
lra::lra(){
    pars.d1     = 0.1;
    pars.d2     = 2.1;
    pars.d3     = 0.9967;
    pars.d5     = 0.2;
    pars.a2     = 0.4;
    pars.c1     = 0.4;
    pars.c0     = 4.0;
    pars.rc     = 7.0;
    pars.rl     = 0.05;
    pars.ver    = 0.9;
    pars.ker    = 0.1;
    pars.ip3    = 0.1;
}

// LRA Parameter Initializer (par structure)
void lra::set_pars(lra_pars parameters){
    pars = parameters;
}

// LRA Parameter Initializer (by PyDictionary)
void lra::set_pars(PyObject* pardict){
    //Assumes that entries of pardict are floats
    pars.d1     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"d1"));
    pars.d2     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"d2"));
    pars.d3     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"d3"));
    pars.d5     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"d5"));
    pars.a2     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"a2"));
    pars.c1     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"c1"));
    pars.c0     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"c0"));
    pars.rc     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"rc"));
    pars.rl     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"rl"));
    pars.ver    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"ver"));
    pars.ker    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"Ker"));
    pars.ip3    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"ip3"));
}

void lra::set_ics(double* ics){
    for(int i=0;i<NEQ;i++) y[i] = ics[i];
}

void lra::set_ics(PyObject* pardict){
    double* ics = (double*)((PyArrayObject*)PyDict_GetItemString(pardict,"ICs"))->data;
    for(int i=0;i<NEQ;i++) y[i] = ics[i];
}

double lra::bias(double t,double dt){
    //CURRENTLY NOT IMPLEMENTED
    //NOTE: t,dt are inherited from the integrator.
    //t in case a function of time is specified.
    //dt in order to deduce depending on time, the counter in case a vector signal
    //is passed.
    return 0;
}

void lra::ode(double t,double dt){
    double I = pars.ip3 + bias(t,dt);
    dy[0] = (pars.rc*pow(HILL(I,pars.d1,1)*HILL(y[0],pars.d5,1)*y[1],3)+pars.rl)*(pars.c0-(1+pars.c1)*y[0])-pars.ver*HILL(y[0],pars.ker,2);
    dy[1] = pars.a2*pars.d2*(I+pars.d1)/(I+pars.d3)*(1.00-y[1])-pars.a2*y[0]*y[1];
}

int lra::ode_gsl(double t, const double y_[], double dy_[]){
        int i;

    // Assign to internal variable
    for(i=0;i<NEQ;i++) {
        y[i] = y_[i];
    }

    // Update internal state and compute derivatives
    lra::ode(t,0);

    // Assign to solver output
    for(i=0;i<NEQ;i++) {
        dy_[i] = dy[i];
    }

    return GSL_SUCCESS;
}

// GSL wrapper (we use the void parameter to pass the pointer to function!)
int lra::ode_gsl_wrapper(double t, const double y_[], double dy_[], void * params){
    assert(params);
    return(static_cast<lra*>(params)->ode_gsl(t,y_,dy_));
}

/*-----------------------------------------------------------------------------
//Output data struct
-----------------------------------------------------------------------------*/
out_lra::out_lra(){};

out_lra::out_lra(long int NUMEL){
    N = NUMEL;
    // Allocate pointers
    ts  = (double*)calloc(N,sizeof(double));
    ca  = (double*)calloc(N,sizeof(double));
    h   = (double*)calloc(N,sizeof(double));
}

PyObject* out_lra::make_PyDict(){
    //WARNING: this method destroy the class that it belongs to
    //Declare list of output arguments
    kvlist oitems;
    PyObject* outdict;
    oitems.push_back(std::make_tuple ("ts",array_double_to_pyobj(ts,N)));
    oitems.push_back(std::make_tuple ("ca",array_double_to_pyobj(ca,N)));
    oitems.push_back(std::make_tuple ("h",array_double_to_pyobj(h,N)));
    //Build PyDict
    outdict = build_PyDict(oitems);
    return outdict;
}

/*-----------------------------------------------------------------------------
CHI model
-----------------------------------------------------------------------------*/
// CHI Model Constructor
chi::chi(){
    pars.vbias  = 0.0;
    pars.vbeta  = 0.0;
    pars.vdelta = 0.5;
    pars.kappad = 1.0;
    pars.kdelta = 0.5;
    pars.v3k    = 2.0;
    pars.kd     = 0.5;
    pars.k3     = 1.0;
    pars.r5p    = 0.0;
    pars.d1     = 0.1;
    pars.d2     = 2.1;
    pars.d3     = 0.9967;
    pars.d5     = 0.2;
    pars.a2     = 0.4;
    pars.c1     = 0.4;
    pars.c0     = 4.0;
    pars.rc     = 7.0;
    pars.rl     = 0.05;
    pars.ver    = 0.9;
    pars.ker    = 0.1;
}

// CHI Parameter Initializer (par structure)
void chi::set_pars(chi_pars parameters){
    pars = parameters;
}

// CHI Parameter Initializer (by PyDictionary)
void chi::set_pars(PyObject* pardict){
    //Assumes that entries of pardict are floats
    pars.vbias  = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"vbias"));
    pars.vbeta  = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"vbeta"));
    pars.vdelta = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"vdelta"));
    pars.kappad = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"kappad"));
    pars.kdelta = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"Kdelta"));
    pars.v3k    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"v3k"));
    pars.kd     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"Kd"));
    pars.k3     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"K3"));
    pars.r5p    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"r5p"));
    pars.d1     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"d1"));
    pars.d2     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"d2"));
    pars.d3     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"d3"));
    pars.d5     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"d5"));
    pars.a2     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"a2"));
    pars.c1     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"c1"));
    pars.c0     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"c0"));
    pars.rc     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"rc"));
    pars.rl     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"rl"));
    pars.ver    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"ver"));
    pars.ker    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"Ker"));
}

void chi::set_ics(double* ics){
    for(int i=0;i<NEQ;i++) y[i] = ics[i];
}

void chi::set_ics(PyObject* pardict){
    double* ics = (double*)((PyArrayObject*)PyDict_GetItemString(pardict,"ICs"))->data;
    for(int i=0;i<NEQ;i++) y[i] = ics[i];
}

double chi::bias(double t,double dt){
    //CURRENTLY NOT IMPLEMENTED
    //NOTE: t,dt are inherited from the integrator.
    //t in case a function of time is specified.
    //dt in order to deduce depending on time, the counter in case a vector signal
    //is passed.
    return 0;
}

void chi::ode(double t,double dt){
    dy[0] = bias(t,dt) + pars.vbeta + pars.vdelta/HILL(pars.kappad,y[0],1)*HILL(y[1],pars.kdelta,2)-pars.v3k*HILL(y[1],pars.kd,4)*HILL(y[0],pars.k3,1)-pars.r5p*y[0];
    dy[1] = (pars.rc*pow(HILL(y[0],pars.d1,1)*HILL(y[1],pars.d5,1)*y[2],3)+pars.rl)*(pars.c0-(1+pars.c1)*y[1])-pars.ver*HILL(y[1],pars.ker,2);
    dy[2] = pars.a2*pars.d2*(y[0]+pars.d1)/(y[0]+pars.d3)*(1.00-y[2])-pars.a2*y[1]*y[2];
}

int chi::ode_gsl(double t, const double y_[], double dy_[]){
        int i;

    // Assign to internal variable
    for(i=0;i<NEQ;i++) {
        y[i] = y_[i];
    }

    // Update internal state and compute derivatives
    chi::ode(t,0);

    // Assign to solver output
    for(i=0;i<NEQ;i++) {
        dy_[i] = dy[i];
    }

    return GSL_SUCCESS;
}

// GSL wrapper (we use the void parameter to pass the pointer to function!)
int chi::ode_gsl_wrapper(double t, const double y_[], double dy_[], void * params){
    assert(params);
    return(static_cast<chi*>(params)->ode_gsl(t,y_,dy_));
}

/*-----------------------------------------------------------------------------
//Output data struct
-----------------------------------------------------------------------------*/
out_chi::out_chi(){};

out_chi::out_chi(long int NUMEL){
    N = NUMEL;
    // Allocate pointers
    ts  = (double*)calloc(N,sizeof(double));
    ip3 = (double*)calloc(N,sizeof(double));
    ca  = (double*)calloc(N,sizeof(double));
    h   = (double*)calloc(N,sizeof(double));
}

PyObject* out_chi::make_PyDict(){
    //WARNING: this method destroy the class that it belongs to
    //Declare list of output arguments
    kvlist oitems;
    PyObject* outdict;
    oitems.push_back(std::make_tuple ("ts",array_double_to_pyobj(ts,N)));
    oitems.push_back(std::make_tuple ("ip3",array_double_to_pyobj(ip3,N)));
    oitems.push_back(std::make_tuple ("ca",array_double_to_pyobj(ca,N)));
    oitems.push_back(std::make_tuple ("h",array_double_to_pyobj(h,N)));
    //Build PyDict
    outdict = build_PyDict(oitems);
    return outdict;
}

/*-----------------------------------------------------------------------------
LRA SIMULATOR
-----------------------------------------------------------------------------*/
// GSL wrapper
out_lra lra_simulator_gsl(PyObject* pars,PyObject* solver_options){
    //Create chi astrocyte
    lra astro;
    //Set parameters
    astro.set_pars(pars);
    //Set ICs
    astro.set_ics(pars);
    // Solver options
    set_solver_gsl opts(solver_options);

    // Instantiate GSL solver
    gsl_odeiv2_system sys = {&astro.ode_gsl_wrapper, nullptr, astro.NEQ, reinterpret_cast<void *>(std::addressof(astro))};
    gsl_odeiv2_driver * d;
    d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_msadams, opts.dt, opts.atol, opts.rtol);

    double ti, t = opts.t0;
    int status; // Integrator status

    //Prepare output
    out_lra out(opts.nstep); //Creates output of the simulator (allowed in C++ and C99)

    int i = 0,j;
    do{
        ti = t + opts.dt; // t is automatically updated by the driver
        // Produce output (since we are running in a do cycle, the output must be assigned
        // before the next integration step)
        out.ts[i]  = ti;
        out.ca[i]  = astro.y[0];
        out.h[i]   = astro.y[1];
        // Integration step
        status = gsl_odeiv2_driver_apply (d, &t, ti, astro.y);
        // Error checking
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            break;
        }
        i++;
    }while(i < opts.nstep);

    // Save only last point of solution
    // Free memory
    gsl_odeiv2_driver_free (d);

    return out;
}

/*-----------------------------------------------------------------------------
CHI SIMULATOR
-----------------------------------------------------------------------------*/
// GSL wrapper
out_chi chi_simulator_gsl(PyObject* pars,PyObject* solver_options){
    //Create chi astrocyte
    chi astro;
    //Set parameters
    astro.set_pars(pars);
    //Set ICs
    astro.set_ics(pars);
    // Solver options
    set_solver_gsl opts(solver_options);

    // Instantiate GSL solver
    gsl_odeiv2_system sys = {&astro.ode_gsl_wrapper, nullptr, astro.NEQ, reinterpret_cast<void *>(std::addressof(astro))};
    gsl_odeiv2_driver * d;
    d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_msadams, opts.dt, opts.atol, opts.rtol);

    double ti, t = opts.t0;
    int status; // Integrator status

    //Prepare output
    out_chi out(opts.nstep); //Creates output of the simulator (allowed in C++ and C99)

    int i = 0,j;
    do{
        ti = t + opts.dt; // t is automatically updated by the driver
        // Produce output (since we are running in a do cycle, the output must be assigned
        // before the next integration step)
        out.ts[i]  = ti;
        out.ip3[i] = astro.y[0];
        out.ca[i]  = astro.y[1];
        out.h[i]   = astro.y[2];
        // Integration step
        status = gsl_odeiv2_driver_apply (d, &t, ti, astro.y);
        // Error checking
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            break;
        }
        i++;
    }while(i < opts.nstep);

    // Save only last point of solution
    // Free memory
    gsl_odeiv2_driver_free (d);

    return out;
}

/*-----------------------------------------------------------------------------
G-CHI model
-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------
Definitions
-----------------------------------------------------------------------------*/
// G-CHI Model Constructor
gchi::gchi(){
    pars.yrel = 0.001;
    pars.vp   = 0.3;
    pars.omp  = 1.8;
    pars.zeta = 0.0;
    pars.kkc  = 0.5;
    pars.vbias  = 0.0;
    pars.vbeta  = 0.0;
    pars.vdelta = 0.5;
    pars.kappad = 1.0;
    pars.kdelta = 0.5;
    pars.v3k    = 2.0;
    pars.kd     = 0.5;
    pars.k3     = 1.0;
    pars.r5p    = 0.0;
    pars.d1     = 0.1;
    pars.d2     = 2.1;
    pars.d3     = 0.9967;
    pars.d5     = 0.2;
    pars.a2     = 0.4;
    pars.c1     = 0.4;
    pars.c0     = 4.0;
    pars.rc     = 7.0;
    pars.rl     = 0.05;
    pars.ver    = 0.9;
    pars.ker    = 0.1;
    pars.stimulus = (double*)calloc(3, sizeof(double));
}

// CHI Parameter Initializer (par structure)
void gchi::set_pars(gchi_pars parameters){
    pars = parameters;
}

// CHI Parameter Initializer (by PyDictionary)
void gchi::set_pars(PyObject* pardict){
    //Assumes that entries of pardict are floats
    pars.yrel   = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"yrel"));
    pars.vp     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"Op"));
    pars.omp    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"OmegaP"));
    pars.zeta   = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"zeta"));
    pars.kkc    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"Kkc"));
    pars.vbias  = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"vbias"));
    pars.vbeta  = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"vbeta"));
    pars.vdelta = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"vdelta"));
    pars.kappad = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"kappad"));
    pars.kdelta = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"Kdelta"));
    pars.v3k    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"v3k"));
    pars.kd     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"Kd"));
    pars.k3     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"K3"));
    pars.r5p    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"r5p"));
    pars.d1     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"d1"));
    pars.d2     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"d2"));
    pars.d3     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"d3"));
    pars.d5     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"d5"));
    pars.a2     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"a2"));
    pars.c1     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"c1"));
    pars.c0     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"c0"));
    pars.rc     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"rc"));
    pars.rl     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"rl"));
    pars.ver    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"ver"));
    pars.ker    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"Ker"));
    // Bias
    pars.stimulus    = (double*)calloc(3, sizeof(double));
    pars.stimulus[0] = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"T"));
    pars.stimulus[1] = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"pw"));
    pars.stimulus[2] = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"yb"));
}

void gchi::set_ics(double* ics){
    for(int i=0;i<NEQ;i++) y[i] = ics[i];
}

void gchi::set_ics(PyObject* pardict){
    double* ics = (double*)((PyArrayObject*)PyDict_GetItemString(pardict,"ICs"))->data;
    for(int i=0;i<NEQ;i++) y[i] = ics[i];
}

double gchi::bias(double t,double dt){
    //CURRENTLY NOT IMPLEMENTED
    //NOTE: t,dt are inherited from the integrator.
    //t in case a function of time is specified.
    //dt in order to deduce depending on time, the counter in case a vector signal
    //is passed.
    // signal = [T,d,I,ybias]
    if(fmod(t,pars.stimulus[0])<pars.stimulus[1]){
        return pars.stimulus[2];
    }
    else{
        return 0.0;
    }
}

void gchi::ode(double t,double dt){
    dy[0] = pars.vbias + pars.vp*(pars.yrel + bias(t,dt))*(1-y[0]) - pars.omp*(1.0+pars.zeta*HILL(y[2],pars.kkc,1))*y[0];
    dy[1] = pars.vbeta*y[0] + pars.vdelta/HILL(pars.kappad,y[1],1)*HILL(y[2],pars.kdelta,2)-pars.v3k*HILL(y[2],pars.kd,4)*HILL(y[1],pars.k3,1)-pars.r5p*y[1];
    dy[2] = (pars.rc*pow(HILL(y[1],pars.d1,1)*HILL(y[2],pars.d5,1)*y[3],3)+pars.rl)*(pars.c0-(1+pars.c1)*y[2])-pars.ver*HILL(y[2],pars.ker,2);
    dy[3] = pars.a2*pars.d2*(y[1]+pars.d1)/(y[1]+pars.d3)*(1.00-y[3])-pars.a2*y[2]*y[3];
}

int gchi::ode_gsl(double t, const double y_[], double dy_[]){

    int i;

    // Assign to internal variable
    for(i=0;i<NEQ;i++) {
        y[i] = y_[i];
    }

    // Update internal state and compute derivatives
    gchi::ode(t,0);

    // Assign to solver output
    for(i=0;i<NEQ;i++) {
        dy_[i] = dy[i];
    }

    return GSL_SUCCESS;
}

// GSL wrapper needed to cast as pointer-to-function a pointer to member function.
int gchi::ode_gsl_wrapper(double t, const double y_[], double dy_[], void * params){
    assert(params);
    return(static_cast<gchi*>(params)->ode_gsl(t,y_,dy_));
}

/*-----------------------------------------------------------------------------
//Output data struct
-----------------------------------------------------------------------------*/
out_gchi::out_gchi(){};

out_gchi::out_gchi(long int NUMEL){
    N = NUMEL;
    // Allocate pointers
    ts  = (double*)calloc(N,sizeof(double));
    rec = (double*)calloc(N,sizeof(double));
    ip3 = (double*)calloc(N,sizeof(double));
    ca  = (double*)calloc(N,sizeof(double));
    h   = (double*)calloc(N,sizeof(double));
}

PyObject* out_gchi::make_PyDict(){
    init_numpy();
    //WARNING: this method destroy the class that it belongs to
    //Declare list of output arguments
    kvlist oitems;
    PyObject* outdict;
    oitems.push_back(std::make_tuple ("ts",array_double_to_pyobj(ts,N)));
    oitems.push_back(std::make_tuple ("rec",array_double_to_pyobj(rec,N)));
    oitems.push_back(std::make_tuple ("ip3",array_double_to_pyobj(ip3,N)));
    oitems.push_back(std::make_tuple ("ca",array_double_to_pyobj(ca,N)));
    oitems.push_back(std::make_tuple ("h",array_double_to_pyobj(h,N)));
    //Build PyDict
    outdict = build_PyDict(oitems);
    //out_mhvastro_e::~out_mhvastro_e();  //MIGHT GIVE PROBLEMS
//    delete this;
    return outdict;
}

/*-----------------------------------------------------------------------------
G-CHI SIMULATOR
-----------------------------------------------------------------------------*/
// GSL wrapper
out_gchi gchi_simulator_gsl(PyObject* pars, PyObject* solver_options){
    //Create chi astrocyte
    gchi astro;
    //Set parameters
    astro.set_pars(pars);
    //Set ICs
    astro.set_ics(pars);
    // Solver options
    set_solver_gsl opts(solver_options);

    // Instantiate GSL solver
    gsl_odeiv2_system sys = {&astro.ode_gsl_wrapper, nullptr, astro.NEQ, reinterpret_cast<void *>(std::addressof(astro))};
    gsl_odeiv2_driver * d;
    d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_msadams, opts.dt, opts.atol, opts.rtol);


    double ti, t = opts.t0;
    int status; // Integrator status

    //Prepare output
    out_gchi out(opts.nstep); //Creates output of the simulator (allowed in C++ and C99)

    int i = 0,j;
    do{
        ti = t + opts.dt; // t is automatically updated by the driver
        // Produce output (since we are running in a do cycle, the output must be assigned
        // before the next integration step)
        out.ts[i]  = ti;
        out.rec[i] = astro.y[0];
        out.ip3[i] = astro.y[1];
        out.ca[i]  = astro.y[2];
        out.h[i]   = astro.y[3];
        // Integration step
        status = gsl_odeiv2_driver_apply (d, &t, ti, astro.y);
        // Error checking
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            break;
        }
        i++;
    }while(i < opts.nstep);

    // Save only last point of solution
    // Free memory
    gsl_odeiv2_driver_free (d);

    return out;
}

/*-----------------------------------------------------------------------------
EXTENDED G-CHI model
-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------
Definitions
-----------------------------------------------------------------------------*/
// eG-CHI Model Constructor
egchi::egchi(){
    pars.yrel = 0.001;
    pars.vp   = 0.3;
    pars.omp  = 1.8;
//    pars.vkp  = 0.5;
//    pars.rtot = 5.0;
//    pars.vk   = 1.0;
//    pars.kkd  = 1.0;
    pars.vkd  = 0.5;
    pars.kkc  = 0.5;
    pars.vk   = 1.0;
    pars.omkd = 2.5;

    pars.vd   = 0.5;
    pars.kdc  = 0.3;
    pars.kdd  = 1.0;
    pars.omd  = 0.1;
//    pars.alpha= 0.0;
//    pars.beta = 0.0;
    pars.vbias  = 0.0;
    pars.vbeta  = 0.0;
    pars.vdelta = 0.5;
    pars.kappad = 1.0;
    pars.kdelta = 0.5;
    pars.v3k    = 2.0;
    pars.kd     = 0.5;
    pars.k3     = 1.0;
    pars.r5p    = 0.0;
    pars.d1     = 0.1;
    pars.d2     = 2.1;
    pars.d3     = 0.9967;
    pars.d5     = 0.2;
    pars.a2     = 0.4;
    pars.c1     = 0.4;
    pars.c0     = 4.0;
    pars.rc     = 7.0;
    pars.rl     = 0.05;
    pars.ver    = 0.9;
    pars.ker    = 0.1;
}

// CHI Parameter Initializer (par structure)
void egchi::set_pars(egchi_pars parameters){
    pars = parameters;
}

// CHI Parameter Initializer (by PyDictionary)
void egchi::set_pars(PyObject* pardict){
    //Assumes that entries of pardict are floats
    pars.yrel   = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"yrel"));
    pars.vp     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"Op"));
    pars.omp    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"OmegaP"));
//    pars.vkp    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"vkp"));
//    pars.rtot   = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"rtot"));
    pars.vkd    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"vkd"));
    pars.kkc    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"Kkc"));
    pars.vk     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"vk"));
    pars.omkd   = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"OmegaKD"));
//    pars.kkd    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"Kkd"));
    pars.vd     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"vd"));
    pars.kdc    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"Kdc"));
    pars.kdd    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"Kdd"));
    pars.omd    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"OmegaD"));
//    pars.alpha  = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"alpha"));
//    pars.beta   = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"beta"));
    pars.vbias  = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"vbias"));
    pars.vbeta  = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"vbeta"));
    pars.vdelta = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"vdelta"));
    pars.kappad = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"kappad"));
    pars.kdelta = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"Kdelta"));
    pars.v3k    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"v3k"));
    pars.kd     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"Kd"));
    pars.k3     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"K3"));
    pars.r5p    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"r5p"));
    pars.d1     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"d1"));
    pars.d2     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"d2"));
    pars.d3     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"d3"));
    pars.d5     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"d5"));
    pars.a2     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"a2"));
    pars.c1     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"c1"));
    pars.c0     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"c0"));
    pars.rc     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"rc"));
    pars.rl     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"rl"));
    pars.ver    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"ver"));
    pars.ker    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"Ker"));
}

void egchi::set_ics(double* ics){
    for(int i=0;i<NEQ;i++) y[i] = ics[i];
}

void egchi::set_ics(PyObject* pardict){
    double* ics = (double*)((PyArrayObject*)PyDict_GetItemString(pardict,"ICs"))->data;
    for(int i=0;i<NEQ;i++) y[i] = ics[i];
}

double egchi::bias(double t,double dt){
    //CURRENTLY NOT IMPLEMENTED
    //NOTE: t,dt are inherited from the integrator.
    //t in case a function of time is specified.
    //dt in order to deduce depending on time, the counter in case a vector signal
    //is passed.
    return 0.0;
}

//void egchi::ode(double t,double dt){
//    dy[0] = pars.vbias + pars.vp*(pars.yrel + bias(t,dt))*(1-y[0]) - pars.omp*(1.0+pars.vkp*y[5]/pars.omp)*y[0];
//    dy[1] = pars.vbeta*y[0] + pars.vdelta/HILL(pars.kappad,y[1],1)*HILL(y[2],pars.kdelta,2)-pars.v3k*HILL(y[2],pars.kd,4)*HILL(y[1],pars.k3,1)-pars.r5p*y[1];
//    dy[2] = (pars.rc*pow(HILL(y[1],pars.d1,1)*HILL(y[2],pars.d5,1)*y[3],3)+pars.rl)*(pars.c0-(1+pars.c1)*y[2])-pars.ver*HILL(y[2],pars.ker,2);
//    dy[3] = pars.a2*pars.d2*(y[1]+pars.d1)/(y[1]+pars.d3)*(1.00-y[3])-pars.a2*y[2]*y[3];
//    dy[4] = pars.vbeta*y[0] + pars.vdelta/HILL(pars.kappad,y[1],1)*HILL(y[2],pars.kdelta,2) - pars.vk*HILL(y[2],pars.kkc,1)*HILL(y[4],pars.kkd,1) - pars.vd*HILL(y[2],pars.kdc,2)*HILL(y[4],pars.kdd,2) -pars.omd*y[4];
//    dy[5] = pars.vk*HILL(y[2],pars.kkc,1)*HILL(y[4],pars.kkd,1) - pars.vkp*pars.rtot*y[0]*y[5];
//}

void egchi::ode(double t,double dt){
    dy[0] = pars.vbias + pars.vp*(pars.yrel + bias(t,dt))*(1-y[0]) - pars.omp*(1.0+pars.vk*y[5]/pars.omp)*y[0];
    dy[1] = pars.vbeta*y[0] + pars.vdelta/HILL(pars.kappad,y[1],1)*HILL(y[2],pars.kdelta,2)-pars.v3k*HILL(y[2],pars.kd,4)*HILL(y[1],pars.k3,1)-pars.r5p*y[1];
    dy[2] = (pars.rc*pow(HILL(y[1],pars.d1,1)*HILL(y[2],pars.d5,1)*y[3],3)+pars.rl)*(pars.c0-(1+pars.c1)*y[2])-pars.ver*HILL(y[2],pars.ker,2);
    dy[3] = pars.a2*pars.d2*(y[1]+pars.d1)/(y[1]+pars.d3)*(1.00-y[3])-pars.a2*y[2]*y[3];
    dy[4] = pars.vbeta*y[0] + pars.vdelta/HILL(pars.kappad,y[1],1)*HILL(y[2],pars.kdelta,2) - pars.vkd*y[4]*HILL(y[2],pars.kkc,1) - pars.vd*HILL(y[2],pars.kdc,2)*HILL(y[4],pars.kdd,2) -pars.omd*y[4];
    dy[5] = pars.vkd*y[4]*HILL(y[2],pars.kkc,1) - pars.omkd*y[5];
}

void egchi::jacobian(double t, double dt){
    // Row-ordered Jacobian, according to y-ordered variables
    // Gamma
//    dfdy[0] = -(pars.vp*pars.yrel + pars.omp + pars.vkp*y[5]);
//    dfdy[1] = 0.0;
//    dfdy[2] = 0.0;
//    dfdy[3] = 0.0;
//    dfdy[4] = 0.0;
//    dfdy[5] = -pars.vkp*y[0];
//
//    // I
//    dfdy[6] = pars.vbeta;
//    dfdy[7] = -pars.vdelta*HILL(y[2],pars.kdelta,2)*Hx(y[1],pars.kappad,1) - pars.v3k*HILL(y[2],pars.kd,4)*Hx(y[1],pars.k3,1)-pars.r5p;
//    dfdy[8] = pars.vdelta*(1.0-HILL(y[1],pars.kappad,1))*Hx(y[2],pars.kdelta,2) - pars.v3k*Hx(y[2],pars.kd,4)*HILL(y[1],pars.k3,1);
//    dfdy[9] = 0.0;
//    dfdy[10] = 0.0;
//    dfdy[11] = 0.0;
//
//    // C
//    dfdy[12] = 0.0;
//    dfdy[13] = 3.0*pars.rc*pow(HILL(y[1],pars.d1,1)*HILL(y[2],pars.d5,1),2)*pow(y[3],3)*Hx(y[1],pars.d1,1)*HILL(y[2],pars.d5,1)*(pars.c0-(1+pars.c1)*y[2]);
//    dfdy[14] = 3.0*pars.rc*pow(HILL(y[1],pars.d1,1)*HILL(y[2],pars.d5,1),2)*pow(y[3],3)*HILL(y[1],pars.d1,1)*Hx(y[2],pars.d5,1)*(pars.c0-(1+pars.c1)*y[2]) - (pars.rc*pow(HILL(y[1],pars.d1,1)*HILL(y[2],pars.d5,1)*y[3],3)+pars.rl)*(1+pars.c1) - pars.ver*Hx(y[2],pars.ker,2);
//    dfdy[15] = 3.0*pars.rc*pow(HILL(y[1],pars.d1,1)*HILL(y[2],pars.d5,1),3)*pow(y[3],2)*(pars.c0-(1+pars.c1)*y[2]);
//    dfdy[16] = 0.0;
//    dfdy[17] = 0.0;
//
//    // h
//    dfdy[18] = 0.0;
//    dfdy[19] = pars.a2*Q2x(y[1],pars.d1,pars.d2,pars.d3)*(1.00-y[3]);
//    dfdy[20] = -pars.a2*y[3];
//    dfdy[21] = -pars.a2*(Q2(y[1],pars.d1,pars.d2,pars.d3)+y[2]);
//    dfdy[22] = 0.0;
//    dfdy[23] = 0.0;
//
//    // D
//    dfdy[24] = pars.vbeta;
//    dfdy[25] = -pars.vdelta*HILL(y[2],pars.kdelta,2)*Hx(y[1],pars.kappad,1);
//    dfdy[26] = pars.vdelta*(1.0-HILL(y[1],pars.kappad,1))*Hx(y[2],pars.kdelta,2) - pars.vk*Hx(y[2],pars.kkc,1)*HILL(y[4],pars.kkd,1) - pars.vd*Hx(y[2],pars.kdc,2)*HILL(y[4],pars.kdd,2);
//    dfdy[27] = 0.0;
//    dfdy[28] = -pars.vk*HILL(y[2],pars.kkc,1)*Hx(y[4],pars.kkd,1) - pars.vd*HILL(y[2],pars.kdc,2)*Hx(y[4],pars.kdd,2) -pars.omd;
//    dfdy[29] = 0.0;
//
//    // P
//    dfdy[30] = -pars.vkp*pars.rtot*y[5];
//    dfdy[31] = 0.0;
//    dfdy[32] = pars.vk*Hx(y[2],pars.kkc,1)*HILL(y[4],pars.kkd,1);
//    dfdy[33] = 0.0;
//    dfdy[34] = pars.vk*HILL(y[2],pars.kkc,1)*Hx(y[4],pars.kkd,1);
//    dfdy[35] = -pars.vkp*pars.rtot*y[0];
}

int egchi::ode_gsl(double t, const double y_[], double dy_[]){

    int i;

    // Assign to internal variable
    for(i=0;i<NEQ;i++) {
        y[i] = y_[i];
//        printf("y[%d]=%f",i,y[i]);
    }

    // Update internal state and compute derivatives
    egchi::ode(t,0.0);

    // Assign to solver output
    for(i=0;i<NEQ;i++) {
        dy_[i] = dy[i];
    }

    return GSL_SUCCESS;
}

int egchi::jacobian_gsl(double t, const double y_[], double * dfdy_, double dfdt_[]){

    int i,j;

    // Assign to internal variable
    for(i=0;i<NEQ;i++) {
        y[i] = y_[i];
    }

    // Update internal Jacobian
    egchi::jacobian(t,0.0);

    // Recast Jacobian
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy_, NEQ, NEQ);
    gsl_matrix * m = &dfdy_mat.matrix;
    for(i=0;i<NEQ;i++){
        for(j=0;j<NEQ;j++){
            gsl_matrix_set (m, i, j, dfdy[i*NEQ + j]);
        }
        dfdt_[i] = 0.0;
    }

    return GSL_SUCCESS;
}


// GSL wrapper needed to cast as pointer-to-function a pointer to member function.
int egchi::ode_gsl_wrapper(double t, const double y_[], double dy_[], void * params){
    assert(params);
    return(static_cast<egchi*>(params)->ode_gsl(t,y_,dy_));
}

int egchi::jacobian_gsl_wrapper(double t, const double y_[], double * dfdy_, double dfdt_[], void * params){
    assert(params);
    return(static_cast<egchi*>(params)->jacobian_gsl(t,y_,dfdy_,dfdt_));
}

/*-----------------------------------------------------------------------------
//Output data struct
-----------------------------------------------------------------------------*/
out_egchi::out_egchi(){};

out_egchi::out_egchi(long int NUMEL){
    N = NUMEL;
    // Allocate pointers
    ts  = (double*)calloc(N,sizeof(double));
    rec = (double*)calloc(N,sizeof(double));
    ip3 = (double*)calloc(N,sizeof(double));
    ca  = (double*)calloc(N,sizeof(double));
    h   = (double*)calloc(N,sizeof(double));
    dag = (double*)calloc(N,sizeof(double));
    pkc = (double*)calloc(N,sizeof(double));
}

PyObject* out_egchi::make_PyDict(){
    init_numpy();
    //WARNING: this method destroy the class that it belongs to
    //Declare list of output arguments
    kvlist oitems;
    PyObject* outdict;
    oitems.push_back(std::make_tuple ("ts",array_double_to_pyobj(ts,N)));
    oitems.push_back(std::make_tuple ("rec",array_double_to_pyobj(rec,N)));
    oitems.push_back(std::make_tuple ("ip3",array_double_to_pyobj(ip3,N)));
    oitems.push_back(std::make_tuple ("ca",array_double_to_pyobj(ca,N)));
    oitems.push_back(std::make_tuple ("h",array_double_to_pyobj(h,N)));
    oitems.push_back(std::make_tuple ("dag",array_double_to_pyobj(dag,N)));
    oitems.push_back(std::make_tuple ("pkc",array_double_to_pyobj(pkc,N)));
    //Build PyDict
    outdict = build_PyDict(oitems);
    //out_mhvastro_e::~out_mhvastro_e();  //MIGHT GIVE PROBLEMS
//    delete this;
    return outdict;
}

/*-----------------------------------------------------------------------------
ExG-CHI SIMULATOR
-----------------------------------------------------------------------------*/
// GSL wrapper
out_egchi egchi_simulator_gsl(PyObject* pars, PyObject* solver_options){
    //Create chi astrocyte
    egchi astro;
    //Set parameters
    astro.set_pars(pars);
    //Set ICs
    astro.set_ics(pars);
    // Solver options
    set_solver_gsl opts(solver_options);

//    // Instantiate GSL solver
    gsl_odeiv2_system sys = {&astro.ode_gsl_wrapper, nullptr, astro.NEQ, reinterpret_cast<void *>(std::addressof(astro))};
    gsl_odeiv2_driver * d;
    d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_msadams, opts.dt, opts.atol, opts.rtol);

//    gsl_odeiv2_system sys = {&astro.ode_gsl_wrapper, &astro.jacobian_gsl_wrapper, astro.NEQ, reinterpret_cast<void *>(std::addressof(astro))};
//    gsl_odeiv2_driver * d;
//    d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_bsimp, opts.dt, opts.atol, opts.rtol);

    double ti, t = opts.t0;
    int status; // Integrator status

    //Prepare output
    out_egchi out(opts.nstep); //Creates output of the simulator (allowed in C++ and C99)

//    printf("\nCheck-2");
    int i = 0,j;
    do{
        ti = t + opts.dt; // t is automatically updated by the driver
        // Produce output (since we are running in a do cycle, the output must be assigned
        // before the next integration step)
        out.ts[i]  = ti;
        out.rec[i] = astro.y[0];
        out.ip3[i] = astro.y[1];
        out.ca[i]  = astro.y[2];
        out.h[i]   = astro.y[3];
        out.dag[i] = astro.y[4];
        out.pkc[i] = astro.y[5];
        // Integration step
        status = gsl_odeiv2_driver_apply (d, &t, ti, astro.y);
        // Error checking
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            break;
        }
        i++;
    }while(i < opts.nstep);

    // Save only last point of solution
    // Free memory
    gsl_odeiv2_driver_free (d);

    return out;
}