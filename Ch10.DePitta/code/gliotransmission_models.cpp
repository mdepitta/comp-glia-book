// Compile as
// gcc -g -std=c++11 -lpython2.7 -lstdc++ -lm -lgsl -lgslcblas -I/usr/include/python2.7 -I/usr/include/numpy -I/home/maurizio/Editorial.Projects/Chapter.DePitta.Gliotransmission/code/pycustommodules -I/home/maurizio/Editorial.Projects/Chapter.DePitta.Gliotransmission/code/pycustommodules/solvers -I/home/maurizio/Editorial.Projects/Chapter.DePitta.Gliotransmission/code/code /home/maurizio/Editorial.Projects/Chapter.DePitta.Gliotransmission/code/pycustommodules/pycapi_utils.cpp /home/maurizio/Editorial.Projects/Chapter.DePitta.Gliotransmission/code/pycustommodules/solvers/solver_options.cpp /home/maurizio/Editorial.Projects/Chapter.DePitta.Gliotransmission/code/pycustommodules/solvers/stochastic_solvers.cpp /home/maurizio/Editorial.Projects/Chapter.DePitta.Gliotransmission/code/code/astrocyte_models.cpp -o test.o

#include "gliotransmission_models.h"

/*-----------------------------------------------------------------------------
Definitions
-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------
Event-handling functions and auxiliary functions for update/integration
-----------------------------------------------------------------------------*/

void clip_h(double *x, double lims[]){
    // lims is a 3-element double array where the first two elements are low_lim, high_lim and the third is the index of
    // the variable in the 'y' in the system to clip
    int index = (int)lims[2];
    x[index] = std::max(lims[0], std::min(x[index], lims[1]));
}

void c_thr(void **x, void **vals){
  // Implement a threshold crossing to detect GRE
  // x = [*c(t-1), *c(t)]; lims = [&cthr,&time,*event]. Where val is 0 or equal to the time instant
  double x0 = *(double*)x[0];
  double x1 = *(double*)x[1];
  double v0 = *(double*)vals[0]; // cthr
//  double v1 = *(double*)vals[1]; // time
  int *v1 = (int*)vals[1]; // event flag
//  if((x[0]-vals[0]<0.0)&&(x[1]-vals[0]>=0)) vals[2] = vals[1];
  if((x0-v0<0.0)&&(x1-v0>=0.0)) *v1 = 1;
}

void clip_thr(void **x, void **vals){
    // Lumps together clipping and c_thr handling
    // First do the threshold handling
    c_thr(x,vals);
    // Second do the clipping
    // This assumes that vals[2] = *lims;
    clip_h((double*)x[1],(double*)vals[2]);
}

void u0_update(double ts,void **x, void **gs_){
    // Update synaptic release probability in the ASN spiking model
    // ts : time of spike
    // x = [(release_pars*)synapse.pars->u0, (asn_pars*)pars->alpha, (double)dt];
    // gs_ = [*gammas]
    release_pars *synpars = (release_pars*)x[0];
    asn_pars *asnpars = (asn_pars*)x[1];
    double dt  = *(double*)x[2];
    int index = (int)floor(ts/dt);
    double *gammas = (double*)gs_[0];
    synpars->u0 = U0(asnpars->synapse.u0,asnpars->alpha,gammas[index]);
}
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

    N_s = 0;
    N_d = NEQ;

    y    = (double*)calloc(NEQ,sizeof(double));
    dy   = (double*)calloc(NEQ,sizeof(double));
    dfdy = (double*)calloc(NEQ*NEQ,sizeof(double));
}

lra::~lra(){
    free(y);
    free(dy);
    free(dfdy);
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
    noise       = (int)PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"noise"));

    if(noise){
        N_s = 1;
        N_d = NEQ-N_s;
    }
}

lra_pars lra::get_pars(){
    lra_pars p = pars;
    return p;
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
    void** voidPtrs = static_cast<void**>(params);
    return(static_cast<lra*>(voidPtrs[0])->ode_gsl(t,y_,dy_));
}

// Jacobian
void lra::jacobian(double t, double dt){
    double I = pars.ip3 + bias(t,dt);
    //C
    dfdy[0] = -(1+pars.c1)*(pars.rc*pow(HILL(I,pars.d1,1)*HILL(y[0],pars.d5,1)*y[1],3)+pars.rl)
              +3*(pars.c0-(1+pars.c1)*y[0])*pars.rc*pow(HILL(I,pars.d1,1)*y[1],3)*pow(HILL(y[0],pars.d5,1),2)*Hx(y[0],pars.d5,1)
              -pars.ver*Hx(y[0],pars.ker,2);
    dfdy[1] = 3*(pars.c0-(1+pars.c1)*y[0])*pars.rc*pow(HILL(I,pars.d1,1)*HILL(y[0],pars.d5,1),3)*pow(y[1],2);
    //h
    dfdy[2] = -pars.a2*y[1];
    dfdy[3] = -pars.a2*pars.d2*(I+pars.d1)/(I+pars.d3) - pars.a2*y[0];
}

int lra::jacobian_gsl(double t, const double y_[], double * dfdy_, double dfdt_[]){
    int i,j;

    // Assign to internal variable
    for(i=0;i<NEQ;i++) {
        y[i] = y_[i];
    }

    // Update internal Jacobian
    lra::jacobian(t,0.0);

    // Recast Jacobian
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy_, NEQ, NEQ);
    gsl_matrix * m = &dfdy_mat.matrix;
    for(i=0;i<NEQ;i++){
        for(j=0;j<NEQ;j++){
            gsl_matrix_set (m, i, j, dfdy[i*NEQ + j]);
        }
        dfdt_[i] = dy[i];
    }

    return GSL_SUCCESS;
}

int lra::jacobian_gsl_wrapper(double t, const double y_[], double * dfdy_, double dfdt_[], void * params){
    assert(params);
    void** voidPtrs = static_cast<void**>(params);
    return(static_cast<lra*>(voidPtrs[0])->jacobian_gsl(t,y_,dfdy_,dfdt_));
}

double lra::tauh(double t,double dt){
    double I = pars.ip3 + bias(t,dt);
    return (I+pars.d3)/(pars.d2*pars.a2*(I+pars.d1)+pars.a2*(I+pars.d3)*y[0]);
}

int lra::ode_s(double t,double dt,const double y_[], double dy_[]){

    // Assign to internal variable
    for(int i=0;i<NEQ;i++) {
        y[i] = y_[i];
    }

    // Update internal state and compute derivatives
    lra::ode(t,dt);

    // Assign to solver output (only h equation is stochastic in this formulation)
    dy_[0] = dy[1]*sqrt(lra::tauh(t,dt));
    return 1;
}

int lra::ode_s_wrapper(double t,double dt,const double y_[], double dy_[], void *params){
    assert(params);
    void** voidPtrs = static_cast<void**>(params);
    return(static_cast<lra*>(voidPtrs[0])->ode_s(t,dt,y_,dy_));
}

int lra::ode_sd(double t,double dt,const double y_[], double dy_[]){
    double I = pars.ip3 + bias(t,dt);
    // Assign to internal variable
    for(int i=0;i<NEQ;i++) {
        y[i] = y_[i];
    }

    // Update internal state and compute derivatives
//    lra::ode(t,dt);

    // Assign to solver output (only h equation is stochastic in this formulation)
    // d(dh/dt)/dh
    dy_[0] = (-pars.a2*pars.d2*(I+pars.d1)/(I+pars.d3) - pars.a2*y[0])*sqrt(lra::tauh(t,dt));
    return 1;
}

int lra::ode_sd_wrapper(double t,double dt,const double y_[], double dy_[], void *params){
    assert(params);
    void** voidPtrs = static_cast<void**>(params);
    return(static_cast<lra*>(voidPtrs[0])->ode_sd(t,dt,y_,dy_));
}

// Testing routine (to be compiled as a standalone)
out_lra lra::simulate(){
    y[1] = 0.9;
    // Set up solver
    // 1. Recast ODEs in a way that is accepted by the stochastic integrator
    sde astro_model(N_d,N_s,y,
                    this->ode_gsl_wrapper,
                    this->jacobian_gsl_wrapper,
                    this->ode_s_wrapper,
                    this->ode_sd_wrapper);
    double *y__ = astro_model.y; // Access the internal variable of the integrator to do the clipping
                                 // NOTE: I did not find so far a better solution to work directly with this->y.
    // 2. Define integrator instance
    de_solver integrator(&astro_model);
    // 3. Solver options
    set_solver_gsl opts;
    // 4. Allocate solution
    out_lra out(opts.nstep); //Creates output of the simulator (allowed in C++ and C99)
    double **sol = new double*[NEQ];
    sol[0] = out.ca;
    sol[1] = out.h;

    // Integrate
    double lims[3] = {0.,1.,1};
    integrator.integrate(opts,static_cast<void*>(this),&clip_h,y__,lims,out.ts,sol);

    // Free memory allocated for sol
    for(int i=0;i<NEQ;i++) delete sol[i];
    delete[] sol;
}

out_lra lra::simulate(PyObject* pars,PyObject* solver_options){
    //Set parameters
    lra::set_pars(pars);
    //Set ICs
    lra::set_ics(pars);

    // Set up solver
    // 1. Recast ODEs in a way that is accepted by the stochastic integrator
    sde astro_model(N_d,N_s,y,
                    this->ode_gsl_wrapper,
                    this->jacobian_gsl_wrapper,
                    this->ode_s_wrapper,
                    this->ode_sd_wrapper);
    double *y__ = astro_model.y; // Access the internal variable of the integrator to do the clipping
                                 // NOTE: I did not find so far a better solution to work directly with this->y.
    // 2. Define integrator instance
    de_solver integrator(&astro_model);
    // 3. Solver options
    set_solver_gsl opts(solver_options);
    // 4. Allocate solution
    out_lra out(opts.nstep); //Creates output of the simulator (allowed in C++ and C99)
    double **sol = new double*[NEQ];
    sol[0] = out.ca;
    sol[1] = out.h;

    // Integrate
    double lims[3] = {0.,1.,1};
    integrator.integrate(opts,static_cast<void*>(this),&clip_h,y__,lims,out.ts,sol);

    // Free memory allocated for sol
    delete[] sol; // Delete 'sol' but not the blocks of memory pointed by elements in sol since those are in out.ca, and out.h

    return out;
}

/*-----------------------------------------------------------------------------
//Output data struct
-----------------------------------------------------------------------------*/
out_lra::out_lra(long int NUMEL){
    N = NUMEL;
    // Allocate pointers
    ts  = (double*)calloc(N,sizeof(double));
    ca  = (double*)calloc(N,sizeof(double));
    h   = (double*)calloc(N,sizeof(double));
}

PyObject* out_lra::make_PyDict(){
    init_numpy();
    //WARNING: this method destroy the class that it belongs to
    //Declare list of output arguments
    kvlist oitems;
    PyObject* outdict;
    oitems.push_back(std::make_tuple ("t",array_double_to_pyobj(ts,N)));
    oitems.push_back(std::make_tuple ("ca",array_double_to_pyobj(ca,N)));
    oitems.push_back(std::make_tuple ("h",array_double_to_pyobj(h,N)));
    //Build PyDict
    outdict = build_PyDict(oitems);
    return outdict;
}

/*-----------------------------------------------------------------------------
//Exocytosis models
-----------------------------------------------------------------------------*/
release::release(int Nv_,int Nsyn_, long int Nspk_){
    NVAR = (size_t)Nv_;
    N_syn = (size_t)Nsyn_;
    N_spk = Nspk_;
    int i;

    // Allocate last_spike where each element stands for the last spike of each synapse
    last_spike = (double*)calloc(N_syn,sizeof(double));
    // Container for index of last spike index
    spike_index = (long int*)calloc(N_syn,sizeof(long int));
    // Variables
    if(NVAR>1) u = (double*)calloc(N_spk,sizeof(double));
    x   = (double*)calloc(N_spk,sizeof(double));
    y   = (double*)calloc(N_spk,sizeof(double));
    allocate_memory = 1;

    // Set default parameters
    pars.u0   = 0.5;
    pars.taud = 0.5;
    pars.tauf = 1.0;
    pars.taui = 0.02;
    pars.psc0 = 30;  // [pA]

    // Initial conditions may be thoughts as if a spike occurred at t=0 (although one must take into account that
    // only u(0) and y(0) can be saved as they are but not x(0) since for this latter we actually need x(0-). So special
    // handling of this variable's ICs is required in ::simulate
    for(i=0;i<N_syn;i++) last_spike[i] = 0;
}

release::release(int Nv_,int Nsyn_,long int Nspk_,PyObject* pardict){
    int i;
    NVAR = (size_t)Nv_;
    N_syn = (size_t)Nsyn_;
    N_spk = Nspk_;

    //Assumes that entries of pardict are floats
    pars.u0   = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"u0"));
    pars.taud = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"taud"));
    pars.tauf = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"tauf"));
    pars.taui = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"taui"));
    pars.psc0 = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"psc0"));

    // Allocate last_spike where each element stands for the last spike of each synapse
    last_spike = (double*)calloc(N_syn,sizeof(double));
    // Container for index of last spike index
    spike_index = (long int*)calloc(N_syn,sizeof(long int));
    // Variables
    // NOTE it ignores NVAR>1 if tauf is <=0.01 (fast facilitation) and set NVAR to 1
    if((NVAR>1)&&(pars.tauf>0.01))
        u = (double*)calloc(N_spk,sizeof(double));
    else
        NVAR = 1;
    x = (double*)calloc(N_spk,sizeof(double));
    y = (double*)calloc(N_spk,sizeof(double));
    allocate_memory = 1;

    // See ::release() for explanation
    for(i=0;i<N_syn;i++) last_spike[i] = 0;
}

release::~release(){
    if(allocate_memory){
        free(last_spike);
        free(spike_index);
        free(x);
        free(y);
        if(NVAR>1) free(u);
    }
    if(create_ics) free(ics);
}

void release::set_pars(release_pars parameters){
    pars = parameters;
}

void release::set_pars(PyObject *pardict){
    //Assumes that entries of pardict are floats
    pars.u0   = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"u0"));
    pars.taud = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"taud"));
    pars.tauf = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"tauf"));
    pars.taui = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"taui"));
    pars.psc0 = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"psc0"));
}

void release::set_ics(double *ics_){
    // Create ICs
    ics = (double*)calloc((NVAR+1)*N_syn,sizeof(double));
    for(int i=0;i<N_syn;i++){
        ics[i*(NVAR+1)] = ics_[i*(NVAR+1)];                //x
        ics[i*(NVAR+1)+1] = ics_[i*(NVAR+1)+1];            //y
        if(NVAR>1) ics[i*(NVAR+1)+2] = ics_[i*(NVAR+1)+2]; //u
    }
    create_ics = 1;
}

void release::set_ics(PyObject *pardict){
    double *ics_ = (double*)((PyArrayObject*)PyDict_GetItemString(pardict,"ICs"))->data;
    release::set_ics(ics_);
}

void release::get_dims(){
    printf("\nNVAR=%d\tN_syn=%d\tN_spk=%d",NVAR,N_syn,N_spk);
}

void release::x_sol(double x0_,double *x_,double dt_){
    *x_ = 1.0 - exp(-dt_/pars.taud) + x0_*exp(-dt_/pars.taud);
}

void release::uy_sol(double u0_,double *u_,double dt_){
    *u_ = u0_*exp(-dt_/pars.tauf);
}

void release::update(double* state, double ts, int isyn){
    // State is a variable vector that contains (NVAR=1) [x,psc,U0]; (NVAR=2) [x,psc,u]
    // Retrieve values of variables from the last spike
    double xv = state[0];
    double yv = state[1];
    double uv = state[2];
    double un = pars.u0;
    double dt = ts - last_spike[isyn];

    // Update variables
    x[counter] = (1.0-exp(-dt/pars.taud)) + xv*(1.0-uv)*exp(-dt/pars.taud);
    if(NVAR>1){
        u[counter] = (1-pars.u0)*uv*exp(-dt/pars.tauf) + pars.u0;
        un = u[counter];
    }
    y[counter] = yv*exp(-dt/pars.taui) + pars.psc0*x[counter]*un;
}

out_release release::simulate(double *spikes){
    double vars[3];
    double *state = vars;  // Pointer to state variable vector
    int isyn;  // Synaptic index
    int i;     // Generic counter

    // Actual simulation
    for(counter=0;counter<N_spk;counter++){
        isyn = (int)spikes[N_spk+counter];
        // Assign xn, un, yn with the first spike starting from [1,0,0]
        // The first spike to a synapse is detected if the last_spike is negative
        //u
        if(NVAR>1){
            if(last_spike[isyn]<=0)
                state[2] = ics[2+isyn*(NVAR+1)];
            else
                state[2] = u[spike_index[isyn]]; // Last u_n
        }
        else
            state[2] = pars.u0; // (1-eq model)
        //x,y
        if(last_spike[isyn]<=0){
            state[0] = ics[isyn*(NVAR+1)]/(1.0-state[2]);
            state[1] = ics[1+isyn*(NVAR+1)];
            }
        else{
            state[0] = x[spike_index[isyn]]; // Last x_n
            state[1] = y[spike_index[isyn]]; // Last y_n
        }
        // Effective update
        release::update(state,spikes[counter],isyn);
        // Update last spike and spike_index
        last_spike[isyn]  = spikes[counter];
        spike_index[isyn] = counter;
    }

    // Generate output
    out_release sol(N_spk,NVAR);
    for(i=0;i<N_spk;i++){
        sol.x[i] = x[i];
        sol.y[i] = y[i];
    }
    if(NVAR>1)
        for(i=0;i<N_spk;i++) sol.u[i] = u[i];

    return sol;
}

out_release release::simulate(double *spikes, void (*fpa)(double,void**,void**),void **f_var, void **f_par){
    double vars[3];
    double *state = vars;  // Pointer to state variable vector
    int isyn;  // Synaptic index
    int i;     // Generic counter

    // Actual simulation
    for(counter=0;counter<N_spk;counter++){
        isyn = (int)spikes[N_spk+counter];
        // Update u0 if needed
        if(fpa!=nullptr) fpa(spikes[counter],f_var,f_par);
        // Assign xn, un, in with the first spike starting from [1,0,0]
        // The first spike to a synapse is detected if the last_spike is negative
        //u
        if(NVAR>1){
            if(last_spike[isyn]<=0)
                state[2] = ics[2+isyn*(NVAR+1)];
            else
                state[2] = u[spike_index[isyn]]; // Last u_n
        }
        else
            state[2] = pars.u0; // (1-eq model)
        //x,y
        if(last_spike[isyn]<=0){
            // Because the update formula works on x(t^-) we assume that at x(0)=x(0^-)(1-u(0+)) --> x(0-)=x(0+)/(1-u(0+))
            state[0] = ics[isyn*(NVAR+1)]/(1-state[2]);
            state[1] = ics[1+isyn*(NVAR+1)];
            }
        else{
            state[0] = x[spike_index[isyn]]; // Last x_n
            state[1] = y[spike_index[isyn]]; // Last y_n
        }
        // Effective update
        release::update(state,spikes[counter],isyn);
        // Update last spike and spike_index
        last_spike[isyn]  = spikes[counter];
        spike_index[isyn] = counter;
    }

    // Generate output
    out_release sol(N_spk,NVAR);
    for(i=0;i<N_spk;i++){
        sol.x[i] = x[i];
        sol.y[i] = y[i];
    }
    if(NVAR>1)
        for(i=0;i<N_spk;i++) sol.u[i] = u[i];

    return sol;
}

/*-----------------------------------------------------------------------------
//Exocytosis models -- output class
-----------------------------------------------------------------------------*/
out_release::out_release(long int Ns_, int Nv_){
    N_spk = Ns_;
    NVAR = (size_t)Nv_;
    x    = (double*)calloc(N_spk,sizeof(double));
    y    = (double*)calloc(N_spk,sizeof(double));
    if(NVAR>1) u = (double*)calloc(N_spk,sizeof(double));
}

out_release::out_release(long int NUMEL, int Nv_, int Nout_){
    N_spk = NUMEL;  // re-use N_spk to keep info on number of data points
    NVAR  = (size_t)Nv_;
    N_out = Nout_;
    assert((N_out>0)&&(N_out<=NVAR+1)); // Some basic error check
    ts = (double*)calloc(N_spk,sizeof(double));
    x  = (double*)calloc(N_spk,sizeof(double));
    if(N_out>1) y  = (double*)calloc(N_spk,sizeof(double));
    if((NVAR>1)&&(N_out>NVAR)) u = (double*)calloc(N_spk,sizeof(double));
    mean_field = 1;  // Set the MF flag to 1 in order to provide correct output dictionary
}

void out_release::get_dims(){
    printf("\nNVAR=%d\tN_spk=%d",NVAR,N_spk);
}

PyObject* out_release::make_PyDict(){
    init_numpy();
    //WARNING: this method destroy the class that it belongs to
    //Declare list of output arguments
    kvlist oitems;
    PyObject* outdict;
    if(mean_field){
        oitems.push_back(std::make_tuple ("t",array_double_to_pyobj(ts,N_spk)));
        oitems.push_back(std::make_tuple ("x",array_double_to_pyobj(x,N_spk)));
        if(N_out>1) oitems.push_back(std::make_tuple ("y",array_double_to_pyobj(y,N_spk)));
        if((NVAR>1)&&(N_out>NVAR)) oitems.push_back(std::make_tuple ("u",array_double_to_pyobj(u,N_spk)));
    }
    else{
        // DEBUG
        oitems.push_back(std::make_tuple ("x",array_double_to_pyobj(x,N_spk)));
        oitems.push_back(std::make_tuple ("y",array_double_to_pyobj(y,N_spk)));
        if(NVAR>1) oitems.push_back(std::make_tuple ("u",array_double_to_pyobj(u,N_spk)));
    }
    //Build PyDict
    outdict = build_PyDict(oitems);
    return outdict;
}

//*-----------------------------------------------------------------------------
// Average model
//-----------------------------------------------------------------------------*/
release_ave::release_ave(int Nv_){
    NVAR = (size_t)Nv_;
    NEQ  = NVAR+1;

    // Set default parameters
    pars.u0   = 0.5;
    pars.taud = 0.5;
    pars.tauf = 1.0;
    pars.taui = 0.02;
    pars.psc0 = 30;  // [pA]

    // Variables
    y    = (double*)calloc(NEQ,sizeof(double));
    dy   = (double*)calloc(NEQ,sizeof(double));
    dfdy = (double*)calloc(NEQ*NEQ,sizeof(double));
}

release_ave::~release_ave(){
    free(y);
    free(dy);
    free(dfdy);
}

void release_ave::set_pars(PyObject* pardict){
    //Assumes that entries of pardict are floats
    pars.u0   = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"u0"));
    pars.taud = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"taud"));
    pars.tauf = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"tauf"));
    pars.taui = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"taui"));
    pars.psc0 = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"psc0"));
}

void release_ave::set_ics(double *ics){
    for(int i=0;i<NEQ;i++) y[i] = ics[i]; // ICs must be provided in the order of x,y and u (if this latter is given)
}

void release_ave::set_ics(PyObject *pardict){
    double *ics = (double*)((PyArrayObject*)PyDict_GetItemString(pardict,"ICs"))->data;
    release_ave::set_ics(ics);
}

double release_ave::bias(double t,double dt){
    //CURRENTLY NOT IMPLEMENTED
    //NOTE: t,dt are inherited from the integrator.
    //t in case a function of time is specified.
    //dt in order to deduce depending on time, the counter in case a vector signal
    //is passed.
    return 0;
}

void release_ave::ode(double t,double dt){
    double u;
    double f = nu + release_ave::bias(t,dt);
    if(NVAR>1){
        dy[2] = (pars.u0-y[2])/pars.tauf + pars.u0*(1-y[2])*f;
        u = y[2];
    }
    else{
        u = pars.u0;
    }
    dy[0] = (1-y[0])/pars.taud - u*y[0]*f;
    dy[1] = -y[1]/pars.taui + pars.psc0*u*y[0]*f;
}

// GSL wrapper
int release_ave::ode_gsl(double t, const double y_[], double dy_[]){
    int i;

    // Assign to internal variable
    for(i=0;i<NEQ;i++) y[i] = y_[i];

    // Update internal state and compute derivatives
    release_ave::ode(t,0);

    // Assign to solver output
    for(i=0;i<NEQ;i++) dy_[i] = dy[i];

    return GSL_SUCCESS;
}

int release_ave::ode_gsl_wrapper(double t, const double y_[], double dy_[], void * params){
    assert(params);
    void** voidPtrs = static_cast<void**>(params);
    return(static_cast<release_ave*>(voidPtrs[0])->ode_gsl(t,y_,dy_));
}

// Jacobian
void release_ave::jacobian(double t, double dt){
    double u;
    double f = nu + release_ave::bias(t,dt);

    if (NVAR<2){
        u = pars.u0;
        //x
        dfdy[0] = -1.0/pars.taud - u*f;
        dfdy[1] = 0.0;
        //y
        dfdy[2] = pars.psc0*u*f;
        dfdy[3] = -1.0/pars.taui;
    }
    else{
        u = y[2];
        //x
        dfdy[0] = -1.0/pars.taud - u*f;
        dfdy[1] = 0.0;
        dfdy[2] = -y[0]*f;
        //y
        dfdy[3] = pars.psc0*u*f;
        dfdy[4] = -1.0/pars.taui;
        dfdy[5] = pars.psc0*y[0]*f;
        //u
        dfdy[6] = 0.0;
        dfdy[7] = 0.0;
        dfdy[8] = 1.0/pars.tauf - pars.u0*f;
    }
}

int release_ave::jacobian_gsl(double t, const double y_[], double * dfdy_, double dfdt_[]){
    int i,j;

    // Assign to internal variable
    for(i=0;i<NEQ;i++) {
        y[i] = y_[i];
    }

    // Update internal Jacobian
    release_ave::jacobian(t,0.0);

    // Recast Jacobian
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy_, NEQ, NEQ);
    gsl_matrix * m = &dfdy_mat.matrix;
    for(i=0;i<NEQ;i++){
        for(j=0;j<NEQ;j++){
            gsl_matrix_set (m, i, j, dfdy[i*NEQ + j]);
        }
        dfdt_[i] = dy[i];
    }

    return GSL_SUCCESS;
}

int release_ave::jacobian_gsl_wrapper(double t, const double y_[], double * dfdy_, double dfdt_[], void * params){
    assert(params);
    void** voidPtrs = static_cast<void**>(params);
    return(static_cast<release_ave*>(voidPtrs[0])->jacobian_gsl(t,y_,dfdy_,dfdt_));
}

out_release release_ave::simulate(double rate_,PyObject *pardict,PyObject* solver_options){
    //Set parameters
    release_ave::set_pars(pardict);
    //Set ICs
    release_ave::set_ics(pardict);
    //Assign rate
    nu = rate_;

    // Set up solver
    // 1. Recast ODEs in a way that is accepted by the stochastic integrator
    sde release_mf(NEQ,0,y,
                   this->ode_gsl_wrapper,
                   this->jacobian_gsl_wrapper,
                   nullptr,
                   nullptr);
    // 2. Define integrator instance
    de_solver integrator(&release_mf);
    // 3. Solver options
    set_solver_gsl opts(solver_options);
    // 4. Allocate solution
    out_release out((long int)opts.nstep,NVAR,NVAR+1); //Creates output of the simulator (allowed in C++ and C99)
    double **sol = new double*[NEQ];
    sol[0] = out.x;
    sol[1] = out.y;
    if(NVAR>1) sol[2] = out.u;

    // Integrate
    integrator.integrate(opts,static_cast<void*>(this),nullptr,y,y,out.ts,sol);

    // Free memory allocated for sol
    delete[] sol; // Delete 'sol' but not the blocks of memory pointed by elements in sol since those are in out.ca, and out.h

    return out;
}

//*-----------------------------------------------------------------------------
// Calcium-dependent gliotransmitter exocytosis (output class)
//-----------------------------------------------------------------------------*/
PyObject* out_gtrelease::make_PyDict(){
    init_numpy();
    //WARNING: this method destroy the class that it belongs to
    //Declare list of output arguments
    kvlist oitems;
    PyObject* outdict;
    oitems.push_back(std::make_tuple ("t",array_double_to_pyobj(calcium.ts,N)));
    oitems.push_back(std::make_tuple ("ca",array_double_to_pyobj(calcium.ca,N)));
    oitems.push_back(std::make_tuple ("h",array_double_to_pyobj(calcium.h,N)));
    oitems.push_back(std::make_tuple ("xa",array_double_to_pyobj(x,N)));
    //Build PyDict
    outdict = build_PyDict(oitems);
    return outdict;
}

//*-----------------------------------------------------------------------------
// Calcium-dependent gliotransmitter exocytosis
//-----------------------------------------------------------------------------*/
void gtrelease::set_pars(PyObject* pardict){
    // Set the calcium model
    calcium.set_pars(pardict);
    pars.calcium = calcium.get_pars();
    pars.cthr    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"cthr"));
    pars.ua      = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"ua"));
    pars.taua    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"taua"));
}

void gtrelease::set_ics(double *ics){
    calcium.set_ics(ics);
    x = ics[2]/(1-pars.ua);
}

void gtrelease::set_ics(PyObject* pardict){
    double *ics = (double*)((PyArrayObject*)PyDict_GetItemString(pardict,"ICs"))->data;
    gtrelease::set_ics(ics);
}

void gtrelease::update(double t){
    r = 0; //No release upon entering the routine
    // Receive variable x state and c variable
    double dt = t - last_spike;
    // Update variables
    r = pars.ua*x; // Release
    x = (1.0-exp(-dt/pars.taua)) + x*(1.0-pars.ua)*exp(-dt/pars.taua);
}

// Simulator
out_gtrelease gtrelease::simulate(PyObject* pardict,PyObject* solver_options){
    //Set parameters
    gtrelease::set_pars(pardict);
    // ICs
    gtrelease::set_ics(pardict);

    // Set up solver
    // 1. Recast ODEs in a way that is accepted by the stochastic integrator
    sde ca_gtrel(calcium.N_d,calcium.N_s,calcium.y,
                 this->calcium.ode_gsl_wrapper,
                 this->calcium.jacobian_gsl_wrapper,
                 this->calcium.ode_s_wrapper,
                 this->calcium.ode_sd_wrapper);
    // 2. Define integrator instance
    de_solver integrator(&ca_gtrel);
    // 3. Solver options
    set_solver_gsl opts(solver_options);
    // 4. Allocate solution
    out_gtrelease out(opts.nstep); //Creates output of the simulator (allowed in C++ and C99)
    double **sol = new double*[calcium.NEQ];
    sol[0] = out.calcium.ca;
    sol[1] = out.calcium.h;

    // Integrate
    // Setup auxiliary function for event detection (i.e. gt. threshold crossing)
    double lims[3] = {0.,1.,1}; // You must pass as last element also the index of the variable in 'y' to be clipped
    double **vars = new double*[2];
    vars[0]  = &integrator.y_old[0]; // Pointer to y_old in the integrator
    vars[1]  = ca_gtrel.y;           // Pointer to y in the integrator-compatible model (we need this y instead of integrator.sys->y because we need to change it by clipping)
    void **p = (void**)calloc(3,sizeof(void*));
    p[0]     = &pars.cthr;    //double*
    p[1]     = &integrator.event; //int*
    p[2]     = lims;
    // Integrator internal counters
    int i, nsteps=0, counter=0;
    const int NSTEP = opts.nstep;

    do{
        for(int i=0;i<calcium.NEQ;i++) sol[i] += nsteps; // Move counter position
        nsteps = integrator.integrate_event(opts,static_cast<void*>(this),&clip_thr,reinterpret_cast<void**>(vars),p,out.calcium.ts+counter,sol);
        // Exit from the integration routine if it crossed the threshold
        counter += nsteps;
        if(integrator.event>0){   // GREs are only detected when c_thr has been crossed. The above simulator may also terminate without any crossing, in this case then no update is performed.
            gtrelease::update(integrator.t);
            out.x[counter-1] = x; // Update x. It is counter-1 because the event break in the integrator is after i++ (in a do--while scheme)
            last_spike = integrator.t;

        }
        opts.t0 = integrator.t;
        opts.nstep -= nsteps; // Remaining number of steps to compute in twin (since the integrator starts from at each run 0)
        integrator.event = 0; // Reset event flag
    }while(counter < NSTEP-1);

    // Free allocated memory for auxiliary function
    delete[] vars; // Delete just vars but not the elements of vars as those will be removed by destructors of classes pointed by those elements
    free(p);
    // Free memory allocated for sol
    delete[] sol; // Delete 'sol' but not the blocks of memory pointed by elements in sol since those are in out.ca, and out.h

    // Return solution
    return out;
}

/*-----------------------------------------------------------------------------
 Tripartite synapse
-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------
 Output class
-----------------------------------------------------------------------------*/
PyObject* out_asn::make_PyDict(){
    init_numpy();
    PyObject* sol;
    PyObject* syn = synapse.make_PyDict();
    kvlist oitems;
    PyObject* outdict;
    oitems.push_back(std::make_tuple ("t",array_double_to_pyobj(ts,N)));
    oitems.push_back(std::make_tuple ("ga",array_double_to_pyobj(ga,N)));
    oitems.push_back(std::make_tuple ("gammas",array_double_to_pyobj(gammas,N)));
    if(mean_field){
        // Mean-field model
        oitems.push_back(std::make_tuple ("xa",array_double_to_pyobj(xa,N)));
        // Build PyDict
        sol = build_PyDict(oitems);
        PyDict_Merge(syn,sol,1);
        sol = PyDict_Copy(syn);
    }
    else{
        // Spiking model
        sol = gliot.make_PyDict();
        // Build PyDict
        outdict = build_PyDict(oitems);
        PyDict_Merge(sol,syn,1);     // Merge solution
        PyDict_Merge(sol,outdict,1); // Merge solution
        // Do some polishing
        if(NEQ<3){
            // We do not need ca and h when NEQ<3 (i.e. only ECS variables, and GRE given)
            PyDict_DelItem(sol,PyString_FromString("ca"));
            PyDict_DelItem(sol,PyString_FromString("h"));
        }
    }
    return sol;
}

/*-----------------------------------------------------------------------------
 Spiking Model class
-----------------------------------------------------------------------------*/
void asn::set_pars(PyObject* pardict){
    gliot.set_pars(pardict);
    synapse.set_pars(pardict);
    // Set parameters in the main class
    pars.gliot = gliot.pars;
    pars.synapse = synapse.pars;

    // ECS parameters
    pars.rhoe   = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"rhoe"));
    pars.Gtot   = 1e3*PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"Gtot")); // Convert Gtot into [uM] from [mM]
    pars.taue   = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"taue"));
    pars.Ou     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"Ou"));
    pars.Ku     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"Ku"));
    pars.js     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"js"));
    pars.Og     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"Og"));
    pars.taug   = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"taug"));
    pars.alpha  = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"alpha"));

    // Adjust system dimensions (for stochastic stimulation)
    if(NEQ>2){
        N_d = 2 + gliot.calcium.N_d;
        N_s = gliot.calcium.N_s;
    }
    else{
        N_d = NEQ;
        N_s = NEQ - N_d;
    }
}

void asn::set_ics(PyObject* pardict){
    int i;
    double *ics, *icr, ic_aux[3] = {0.,0.,0.};
    ics = (double*)((PyArrayObject*)PyDict_GetItemString(pardict,"ICs"))->data;
    icr = (double*)((PyArrayObject*)PyDict_GetItemString(pardict,"ICr"))->data;
    synapse.set_ics(icr);
    for(i=0;i<2;i++) y[i] = ics[i];  // The first two ICs are always ECS variables
    if(NEQ>2){
        ic_aux[0] = ics[2];
        ic_aux[1] = ics[3];
        ic_aux[2] = icr[3];
        gliot.set_ics(ic_aux);
        for(i=0;i<2;i++) y[i+2] = gliot.calcium.y[i]; // ICs is of 4 elements in this case and the last two are for  astrocyte
    }
    else{
        // The last ics in the ICr vector refers to xa
        ic_aux[2] = icr[3];
        gliot.set_ics(ic_aux);
    }
}

void asn::ode(double t,double dt){
    int i;

    // ECS variables G_A and Gamma_S
    dy[0] = -pars.Ou*HILL(y[0],pars.Ku,1) - y[0]/pars.taue;
    dy[1] = pars.Og*y[0]*(1-y[1]) - y[1]/pars.taug;

    if(NEQ>2){
        // Compute ODEs for the calcium part
        gliot.calcium.ode(t,dt);
        for(i=0;i<2;i++) dy[i+2] = gliot.calcium.dy[i];
    }
}

int asn::ode_gsl(double t, const double y_[], double dy_[]){
    int i;

    // Assign to internal variable
    for(i=0;i<2;i++) {
        y[i] = y_[i];
    }
    if(NEQ>2){
        for(i=0;i<2;i++)
            gliot.calcium.y[i] = y_[i+2];
    }

    // Update internal state and compute derivatives
    asn::ode(t,0);

    // Assign to solver output
    for(i=0;i<NEQ;i++) {
        dy_[i] = dy[i];
    }

    return GSL_SUCCESS;
}

int asn::ode_gsl_wrapper(double t, const double y_[], double dy_[], void * params){
    assert(params);
    void** voidPtrs = static_cast<void**>(params);
    return(static_cast<asn*>(voidPtrs[0])->ode_gsl(t,y_,dy_));
}

void asn::jacobian(double t, double dt){
    if(NEQ>2){
        // G_A
        dfdy[0] = -pars.Ou*Hx(y[0],pars.Ku,1) - 1/pars.taue;
        dfdy[1] = 0.0;
        dfdy[2] = 0.0;
        dfdy[3] = 0.0;
        // Gamma_S
        dfdy[4] = pars.Og*(1-y[1]);
        dfdy[5] = -pars.Og*y[0] - 1.0/pars.taug;
        dfdy[6] = 0.0;
        dfdy[7] = 0.0;
        // Astrocyte part
        gliot.calcium.jacobian(t, dt);
        // C
        dfdy[8] = 0.0;
        dfdy[9] = 0.0;
        dfdy[10] = gliot.calcium.dfdy[0];
        dfdy[11] = gliot.calcium.dfdy[1];
        // h
        dfdy[12] = 0.0;
        dfdy[13] = 0.0;
        dfdy[14] = gliot.calcium.dfdy[2];
        dfdy[15] = gliot.calcium.dfdy[3];
    }
    else{
        // G_A
        dy[0] = -pars.Ou*Hx(y[0],pars.Ku,1) - 1.0/pars.taue;
        dy[1] = 0.0;
        // Gamma_S
        dy[2] = pars.Og*(1-y[1]);
        dy[3] = -pars.Og*y[0] - 1.0/pars.taug;
    }
}

int asn::jacobian_gsl(double t, const double y_[], double * dfdy_, double dfdt_[]){
    int i,j;

    // Assign to internal variable
    for(i=0;i<2;i++) {
        y[i] = y_[i];
    }
    if(NEQ>2){
        for(i=0;i<2;i++)
            gliot.calcium.y[i] = y_[i+2];
    }

    // Update internal Jacobian
    asn::jacobian(t,0.0);

    // Recast Jacobian
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy_, NEQ, NEQ);
    gsl_matrix * m = &dfdy_mat.matrix;
    for(i=0;i<NEQ;i++){
        for(j=0;j<NEQ;j++){
            gsl_matrix_set (m, i, j, dfdy[i*NEQ + j]);
        }
        dfdt_[i] = dy[i];
    }

    return GSL_SUCCESS;
}

int asn::jacobian_gsl_wrapper(double t, const double y_[], double * dfdy_, double dfdt_[], void * params){
    assert(params);
    void** voidPtrs = static_cast<void**>(params);
    return(static_cast<asn*>(voidPtrs[0])->jacobian_gsl(t,y_,dfdy_,dfdt_));
}

int asn::ode_s_wrapper(double t,double dt,const double y_[], double dy_[], void *params){
    assert(params);
    void** voidPtrs = static_cast<void**>(params);
    // You need to feed into ode_s a vector of size lra:NEQ with (c,h)
    return(static_cast<asn*>(voidPtrs[0])->gliot.calcium.ode_s(t,dt,y_+2,dy_));
}

int asn::ode_sd_wrapper(double t,double dt,const double y_[], double dy_[], void *params){
    assert(params);
    void** voidPtrs = static_cast<void**>(params);
    // You need to feed into ode_sd a vector of size lra:NEQ with (c,h)
    return(static_cast<asn*>(voidPtrs[0])->gliot.calcium.ode_sd(t,dt,y_+2,dy_));
}

out_asn asn::simulate(PyObject* pardict,PyObject* solver_options,double *spikes, double *gres){
    int i; // Generic counter

    // Set-up model
    asn::set_pars(pardict);
    asn::set_ics(pardict);

    // Set up solver
    // 1. Recast ODEs in a way that is accepted by the stochastic integrator
    double y0[4] = {0.,0.,0.,0.};
    for(i=0;i<2;i++) {
        y0[i] = y[i];
        if(NEQ>2) y0[i+2] = gliot.calcium.y[i]; // CHECK : might SEGFAULT for NEQ<3
    }
    // Solver compatible class instance (N_d and N_s where set by asn::set_pars(PyObj*)
    sde tsn(N_d,N_s,y0,
            this->ode_gsl_wrapper,
            this->jacobian_gsl_wrapper,
            this->ode_s_wrapper,
            this->ode_sd_wrapper);
    // 2. Define integrator instance
    de_solver integrator(&tsn);
    // 3. Solver options
    set_solver_gsl opts(solver_options);
    // 4. Allocate solution
    out_asn out(NEQ,NVAR,opts.nstep,N_spk,N_gre); //Creates output of the simulator (allowed in C++ and C99)
    // Solution pointers
    double **sol = new double*[NEQ];
    sol[0] = out.ga;
    sol[1] = out.gammas;
    if(NEQ>2){
        sol[2] = out.gliot.calcium.ca;
        sol[3] = out.gliot.calcium.h;
    }

    // Integrator internal counters
    int nsteps=0, counter=0;
    const int NSTEP = opts.nstep;

    if(NEQ>2){
        // Integrate
        // Setup auxiliary function for event detection (i.e. gt. threshold crossing)
        double lims[3] = {0.,1.,1};      // h is now in 4th position, but because we are passing y from the 3rd position, then h is the second element in y
        double **vars = new double*[2];
        vars[0]  = &integrator.y_old[2]; // Pointer to "c"in y_old in the integrator
        vars[1]  = tsn.y+2;              // Pointer to "c" in y in the integrator-compatible model (we need this y instead of integrator.sys->y because we need to change it by clipping)
        void **p = (void**)calloc(3,sizeof(void*));
        p[0]     = &pars.gliot.cthr;    //double*
        p[1]     = &integrator.event; //int*
        p[2]     = lims;

        // Actual integration: We integrate up to a C_thr crossing and a GRE. We exit the integrator, update xa, and restart integration
        // with an updated G_A by r_A. The interruption of integration and GRE detecton is by clip_thr, setting integrator.event=1;
        do{
            for(i=0;i<NEQ;i++) sol[i] += nsteps; // Move counter position
            nsteps = integrator.integrate_event(opts,static_cast<void*>(this),&clip_thr,reinterpret_cast<void**>(vars),p,out.ts+counter,sol);
            // Exit from the integration routine if it crossed the threshold
            counter += nsteps;
            if(integrator.event>0){   // GREs are only detected when c_thr has been crossed. The above simulator may also terminate without any crossing, in this case then no update is performed.
                gliot.update(integrator.t);
                out.gliot.x[counter-1] = gliot.x; // Update x. It is counter-1 because the event break in the integrator is after i++ (in a do--while scheme)
                gliot.last_spike = integrator.t;
                // You need to update the ICs on G_A --> I did not find a better way to do this
                integrator.sys->y[0] += pars.rhoe*pars.Gtot*gliot.r;
            }
            opts.t0 = integrator.t;
            opts.nstep -= nsteps; // Remaining number of steps to compute in twin (since the integrator starts from at each run 0)
            integrator.event = 0; // Reset event flag
        }while(counter < NSTEP-1);
        // Free temporarily allocate memory
        delete[] vars; // Delete just vars but not the elements of vars as those will be removed by destructors of classes pointed by those elements
        free(p);
    }
    else{
        // Deals with the case where no calcium dynamics is needed and spikes for GRE are provided
        // In this case since we already have the instants of GREs, we rely on 'explicit' handling of events, interrupting
        // integration every time we reach a new spike
        int gre_index = 0;
        double tfin_;
        do{
            // First see where we are on the given GRE sequence
            if(gre_index<N_gre) {
                // The integrator takes an integer version of time (via counter and opts.nstep) so we need to modify these
                // parameters rather than opts.tfin
                if(gre_index>0)
                    opts.nstep = (long int)((gres[gre_index]-opts.t0)/opts.dt);
                else
                    opts.nstep = (long int)((gres[gre_index]-opts.t0)/opts.dt + 1);
            }
            else{
                opts.nstep = NSTEP-counter;
            }
            // Proceed with normal integration (according to the above scheme for NEQ>2)
            for(int i=0;i<NEQ;i++) sol[i] += nsteps; // Move counter position
            nsteps = integrator.integrate_event(opts,static_cast<void*>(this),nullptr,static_cast<void**>(nullptr),static_cast<void**>(nullptr),out.ts+counter,sol);
            // Exit from the integration routine if it crossed the threshold
            counter += nsteps;
            if(gre_index<N_gre){ // Assumes that all GREs occurs for t<
                gliot.update(integrator.t);
                out.gliot.x[gre_index] = gliot.x; // Update x. In this case out.gliot.x has N_gre elements.
                gliot.last_spike = integrator.t;
                // You need to update the ICs on G_A --> I did not find a better way to do this
                integrator.sys->y[0] += pars.rhoe*pars.Gtot*gliot.r;
            }
            opts.t0 = integrator.t;
            // Update GRE index
            gre_index++;
        }while(counter < NSTEP);
    }

    // Append synaptic dynamics
    // Reallocate new memory
    void **vars_s = (void**)calloc(3,sizeof(void*));
    vars_s[0] = &synapse.pars;
    vars_s[1] = &pars;
    vars_s[2] = &opts.dt;
    double **p_s = new double*[1];
    p_s[0] = out.gammas;

    // Use a version of synapse.simulation that allows for changing u0 in the synaptic model
    out.synapse = synapse.simulate(spikes,&u0_update,vars_s,reinterpret_cast<void**>(p_s));
    // Free allocated memory
    free(vars_s);
    delete[] p_s;
    // Free memory allocated for sol
    delete[] sol; // Delete 'sol' but not the blocks of memory pointed by elements in sol since those are in out.ca, and out.h

    // Return solution
    return out;
}

/*-----------------------------------------------------------------------------
 Mean-field Model class
-----------------------------------------------------------------------------*/
void asn_ave::set_pars(PyObject* pardict){
    synapse.set_pars(pardict);
    // Set parameters in the main class
    pars.synapse = synapse.pars;

    // ECS parameters
    pars.rhoe   = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"rhoe"));
    pars.Gtot   = 1e3*PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"Gtot")); // Convert into [uM]
    pars.taue   = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"taue"));
    pars.Ou     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"Ou"));
    pars.Ku     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"Ku"));
    pars.js     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"js"));
    pars.Og     = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"Og"));
    pars.taug   = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"taug"));
    pars.alpha  = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"alpha"));

    // Set Gt dynamics parameters
    gliot.pars.u0      = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"ua"));
    gliot.pars.taud    = PyFloat_AS_DOUBLE(PyDict_GetItemString(pardict,"taua"));
    gliot.pars.psc0    = 0.0; // This allows some short-circuiting in the equations of release_ave
    // Copy to the internal model parameter (only ua, and taua are necessary) --> this is actually just for completeness
    pars.gliot.ua   = gliot.pars.u0;
    pars.gliot.taua = gliot.pars.taud;
}

void asn_ave::set_ics(PyObject* pardict){
    int i;
    double *ics, *icr, ic_aux[3] = {0.,0.,0.};
    ics = (double*)((PyArrayObject*)PyDict_GetItemString(pardict,"ICs"))->data;
    icr = (double*)((PyArrayObject*)PyDict_GetItemString(pardict,"ICr"))->data;
    ic_aux[0] = icr[3]; // xa(0)
    gliot.set_ics(ic_aux); // feeds [xa(0),0,0]
    for(i=0;i<NVAR+1;i++) ic_aux[i] = icr[i];
    synapse.set_ics(ic_aux); // feeds [xs(0),ys(0),us(0)]

    // Assign ICs
    y[0] = ics[0];       // G_A
    y[1] = ics[1];       // Gamma_S
    y[2] = gliot.y[0];   // X_A
    y[3] = synapse.y[0]; // X_S
    y[4] = synapse.y[1]; // Y_S
    if(NVAR>1) y[5] = synapse.y[2];
}

void asn_ave::ode(double t,double dt){
    double f_s = synapse.nu + synapse.bias(t,dt);
    double f_a = gliot.nu + gliot.bias(t,dt);

    // Update u0
    synapse.pars.u0 = U0(pars.synapse.u0,pars.alpha,y[1]);

    // Compute derivatives of subsystems
    gliot.ode(t,dt);
    synapse.ode(t,dt);

    dy[0] = pars.rhoe*pars.Gtot*pars.gliot.ua*y[2]*f_a - y[0]/pars.taue;
    dy[1] = pars.Og*y[0]*(1-y[1]) - y[1]/pars.taug;
    dy[2] = gliot.dy[0];
    dy[3] = synapse.dy[0];
    dy[4] = synapse.dy[1];
    if(NVAR>1) dy[5] = synapse.dy[2];
}

int asn_ave::ode_gsl(double t, const double y_[], double dy_[]){
    int i;

    // Assign to internal variable
    for(i=0;i<2;i++) {
        y[i] = y_[i];
    }
    gliot.y[0]   = y_[2];
    synapse.y[0] = y_[3];
    synapse.y[1] = y_[4];
    if(NVAR>1) synapse.y[2] = y_[5];

    // Update internal state and compute derivatives
    asn_ave::ode(t,0);

    // Assign to solver output
    for(i=0;i<NEQ;i++) {
        dy_[i] = dy[i];
    }

    return GSL_SUCCESS;
}

int asn_ave::ode_gsl_wrapper(double t, const double y_[], double dy_[], void * params){
    assert(params);
    void** voidPtrs = static_cast<void**>(params);
    return(static_cast<asn_ave*>(voidPtrs[0])->ode_gsl(t,y_,dy_));
}

void asn_ave::jacobian(double t, double dt){
    double f_s = synapse.nu + synapse.bias(t,dt);
    double f_a = gliot.nu + gliot.bias(t,dt);

    // Update u0
    synapse.pars.u0 = U0(pars.synapse.u0,pars.alpha,y[1]);

    // Compute jacobians of subsystems
    gliot.jacobian(t,dt);
    synapse.jacobian(t,dt);

    if(NVAR>1){
        // G_A
        dfdy[0] = -1.0/pars.taue;
        dfdy[1] = 0.0;
        dfdy[2] = pars.rhoe*pars.Gtot*pars.gliot.ua*f_a;
        dfdy[3] = 0.0;
        dfdy[4] = 0.0;
        dfdy[5] = 0.0;
        // Gamma_S
        dfdy[6] = pars.Og*(1-y[1]);
        dfdy[7] = -pars.Og*y[0] - 1.0/pars.taug;
        dfdy[8] = 0.0;
        dfdy[9] = 0.0;
        dfdy[10] = 0.0;
        dfdy[11] = 0.0;
        // X_A
        dfdy[12] = 0.0;
        dfdy[13] = 0.0;
        dfdy[14] = gliot.dfdy[0];
        dfdy[15] = 0.0;
        dfdy[16] = 0.0;
        dfdy[17] = 0.0;
        // X_S
        dfdy[18] = 0.0;
        dfdy[19] = 0.0;
        dfdy[20] = 0.0;
        dfdy[21] = synapse.dfdy[0];
        dfdy[22] = synapse.dfdy[1];
        dfdy[23] = synapse.dfdy[2];
        // Y_S
        dfdy[24] = 0.0;
        dfdy[25] = 0.0;
        dfdy[26] = 0.0;
        dfdy[27] = synapse.dfdy[3];
        dfdy[28] = synapse.dfdy[4];
        dfdy[29] = synapse.dfdy[5];
        // U_S
        dfdy[30] = 0.0;
        dfdy[31] = (1.0/pars.synapse.tauf + (1-y[5])*f_s)*(pars.alpha*pars.synapse.u0);
        dfdy[32] = 0.0;
        dfdy[33] = synapse.dfdy[6];
        dfdy[34] = synapse.dfdy[7];
        dfdy[35] = synapse.dfdy[8];
    }
    else{
        dfdy[0] = -1.0/pars.taue;
        dfdy[1] = 0.0;
        dfdy[2] = pars.rhoe*pars.Gtot*pars.gliot.ua*f_a;
        dfdy[3] = 0.0;
        dfdy[4] = 0.0;
        // Gamma_S
        dfdy[5] = pars.Og*(1-y[1]);
        dfdy[6] = -pars.Og*y[0] - 1.0/pars.taug;
        dfdy[7] = 0.0;
        dfdy[8] = 0.0;
        dfdy[9] = 0.0;
        dfdy[10] = 0.0;
        // X_A
        dfdy[11] = 0.0;
        dfdy[12] = 0.0;
        dfdy[13] = gliot.dfdy[0];
        dfdy[14] = 0.0;
        dfdy[15] = 0.0;
        // X_S
        dfdy[16] = 0.0;
        dfdy[17] = 0.0;
        dfdy[18] = 0.0;
        dfdy[19] = synapse.dfdy[0];
        dfdy[20] = synapse.dfdy[1];
        // Y_S
        dfdy[21] = 0.0;
        dfdy[22] = 0.0;
        dfdy[23] = 0.0;
        dfdy[24] = synapse.dfdy[2];
        dfdy[25] = synapse.dfdy[3];
    }
}

int asn_ave::jacobian_gsl(double t, const double y_[], double * dfdy_, double dfdt_[]){
    int i,j;

    // Assign to internal variable
    for(i=0;i<2;i++) {
        y[i] = y_[i];
    }
    gliot.y[0]   = y_[2];
    synapse.y[0] = y_[3];
    synapse.y[1] = y_[4];
    if(NVAR>1) synapse.y[2] = y_[5];

    // Update internal Jacobian
    asn_ave::jacobian(t,0.0);

    // Recast Jacobian
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy_, NEQ, NEQ);
    gsl_matrix * m = &dfdy_mat.matrix;
    for(i=0;i<NEQ;i++){
        for(j=0;j<NEQ;j++){
            gsl_matrix_set (m, i, j, dfdy[i*NEQ + j]);
        }
        dfdt_[i] = dy[i];
    }

    return GSL_SUCCESS;
}

int asn_ave::jacobian_gsl_wrapper(double t, const double y_[], double * dfdy_, double dfdt_[], void * params){
    assert(params);
    void** voidPtrs = static_cast<void**>(params);
    return(static_cast<asn_ave*>(voidPtrs[0])->jacobian_gsl(t,y_,dfdy_,dfdt_));
}

out_asn asn_ave::simulate(double rate_s,double rate_a,PyObject *pardict,PyObject* solver_options){
    asn_ave::set_pars(pardict);
    asn_ave::set_ics(pardict);

    // Set internal rate in synaptic dynamics
    synapse.nu = rate_s;
    // Set internal rate in gt. release
    gliot.nu = rate_a;

    // Set up solver
    // 1. Recast ODEs in a way that is accepted by the stochastic integrator
    sde tsn_mf(NEQ,0,y,
               this->ode_gsl_wrapper,
               this->jacobian_gsl_wrapper,
               nullptr,
               nullptr);
    // 2. Define integrator instance
    de_solver integrator(&tsn_mf);
    // 3. Solver options
    set_solver_gsl opts(solver_options);
    // 4. Allocate solution
    out_asn out(NEQ,NVAR,(long int)opts.nstep);
    double **sol = new double*[NEQ];
    sol[0] = out.ga;
    sol[1] = out.gammas;
    sol[2] = out.xa;
    sol[3] = out.synapse.x;
    sol[4] = out.synapse.y;
    if(NVAR>1) sol[5] = out.synapse.u;

    // Integrate
    integrator.integrate(opts,static_cast<void*>(this),
                         nullptr,static_cast<double*>(nullptr),static_cast<double*>(nullptr),out.ts,sol);

    // Free memory allocated for sol
    delete[] sol; // Delete 'sol' but not the blocks of memory pointed by elements in sol since those are in out.ca, and out.h

    return out;

}

///*-----------------------------------------------------------------------------
//Testing
//-----------------------------------------------------------------------------*/
//int main(){
//    // Define astrocyte model
//    lra astrocyte;
//
//    // Declare output structure
//    out_lra out;
//
//    // Simulator
//    out = astrocyte.simulate();
//
//    return 1;
//}
