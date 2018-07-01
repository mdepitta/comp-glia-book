// gcc -g -std=c++11 -lpython2.7 -lstdc++ -lm -lgsl -lgslcblas -I/usr/include/python2.7 -I/usr/include/numpy -I/home/maurizio/Editorial.Projects/Chapter.DePitta.Gliotransmission/code/pycustommodules -I/home/maurizio/Editorial.Projects/Chapter.DePitta.Gliotransmission/code/pycustommodules/solvers /home/maurizio/Editorial.Projects/Chapter.DePitta.Gliotransmission/code/pycustommodules/pycapi_utils.cpp /home/maurizio/Editorial.Projects/Chapter.DePitta.Gliotransmission/code/pycustommodules/solvers/solver_options.cpp /home/maurizio/Editorial.Projects/Chapter.DePitta.Gliotransmission/code/pycustommodules/solvers/stochastic_solvers.cpp -o test.o

#include "stochastic_solvers.h"

//----------------------------------------------------------------------------------------------------------------------
// Testing model
//----------------------------------------------------------------------------------------------------------------------
class model{
        public:
        int noise = 0;  // noise flag (default no noise)
        size_t NEQ = 2;
        // Variables
        double y[2],dy[2];  //C,h,I state variables; their derivatives; and the old values
        model(){y[0]=1.0; y[1]=0.1;};  //constructor
        ~model(){}; //destructor
        // RHS (deterministic version)
        void ode(double t);      //ODE (used for integration) (standard form)
        int ode_gsl(double t, const double y_[], double dy_[]); // GSL wrapper
        static int ode_gsl_wrapper(double t, const double y_[], double dy_[], void * params);
        // Simulator
        void simulate();
};

void model::ode(double t){
    dy[0] = y[0] - y[1];
    dy[1] = y[1]*y[0]*y[0] - 2.5*y[0]/y[1];
}

int model::ode_gsl(double t, const double y_[], double dy_[]){
    size_t i;
    // DEBUG
//    printf("\nNEQ=%d",NEQ);
    for(i=0;i<NEQ;i++) y[i] = y_[i];
    model::ode(t);
    for(i=0;i<NEQ;i++) dy_[i] = dy[i];
    return GSL_SUCCESS;
}

int model::ode_gsl_wrapper(double t, const double y_[], double dy_[], void * params){
    assert(params);
    void** voidPtrs = static_cast<void**>(params);
    return(static_cast<model*>(voidPtrs[0])->ode_gsl(t,y_,dy_));
}

void model::simulate(){
    sde recast_model(2,0,y,this->ode_gsl_wrapper,nullptr,nullptr,nullptr);

    // 2. Define integrator instance
    de_solver integrator(&recast_model);
    // 3. Solver options
    set_solver_gsl opts;
    // 4. Allocate solution
    // Integrate
//    integrator.integrate("msadams",opts,static_cast<void*>(this));
}

//----------------------------------------------------------------------------------------------------------------------
// Stochastic ODE system definition
//----------------------------------------------------------------------------------------------------------------------
sde::sde(size_t Nd_,size_t Ns_, double *y0,
         int (*dfp)(double, const double*, double*, void*),
         int (*jcb)(double, const double*, double*, double*, void*),
         int (*sfp)(double, double, const double*, double*, void*),
         int (*sdfp)(double, double, const double*, double*, void*)){
    // Assign system dimensions
    N_d = Nd_;
    N_s = Ns_;
    NEQ = N_d + N_s;

    //DEBUG
//    printf("\nSDE: N_d=%d\tN_s=%d\tNEQ=%d",N_d,N_s,NEQ);

    if(N_s>0) noise=1;

    // Assign RHS
    drift = dfp;
    jac   = jcb;
    diffusion            = sfp;
    diffusion_derivative = sdfp;

    // Create internal variables (creates a copy of the original model variables and will work on this copy)
    y     = new double[NEQ];  // NOTE: The first N_d variables are deterministic and the last N_s are stochastic
    ddy   = new double[NEQ];  // NOTE: The first N_d variables are deterministic and the last N_s are stochastic
    sdy   = new double[N_s];
    sdy_d = new double[N_s];

//    // Initialize variables
    for(int i=0;i<NEQ;i++) y[i] = y0[i];
}

sde::~sde(){
    delete[] y;
    delete[] ddy;
    delete[] sdy;
    delete[] sdy_d;
}

int sde::drift_wrapper(double t, const double y_[], double dy_[], void * params){
    assert(params);
    void** voidPtrs = static_cast<void**>(params);
    return(static_cast<sde*>(voidPtrs[1])->drift(t,y_,dy_,params));
}

int sde::jac_wrapper(double t, const double y_[], double *dfdy_, double dfdt_[], void *params){
    assert(params);
    void** voidPtrs = static_cast<void**>(params);
    return(static_cast<sde*>(voidPtrs[1])->jac(t,y_,dfdy_,dfdt_,params));
}

int sde::diffusion_wrapper(double t, double dt, const double y[], double dydt[], void *params){
    assert(params);
    void** voidPtrs = static_cast<void**>(params);
    return(static_cast<sde*>(voidPtrs[1])->diffusion(t,dt,y,dydt,params));
}

int sde::diffusion_derivative_wrapper(double t, double dt, const double y[], double dydt[], void *params){
    assert(params);
    void** voidPtrs = static_cast<void**>(params);
    return(static_cast<sde*>(voidPtrs[1])->diffusion_derivative(t,dt,y,dydt,params));
}

//----------------------------------------------------------------------------------------------------------------------
// Generic stochastic solver
//----------------------------------------------------------------------------------------------------------------------
de_solver::de_solver(sde *sys_){
    sys   = sys_;
    noise = sys->noise;
    N_d = sys->N_d;
    N_s = sys->N_s;
    NEQ = sys->NEQ;
    // DEBUG
    y_old = new (std::nothrow) double [NEQ];

    // DEBUG
//    printf("\nSOLV: N_d=%d\tN_s=%d\tNEQ=%d",N_d,N_s,NEQ);
}

de_solver::~de_solver(){
    delete[] y_old;
}

void de_solver::milstein(double t, double dt, void *params){
    size_t i;
    auto dist = std::bind(std::normal_distribution<double>{0.0, 1.0},
                          std::mt19937(std::random_device{}()));
    double dw = sqrt(dt)*dist(); // Noise

//    // DEBUG
//    printf("\nN_d=%d\tN_s=%d",N_d,N_s);
//    for(i=0;i<NEQ;i++) printf("\nM: sys->y[%d]=%f",i,sys->y[i]);
//    printf("\nM:\t");
//    for(i=0;i<NEQ;i++) printf("\ty_old[%d]=%f",i,y_old[i]);
//    printf("\nMS:\t");
//    for(i=0;i<N_s;i++) printf("\tsdy[%d]=%f\t sdy_d[%d]=%f",i,sys->sdy[i],i,sys->sdy_d);

    // Compute derivatives
    sys->diffusion(t,dt,y_old,sys->sdy,params);
    sys->diffusion_derivative(t,dt,y_old,sys->sdy_d,params);
//    for(i=0;i<N_s;i++) printf("\t-->sdy[%d]=%f\t sdy_d[%d]=%f",i,sys->sdy[i],i,sys->sdy_d);
    for(i=0;i<N_s;i++){
        sys->y[i+N_d] = sys->y[i+N_d] + dw*sys->sdy[i] + 0.5*(dw*dw)*sys->sdy[i]*sys->sdy_d[i];
    }
}

void de_solver::integrate(set_solver_gsl opts,void *params_base,void (*fpa)(double*, double p[]),
                          double *f_var, double *f_par, double *ts,double **sol){
    // ts must be double[opts.npts]
    // sol must be double[NEQ][opts.npts]
    // Standard solver version with basic function

    // Instantiate GSL solver
    void **params = (void**)calloc(2,sizeof(void*));
    params[0] = params_base;
    params[1] = static_cast<void *>(sys);
    const gsl_odeiv2_step_type *method;
    jacPtr jacobian;
    const char *name = (const char*)opts.solver_id;
    // Implicit methods
    if(strcmp(name,"gsl_rk1imp")==0){
        method = gsl_odeiv2_step_rk1imp;
        jacobian = sys->jac_wrapper;
    }
    else if(strcmp(name,"gsl_rk2imp")==0){
        method = gsl_odeiv2_step_rk2imp;
        jacobian = sys->jac_wrapper;
    }
    else if(strcmp(name,"gsl_rk4imp")==0){
        method = gsl_odeiv2_step_rk4imp;
        jacobian = sys->jac_wrapper;
    }
    else if(strcmp(name,"gsl_bsimp")==0){
        method = gsl_odeiv2_step_msbdf;
        jacobian = sys->jac_wrapper;
    }
    else if(strcmp(name,"gsl_msbdf")==0){
        method = gsl_odeiv2_step_msbdf;
        jacobian = sys->jac_wrapper;
    }
    // Explicit methods
    else if(strcmp(name,"gsl_rk4")==0){
        method = gsl_odeiv2_step_rk4;
        jacobian = nullptr;
    }
    else if(strcmp(name,"gsl_rk8pd")==0){
        method = gsl_odeiv2_step_rk8pd;
        jacobian = nullptr;
    }
    else{
        // ms_adams method (no jacobian) (respond to solver_id='gsl_msadams')
        method = gsl_odeiv2_step_msadams;
        jacobian = nullptr;
    }
    gsl_odeiv2_system system = {sys->drift_wrapper, jacobian, sys->NEQ, params};
    gsl_odeiv2_driver * d;
    d = gsl_odeiv2_driver_alloc_y_new (&system, method, opts.dt, opts.atol, opts.rtol);

    double ti;
    t = opts.t0;
    int status; // Integrator status
//    double lims[2] = {0.,1.};

    int i = 0,j;
    do{
        ti = t + opts.dt; // t is automatically updated by the driver
        // Produce output (since we are running in a do cycle, the output must be assigned
        // before the next integration step) -- ts and sol are class public variables
        ts[i] = ti;
        for(j=0;j<NEQ;j++){
            sol[j][i] = sys->y[j];
            y_old[j] = sys->y[j]; // Save old status internally
        }
        // Integration step (deterministic part)
        status = gsl_odeiv2_driver_apply (d, &t, ti, sys->y);

        // Error checking
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            break;
        }
        if(sys->noise){
            de_solver::milstein(t,opts.dt,params); // Add the stochastic part to y[1]
        }
        if(fpa!=nullptr) fpa(f_var,f_par);
        i++;
    }while(i < opts.nstep);

    // Save only last point of solution
    // Free memory
    free(params);
    gsl_odeiv2_driver_free (d);
}

int de_solver::integrate_event(set_solver_gsl opts,void *params_base,void (*fpa)(void**, void**),
                          void **f_var, void **f_par, double *ts,double **sol){
    // ts must be double[opts.npts]
    // sol must be double[NEQ][opts.npts]
    // NOTE: use this method when you need to perform complex manipulations of data, or interrupt cycle through events
    // The pointer to the generic function void (*fpa)(void**, void**) allows this.

    // Instantiate GSL solver
    void **params = (void**)calloc(2,sizeof(void*));
    params[0] = params_base;
    params[1] = static_cast<void *>(sys);
    const gsl_odeiv2_step_type *method;
    jacPtr jacobian;
    const char *name = (const char*)opts.solver_id;
    // Implicit methods
    if(strcmp(name,"gsl_rk1imp")==0){
        method = gsl_odeiv2_step_rk1imp;
        jacobian = sys->jac_wrapper;
    }
    else if(strcmp(name,"gsl_rk2imp")==0){
        method = gsl_odeiv2_step_rk2imp;
        jacobian = sys->jac_wrapper;
    }
    else if(strcmp(name,"gsl_rk4imp")==0){
        method = gsl_odeiv2_step_rk4imp;
        jacobian = sys->jac_wrapper;
    }
    else if(strcmp(name,"gsl_bsimp")==0){
        method = gsl_odeiv2_step_msbdf;
        jacobian = sys->jac_wrapper;
    }
    else if(strcmp(name,"gsl_msbdf")==0){
        method = gsl_odeiv2_step_msbdf;
        jacobian = sys->jac_wrapper;
    }
    // Explicit methods
    else if(strcmp(name,"gsl_rk4")==0){
        method = gsl_odeiv2_step_rk4;
        jacobian = nullptr;
    }
    else if(strcmp(name,"gsl_rk8pd")==0){
        method = gsl_odeiv2_step_rk8pd;
        jacobian = nullptr;
    }
    else{
        // ms_adams method (no jacobian) (respond to solver_id='gsl_msadams')
        method = gsl_odeiv2_step_msadams;
        jacobian = nullptr;
    }
    gsl_odeiv2_system system = {sys->drift_wrapper, jacobian, sys->NEQ, params};
    gsl_odeiv2_driver * d;
    d = gsl_odeiv2_driver_alloc_y_new (&system, method, opts.dt, opts.atol, opts.rtol);

    double ti;
    t = opts.t0;
    int status; // Integrator status

    int i = 0,j;
    // DEBUG
    int k;
    do{
        ti = t + opts.dt; // t is automatically updated by the driver
        // Produce output (since we are running in a do cycle, the output must be assigned
        // before the next integration step) -- ts and sol are class public variables
        ts[i] = ti;
        for(j=0;j<NEQ;j++) {
            sol[j][i] = sys->y[j];
            y_old[j] = sys->y[j];  // It is useful to have an internal copy of the value of y before update
        }

//        // DEBUG
//        printf("\nB:\t");
//        for(k=0;k<NEQ;k++) printf("\ty[%d]=%f",k,sys->y[k]);
        // Integration step (deterministic part)
        status = gsl_odeiv2_driver_apply (d, &t, ti, sys->y);

//        // DEBUG
//        printf("\nA:\t");
//        for(k=0;k<NEQ;k++) printf("\ty[%d]=%f",k,sys->y[k]);

        // Error checking
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            break;
        }
        if(sys->noise){
            de_solver::milstein(t,opts.dt,params); // Add the stochastic part to y[1]
//            // DEBUG
//            printf("\nF:\t");
//            for(k=0;k<NEQ;k++) printf("\ty[%d]=%f",k,sys->y[k]);
        }
        if(fpa!=nullptr) fpa(f_var,f_par);
        i++;
        if(event>0) break;  // Exit from cycle, this is only if (fpa) create an event.
    }while(i < opts.nstep);

    // Save only last point of solution
    // Free memory
    free(params);
    gsl_odeiv2_driver_free (d);

    return i;
}


//----------------------------------------------------------------------------------------------------------------------
// Testing
//----------------------------------------------------------------------------------------------------------------------
//int main(){
//    int i, N=100;
//    sde system(0,1,nullptr,drift,diff,diff_der);
//    system.y[0] = 0.0;
//    double t=0, T = 10, dt=1e-3;
//    // Instantiate integrator
//    de_solver integrator(&system,dt);
//
//    for(t=0;t<T;t+=dt){
//        printf("\nt=%f\ty=%f",t,system.y[0]);
//        integrator.milstein(t);
//}
///*-----------------------------------------------------------------------------
//Testing
//-----------------------------------------------------------------------------*/
//int main(){
//    // Define astrocyte model
//    model dummy;
//
//    // Simulator
//    dummy.simulate();
//
//    return 1;
//}