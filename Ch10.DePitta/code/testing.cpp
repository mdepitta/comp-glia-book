#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <functional> //std::bind
#include <cassert>
#include <cstring>
#include <iterator>   //noise generator
#include <random>     //noise generator

void clip_h(double *x){
  x[1] = std::max(0.0, std::min(x[1], 1.0));
}

class A{
    public:
        double y[2];
        A(){y[0]= 10.; y[1]=20.5;};
        ~A(){};
        int square(double t, double *y_);
        static int square_wrapper(double t, double *y_, void *params);
        void simulate();
};

class B{
    public:
        double *y;
        B(double*, int (*fp)(double, double*, void*));
        ~B();
        int (*B_fun)(double, double *y_, void*);
        static int B_wrapper(double t, double *y_, void* params);
};

class C{
    B *sys;
    public:
        C(B *obj){sys = obj;};
        ~C(){};
        void integrate(void *params, void (*fpa)(double*),double *f_var);
};

/*-------------------------------------------------------------------------
// Implementation
-------------------------------------------------------------------------*/
int A::square(double t, double *y_){
    y[0] = y_[0]; y[1] = y_[1];
    y[0] = y[0] + y[1]*t;
    y[1] = y[0]*y[0] + y[1]*y[1]*t + y[0]*y[1];
    return 1;
}

int A::square_wrapper(double t, double *y_, void *params){
    assert(params);
    void** voidPtrs = static_cast<void**>(params);
    return(static_cast<A*>(voidPtrs[0])->square(t,y_));
}

void A::simulate(){
    B A_solver(y,this->square_wrapper);
    C integrator(&A_solver);
    // The next instance cannot modify y
    integrator.integrate(static_cast<void*>(this),&clip_h,static_cast<double*>(y));
}

B::B(double *x, int (*fp)(double, double*, void*)){
//    y  = new double[2];
    y = x;
//    y[0] = x[0];
//    y[1] = x[1];
    B_fun = fp;
}

B::~B(){
//    delete[] y;
}

int B::B_wrapper(double t, double *y_, void* params){
    assert(params);
    void** voidPtrs = static_cast<void**>(params);
    return(static_cast<B*>(voidPtrs[1])->B_fun(t,y_,params));
}

void C::integrate(void *A_par, void (*fpa)(double*),double *f_var){
    void **params = (void**)calloc(2,sizeof(void*));
    params[0] = A_par;
    params[1] = static_cast<void *>(sys);
    double t=0, dt=0.05, tf=1.0;

    // Some loop over time
    while(t<tf){
        sys->B_fun(t,sys->y,params);
        printf("\ny[0]=%f\ty[1]=%f",sys->y[0],sys->y[1]);
        fpa(f_var); // This is the clipping
        printf("\tCLIPPING: y[1]=%f",sys->y[1]);
        t += dt;
    }
}

//----------------------------------------------------------
// Allocation issue
//----------------------------------------------------------
class A2{
    int N;
    public:
        double *y;
        A2(int N_);
        ~A2(){free(y);};
};

A2::A2(int N_){
    N = N_;
    y = (double*)calloc(N,sizeof(double));
}

class B2{
    int N;
    public:
        A2 obj;
        B2(int N_) : N(N_), obj(N_){printf("\nHello!");};
        ~B2(){};
};

int main(){
    // Testing integration
//    A model;
//    model.simulate();

    // Testing allocation
    int N = 10;
    B2 model(N);
    for(int i=0;i<N;i++) model.obj.y[i] = i;
    for(int i=0;i<N;i++) printf("\ny[%d]=%f",i,model.obj.y[i]);
}