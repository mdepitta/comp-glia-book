from __future__ import division
import numpy as np
from scipy import interpolate
import scipy.signal as sps
import matplotlib.pylab as plt
import gc
import sys, os, getopt
sys.path.append(os.path.join(os.path.expanduser('~'),'Ch6.DePitta.GPCR'))
import pycustommodules.solvers.solver_utils as su
import pycustommodules.save_utils as svu
import astrocyte_models as models

## NOTE: add to PYTHONPATH the directory where PyGMO is installed
sys.path.append('/usr/local/lib/python2.7/site-packages/')
import pygmo as pg

#----------------------------------------------------------------------------------------
# Data preparation
#----------------------------------------------------------------------------------------
def prepare_data(file):
    pass

def problem_dimension(model):
    """
    Return problem dimension (aka number of variables in the optimization problem)

    Input:
    - model  :  {'astro'} | 'ip3' | 'pkc'
    - config :  String    One of the allowed configurations for each model

    Return:
    - dim    : Integer    Number of variables (aka degrees of freedom / dimension) of the problem to be optimized

    v1.0
    Maurizio De Pitta', INRIA Rhone-Alpes, Villeurbanne, November 3, 2017.
    """

    if model=='lra':
        dim = 5
    elif model=='lra2':
        # Assumes d1==d3
        dim = 4
    elif model=='lra_fit':
        dim = 5
    elif model=='chi':
        dim = 7
    elif model=='pkc2':
        dim = 11
    elif 'pkc' in model:
        dim = 13
    return dim

def rescale_var(var_val,normalized=False,bounds=None,scale='lin'):
    """
    Rescale variable (whether it is given as normalized or not)

    Input arguments:
    - var_val    : Float   Variable/Parameter value
    - normalized : {False} | True
    - bounds     : {None} | 2x1 Tuple   Must be provided only if normalized==True
    - scale      : {'lin'} | 'log'      Scaling of the value/parameter (and boundaries accordingly)

    Return:
    - rescaled_value : Float   Rescaled variable value

    v1.0
    Maurizio De Pitta', INRIA Rhone-Alpes, Villeurbanne, Feb 23, 2017.
    """

    # Sanity check
    if normalized:
        if not bounds:
            raise ValueError("Bounds required by normalization (True) were not specified")

    if scale=='lin':
        if normalized:
            rescaled_value = var_val*(bounds[1]-bounds[0]) + bounds[0]
        else:
            rescaled_value = var_val
    else:
        # scale is 'log'
        if normalized:
            rescaled_value = 10.0**(bounds[0]+var_val*(bounds[1]-bounds[0]))
        else:
            rescaled_value = 10.0**var_val

    # Return rescaled value
    return rescaled_value

def rescale_vector(x,model,normalized=False,bounds=None,scaling='lin'):

    # Define some lambdas
    extract_bounds  = lambda bounds,idx_var : bounds if not bounds else bounds[idx_var]
    extract_scaling = lambda scaling,idx_var : 'lin' if not scaling else scaling[idx_var]

    # Sanity check
    if normalized:
        if not bounds:
            raise ValueError("Bounds required by normalization (True) were not specified")

    # Retrieve size of the solution vector
    dim = problem_dimension(model=model)

    # Allocate space
    x_rescaled = np.zeros(dim)
    for i in xrange(dim):
        x_rescaled[i] = rescale_var(x[i],normalized=normalized,bounds=extract_bounds(bounds,i),scale=extract_scaling(scaling,i))

    return x_rescaled

def bounds_for_exploration(bounds_dict,normalized=False):
    """
    Provides bounds for parameters in the format required by PyGMO algorithm for exploration

    Input parameters:
    - bounds     : Dict   (Ordered dictionary) of boundary tuples
    - normalized : {False} | True  Normalized boundaries (overrides the one provided)

    Return:
    - bounds     : list of tuples for bounds of each parameter
    - scaling    : list of strings specifying scaling of each parameter within its bounds

    v1.0
    Maurizio De Pitta', INRIA Rhone-Alpes, Villeurbanne, Feb 23, 2017.
    """

    if normalized:
       for k,b in bounds_dict.iteritems():
           bounds_dict[k] = [0.,1.,b[-1]]

    # Provide (ordered) bounds (implied by bounds_dict data type)
    # NOTE: bounds and scaling are implicitly ordered sequences in this fashion, and this is later used inside the cost
    # function to properly reconvert the parameters to their accepted value
    # bounds = [(float(b[:2])) for b in bounds_dict.values()]
    bounds = [tuple(np.asarray(b[:2],dtype=float)) for b in bounds_dict.values()]
    scaling = [b[-1] for b in bounds_dict.values()]

    return bounds,scaling

def problem_solution(popsol,problem,model='ei-network',normalized=True,bounds=None,scaling=None):
    """
    Return the solution from the optimizer in an orderly (standard) fashion for proper handling. The solution assignment
    MUST reflect the case-by-case separations for variable investigation in the objective function routines ei_cv_function
    and ng_cv_function.

    Input parameters:
    - popsol     : Last-evolved population class with best solution
    - problem    : Problem class fed into the evolution algorithm (i.e. cv_problem in ei_parallel.py or ng_parallel.py)
    - model      : {'ei-network'} | 'ng-network'
    - normalized : {True} | False   Flag specifying whether the problem was normalized or not in the evolution

    Return:
    - sol : Array   Problem-specific ordered solution.
      *if model=='lra' : sol = [fit, d1, d2, d3, d5, a2]

    REMINDER : Any change in the assignment of variables in the objective functions MUST be reflected in this solution
    routine as well.

    v1.0
    Maurizio De Pitta', INRIA Rhone-Alpes, Villeurbanne, November 3, 2017.
    """

    if model=='lra':
        N_var = problem_dimension(model='lra') # Number of variables to include in the solution
        N_par = 0 # Number of parameters to include in the solution
        sol = np.zeros(N_par+N_var+1,dtype=float) # [fit, d1, d2, d3, d5, a2]
        # The following is valid for all configurations
        # Retrieve original CV for NG network
        sol[0] = popsol.population.champion.f[0]
        # Retrieve rescaled values for solution
        x_rescaled = rescale_vector(popsol.population.champion.x,model=model,normalized=normalized,bounds=bounds,scaling=scaling)
        for i in xrange(N_var):
            sol[i+1] = x_rescaled[i]

    return sol

# -------------------------------------------------------------------------------------------------
# Cost functions for optimization problems
# -------------------------------------------------------------------------------------------------
def fit_popen(x,
              p_ca,p_ip3,
              normalized=False,
              bounds=None,
              scaling=None):

    # This function minimizes p_open vs. Ca2+ and IP3 from Ramos-Franco et al.(Biophys. J. 2000)
    # Rescaled variables properly, according to normalization and config flags
    x_rescaled = rescale_vector(x,model='lra',normalized=normalized,bounds=bounds,scaling=scaling)
    d1 = x_rescaled[0]
    d2 = x_rescaled[1]
    d3 = x_rescaled[2]
    d5 = x_rescaled[3]
    a2 = x_rescaled[4]

    # print d1,d2,d3,d5,a2

    astro = models.Astrocyte(model='lra',
                             d1=d1,d2=d2,d3=d3,d5=d5,a2=a2)
    # Retrieve Ca2+ points
    x_ca = p_ca[0]
    x_ip3 = p_ip3[0]

    # # p_open vs. IP3 was obtained with 0.25 uM free Ca2+
    astro.popen(0.25,x_ip3)
    po_ip3 = astro.po
    # print po_ip3

    # p_open vs. Ca2+ was obtained with 1 uM IP3
    astro.popen(x_ca,1.0)
    po_ca = astro.po
    # print po_ca

    # LMSE of the two curves
    fc = np.sum((po_ip3-p_ip3[1])**2) + np.sum((po_ca-p_ca[1])**2)
    # fc = np.sum((po_ip3-p_ip3[1])**2)
    # fc = np.sum((po_ca - p_ca[1]) ** 2)

    del astro  # Clean object
    return fc

def fit_popen2(x,
               p_ca, p_ip3,
               normalized=False,
               bounds=None,
               scaling=None):
    # This function minimizes p_open vs. Ca2+ and IP3 from Ramos-Franco et al.(Biophys. J. 2000)
    # Rescaled variables properly, according to normalization and config flags
    x_rescaled = rescale_vector(x, model='lra2', normalized=normalized, bounds=bounds, scaling=scaling)
    d1 = x_rescaled[0]
    d2 = x_rescaled[1]
    d5 = x_rescaled[2]
    a2 = x_rescaled[3]

    # print d1,d2,d3,d5,a2

    astro = models.Astrocyte(model='lra2',
                             d1=d1, d2=d2, d5=d5, a2=a2)
    # Retrieve Ca2+ points
    x_ca = p_ca[0]
    x_ip3 = p_ip3[0]

    # # p_open vs. IP3 was obtained with 0.25 uM free Ca2+
    astro.popen(0.25, x_ip3)
    po_ip3 = astro.po
    # print po_ip3

    # p_open vs. Ca2+ was obtained with 1 uM IP3
    astro.popen(x_ca, 1.0)
    po_ca = astro.po
    # print po_ca

    # LMSE of the two curves
    fc = np.sum((po_ip3 - p_ip3[1]) ** 2) + np.sum((po_ca - p_ca[1]) ** 2)
    # fc = np.sum((po_ip3-p_ip3[1])**2)
    # fc = np.sum((po_ca - p_ca[1]) ** 2)

    del astro  # Clean object
    return fc

def fit_ca(x,
           params,
           ca_exp,
           model='chi',
           normalized=False,
           bounds=None,
           scaling=None):

    # This function minimizes the SE of normalized calcium traces and solution from integration of the chi model
    # Rescaled variables properly, according to normalization and config flags

    # Retrieve fixed parameters
    d1 = params[0]
    d2 = params[1]
    d3 = params[2]
    d5 = params[3]
    a2 = params[4]
    c0 = params[5]

    # Rescale from normalized to regular values
    x_rescaled = rescale_vector(x,model=model,normalized=normalized,bounds=bounds,scaling=scaling)
    if model=='lra_fit':
        # Retrieve independent variables
        rc     = x_rescaled[0]
        ver    = x_rescaled[1]
        ip3    = x_rescaled[2]
        C0     = x_rescaled[3]
        h0     = x_rescaled[4]

        # Set new ICs (ca_trace is a handle to a spline object)
        ICs = np.asarray([C0, h0], dtype=float)

        # Generate astrocyte model
        astro = models.Astrocyte(model='lra',
                                 c1=0.5, rl=0.1, c0 = c0, Ker=0.1,
                                 d1=d1, d2=d2, d3=d3, d5=d5, a2=a2, ICs=ICs,
                                 rc=rc, ver=ver, ip3=ip3)
    elif model=='chi':
        # Retrieve case-specific fixed parameters
        rc  = params[6]
        ver = params[7]
        iset = int(params[8])

        # Assign rescaled variables
        vbeta  = x_rescaled[0]
        vdelta = x_rescaled[1]
        v3k    = x_rescaled[2]
        r5p    = x_rescaled[3]
        C0     = x_rescaled[4]
        h0     = x_rescaled[5]
        I0     = x_rescaled[6]

        # Set new ICs (ca_trace is a handle to a spline object)
        ICs = np.asarray([I0,C0,h0], dtype=float)
        # Generate astrocyte model
        astro = models.Astrocyte(model='chi',
                                 c1=0.5, rl=0.1, c0=c0, Ker=0.1, rc=rc, ver=ver,
                                 d1=d1, d2=d2, d3=d3, d5=d5, a2=a2, ICs=ICs,
                                 vbeta=vbeta,vdelta=vdelta, v3k=v3k, r5p=r5p)
    elif (model!='pkc2') and ('pkc' in model):
        rc     = params[6]
        ver    = params[7]
        r5p    = params[8]
        rtot   = params[9]

        yrel   = x_rescaled[0]
        vbeta  = x_rescaled[1]
        vdelta = x_rescaled[2]
        v3k = x_rescaled[3]
        vkp = x_rescaled[4]
        vk  = x_rescaled[5]
        vd  = x_rescaled[6]
        # ICs
        C0  = x_rescaled[7]
        h0  = x_rescaled[8]
        I0  = x_rescaled[9]
        g0  = x_rescaled[10]
        D0  = x_rescaled[11]
        P0  = x_rescaled[12]

        # Set new ICs (ca_trace is a handle to a spline object)
        ICs = np.asarray([g0, I0, C0, h0, D0, P0], dtype=float)
        # Generate astrocyte model
        astro = models.Astrocyte(model='pkc',
                                 c1=0.5, rl=0.1, c0=c0, Ker=0.1, rc=rc, ver=ver,
                                 d1=d1, d2=d2, d3=d3, d5=d5, a2=a2,
                                 vbeta=vbeta, vdelta=vdelta, v3k=v3k, r5p=r5p,yrel=yrel,
                                 vkp=vkp,vk=vk,vd=vd,rtot=rtot,
                                 ICs=ICs)
    elif model=='pkc2':
        rc     = params[6]
        ver    = params[7]
        vbeta  = params[8]
        vdelta = params[9]
        v3k    = params[10]
        r5p    = params[11]
        rtot   = params[12]

        yrel   = x_rescaled[0]
        Kkc    = x_rescaled[1]
        vkp = x_rescaled[2]
        vk  = x_rescaled[3]
        vd  = x_rescaled[4]
        # ICs
        C0  = x_rescaled[5]
        h0  = x_rescaled[6]
        I0  = x_rescaled[7]
        g0  = x_rescaled[8]
        D0  = x_rescaled[9]
        P0  = x_rescaled[10]

        # Set new ICs (ca_trace is a handle to a spline object)
        ICs = np.asarray([g0, I0, C0, h0, D0, P0], dtype=float)
        # Generate astrocyte model
        astro = models.Astrocyte(model='pkc',
                                 c1=0.5, rl=0.01, c0=c0, Ker=0.1, rc=rc, ver=ver,
                                 d1=d1, d2=d2, d3=d3, d5=d5, a2=a2,
                                 vbeta=vbeta, vdelta=vdelta, v3k=v3k, r5p=r5p,yrel=yrel,
                                 Kkc=Kkc,vkp=vkp,vk=vk,vd=vd,rtot=rtot,
                                 ICs=ICs)


    # Integrate
    if model in ['lra_fit','chi']:
        tfin = ca_exp['time'][0][-1]
    elif 'pkc' in model:
        tfin = ca_exp['time'][-1]

    # 'lra' and 'chi'
    ATOL = 1e-8
    RTOL = 1e-6
    # 'pkc'
    ATOL = 1e-4
    RTOL = 1e-4

    options = su.solver_opts(t0=0.0, tfin=tfin, dt=1e-2, atol=ATOL, rtol=RTOL, method="gsl_msadams")
    astro.integrate(algparams=options,normalized=True)
    # Generate spline from solution for later comparison with experimental data
    ca_sol = interpolate.UnivariateSpline(astro.sol['ts'], astro.sol['ca'], k=3, s=0)
    pkc_sol = interpolate.UnivariateSpline(astro.sol['ts'], astro.sol['pkc'], k=3, s=0)

    if model=='lra_fit':
        # Compute max/min
        # _, md = relextrema(ca_exp['ca'])
        mns, mxs = relextrema(astro.sol['ca'])
        index = np.sort(np.r_[ca_exp['mean_imin'],ca_exp['mean_imax']])
        # Fitness on max/min + number of oscillations
        fc = np.sum((ca_sol(ca_exp['time'][0][index])-ca_exp['mean_smoothed'][index])**2)
        fc += 1000.*((np.size(mxs)-np.size(ca_exp['mean_imax']))**2 + (np.size(mns)-np.size(ca_exp['mean_imin']))**2)
        # fc += 1000. * ((np.size(mxs) - np.size(ca_exp['mean_imax'])) ** 2)
        # Add mean oscillation amplitude
        # fc += 1500.*(1.-(np.mean(ca_exp['mean_smoothed'][ca_exp['mean_imax']]) - np.mean(ca_exp['mean_smoothed'][ca_exp['mean_imin']]))/(np.mean(astro.sol['ca'][mxs])-np.mean(astro.sol['ca'][mns])))**2
        # Add distance between oscillation amplitudes
        if (np.size(mxs)==np.size(ca_exp['mean_imax'])) and (np.size(mns)==np.size(ca_exp['mean_imin'])):
            fc += 2000.*np.sum((osc_amplitude(ca_exp['mean_smoothed'],ca_exp['mean_imin'],ca_exp['mean_imax'])-osc_amplitude(astro.sol['ca'],mns,mxs))**2)
        else:
            fc *= 1e6
        # print fc
    elif model=='chi':
        # (Same as 'lra_fit' but considering iset)
        mns, mxs = relextrema(astro.sol['ca'])
        index = np.sort(np.r_[ca_exp['imin'][iset],ca_exp['imax'][iset]])
        # Fitness on max/min + number of oscillations
        fc = np.sum((ca_sol(ca_exp['time'][iset][index])-ca_exp['smoothed'][iset][index])**2)
        fc += 1000.*((np.size(mxs)-np.size(ca_exp['imax'][iset]))**2 + (np.size(mns)-np.size(ca_exp['imin'][iset]))**2)
        # fc += 1000. * ((np.size(mxs) - np.size(ca_exp['imax'][iset])) ** 2)
        # Add mean oscillation amplitude
        # fc += 1500.*(1.-(np.mean(ca_exp['smoothed'][ca_exp['imax'][iset]]) - np.mean(ca_exp['smoothed'][ca_exp['imin'][iset]]))/(np.mean(astro.sol['ca'][mxs])-np.mean(astro.sol['ca'][mns])))**2
        # Add distance between oscillation amplitudes and phase difference
        if (np.size(mxs)==np.size(ca_exp['imax'][iset])) and (np.size(mns)==np.size(ca_exp['imin'][iset])):
            fc += 2000.*np.sum((osc_amplitude(ca_exp['smoothed'][iset],ca_exp['imin'][iset],ca_exp['imax'][iset])-osc_amplitude(astro.sol['ca'],mns,mxs))**2)
            fc += 2000.*np.sum(np.abs(ca_exp['time'][iset][ca_exp['imax'][iset]]- astro.sol['ts'][mxs]))
        else:
            fc *= 1e6

    del astro  # Clean object
    return fc

# -------------------------------------------------------------------------------------------------
# Problem classes
# -------------------------------------------------------------------------------------------------
class fit_problem:
    # def __init__(self, model,config,pars,**kwargs):
    def __init__(self, model, pars):
        # First we call the constructor of the base class telling
        # essentially to PyGMO what kind of problem to expect (1 objective, 0 contraints etc.)
        self.__model = model
        self.__normalized = True  # Parameters are in the unitary hypercube of self.__dim dimensions

        if model in ['lra','lra2']:
            # Parameter assignment
            self.p_ca = pars[0]
            self.p_ip3 = pars[1]
        elif model in ['lra_fit','chi']:
            self.params = pars[0]
            self.ca_trace = pars[1]

        # System dimension
        self.__dim = problem_dimension(model)

        # Get boundary dictionary
        bounds_dict = models.model_bounds(model=model)
        # Save effective bounds and scaling
        self.bounds, self.scaling = bounds_for_exploration(bounds_dict, normalized=False)


        # Generate effective bounds for exploration (normalized)
    def get_bounds(self):
        lb = tuple([0.0] * self.__dim)
        ub = tuple([1.0] * self.__dim)
        return (lb, ub)

    # Reimplementation of the virtual method that defines the objective function
    def fitness(self, x):
        # In-Class Wrapper of objective function
        if self.__model == 'lra':
            fobj = fit_popen(x,
                             self.p_ca, self.p_ip3,
                             normalized=self.__normalized,
                             bounds=self.bounds,
                             scaling=self.scaling)
        elif self.__model == 'lra2':
            fobj = fit_popen2(x,
                             self.p_ca, self.p_ip3,
                             normalized=self.__normalized,
                             bounds=self.bounds,
                             scaling=self.scaling)
        elif self.__model in ['lra_fit','chi']:
            fobj = fit_ca(x,
                          self.params,
                          self.ca_trace,
                          model=self.__model,
                          normalized=self.__normalized,
                          bounds=self.bounds,
                          scaling=self.scaling)

        if 'pkc' not in self.__model:
            # Return a tuple with one element only. In PyGMO the objective functions must return tuples.
            return (fobj,)
        else:
            # 'PKC' model returns already a tuple
            return fobj

    def get_nobj(self):
        # Specify number of objective functions to maximize
        if 'pkc' not in self.__model:
            return 1
        else:
            return 6
# -------------------------------------------------------------------------------------------------
# Effective optimization problems
# -------------------------------------------------------------------------------------------------
def popen_optimize(method='de-simple'):
    # Load experimental data
    p_open = svu.loaddata('../data/p_open.pkl')[0]
    # Get boundary dictionary
    # bounds_dict = models.model_bounds(model='lra')
    bounds_dict = models.model_bounds(model='lra2')
    # Save effective bounds and scaling
    bounds, scaling = bounds_for_exploration(bounds_dict, normalized=False)
    args = (p_open['ca'],p_open['ip3'])

    #-----------------------------------------------------------------------------------------------------------
    # Effective evolution
    #-----------------------------------------------------------------------------------------------------------
    # Size parameters
    # N_var = problem_dimension(model='lra')
    N_var = problem_dimension(model='lra2')
    N_individuals  = 100
    N_generations = 500*N_var
    N_evolutions  = 50

    # Tolerances
    FTOL = 1e-3
    XTOL = 1e-3

    # prob = pg.problem(fit_problem('lra',args))
    prob = pg.problem(fit_problem('lra2',args))

    # Optimization
    # NOTE: We keep each algorithm self contained here because different setting of parameters for evolution are
    # needed depending on the chose algorithm
    if method=='de-simple':
        ##-----------------
        # # Single node
        ##-----------------
        # algo = pg.algorithm(pg.de(gen=N_generations,F=0.75,CR=0.9,variant=6,ftol=FTOL,xtol=XTOL))
        # algo = pg.algorithm(pg.de(gen=N_generations, ftol=1e-12, tol=1e-6))
        # algo = pg.algorithm(pg.pso(gen=N_generations))
        algo = pg.algorithm(pg.bee_colony(gen=N_generations))
        algo.set_verbosity(100)
        pop = pg.population(prob, size=N_individuals)
        pop = algo.evolve(pop)

        best_fitness = pop.get_f()[pop.best_idx()]
        print(best_fitness)

        x_rescaled = rescale_vector(pop.champion_x,model='lra2',normalized=True,bounds=bounds,scaling=scaling)
        astro = models.Astrocyte(model='lra2',
                          d1=x_rescaled[0],
                          d2=x_rescaled[1],
                          d5=x_rescaled[2],
                          a2=x_rescaled[3])
        astro.popen(p_open['ca'][0],1.0)
        po_ca = astro.po
        astro.popen(0.25, p_open['ip3'][0])
        po_ip3 = astro.po
        print x_rescaled

        svu.savedata([x_rescaled],'../data/fit_po.pkl')
        plt.semilogx(p_open['ca'][0],p_open['ca'][1],'k-',p_open['ca'][0],po_ca,'ro',
                     p_open['ip3'][0],p_open['ip3'][1],'m-',p_open['ip3'][0],po_ip3,'bo')
        plt.show()

def lra_optimize(method='pso',N_ind=100,N_gen=30,
                 parallel=False,
                 show_results=False,
                 dir='../data/',file_name='lra_fit'):
    # Load experimental data
    ca_data = svu.loaddata('../data/herzog_data.pkl')[0]
    po_data = svu.loaddata('../data/fit_po.pkl')[0]

    # Get boundary dictionary
    bounds_dict = models.model_bounds(model='lra_fit')
    # Save effective bounds and scaling
    bounds, scaling = bounds_for_exploration(bounds_dict, normalized=False)

    # Prepare model parameters
    c0 = 5.0
    params = np.asarray(np.r_[po_data[[0,1,0,2,3]],c0])
    args = (params, ca_data)

    #-----------------------------------------------------------------------------------------------------------
    # Effective evolution
    #-----------------------------------------------------------------------------------------------------------
    # Size parameters
    N_individuals = N_ind
    N_var = problem_dimension(model='lra_fit')
    N_generations = N_gen * N_var

    # Tolerances
    FTOL = 1e-6
    XTOL = 1e-6

    # prob = pg.problem(fit_problem('lra',args))
    prob = pg.problem(fit_problem('lra_fit',args))

    # Optimization
    # NOTE: We keep each algorithm self contained here because different setting of parameters for evolution are
    # needed depending on the chose algorithm
    if method=='pso':
        algo = pg.algorithm(pg.pso(gen=N_generations,variant = 5, neighb_type = 4))
    elif method=='bee':
        algo = pg.algorithm(pg.bee_colony(gen=N_generations,limit=1))
    elif method=='de':
        algo = pg.algorithm(pg.de(gen=N_generations, ftol=FTOL, tol=XTOL))

    if not parallel:
        # Single-node optimation to run on local machine / laptop
        verbosity = 20

        algo.set_verbosity(verbosity)
        pop = pg.population(prob, size=N_individuals)
        pop = algo.evolve(pop)

        best_fitness = pop.get_f()[pop.best_idx()]
        print(best_fitness)
        x_rescaled = rescale_vector(pop.champion_x,model='lra_fit',normalized=True,bounds=bounds,scaling=scaling)

        # Show results of fit
        if show_results:
            astro = models.Astrocyte(model='lra',
                                     d1=po_data[0], d2=po_data[1], d3=po_data[0], d5=po_data[2], a2=po_data[3],
                                     c0=x_rescaled[5], c1=0.5, rc=x_rescaled[0], rl=0.1, ver=x_rescaled[1], Ker=0.1,
                                     ip3=x_rescaled[2],
                                     ICs=np.asarray([x_rescaled[3], x_rescaled[4]]))
            options = su.solver_opts(t0=ca_data['time'][0][0], tfin=ca_data['time'][0][-1], dt=1e-4, atol=1e-8, rtol=1e-6, method="gsl_msadams")
            astro.integrate(algparams=options, normalized=True)

            ca_trace = ca_data['mean_smoothed']
            plt.plot(ca_data['time'][0],ca_trace,'k-',
                     astro.sol['ts'],astro.sol['ca'],'r-')
            plt.show()

    else:
        # Parallel version to run on cluster
        N_evolutions = 100
        N_islands = 10
        verbosity = 100
        # Initiate optimization
        archi = pg.archipelago(n=N_islands, algo=algo, prob=prob, pop_size=N_individuals)
        archi.evolve(N_evolutions)
        archi.wait()
        imin = np.argmin(archi.get_champions_f())
        x_rescaled = rescale_vector(archi.get_champions_x()[imin],model='chi',normalized=True,bounds=bounds,scaling=scaling)
        print archi.get_champions_f()[imin]
        print x_rescaled
        # svu.savedata([x_rescaled],'../data/ca_trace_0.pkl')

    # Save data from fit simulation
    svu.savedata([x_rescaled],dir+file_name+'.pkl')

def chi_optimize(method='pso',
                 parallel=False,
                 N_ind=100,N_gen=30,
                 iset=0,
                 show_results=False,
                 dir='../data/', file_name='chi_fit'):
    # Load experimental data
    N_set = 2
    ca_data = svu.loaddata('../data/herzog_data.pkl')[0]
    po_data = svu.loaddata('../data/fit_po.pkl')[0]
    lr_data = svu.loaddata('../data/fit_lra.pkl')[0]

    # Get boundary dictionary
    bounds_dict = models.model_bounds(model='chi')
    # Save effective bounds and scaling
    bounds, scaling = bounds_for_exploration(bounds_dict, normalized=False)

    # Prepare model parameters
    c0 = 5.0
    # c0 = 15.0
    params = np.asarray(np.r_[po_data[[0,1,0,2,3]], lr_data[[0,1]], c0, iset])
    args = (params,ca_data)

    #-----------------------------------------------------------------------------------------------------------
    # Effective evolution
    #-----------------------------------------------------------------------------------------------------------
    # Size parameters
    # N_var = problem_dimension(model='lra')
    N_individuals = N_ind
    N_var = problem_dimension(model='chi')
    N_generations = N_gen * N_var

    # Tolerances
    FTOL = 1e-6
    XTOL = 1e-6

    # prob = pg.problem(fit_problem('lra',args))
    prob = pg.problem(fit_problem('chi',args))

    # Optimization
    # NOTE: We keep each algorithm self contained here because different setting of parameters for evolution are
    # needed depending on the chose algorithm
    if method=='pso':
        algo = pg.algorithm(pg.pso(gen=N_generations,variant = 5, neighb_type = 4, max_vel=0.8))
    elif method=='bee':
        algo = pg.algorithm(pg.bee_colony(gen=N_generations,limit=2))
    elif method=='de':
        algo = pg.algorithm(pg.de(gen=N_generations, ftol=FTOL, tol=XTOL))
        # Single-node optimation to run on local machine / laptop
        # N_individuals = 100
        # N_generations = 40 * N_var
        # verbosity = 20

    if not parallel:
        verbosity = 20
        algo.set_verbosity(verbosity)
        pop = pg.population(prob, size=N_individuals)
        pop = algo.evolve(pop)

        best_fitness = pop.get_f()[pop.best_idx()]
        print(best_fitness)

        x_rescaled = rescale_vector(pop.champion_x,model='chi',normalized=True,bounds=bounds,scaling=scaling)

        # Show results of fit
        if show_results:
            astro = models.Astrocyte(model='chi',
                                     d1=po_data[0], d2=po_data[1], d3=po_data[0], d5=po_data[2], a2=po_data[3],
                                     c0=c0, c1=0.5, rl=0.1, Ker=0.1, rc=lr_data[0], ver=lr_data[1],
                                     vbeta=x_rescaled[0], vdelta=x_rescaled[1], v3k=x_rescaled[2], r5p=x_rescaled[3],
                                     ICs=np.asarray([x_rescaled[6], x_rescaled[4], x_rescaled[5]]))

            options = su.solver_opts(t0=ca_data['time'][iset][0], tfin=ca_data['time'][iset][-1], dt=1e-4, atol=1e-8, rtol=1e-6, method="gsl_msadams")
            astro.integrate(algparams=options, normalized=True)

            ca_trace = ca_data['smoothed'][iset]
            plt.plot(ca_data['time'][0],ca_trace,'k-',
                     astro.sol['ts'],astro.sol['ca'],'r-')
            plt.show()

    else :
        # Parallel version to run on cluster
        N_evolutions = 100
        N_islands = 10
        # Initiate optimization
        algo = pg.algorithm(pg.pso(gen=N_generations,variant = 5, neighb_type = 4, max_vel=0.8))
        archi = pg.archipelago(n=N_islands, algo=algo, prob=prob, pop_size=N_individuals)
        archi.evolve(N_evolutions)
        archi.wait()
        imin = np.argmin(archi.get_champions_f())
        x_rescaled = rescale_vector(archi.get_champions_x()[imin], model='chi', normalized=True, bounds=bounds,
                                    scaling=scaling)
        print archi.get_champions_f()[imin]
        print x_rescaled

    svu.savedata([x_rescaled],dir+file_name+'.pkl')

#-----------------------------------------------------------------------------------------------------------------------
# Attempt to fit by scipy
#-----------------------------------------------------------------------------------------------------------------------
def relextrema(data,thr=20):
    """
    Retrieve relative max and minima from calcium data. Works only when the oscillations are regular.
    :param data:
    :return:
    """
    mx = sps.argrelextrema(data,np.greater)[0]
    mn = sps.argrelextrema(data,np.less)[0]
    if (np.size(mx)>0) and (np.size(mn)>0):
        imx = [mx[0]]
        imn = [mn[0]]
    else :
        return [], []
    # Compute relative maxima
    df = 0  # difference between point positions
    for i in xrange(1,np.size(mx)):
        if np.diff(mx[[df,i]]) >= thr:
            df = i
            imx.append(mx[i])
        else :
            if data[mx[i]] > data[mx[i-1]]:
                del imx[-1]
                imx.append(mx[i])
    # Compute relative minima
    df = 0  # difference between point positions
    for i in xrange(1,np.size(mn)):
        if np.diff(mn[[df,i]]) >= thr:
            df = i
            imn.append(mn[i])
        else :
            if data[mn[i]] < data[mn[i-1]]:
                del imn[-1]
                imn.append(mn[i])
    return imn, imx

def osc_amplitude(ca_data,imin,imax):
    """
    Estimates amplitude of individual oscillations
    :param ca_data:
    :param imin:
    :param imax:
    :return:
    """
    sx = np.size(imax)
    sm = np.size(imin)
    ds = np.abs(sm-sx)
    if sx > sm :
        dist = ca_data[imax[ds:]]-ca_data[imin]
    elif sx < sm :
        dist = ca_data[imax]-ca_data[imin[ds:]]
    else :
        dist = ca_data[imax]-ca_data[imin]
    return np.asarray(dist, dtype=float)

if __name__ == "__main__":
    #----------------------------------------------------------------
    # Open probability
    #----------------------------------------------------------------
    # popen_optimize(method='de-simple')

    # Handling of Input arguments passed from command line
    if np.size(sys.argv)>1:
        # Arguments passed by command line
        try:
            opts,args = getopt.getopt(sys.argv[1:],"hm:a:p:i:g:s:d:f:",
                                      ["--model=","--algorithm=","--parallel=","--individuals=","--generations=","--dset=","--save-dir=","--file-name="])
            print opts, np.size(opts)
            if np.size(opts)==0: raise getopt.GetoptError
            for opt,arg in opts:
                if opt=='-h':
                    print 'Usage: fit_data -m <model> -a <method> -p <parallel-flag> -i <individuals> -g <generations> -s <data_set> -d <save-dir> -f <file-name>'
                    sys.exit()
                elif opt in ("-m","--model"):
                    model = arg
                elif opt in ("-a","--algorithm"):
                    method = arg
                elif opt in ("-p", "--parallel"):
                    parallel = bool(int(arg))
                elif opt in ("-i","--individuals"):
                    N_ind = int(arg)
                elif opt in ("-g","--generations"):
                    N_gen = int(arg)
                elif opt in ("-s","--dset"):
                    iset = int(arg)
                elif opt in ("-d","--save-dir"):
                    save_dir = arg
                elif opt in ("-f","--file-name"):
                    file_name = arg
        except getopt.GetoptError, exc:
            print exc.msg
            print 'Usage: fit_data -m <model> -a <method> -p <parallel-flag> -i <individuals> -g <generations> -s <data_set> -d <save-dir> -f <file-name>'
            sys.exit()
    else:
        # No input arguments are given in this case (default)
        print '(default run): fit_data.py -m \'lra\' -a \'pso\' -p 0 -i 100 -g 30 -s 0 -d \'../\' -f \'lra_fit_0\' '
        # model = 'lra_fit'
        # method = 'pso'
        # parallel = False
        # N_ind = 100
        # N_gen = 30
        # iset = 0
        # save_dir = '../data/' # Alternatively (on RCC): '/home/mdepitta/Dropbox/Stability_Analysis_network/'
        # file_name = 'lra_fit' + '_' + method

        model = 'pkc2'
        # method = 'pso'
        method = 'moead'
        # method = 'nsga2'
        parallel = False
        N_ind = 21
        N_gen = 3
        iset = 1
        save_dir = '../data/' # Alternatively (on RCC): '/home/mdepitta/Dropbox/Stability_Analysis_network/'
        file_name = model+ '_fit' + '_' + method

    gc.collect()
    if model=='lra_fit':
        #----------------------------------------------------------------
        # LRA optimization on mean set parameters
        #----------------------------------------------------------------
        lra_optimize(method=method, parallel=parallel,
                     N_ind=N_ind, N_gen=N_gen,
                     show_results=False, dir=save_dir, file_name=file_name+ '_' + method)

    elif model=='chi':
        #----------------------------------------------------------------
        # ChI optimization
        #----------------------------------------------------------------
        chi_optimize(method=method, parallel=parallel,
                     N_ind=N_ind, N_gen=N_gen,
                     iset=iset,
                     show_results=False,
                     dir=save_dir, file_name=file_name+ '_' + method)
        # chi_optimize(method='pso-parallel')
    gc.collect()