'''
astrocyte_models.py

Library of astrocyte models.
- chi_model       : CHI model simulator;
- gchi_model      : GCHI model simulator;
- gchi_norm_model : GCHI model simulator w/ normalized parameters (for CMAES)

Maurizio De Pitta', The University of Chicago, Feb 27th, 2015.
'''
from __future__ import division
import numpy as np
import copy as cp
from scipy import *
import weave
from weave import converters
import sys, os
base_dir = '/home/maurizio/Ch10.DePitta'
sys.path.append(os.path.join(os.path.expanduser('~'),base_dir))
import pycustommodules.general_utils as gu
import pycustommodules.solvers.solver_utils as su
import pycustommodules.spk_generators as spg

import matplotlib.pylab as plt

###############################################################################
# General functions
###############################################################################
normalize = lambda x : (x - x.min())/(x.max()-x.min())
unnormalize = lambda xn, mm : mm[0] + xn*(np.diff(mm)) # mm must be in [min, max]

# Some lambdas to reconstruct solutions
xsol = lambda xn, dt, tau: 1 + (xn - 1) * np.exp(-dt / tau)
usol = lambda un, dt, tau: un * np.exp(-dt / tau)

def Hill(x,k,n):
    return np.asarray(x**n/(x**n+k**n), dtype=float)

###############################################################################
# Parameters definition
###############################################################################
def lra_parameters(**kwargs):
    """
    Parameters for the Li-Rinzel astrocyte (Li and Rinzel, JTB 1994).

    Input:
    - noise : {False} | True  : tell to the simulator whether or not to simulate the stochastic LR model (Default: False)

    Return: Dictionary of model parameters.

    Maurizio De Pitta', Basque Center of Applied Mathematics, June 8, 2017.
    """
    pars = {'d1'    : 0.1,
            'd2'    : 2.1,
            'd3'    : 0.9967,
            'd5'    : 0.2,
            'a2'    : 0.4,
            'c1'    : 0.4,
            'c0'    : 4,
            'rc'    : 7,
            'rl'    : 0.05,
            'ver'   : 0.9,
            'Ker'   : 0.1,
            'ip3'   : 0.1,
            'noise' : False,
            'ICs'   : [0.05,0.99]
            }

    ## User-defined parameters
    pars = gu.varargin(pars,**kwargs)
    # Parameters must be floats
    for k,item in pars.iteritems():
        if isscalar(item):
            pars[k] = float(item)
        elif isinstance(item, bool):
            pars[k] = np.intc(item)
        else:
            pars[k] = array(item,dtype=float)
    return pars

def exocytosis_parameters(**kwargs):
    """
    Parameters for the Tsodyks-Markram model of exocytosis and synaptic currents (1- or 2-equations) (Tsodyks, Les Houches  2005).

    Input:
    - ICs : Must be an array-like of three elements for (x,y,u).

    Return: Dictionary of model parameters.

    Maurizio De Pitta', Basque Center for Applied Mathematics, June 11, 2018.
    """
    pars = {'u0'   : 0.5,     # Basal release probability
            'taud' : 0.5,     # Depression time constant
            'tauf' : 1.0,     # Facilitation time constant
            'taui' : 0.02,    # Time decay of AMPA-like currents
            'psc0' : 30,      # Max PSC in [pA]
            'ICs'  : [1.,0.,0.]  #x,y,u
            }
    ## User-defined parameters
    pars = gu.varargin(pars,**kwargs)
    for k,item in pars.iteritems():
        if isscalar(item):
            pars[k] = float(item)
        else:
            pars[k] = array(item,dtype=float)
    return pars

def stimulus_params(**kwargs):
    """
    Build pars vector with stimulus and system details (used in gtrs class for calls to models including exocytosis)

    Input:
    **kwargs:
        - rate         : Array-like or scalar for rate of presynaptic APs
        - twin         : Array-like [t0,tf] for duration of presynaptic stimulation
        - stimulus     : 'poisson' | 'periodic' | 'fixed' (requires specification of spikes_pre)
        - spikes_pre   : {None} | array-like with AP instants to pass to the simulator (if stimulus='fixed')
        - Nsyn         : Number of synapses (>0)
        - Nvar         : Number of variables in the synaptic model (1 or 2)
        - rate_gtr     : Array-like or scalar for rate of GREs
        - twin_gtr     : Array-like [t0,tf] for duration of gliotransmitter release
        - stimulus_gtr : 'poisson' | 'periodic' | 'fixed' (requires specification of pre_gtr)
        - pre_gtr      : {None} | array-like with AP instants to pass to the simulator (if stimulus_gtr='fixed')

    Return: Dictionary of model parameters.

    Maurizio De Pitta', Basque Center for Applied Mathematics, June 12, 2018.
    """
    pars = {'rate'      : 1.,
            'twin'      : [0.,10.],
            'stimulus'  : 'poisson',
            'spikes_pre': None,
            'Nsyn'      : 1,
            'Nvar'      : 2,
            # GRE parameters
            'rate_gtr'  : None,
            'twin_gtr'  : None,
            'stimulus_gtr' : None,
            'pre_gtr'      : None
             }

    ## User-defined parameters
    pars = gu.varargin(pars,**kwargs)
    # Parameters must be floats
    for k,item in pars.iteritems():
        if k in ['Nsyn','Nvar']:
            pars[k] = int(item)
        elif k=='rate':
            if np.isscalar(item):
                pars[k] = float(item)
            else:
                pars[k] = np.array(item,dtype=float)
        elif k=='twin':
            if not np.isscalar(pars['rate']):
                assert len(pars['rate'])==len(item), "twin must be of the same size of rate"
            pars[k] = np.array(item, dtype=float)
        elif (k=='rate_gtr') and (item!=None):
            if np.isscalar(item):
                pars[k] = float(item)
            else:
                pars[k] = np.array(item,dtype=float)
        elif (k=='twin_gtr') and (item!=None):
            if (not pars['rate_gtr']) and (not np.isscalar(pars['rate_gtr'])):
                assert len(pars['rate_gtr'])==len(item), "twin_gtr must be of the same size of rate_gtr"
            pars[k] = np.array(item, dtype=float)
    return pars

def gtrelease_parameters(**kwargs):
    """
    Calcium-dependent gliotransmitter release.

    Input:
    **kwargs:
      - lra_parameters kwargs
      - cthr : calcium threshold for exocytosis [uM]
      - ua   : basal Gt. release probability
      - taua : recovery time constant [s]

    Return:
    - pars : Parameter dictionary

    Maurizio De Pitta', Basque Center of Applied Mathematics, June 17, 2018.
    """
    pars = {'cthr' : 0.5,
            'ua'   : 0.6,
            'taua' : 1./0.6}
    pars = gu.merge_dicts(pars,lra_parameters())
    pars['ICs'] = np.r_[pars['ICs'],1.0]
    ## User-defined parameters
    pars = gu.varargin(pars, **kwargs)
    for k,item in pars.iteritems():
        if isscalar(item):
            pars[k] = float(item)
        else:
            pars[k] = array(item,dtype=float)
    return pars

def asn_parameters(model='spk',**kwargs):
    """
    Parameter for the model of a gliotransmitter-modulated synaptic release

    Input:
    - model:
    **kwargs:
      - exocytosis_parameters **kwargs
      - gtrelease_parameters **kwargs
      - rhoe  : Volume ratio
      - Gtot  : Total Gt. vesicular concentration [mM] [!]
      - Ou    : Max uptake rate [uM/s]
      - Ku    : Transporter affinity [uM]
      - Og    : Presynaptic receptor binding rate [1/(uM.s)]
      - taug  : Presynaotic receptor decay time constant [s]
      - alpha : Effect parameter

    The value of 'js' will be computed internally and automatically included in the returned dictionary.

    Return: Dictionary of model parameters.

    Maurizio De Pitta', Basque Center for Applied Mathematics, June 18, 2017.
    """

    pars = {'rhoe' : 6.5e-4,
            'Ou'   : 0.,
            'Ku'   : 100.,
            'taue' : 1./60,
            'Gtot' : 200., # MUST BE in [mM]
            'Og'   : 1.5,
            'taug' : 30.,
            'alpha': 0.5
            }
    pars = gu.merge_dicts(pars, gtrelease_parameters(),exocytosis_parameters())
    pars['ICs'] = np.asarray([0.,0.,0.05,0.99])  # [G_A,\Gamma_S,c,h]
    pars['ICr'] = np.asarray([1,0.,0.,1.])  # [x_S,y_S,u_S,x_A]
    ## User-defined parameters
    pars = gu.varargin(pars, **kwargs)
    ## Takes only the first two elements of ICs in the MF model
    if model=='ave':
        pars['ICs'] = pars['ICs'][:2]
    if 'js' in kwargs:
        pars['js'] = kwargs['js']
    else:
        pars['js'] = pars['rhoe']*pars['Og']*1e3*pars['Gtot']*pars['taue']
    for k,item in pars.iteritems():
        if isscalar(item):
            pars[k] = float(item)
        else:
            pars[k] = array(item,dtype=float)
    # pars['Gtot'] *= 1e3 # Convert to [uM]
    return pars

###############################################################################
# Simulators
###############################################################################
def exocytosis_model(rate,twin,N_syn,
                     pars,N_var,
                     recompile=0,
                     model='spiking',
                     solver_options=None,
                     stimulus='poisson',spikes_pre=None):
    """
    Simulator of the Tsodyks-Markram model (Tsodyks, Les Houches 2005)

    Input:
     - rate           : Array-like or scalar for rate of presynaptic APs
     - twin           : Array-like [t0,tf] for duration of presynaptic stimulation
     - N_syn          : Number of independent synapses / release sites to simulate
     - pars           : exocytosis_parameters
     - N_var          : Number of variables for the exocytosis model (1 | 2)
     - recompile      : 0 | 1 (debug flag)
     - model          : 'spiking' | 'average' (mean-field)
     - solver_options : solver options (if model=='average')
     - stimulus       : 'poisson' | 'periodic' | 'fixed' (requires specification of spikes_pre)
     - spikes_pre     : {None} | array-like with AP instants to pass to the simulator (if stimulus='fixed')

    Return: Solution dictionary.

    Maurizio De Pitta', Basque Center of Applied Mathematics, June 11, 2018.
    """

    # First check that N_var is compatible
    assert N_var<3 and N_var>=1, "Number of variables of exocytosis model must be either 1 (x) or 2 (x,u)"

    # Assures that twin is Numpy array for later handling
    twin  = np.asarray(twin,dtype=float)

    # Also convert make sure to recast N_eq in a way that is suitable for C
    N_var = np.intc(N_var)
    if not N_syn: N_syn = 0  # Dummy value of N_syn in the case of the mean-field model
    N_syn = np.intc(N_syn)

    if model=='spiking':
        # Create input_spikes
        if twin.size == 2:
            # One rate or multiple ones in the same interval
            spikes = spg.input_spikes(N_syn, twin[1], rate, 0, stimulus=stimulus, spikes_pre=spikes_pre)
        else:
            # Multiple rates in different intervals
            spikes = spg.input_spikes(N_syn, twin, rate, 0, stimulus=stimulus, spikes_pre=spikes_pre)
            twin = np.r_[0, np.sum(twin)]
        N_spk = int(np.shape(spikes)[1])

        # Check that ICs are of size N_var x N_syn
        # NOTE: you will have to pass the whole vector of ICs to compute solutions from last point
        if pars['ICs'].size != (N_var+1)*N_syn:
            pars['ICs'] = np.tile(pars['ICs'][:N_var+1],(1,N_syn))[0]

        support_code = """
                     #include "gliotransmission_models.h"
                     """
        source_files = [os.path.join(os.path.expanduser('~'), base_dir + '/pycustommodules/pycapi_utils.cpp'),
                        os.path.join(os.path.expanduser('~'), base_dir + '/pycustommodules/solvers/solver_options.cpp'),
                        os.path.join(os.path.expanduser('~'), base_dir + '/pycustommodules/solvers/stochastic_solvers.cpp'),
                        os.path.join(os.path.expanduser('~'), base_dir + '/code/gliotransmission_models.cpp')]
        code = """
               // Version
               double version = 1.0;
    
               // Define astrocyte model
               release synapse(N_var,N_syn,N_spk,pars);
               synapse.set_ics(pars);
    
               // Declare output structure
               out_release out(N_spk,N_var);
    
               // Simulator
               out = synapse.simulate(spikes.data());
    
               //Output 
               return_val = out.make_PyDict();
               """
        libs = ['gsl', 'gslcblas', 'm']
        dirs = [os.path.join(os.path.expanduser('~'), base_dir + '/code/'),
                os.path.join(os.path.expanduser('~'), base_dir + '/pycustommodules'),
                os.path.join(os.path.expanduser('~'), base_dir + '/pycustommodules/solvers')]
        vars = ['pars', 'spikes','N_spk','N_var','N_syn']
        otm = weave.inline(code,
                           vars,
                           support_code=support_code,
                           sources=source_files,
                           libraries=libs,
                           library_dirs=dirs,
                           include_dirs=dirs,
                           runtime_library_dirs=dirs,
                           type_converters=converters.blitz,
                           compiler='gcc',
                           extra_compile_args=['-std=c++11'],
                           force=recompile)
        # Post-stimulus processing
        otm['spk'] = spikes[0]     # Spike instants
        otm['is']  = spikes[-1]    # Synapse indexes in the spike train
        otm['ICs'] = pars['ICs']
        if N_var>1:
            u_ = otm['u']
        else:
            u_ = None
        otm['LCs'] = last_point(pars,otm['spk'],twin[-1],otm['x'],otm['y'],otm['is'],uval=u_)
    elif model=='average':
        # Check that rate is scalar
        assert isscalar(rate), "Mean-field rate must be scalar"

        support_code = """
                     #include "gliotransmission_models.h"
                     """
        source_files = [os.path.join(os.path.expanduser('~'), base_dir + '/pycustommodules/pycapi_utils.cpp'),
                        os.path.join(os.path.expanduser('~'), base_dir + '/pycustommodules/solvers/solver_options.cpp'),
                        os.path.join(os.path.expanduser('~'), base_dir + '/pycustommodules/solvers/stochastic_solvers.cpp'),
                        os.path.join(os.path.expanduser('~'), base_dir + '/code/gliotransmission_models.cpp')]
        code = """
               // Version
               double version = 0.0;

               // Define astrocyte model
               release_ave synapse(N_var);

               // Declare output structure
               out_release out;

               // Simulator
               out = synapse.simulate(rate,pars,solver_options);

               //Output 
               return_val = out.make_PyDict();
               """
        libs = ['gsl', 'gslcblas', 'm']
        dirs = [os.path.join(os.path.expanduser('~'), base_dir + '/code/'),
                os.path.join(os.path.expanduser('~'), base_dir + '/pycustommodules'),
                os.path.join(os.path.expanduser('~'), base_dir + '/pycustommodules/solvers')]
        vars = ['rate', 'pars', 'solver_options', 'N_var']
        otm = weave.inline(code,
                           vars,
                           support_code=support_code,
                           sources=source_files,
                           libraries=libs,
                           library_dirs=dirs,
                           include_dirs=dirs,
                           runtime_library_dirs=dirs,
                           type_converters=converters.blitz,
                           compiler='gcc',
                           extra_compile_args=['-std=c++11'],
                           force=recompile)

    # Post-stimulus processing
    otm['twin'] = twin  # Simulate interval
    # Add released resources
    if (N_var < 2):
        otm['r'] = pars['u0'] * otm['x']
    else:
        otm['r'] = np.multiply(otm['x'], otm['u'])
    return otm

def lr_model(pars, solver_options,
             recompile=0):
    """
    Simulator for the plain Li-Rinzel model (Li and Rinzel, JTB 1994) and stochastic Li-Rinzel model (Shuai and Jung, Biophys. J. 2002)

    Input:
    - pars           : lra_parameters
    - solver_options : solver options
    - recompile      : 0 | 1 (debug routine)

    Return: Solution dictionary.

    Maurizio De Pitta', Basque Center of Applied Mathematics, June 11, 2018.
    """
    support_code = """
                 #include "gliotransmission_models.h"
                 """
    source_files = [os.path.join(os.path.expanduser('~'), base_dir+'/pycustommodules/pycapi_utils.cpp'),
                    os.path.join(os.path.expanduser('~'), base_dir+'/pycustommodules/solvers/solver_options.cpp'),
                    os.path.join(os.path.expanduser('~'), base_dir+'/pycustommodules/solvers/stochastic_solvers.cpp'),
                    os.path.join(os.path.expanduser('~'), base_dir+'/code/gliotransmission_models.cpp')]
    code = """
           // Version
           double version = 0.0;
           
           // Define astrocyte model
           lra astrocyte;
           
           // Declare output structure
           out_lra out;

           // Simulator
           out = astrocyte.simulate(pars,solver_options);

           //Output 
           return_val = out.make_PyDict();
           """
    libs = ['gsl', 'gslcblas', 'm']
    dirs = [os.path.join(os.path.expanduser('~'), base_dir+'/code/'),
            os.path.join(os.path.expanduser('~'), base_dir+'/pycustommodules'),
            os.path.join(os.path.expanduser('~'), base_dir+'/pycustommodules/solvers')]
    vars = ['pars', 'solver_options']
    olr = weave.inline(code,
                        vars,
                        support_code=support_code,
                        sources=source_files,
                        libraries=libs,
                        library_dirs=dirs,
                        include_dirs=dirs,
                        runtime_library_dirs=dirs,
                        type_converters=converters.blitz,
                        compiler='gcc',
                        extra_compile_args=['-std=c++11'],
                        force=recompile)
    # Post-stimulus processing
    return olr

def ca_gtrel_model(pars, solver_options,
                   recompile=0):
    """
    Simulator calcium-dependent exocytosis, i.e. Li-Rinzel model with 1-var TM model.

    Input:
    - pars           : gtrelease_parameters
    - solver_options : solver options
    - recompile      : 0 | 1 (debug routine)

    Return: Solution dictionary.

    Maurizio De Pitta', Basque Center of Applied Mathematics, June 11, 2018.
    """
    support_code = """
                 #include "gliotransmission_models.h"
                 """
    source_files = [os.path.join(os.path.expanduser('~'), base_dir + '/pycustommodules/pycapi_utils.cpp'),
                    os.path.join(os.path.expanduser('~'), base_dir + '/pycustommodules/solvers/solver_options.cpp'),
                    os.path.join(os.path.expanduser('~'), base_dir + '/pycustommodules/solvers/stochastic_solvers.cpp'),
                    os.path.join(os.path.expanduser('~'), base_dir + '/code/gliotransmission_models.cpp')]
    code = """
           // Version
           double version = 0.0;

           // Define astrocyte model
           gtrelease gliot;

           // Declare output structure
           out_gtrelease out;

           // Simulator
           out = gliot.simulate(pars,solver_options);

           //Output 
           return_val = out.make_PyDict();
           """
    libs = ['gsl', 'gslcblas', 'm']
    dirs = [os.path.join(os.path.expanduser('~'), base_dir + '/code/'),
            os.path.join(os.path.expanduser('~'), base_dir + '/pycustommodules'),
            os.path.join(os.path.expanduser('~'), base_dir + '/pycustommodules/solvers')]
    vars = ['pars', 'solver_options']
    gtr = weave.inline(code,
                       vars,
                       support_code=support_code,
                       sources=source_files,
                       libraries=libs,
                       library_dirs=dirs,
                       include_dirs=dirs,
                       runtime_library_dirs=dirs,
                       type_converters=converters.blitz,
                       compiler='gcc',
                       extra_compile_args=['-std=c++11'],
                       force=recompile)
    # Post-stimulus processing
    gtr['twin'] = np.asarray([solver_options['t0'],solver_options['tfin']],dtype=float)
    gtr['twin_gtr'] = np.asarray([solver_options['t0'], solver_options['tfin']], dtype=float) # Time window used in the reconstruction
    # Clean 'x' and provide GRE vector
    i_gre = gtr['xa']>0
    gtr['xa']   = gtr['xa'][i_gre]
    gtr['gre'] = gtr['t'][i_gre]
    # A vector with all the indexes of Gt. CONVENTION: we use negative indexes for astrocytic release. Only one release
    # site in this implementation, i.e. index -1
    gtr['ig']   = -1*np.ones(len(gtr['gre']))
    # Add released GTRs
    gtr['ra'] = pars['ua']*gtr['xa']

    return gtr

def asn_model(pars, solver_options,
              N_syn,N_var,rate_syn,twin_syn,
              model='spiking',
              stimulus_syn='poisson',pre_syn=None,
              gtr=None,
              rate_gtr=None,twin_gtr=None,
              stimulus_gtr=None,pre_gtr=None,
              recompile=0):
    """
    Simulator for the plain Li-Rinzel model (Li and Rinzel, JTB 1994) and stochastic Li-Rinzel model (Shuai and Jung, Biophys. J. 2002)

    Input:
    - pars          : asn_parameters
    - solver_options: solver options
    - N_syn         : Number of independent synapses regulated by the same astrocytic domain
    - N_var         : Number of variables in the synaptic model
    - rate_syn      : Array-like or scalar for rate of presynaptic APs
    - twin_syn      : Array-like [t0,tf] for duration of presynaptic stimulation
    - model         : 'spiking' |
    - stimulus_syn  : 'poisson' | 'periodic' | 'fixed' (requires specification of pre_syn)
    - pre_syn       : {None} | array-like with AP instants to pass to the simulator (if stimulus_pre='fixed')
    - gtr           : None | 'ca_dep' If specified uses Li-Rinzel model to find GREs (not properly tested!)
    - rate_gtr      : Array-like or scalar for rate of GREs
    - twin_gtr      : Array-like [t0,tf] for duration of gliotransmitter release
    - stimulus_gtr  : 'poisson' | 'periodic' | 'fixed' (requires specification of pre_gtr)
    - pre_gtr       : {None} | array-like with AP instants to pass to the simulator (if stimulus_gtr='fixed')
    - recompile     : 0 | 1 (debug flag)

    Return: Solution dictionary (model specific).

    Maurizio De Pitta', Basque Center of Applied Mathematics, June 21, 2018.
    """

    # First check that N_var is compatible
    assert N_var<3 and N_var>=1, "Number of variables of exocytosis model must be either 1 (x) or 2 (x,u)"

    # Assures that twin is Numpy array for later handling
    twin_syn  = np.asarray(twin_syn,dtype=float)
    twin_gtr  = np.asarray(twin_gtr,dtype=float)

    # Also convert make sure to recast N_eq in a way that is suitable for C
    N_var = int(N_var)
    if not N_syn: N_syn = 0  # Dummy value of N_syn in the case of the mean-field model
    N_syn = int(N_syn)

    if model=='spiking':
        twin = {}
        # Create input_spikes
        if twin_syn.size == 2:
            # One rate or multiple ones in the same interval
            spikes = spg.input_spikes(N_syn, twin_syn[1], rate_syn, 0, stimulus=stimulus_syn, spikes_pre=pre_syn)
        else:
            # Multiple rates in different intervals
            spikes = spg.input_spikes(N_syn, twin_syn, rate_syn, 0, stimulus=stimulus_syn, spikes_pre=pre_syn)
        twin['syn'] = np.r_[0, np.sum(twin_syn)]
        # To simplify simulations and avoid complex event handling, we smear spike instants in bins of 'dt'
        spikes[0] = np.round(spikes[0],int(np.abs(np.floor(np.log10(solver_options['dt']))))) # Rounds at 10^(-x) (where x=0,1,...)

        # Handle the case where GREs are specified
        # Setting of twin_gtr should be either identical to twin_syn or be such that twin['gtr]==twin['syn']
        # Currently issue a warning
        gres = np.empty((2,0))
        if (rate_gtr!=None) or (stimulus_gtr!= None):
            # Handle the case of a fixed stimulus with respect to a Poisson one
            if stimulus_gtr!='fixed':
                assert (twin_gtr!=None), "twin_gtr not specified"
                assert (stimulus_gtr != None), "stimulus_gtr not specified"
            else:
                assert (pre_gtr!=None), "pre_gtr not specified" # NOTE: will SEGFAULT if pre_gtr contains values >twin[1]
            if stimulus_gtr!='fixed':
                if twin_gtr.size == 2:
                    # One rate or multiple ones in the same interval
                    gres = spg.input_spikes(1, twin_gtr[1], rate_gtr, 0, stimulus=stimulus_gtr, spikes_pre=pre_gtr)
                else:
                    # Multiple rates in different intervals
                    gres = spg.input_spikes(1, twin_gtr, rate_gtr, 0, stimulus=stimulus_gtr, spikes_pre=pre_gtr)
                twin['gtr'] = np.r_[0, np.sum(twin_gtr)]
                if (twin['gtr'][1]!=twin['syn'][1]): print "WARNING: twin['gtr'][1] != twin['syn'][1]"
            else:
                gres = spg.input_spikes(1, [], [], 0, stimulus=stimulus_gtr, spikes_pre=pre_gtr)
                twin['gtr'] = np.asarray([solver_options['t0'], solver_options['tfin']], dtype=float)
            # To simplify simulations and avoid complex event handling, we smear GRE instants in bins of 'dt'
            gres[0] = np.round(gres[0], int(np.abs(np.floor(np.log10(solver_options['dt'])))))  # Rounds at 10^(-x) (where x=0,1,...)
        else:
            gtr = 'ca_dep'
            # In the case of calcium-dependent GTR the time window to perform reconstruction (if requested) will be the
            # whole integration window
            twin['gtr'] = np.asarray([solver_options['t0'],solver_options['tfin']],dtype=float)

        # Number of equations
        if gtr=='ca_dep':
            NEQ = 4
        else:
            NEQ = 2

        # Compute size of spike vectors
        N_spk = int(np.shape(spikes)[1])
        N_gre = int(np.shape(gres)[1])

        # Check on ICs based on NEQ
        if gtr=='ca_dep':
            assert len(pars['ICs'])==NEQ, "ICs with calcium-dependent GTR must be of size 4"
        else:
            if len(pars['ICs'])!=NEQ:         # This check allows to pass ICs from a previous simulation and suppose to continue with a different model
                pars['ICs'] = pars['ICs'][:2] # Takes only the last two elements in the original ICs

        # Check that ICr are of size N_var x N_syn
        # NOTE: you will have to pass the whole vector of ICs to compute solutions from last point
        if (pars['ICr'].size != (N_var+1)*N_syn + 1):
            pars['ICr'] = np.r_[np.tile(pars['ICr'][:N_var+1],(1,N_syn))[0],pars['ICr'][-1]]

        # C-kernel
        support_code = """
                     #include "gliotransmission_models.h"
                     """
        source_files = [os.path.join(os.path.expanduser('~'), base_dir + '/pycustommodules/pycapi_utils.cpp'),
                        os.path.join(os.path.expanduser('~'), base_dir + '/pycustommodules/solvers/solver_options.cpp'),
                        os.path.join(os.path.expanduser('~'), base_dir + '/pycustommodules/solvers/stochastic_solvers.cpp'),
                        os.path.join(os.path.expanduser('~'), base_dir + '/code/gliotransmission_models.cpp')]
        code = """
               // Version
               double version = 0.0;
    
               // Define astrocyte model
               asn tsn(NEQ,N_var,N_syn,N_spk,N_gre);
    
               // Declare output structure
               out_asn out;
    
               // Simulator
               out = tsn.simulate(pars,solver_options,spikes.data(),gres.data());
    
               //Output 
               return_val = out.make_PyDict();
               """
        libs = ['gsl', 'gslcblas', 'm']
        dirs = [os.path.join(os.path.expanduser('~'), base_dir + '/code/'),
                os.path.join(os.path.expanduser('~'), base_dir + '/pycustommodules'),
                os.path.join(os.path.expanduser('~'), base_dir + '/pycustommodules/solvers')]
        vars = ['pars','solver_options',
                'spikes','N_spk',
                'gres','N_gre',
                'NEQ','N_var','N_syn',]
        asn = weave.inline(code,
                           vars,
                           support_code=support_code,
                           sources=source_files,
                           libraries=libs,
                           library_dirs=dirs,
                           include_dirs=dirs,
                           runtime_library_dirs=dirs,
                           type_converters=converters.blitz,
                           compiler='gcc',
                           extra_compile_args=['-std=c++11'],
                           force=recompile)
        # Process synaptic stimuli
        asn['spk']  = spikes[0]     # Spike instants
        asn['is']   = spikes[-1]    # Synapse indexes in the spike train
        if (N_var < 2):
            asn['r'] = pars['u0']*asn['x']
        else:
            asn['r'] = np.multiply(asn['x'],asn['u'])
        if not stimulus_gtr:
            # No specified spikes are fed into the GTR model, so this means that we are using the Ca2+-dependent model
            # of GTR
            # Clean 'x' and provide GRE vector
            i_gre = asn['xa'] > 0
            asn['xa'] = asn['xa'][i_gre]
            asn['gre'] = asn['t'][i_gre]
        else:
            asn['gre'] = gres[0]
        # A vector with all the indexes of Gt. CONVENTION: we use negative indexes for astrocytic release. Only one release
        # site in this implementation, i.e. index -1
        asn['ig'] = -1*np.ones(len(asn['gre']))
        # Add released Gt.
        asn['ra'] = pars['ua']*asn['xa']

        # Append Last point
        if N_var>1:
            u_ = asn['u']
        else:
            u_ = None
        # Check that spk stimulus is not empty. If so then Last point is the same initial point (needed so far to use
        # last_point method with asn['is']=[])
        if (asn['is'].size>0) :
            LCr_syn = last_point(pars,asn['spk'],twin['syn'][1],asn['x'],asn['y'],asn['is'],uval=u_,gtr=False)
        else:
            LCr_syn = pars['ICr'][:-1]
        if (asn['ig'].size>0) :
            LCr_gtr = last_point(pars,asn['gre'],twin['gtr'][1],asn['xa'],asn['xa'],asn['ig'],uval=None,gtr=True)
        else:
            LCr_gtr = pars['ICr'][-1]
        asn['LCr'] = np.r_[LCr_syn,LCr_gtr]
    elif model=='average':
        # Check that rate is scalar
        assert isscalar(rate_syn), "Mean-field synaptic rate must be a non-negative scalar"
        assert isscalar(rate_gtr), "Mean-field gliotransmitter release rate must be a non-negative scalar"
        support_code = """
                     #include "gliotransmission_models.h"
                     """
        source_files = [os.path.join(os.path.expanduser('~'), base_dir + '/pycustommodules/pycapi_utils.cpp'),
                        os.path.join(os.path.expanduser('~'), base_dir + '/pycustommodules/solvers/solver_options.cpp'),
                        os.path.join(os.path.expanduser('~'),
                                     base_dir + '/pycustommodules/solvers/stochastic_solvers.cpp'),
                        os.path.join(os.path.expanduser('~'), base_dir + '/code/gliotransmission_models.cpp')]
        code = """
               // Version
               double version = 0.0;

               // Define astrocyte model
               asn_ave tsn(N_var);

               // Declare output structure
               out_asn out;

               // Simulator
               out = tsn.simulate(rate_syn,rate_gtr,pars,solver_options);

               //Output 
               return_val = out.make_PyDict();
               """
        libs = ['gsl', 'gslcblas', 'm']
        dirs = [os.path.join(os.path.expanduser('~'), base_dir + '/code/'),
                os.path.join(os.path.expanduser('~'), base_dir + '/pycustommodules'),
                os.path.join(os.path.expanduser('~'), base_dir + '/pycustommodules/solvers')]
        vars = ['pars', 'solver_options','N_var',
                'rate_syn', 'rate_gtr']
        asn = weave.inline(code,
                           vars,
                           support_code=support_code,
                           sources=source_files,
                           libraries=libs,
                           library_dirs=dirs,
                           include_dirs=dirs,
                           runtime_library_dirs=dirs,
                           type_converters=converters.blitz,
                           compiler='gcc',
                           extra_compile_args=['-std=c++11'],
                           force=recompile)

    # Post processing of twin
    asn['twin'] = np.asarray([solver_options['t0'], solver_options['tfin']], dtype=float)
    if model=='spiking':
        asn['twin_syn'] = twin['syn']
        asn['twin_gtr'] = twin['gtr']

    # Provide ICs (needed for reconstruction)
    asn['ICs'] = pars['ICs']
    asn['ICr'] = pars['ICr']

    return asn

#-------------------------------------------------------------------------------------------------
# Utils to handle solution
#-------------------------------------------------------------------------------------------------
def last_point(pars,spikes,tfin,xval,yval,idx_,uval=None,gtr=False):
    """
    Estimate the last point of continuous solution.

    Input:
    - pars    :  Dictionary of model parameters (must contain 'u0', 'taud', 'tauf', 'taui')
    - spikes  :  Spike train
    - tfin    :  Last instant of simulation
    - xval    :  x(t_i^-) value array
    - yval    :  y(t_i^+) value array
    - idx     :  Indexes of synapses
    - uval    :  {0.5} | u(t_i^+) or 'u0' value

    Return:
    - LCs     : array of Last points in time for provided input.

    Maurizio De Pitta', Basque Center of Applied Mathematics, June 25, 2018.
    """

    # Make a temporary copy of synaptic indexes (needed to avoid modification of the index internally when treating the
    # the case for gtr=True
    idx = cp.copy(idx_)

    if not gtr:
        Nsyn = int(max(idx))+1
    else:
        Nsyn = 1
        idx += 1 # Make sure that indexes of synapses are 0
    idx = np.asarray(idx, dtype=int)
    x_last,y_last,dt = np.zeros(Nsyn),np.zeros(Nsyn),np.zeros(Nsyn)
    if hasattr(uval,'__len__') or (uval!=None):
        u_last = np.zeros(Nsyn)
    i,j = 0,len(spikes)-1
    # Find the last instant of spike for each synapse
    while (i<Nsyn) and (j>=0):
        if dt[idx[j]]==0:
            dt[idx[j]] = spikes[j]
            if hasattr(uval,'__len__') or (uval!=None):
                u_last[idx[j]] = uval[j]
                x_last[idx[j]] = xval[j]*(1-uval[j])
            else:
                if not gtr:
                    x_last[idx[j]] = xval[j]*(1-pars['u0'])
                else:
                    x_last[idx[j]] = xval[j]*(1-pars['ua'])
            if not gtr: y_last[idx[j]] = yval[j]
            i += 1 # as soon as i==Nsyn we found the last spike for all synapses
        j -= 1
    # Retrieve the interval from last spike to the end of simulation
    dt = tfin-dt

    # Compute last value
    if not gtr:
        x_last = xsol(x_last,dt,pars['taud'])
    else:
        x_last = xsol(x_last, dt, pars['taua'])

    if not gtr:
        y_last = usol(y_last, dt, pars['taui'])
        if hasattr(uval,'__len__') or (uval!=None):
            u_last = usol(u_last, dt, pars['tauf'])

    # Recast last values in the order of ICs
    if not gtr:
        if hasattr(uval,'__len__') or (uval!=None):
            LCs = (np.vstack((x_last, y_last, u_last))).flatten(order='F')
        else:
            LCs = (np.vstack((x_last,y_last))).flatten(order='F')
    else:
        LCs = (np.asarray(x_last)).flatten(order='F')

    return LCs

def reconstruct_solution(spikes,sol,uval,twin,ics,tau,variable,**kwargs):
    """
    Reconstruct solution from spiking solutions.

    Input:
    - spikes : AP/GRE train
    - sol  : Solution points (aka xn, un, yn)
    - uval : value for un (if needed)
    - twin : time interval within which to reconstruct the solution
    - ics  : initial conditions for given sol vector
    - tau  : decay time for sol
    - variable : 'x' | 'u' | 'y'
    **kwargs:
      - dt  : time bin for reconstructed solution

    Return:
    - solution : a list of 2xT/dt arrays where the first row is a time vector, and the second row is the reconstructed
                 solution for the specified variable.

    Maurizio De Pitta', Basque Center of Applied Mathematics, June 13, 2018.
    """

    # Model parameters
    pars = {'dt' : 1e-3}
    pars = gu.varargin(pars,**kwargs)

    # Generate time vector
    time = np.arange(twin[0],twin[-1],pars['dt'])
    time = np.sort(np.r_[time,spikes])
    # Generate spike vector
    tspk = np.copy(time)
    for i in range(1,len(spikes)):
        tspk[np.where(np.logical_and(time>=spikes[i-1],time<spikes[i]))[0]] = spikes[i-1]
    tspk[np.where(time >= spikes[len(spikes)-1])[0]] = spikes[len(spikes)-1]
    tspk[np.where(time < spikes[0])[0]] = 0
    # Generate general solution vector
    vsol = np.ones(time.size)
    if (variable=='x') and isscalar(uval):
        uval = uval * np.ones(sol.size)
    if variable=='x':
        for i in range(1, len(spikes)):
            # x must be given at x(t_i^+) according to xsol
            vsol[np.where(np.logical_and(time >= spikes[i - 1], time < spikes[i]))[0]] = sol[i-1]*(1-uval[i-1])
        vsol[np.where(time >= spikes[len(spikes) - 1])[0]] = sol[len(spikes) - 1]*(1-uval[len(spikes)-1])
    else:
        for i in range(1, len(spikes)):
            vsol[np.where(np.logical_and(time >= spikes[i - 1], time < spikes[i]))[0]] = sol[i-1]
        vsol[np.where(time >= spikes[len(spikes) - 1])[0]] = sol[len(spikes) - 1]
    vsol[np.where(time < spikes[0])[0]] = ics
    # Compute effective solution
    solution = np.zeros((2, time.size))
    solution[0] = time

    if variable=='x':
        # Assumes that the first ICs is x(0)
        solution[1] = xsol(vsol,time-tspk,tau)
    else:
        solution[1] = usol(vsol,time-tspk,tau)

    return solution

#-------------------------------------------------------------------------------------------------
# Gliotransmitter-regulated synapse (GTRS) class
#-------------------------------------------------------------------------------------------------
"""
GTRS is the class that allows simulating all the above models. It is essentially a wrapper of the above routines with 
some additional features such as reconstruction of solutions by default for all synapses and the astrocyte. Class
definition is as following:

    gtrs model_object(model=...)
    
where model is 
    
    -'lra'     : Li-Rinzel astrocyte;
    -'exo_spk' : TM spiking model;
    -'exo_ave' : TM mean field model;
    -'lra_exo' : Calcium-dependent release model;
    -'asn_spk' : Tripartite synapse spiking model;
    -'asn_ave' : Tripartite synapse mean field model.
    
Methods:
    self.stimulation : wrapper of stimulation_pars;
    self.simulate    : simulate the model
    self.reconstruct : reconstruct time solution from spiking solution.

Maurizio De Pitta', Basque Center of Applied Mathematics, June 20, 2018.
"""
class gtrs(object):
    def __init__(self,model,**kwargs):
        self._model = model
        if model=='lra':
            self.pars = lra_parameters(**kwargs)
        elif (model!='lra_exo') and 'exo' in model:
            self.pars = exocytosis_parameters(**kwargs)
        elif model=='lra_exo':
            self.pars = gtrelease_parameters(**kwargs)
        elif 'asn' in model:
            if model=='asn_ave':
                self.pars = asn_parameters('ave',**kwargs)
            else:
                self.pars = asn_parameters(**kwargs)

        # Assign initial conditions
        self.ICs = self.pars['ICs']
        if 'asn' in model:
            self.ICr = self.pars['ICr']

    def stimulation(self,**kwargs):
        self.stim = stimulus_params(**kwargs)

    def simulate(self,
                 algparams=None,
                 reconstruct=False,
                 recompile=0,
                  **kwargs):
        # Solver options handling
        self._algparams = algparams # These are the solver_options
        if (not self._algparams):
            self._algparams = su.solver_opts(method='gsl_rk2imp',**kwargs)
        # Make sure to pass last ICs (necessary if we start simulation from last step)
        self.pars['ICs'] = self.ICs
        if 'asn' in self._model:
            self.pars['ICr'] = self.ICr
        # Handling of different cases
        if self._model=='lra':
            assert self._algparams!=None, "Solver options must be specified"
            self.sol = lr_model(self.pars,self._algparams,recompile=recompile)
        elif self._model=='exo_spk':
            try:
                self.sol = exocytosis_model(self.stim['rate'],self.stim['twin'],
                                            self.stim['Nsyn'],self.pars,self.stim['Nvar'],
                                            model='spiking',
                                            solver_options=None,
                                            stimulus=self.stim['stimulus'],spikes_pre=self.stim['spikes_pre'],
                                            recompile=recompile)
                if reconstruct:
                    # Provide a default reconstruction of all spiking solutions
                    self.ICs = self.sol['ICs']
                    self.reconstruct(synapse_index=np.arange(self.stim['Nsyn']),var='all')
            except AttributeError:
                print "Stimulation was not set"
        elif self._model=='exo_ave':
            assert self._algparams!=None, "Solver options must be specified"
            try:
                self.sol = exocytosis_model(self.stim['rate'],self.stim['twin'],
                                            None,self.pars,self.stim['Nvar'],
                                            model='average',
                                            solver_options=algparams,
                                            recompile=recompile)
            except AttributeError:
                print "Stimulation was not set"
        elif self._model=='lra_exo':
            assert self._algparams!=None, "Solver options must be specified"
            self.sol = ca_gtrel_model(self.pars, self._algparams, recompile=recompile)
            if reconstruct:
                self.reconstruct(synapse_index=[-1], var='x')  # Negative indexes are used for astrocytic release
        elif self._model=='asn_spk':
            assert self._algparams!=None, "Solver options must be specified"
            try:
                self.sol = asn_model(self.pars,algparams,
                                     self.stim['Nsyn'],self.stim['Nvar'],
                                     self.stim['rate'],self.stim['twin'],
                                     model='spiking',gtr=None,
                                     stimulus_syn=self.stim['stimulus'],pre_syn=self.stim['spikes_pre'],
                                     rate_gtr=self.stim['rate_gtr'],twin_gtr=self.stim['twin_gtr'],
                                     stimulus_gtr=self.stim['stimulus_gtr'],pre_gtr=self.stim['pre_gtr'],
                                     recompile=recompile)
                if reconstruct:
                    # Update ICr with the one produced inside the simulator
                    self.ICr = self.sol['ICr']
                    # Provide a default reconstruction of all spiking solutions
                    if self.stim['Nvar'] > 1:
                        self.reconstruct(synapse_index=np.arange(self.stim['Nsyn']),var='all') # Reconstruct synaptic inputs
                    else:
                        print "WARNING: 'reconstruct' does not support 1-var synaptic model in the presence of GTR. Producing only xa."
                    self.reconstruct(synapse_index=[-1], var='x') # Reconstruct GTR
            except AttributeError:
                print "Stimulation was not set"
        elif self._model=='asn_ave':
            try:
                self.sol = asn_model(self.pars,algparams,
                                     self.stim['Nsyn'],self.stim['Nvar'],
                                     self.stim['rate'],self.stim['twin'],
                                     model='average',gtr=None,
                                     rate_gtr=self.stim['rate_gtr'],
                                     recompile=recompile)
            except AttributeError:
                print "Stimulation was not set"

    def reconstruct(self,synapse_index=0,var='all',**kwargs):
        """
        Reconstruct specific solution from spiking model.

        Input:
        synapse_index : Scalar or array of indexes of synapses
        var           : which variable to reconstruct

        Return:
        sol           : Modified solution dictionary with 'ts','xs','us', 'ys' additional entries
        """
        try:
            # Fix dt for reconstruction
            DT = 1e-3

            # First check that synapse_index is not out of range
            synapse_index = np.asarray(np.atleast_1d(synapse_index),dtype=int)
            # assert (synapse_index.min() >=0)&(synapse_index.max()<=self.stim['Nsyn']), "Synapse index out of range [0,N_syn]"

            # Basic handling of ICs (since ICs for the discrete compartments are handled differently in 'asn_spk'
            if self._model=='asn_spk':
                ICs = self.ICr
            else:
                ICs = self.ICs

            for i,si in enumerate(synapse_index):
                # Generate a 'tag' to append to dict keys to distinguish between synaptic ('s') and 'astrocytic ('a') release
                if si>=0:
                    tag = 's'
                    twin = self.sol['twin']
                else:
                    tag = 'g'
                    twin = self.sol['twin_gtr']
                # Retrieve spike indexes
                index = self.sol['i'+tag]==si  # Spike indexes per given synapse
                # Reconstruction is performed at fixed dt for any solution
                if var in ['all','x']:
                    if si>=0:
                        if self.stim['Nvar']>1:
                            u_ = self.sol['u'][index]
                        else:
                            u_ = self.pars['u0'] # Will not work in the case of GTR
                        x = reconstruct_solution(self.sol['spk'][index], self.sol['x'][index], u_,
                                                 twin, ICs[si*(self.stim['Nvar']+1)],
                                                 self.pars['taud'],variable='x',dt=DT)
                    else:
                        x = reconstruct_solution(self.sol['gre'][index], self.sol['xa'][index], self.pars['ua'],
                                                 twin, ICs[-1],
                                                 self.pars['taua'], variable='x', dt=DT)
                    if (i==0):
                        if 't'+tag not in self.sol: self.sol['t'+tag] = [x[0]]
                        self.sol['x'+tag] = [x[-1]]
                    else:
                        if len(self.sol['t'+tag]) <= i: self.sol['t'+tag].append(x[0])
                        self.sol['x'+tag].append(x[-1])
                if (hasattr(self,'stim'))and(self.stim['Nvar']>1)and(var in ['all','u']):
                    u = reconstruct_solution(self.sol['spk'][index], self.sol['u'][index], [],
                                             twin, ICs[si*(self.stim['Nvar']+1)+2],
                                             self.pars['tauf'], variable='u',dt=DT)
                    if (i==0):
                        if 't'+tag not in self.sol: self.sol['t'+tag] = [u[0]]
                        self.sol['u'+tag] = [u[-1]]
                    else:
                        if len(self.sol['t'+tag]) <= i: self.sol['t'+tag].append(u[0])
                        self.sol['u'+tag].append(u[-1])
                if var in ['all','y']:
                    y = reconstruct_solution(self.sol['spk'][index],self.sol['y'][index], [],
                                             twin, ICs[si*(self.stim['Nvar']+1)+1],
                                             self.pars['taui'], variable='y',dt=DT)
                    if (i == 0):
                        if 't'+tag not in self.sol: self.sol['t'+tag] = [y[0]]
                        self.sol['y'+tag] = [y[-1]]
                    else:
                        if len(self.sol['t'+tag]) <= i: self.sol['t'+tag].append(y[0])
                        self.sol['y'+tag].append(y[-1])
        except AttributeError:
            print "Spiking solution does not exist"

    def __del__(self):
        del self.pars
        del self.ICs
        if hasattr(self,'stim'):
            del self.stim
        if hasattr(self,'ICr'):
            del self.ICr
        if hasattr(self,'sol'):
            del self.sol

if __name__ == "__main__":
    # -------------------------------------------------------------------------------------------------
    # Testing Astrocyte class
    # -------------------------------------------------------------------------------------------------
    ####################################################################################################################
    # #LRA model
    ####################################################################################################################
    # options = su.solver_opts(t0=0.0, tfin=10.0, dt=1e-2, atol=1e-8, rtol=1e-8, method="gsl_rk2imp")
    # astro = gtrs(model='lra',ip3=0.4,noise=True)
    # astro.simulate(algparams=options,recompile=1)
    # plt.plot(astro.sol['t'], astro.sol['ca'], 'y-', astro.sol['t'], astro.sol['h'], 'm-')

    # ####################################################################################################################
    # # # Exocytosis (TM model)
    # ####################################################################################################################
    # exo = gtrs(model='exo_spk',u0=0.2, ICs=[1.0,0.,0.0])
    # # exo.stimulation(rate=5,twin=[0.,2.],Nsyn=3,Nvar=1)
    # # exo.stimulation(rate=5,twin=[0.,2.],Nsyn=3,Nvar=1)
    # exo.stimulation(rate=5, twin=[0., 2.], Nsyn=3, Nvar=2, stimulus='poisson')
    # exo.simulate(reconstruct=True,recompile=0)
    # # print exo.sol['LCs']
    # # # 1-VAR
    # # index = exo.sol['is'] == 0
    # # plt.plot(exo.sol['ts'][0], exo.sol['xs'][0], 'y-', exo.sol['spk'][index], exo.sol['x'][index], 'yo')
    # # index = exo.sol['is']==1
    # # plt.plot(exo.sol['ts'][1], exo.sol['xs'][1], 'y-', exo.sol['spk'][index], exo.sol['x'][index], 'yo')
    # # index = exo.sol['is']==2
    # # plt.plot(exo.sol['ts'][2], exo.sol['xs'][2], 'y-', exo.sol['spk'][index], exo.sol['x'][index], 'yo')
    #
    # # 2-VAR
    # index = exo.sol['is']==0
    # plt.plot(exo.sol['ts'][0], exo.sol['xs'][0], 'y-', exo.sol['spk'][index], exo.sol['x'][index], 'yo',
    #          exo.sol['ts'][0], exo.sol['us'][0], 'm-', exo.sol['spk'][index], exo.sol['u'][index], 'mo')
    # index = exo.sol['is']==1
    # plt.plot(exo.sol['ts'][1], exo.sol['xs'][1], 'y-', exo.sol['spk'][index], exo.sol['x'][index], 'yo',
    #          exo.sol['ts'][1], exo.sol['us'][1], 'm--' ,exo.sol['spk'][index], exo.sol['u'][index], 'mo')
    # index = exo.sol['is']==2
    # plt.plot(exo.sol['ts'][2], exo.sol['xs'][2], 'y-', exo.sol['spk'][index], exo.sol['x'][index], 'yo',
    #          exo.sol['ts'][2], exo.sol['us'][2], 'm--' ,exo.sol['spk'][index], exo.sol['u'][index], 'mo')

    ####################################################################################################################
    # # Exocytosis (MF model)
    ####################################################################################################################
    # exo = gtrs(model='exo_ave')
    # exo.stimulation(rate=3,twin=[0.,1.0],Nvar=1)
    # options = su.solver_opts(t0=exo.stim['twin'][0], tfin=exo.stim['twin'][1], dt=1e-4, atol=1e-8, rtol=1e-8, method="gsl_rk2imp")
    # exo.simulate(algparams=options,recompile=0)
    # if exo.stim['Nvar']>1:
    #     plt.plot(exo.sol['t'], exo.sol['x'], 'y', exo.sol['t'], exo.sol['u'],'m')
    # else:
    #     plt.plot(exo.sol['t'], exo.sol['x'], 'y')

    ####################################################################################################################
    # # #Calcium-dependent gliotransmitter exocytosis model
    ####################################################################################################################
    # options = su.solver_opts(t0=0.0, tfin=60.0, dt=1e-3, atol=1e-8, rtol=1e-8, method="gsl_rk2imp")
    # astro = gtrs(model='lra_exo',ip3=0.4,noise=False,ICs=[0.05,0.9,1.0])
    # astro.simulate(algparams=options,reconstruct=True,recompile=0)
    # # astro.simulate(algparams=options,reconstruct=True,recompile=1)
    # # You must take the first element regardless
    # plt.plot(astro.sol['tg'][0], astro.sol['xg'][0], 'y-',astro.sol['gre'],astro.sol['xa'], 'yo',
    #          astro.sol['t'], astro.sol['ca'],'k-',astro.sol['t'], astro.sol['h'],'m-') # You must take the first element regardless

    ####################################################################################################################
    # Tripartite synapse model (with calcium-dependent GTR)
    ####################################################################################################################
    options = su.solver_opts(t0=0.0, tfin=20.0, dt=1e-3, atol=1e-8, rtol=1e-8, method="gsl_rk2imp")
    # tsn = gtrs(model='asn_spk',ip3=0.4,alpha=0.0,noise=False)
    # # # Case of NEQ>2 (i.e. calcium-dependent GTR)
    # # tsn.stimulation(rate=5, twin=[0., 2.], Nsyn=1, Nvar=2)
    # # tsn.simulate(algparams=options,reconstruct=True,recompile=0)
    # # # Case of NEQ<3 (GREs given)
    # # tsn.stimulation(rate=5, twin=[0., 10.], Nsyn=1, Nvar=2,
    # #                 rate_gtr=1, twin_gtr=[0.,10.],stimulus_gtr='poisson')
    # # tsn.simulate(algparams=options,reconstruct=False,recompile=1)
    #
    # # General Output
    # plt.plot(tsn.sol['t'], tsn.sol['ga'], 'y-', tsn.sol['t'], tsn.sol['gammas'], 'm-')
    # if not tsn.stim['rate_gtr']: plt.plot(tsn.sol['t'], tsn.sol['ca'], 'k-', tsn.sol['t'],tsn.sol['h'],'b-')

    # ## Checking on the effect of GTR on u0 using a periodic/poisson stimulus
    # tsn = gtrs(model='asn_spk',ip3=0.4,alpha=1.0,u0=0.1,taud=0.5,tauf=1.,Og=1.5,
    #            noise=False)
    # tsn.stimulation(rate=5, twin=[0., 10.], stimulus='poisson',Nsyn=1, Nvar=2,
    #                 rate_gtr=1, twin_gtr=[0.,10.],stimulus_gtr='fixed',pre_gtr=[4.,6.])
    # # tsn.stimulation(rate=1, twin=[0., 10.], stimulus='periodic', Nsyn=1, Nvar=2)
    # tsn.simulate(algparams=options,reconstruct=True,recompile=0)
    # # Comparison between gammas and y
    # plt.plot(tsn.sol['t'], tsn.sol['gammas'], 'm-', tsn.sol['ts'][0],tsn.sol['ys'][0],'k-')

    # ## Checking on the effect of GTR on u0 using a periodic/poisson stimulus
    # tsn = gtrs(model='asn_spk',noise=False,ua=0.)
    # tsn.stimulation(rate=1, twin=[0., 20.], stimulus='periodic',Nsyn=1, Nvar=2, rate_gtr=0.1, twin_gtr=[0.,20.], stimulus_gtr='periodic')
    # # tsn.stimulation(rate=1, twin=[0., 10.], stimulus='periodic', Nsyn=1, Nvar=2)
    # tsn.simulate(algparams=options,reconstruct=False,recompile=0)
    # # Comparison between gammas and y
    # plt.plot(tsn.sol['t'], tsn.sol['gammas'], 'm-', tsn.sol['t'],tsn.sol['ga'],'k-')

    # Nsyn = 10
    # options = su.solver_opts(t0=0.0, tfin=20.0, dt=1e-3, atol=1e-8, rtol=1e-8, method="gsl_rk2imp")
    # tsn = gtrs(model='asn_spk',Ou=1e2, Ku=0.1, taue=1./30.,Og=0.1,taug=30.,alpha=1.0,u0=0.05,taud=0.3,tauf=1.0)
    # tsn.stimulation(rate=1, twin=[0., 20.], stimulus='periodic',Nsyn=Nsyn, Nvar=2,
    #                 stimulus_gtr='fixed',pre_gtr=[7.0])
    # tsn.simulate(algparams=options,reconstruct=True,recompile=0)
    # # plt.plot(tsn.sol['t'], tsn.sol['gammas'], 'm-', tsn.sol['ts'][0], tsn.sol['ys'][0], 'k-')
    # # for i in xrange(Nsyn):
    # #     plt.plot(tsn.sol['ts'][i], tsn.sol['us'][i])
    # plt.plot(tsn.sol['spk'], tsn.sol['x'], 'ko')
    # plt.plot(tsn.sol['spk'], tsn.sol['u'], 'ro')

    # # Check on synaptic input reconstruction (single variable)
    # index = tsn.sol['is']==0
    # print tsn.sol['x'][index]
    # plt.plot(tsn.sol['ts'][0], tsn.sol['xs'][0], 'k-',tsn.sol['spk'][index],tsn.sol['x'][index],'ko')
    # # plt.plot(tsn.sol['ts'][0], tsn.sol['ys'][0], 'm-',tsn.sol['spk'][index],tsn.sol['y'][index],'yo')
    # plt.plot(tsn.sol['ts'][0], tsn.sol['us'][0], 'k--',tsn.sol['spk'][index],tsn.sol['u'][index],'ks')
    # index = tsn.sol['is']==1
    # print tsn.sol['x'][index]
    # plt.plot(tsn.sol['ts'][1], tsn.sol['xs'][1], 'r-',tsn.sol['spk'][index],tsn.sol['x'][index],'ro')
    # plt.plot(tsn.sol['ts'][1], tsn.sol['us'][1], 'r--', tsn.sol['spk'][index], tsn.sol['u'][index],'rs')

    # # Check on reconstruction of xa (even in the case of a single variable such as xa, the reconstructed variable must be indexed)
    # plt.plot(tsn.sol['tg'][0], tsn.sol['xg'][0],'y-',tsn.sol['gre'], tsn.sol['xa'], 'yo')
    # if not tsn.stim['rate_gtr']: plt.plot(tsn.sol['t'], tsn.sol['ca'], 'k-', tsn.sol['t'],tsn.sol['h'],'b-')

    ####################################################################################################################
    # MF Tripartite synapse model
    ####################################################################################################################
    options = su.solver_opts(t0=0.0, tfin=30.0, dt=1e-3, atol=1e-8, rtol=1e-8, method="gsl_bsimp")
    tsn = gtrs(model='asn_ave',u0=0.6,alpha=0.0)
    tsn.stimulation(rate=5, twin=[0., 30.], Nvar=2,rate_gtr=2)
    tsn.simulate(algparams=options,recompile=1)

    # plt.plot(tsn.sol['t'],tsn.sol['xa'],'y-',tsn.sol['t'],tsn.sol['gammas'],'k-',tsn.sol['t'],tsn.sol['x'],'r-')
    plt.plot(tsn.sol['t'], tsn.sol['ga'], 'k-')
    # plt.plot(tsn.sol['t'],tsn.sol['xa'],'y-',tsn.sol['t'],tsn.sol['gammas'],'k-',tsn.sol['t'],tsn.sol['x'],'r-')
    # if tsn.stim['Nvar']>1:
    #     plt.plot(tsn.sol['t'],tsn.sol['u'],'b-')

    # # Show plot
    plt.show()