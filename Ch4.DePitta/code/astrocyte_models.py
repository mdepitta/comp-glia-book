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
import collections # OrderedDict module
import weave
from weave import converters
import sys, os
base_dir = '/Ch4.DePitta'
sys.path.append(os.path.join(os.path.expanduser('~'),base_dir))
import pycustommodules.general_utils as gu
import pycustommodules.solvers.solver_utils as su
import pycustommodules.save_utils as svu

import matplotlib.pylab as plt

###############################################################################
# General functions
###############################################################################
normalize = lambda x : (x - x.min())/(x.max()-x.min())
unnormalize = lambda xn, mm : mm[0] + xn*(np.diff(mm)) # mm must be in [min, max]

def Hill(x,k,n):
    return np.asarray(x**n/(x**n+k**n), dtype=float)

###############################################################################
# Parameters definition
###############################################################################
def lra_parameters(**kwargs):
    """
    Parameters for the Li-Rinzel astrocyte (Li and Rinzel, JTB 1994).

    Maurizio De Pitta', INRIA Rhone-Alpes, November 3rd, 2017.
    """
    pars = {'d1' : 0.1,
            'd2' : 2.1,
            'd3' : 0.9967,
            'd5' : 0.2,
            'a2' : 0.4,
            'c1' : 0.4,
            'c0' : 4,
            'rc' : 7,
            'rl' : 0.05,
            'ver': 0.9,
            'Ker': 0.1,
            'ip3': 0.1,
            'ICs': [0.05,0.99]
            }

    ## User-defined parameters
    pars = gu.varargin(pars,**kwargs)
    # Parameters must be floats
    for k,item in pars.iteritems():
        if isscalar(item):
            pars[k] = float(item)
        else:
            pars[k] = array(item,dtype=float)
    return pars

def chi_parameters(**kwargs):
    """
    Parameters for the plain ChI model (DePitta' et al., JOBP 2009).
    
    Maurizio De Pitta', The University of Chicago, February 28th, 2015.
    """
    pars_lra = lra_parameters()
    pars_lra.pop('ip3',None)
    pars = {'vbias': 0,            
            'vbeta': 3,
            'vdelta': 0.5,            
            'kappad': 1,
            'Kdelta': 0.5,
            'v3k': 2,
            'Kd': 0.5,
            'K3': 1,
            'r5p': 0,
            'ICs': [0.05, 0.05, 0.99]
            }

    # Merge the two parameter dictionaries
    pars = gu.merge_dicts(pars_lra,pars)
    ## User-defined parameters
    pars = gu.varargin(pars,**kwargs)
    # Parameters must be floats
    for k,item in pars.iteritems():
        if isscalar(item):
            pars[k] = float(item)
        else:
            pars[k] = array(item,dtype=float)
    return pars

def gchi_parameters(**kwargs):
    """
    Parameters for the G-ChI model (DePitta' et al., JOBP 2009).
    
    Maurizio De Pitta', The University of Chicago, February 28th, 2015.
    """
    # TODO: should include here vbeta and leave only vbias in chi_parameters but this requires
    # implementing a bias in the astrocyte model -- right now we use vbeta for constant IP3 production
    # in the ChI model
    pars_chi = chi_parameters()
    pars = {'yrel': 0.02,
            'Op': 0.3,
            'OmegaP': 1.8,
            'zeta': 0,            
            'Kkc': 0.5,
            'ICs': [0.01, 0.05, 0.05, 0.99],
            # Bias / Exogenous stimulation (default: no bias)
            'T'   : 0.,   # Period
            'pw'  : 0.,   # pulse-width
            'yb'  : 0.    # amplitude
            }
    # Merge the two parameter dictionaries
    pars = gu.merge_dicts(pars_chi,pars)
    ## User-defined parameters
    pars = gu.varargin(pars,**kwargs)
    # Parameters must be floats
    for k,item in pars.iteritems():
        if isscalar(item):
            pars[k] = float(item)
        else:
            pars[k] = array(item,dtype=float)
    return pars

def egchi_parameters(**kwargs):
    """
    Parameters for the extended G-ChI model that includes PKC, DAG and AA

    Maurizio De Pitta', INRIA Rhone-Alpes, November 27th, 2017.
    """

    pars_gchi = gchi_parameters(Kkc=0.6)
    del pars_gchi['zeta']
    del pars_gchi['T']
    del pars_gchi['pw']
    del pars_gchi['yb']
    pars = {'vkd'    : 0.5,
            'vk'     : 1.0,
            'OmegaKD': 2.5,
            'vd'     : 1.5,
            'Kdc'    : 0.3,
            'Kdd'    : 0.05,
            'OmegaD' : 0.1,
            'ICs'    : [0.01, 0.05, 0.05, 0.99, 0.05, 0.0] # [ago, I, C, h, D, P]
            }

    # Merge the two parameter dictionaries
    pars = gu.merge_dicts(pars_gchi,pars)
    ## User-defined parameters
    pars = gu.varargin(pars,**kwargs)
    # Parameters must be floats
    for k,item in pars.iteritems():
        if isscalar(item):
            pars[k] = float(item)
        else:
            pars[k] = array(item,dtype=float)
    return pars

def model_bounds(model='lra',**kwargs):
    # Generate empty Ordered Dictionary for boundaries. The order of keys MUST mirror the order of variables 'x' in the
    # cost function
    bounds_dict = collections.OrderedDict()
    if model=='lra':
        bounds_dict['d1'] = [0.1, 10.0, 'log']
        bounds_dict['d2'] = [0.1, 10.0, 'log']
        bounds_dict['d3'] = [0.1, 10.0, 'log']
        bounds_dict['d5'] = [0.1, 10.0, 'log']
        bounds_dict['a2'] = [0.1, 5.0, 'log']
    elif model=='lra2':
        bounds_dict['d1'] = [0.1, 0.5, 'lin']
        bounds_dict['d2'] = [1.0, 4.5, 'lin']
        bounds_dict['d5'] = [0.05,0.5, 'log']
        bounds_dict['a2'] = [0.1, 0.5, 'lin']
    elif model=='lra_fit':
        bounds_dict['rc']  = [2.0, 20.0, 'lin']
        bounds_dict['ver'] = [2.0, 20.0, 'lin']
        bounds_dict['ip3'] = [0.05, 0.5, 'lin']
        bounds_dict['C0']  = [0.05, 0.5, 'lin']
        bounds_dict['h0']  = [0.0,  1.0, 'lin']
        bounds_dict['c0']  = [4.0, 12.0, 'lin']
    elif model=='chi':
        bounds_dict['vbeta']  = [0.001, 5.0, 'log']
        bounds_dict['vdelta'] = [0.001, 0.5, 'log']
        bounds_dict['v3k']    = [0.1,  5.0, 'log']
        bounds_dict['r5p']    = [0.1,  1.0, 'lin']
        bounds_dict['C0']     = [0.05, 5.0, 'log']
        bounds_dict['h0']     = [0.0,  1.0, 'lin']
        bounds_dict['I0']     = [0.05, 5.0, 'log']

    # Customer defined bounds
    bounds_dict = gu.varargin(bounds_dict,**kwargs)

    # Treatment of lin/log scaling parameters
    for b in bounds_dict.values():
        if b[-1]=='log': b[:2] = np.log10(b[:2])

    return bounds_dict

###############################################################################
# Simulators
###############################################################################
def lr_model(pars, solver_options):
    """
    Simulator for the plain Li-Rinzel model (Li and Rinzel, JTB 1994).

    Maurizio De Pitta', INRIA Rhone-Alpes, February 28th, 2015.
    """
    support_code = """
                 #include "astrocyte_models.h"
                 """
    source_files = [os.path.join(os.path.expanduser('~'), base_dir+'/pycustommodules/pycapi_utils.cpp'),
                    os.path.join(os.path.expanduser('~'), base_dir+'/pycustommodules/solvers/solver_options.cpp'),
                    os.path.join(os.path.expanduser('~'), base_dir+'/code/astrocyte_models.cpp')]
    code = """
           // Version
           double version = 1.0; 
           // Declare output structure
           out_lra out;

           // Simulator
           out = lra_simulator_gsl(pars,solver_options);

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
                        force=0)
    # Post-stimulus processing
    return olr

def chi_model(pars,solver_options):
    """
    Simulator for the plain CHI model (DePitta' et al., JOBP 2009).
    
    Maurizio De Pitta', The University of Chicago, February 28th, 2015.
    """
    support_code="""
                 #include "astrocyte_models.h"
                 """
    source_files = [os.path.join(os.path.expanduser('~'), base_dir+'/pycustommodules/pycapi_utils.cpp'),
                    os.path.join(os.path.expanduser('~'), base_dir+'/pycustommodules/solvers/solver_options.cpp'),
                    os.path.join(os.path.expanduser('~'), base_dir+'/code/astrocyte_models.cpp')]
    code = """
           // Version
           double version = 1.0; 
           // Declare output structure
           out_chi out;

           // Simulator
           out = chi_simulator_gsl(pars,solver_options);

           //Output 
           return_val = out.make_PyDict();
           """
    libs=['gsl','gslcblas','m']
    dirs=[os.path.join(os.path.expanduser('~'),base_dir+'/code/'),
          os.path.join(os.path.expanduser('~'),base_dir+'/pycustommodules'),
          os.path.join(os.path.expanduser('~'),base_dir+'/pycustommodules/solvers')]
    vars = ['pars','solver_options']      
    ochi = weave.inline(code,
                      vars,
                      support_code=support_code,
                      sources = source_files,
                      libraries=libs,
                      library_dirs=dirs,
                      include_dirs=dirs,
                      runtime_library_dirs=dirs,
                      type_converters=converters.blitz,
                      compiler='gcc',
                      extra_compile_args=['-std=c++11'],
                      force=0)
    # Post-stimulus processing
    return ochi
    
def gchi_model(pars,solver_options):    
    """
    Simulator for the G-CHI model (DePitta' et al., JOBP 2009).
    
    Maurizio De Pitta', The University of Chicago, March 1st, 2015.
    """
    support_code="""
                 #include "astrocyte_models.h"
                 """
    source_files = [os.path.join(os.path.expanduser('~'),base_dir+'/pycustommodules/pycapi_utils.cpp'),
                    os.path.join(os.path.expanduser('~'),base_dir+'/pycustommodules/solvers/solver_options.cpp'),
                    os.path.join(os.path.expanduser('~'),base_dir+'/code/astrocyte_models.cpp')]
    code = """
           // Version
           double version = 1.0;         
           // Declare output structure
           out_gchi out;

           // Simulator
           out = gchi_simulator_gsl(pars,solver_options);

           //Output
           return_val = out.make_PyDict();
           """
    libs=['gsl','gslcblas','m']
    dirs=[os.path.join(os.path.expanduser('~'),base_dir+'/code/'),
          os.path.join(os.path.expanduser('~'),base_dir+'/pycustommodules'),
          os.path.join(os.path.expanduser('~'),base_dir+'/pycustommodules/solvers')]
    vars = ['pars','solver_options']      
    ogchi = weave.inline(code,
                      vars,
                      support_code=support_code,
                      sources = source_files,
                      libraries=libs,
                      library_dirs=dirs,
                      include_dirs=dirs,
                      runtime_library_dirs=dirs,
                      type_converters=converters.blitz,
                      compiler='gcc',
                      extra_compile_args=['-std=c++11'],
                      force=0)
    # Post-stimulus processing
    return ogchi

def egchi_model(pars, solver_options):
    """
    Simulator for the Extended G-CHI model (DePitta' et al., 2018).

    Maurizio De Pitta', INRIA Rhone-Alpes, November 28, 2017.
    """
    support_code = """
                 #include "astrocyte_models.h"
                 """
    source_files = [os.path.join(os.path.expanduser('~'),base_dir+'/pycustommodules/pycapi_utils.cpp'),
                    os.path.join(os.path.expanduser('~'),base_dir+'/pycustommodules/solvers/solver_options.cpp'),
                    os.path.join(os.path.expanduser('~'),base_dir+'/code/astrocyte_models.cpp')]
    code = """
           // Version
           double version = 1.5;         
           // Declare output structure
           out_egchi out;

           // Simulator
           out = egchi_simulator_gsl(pars,solver_options);

           //Output
           return_val = out.make_PyDict();
           """
    libs = ['gsl', 'gslcblas', 'm']
    dirs = [os.path.join(os.path.expanduser('~'),base_dir+'/code/'),
            os.path.join(os.path.expanduser('~'),base_dir+'/pycustommodules'),
            os.path.join(os.path.expanduser('~'),base_dir+'/pycustommodules/solvers')]
    vars = ['pars', 'solver_options']
    oegchi = weave.inline(code,
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
                         force=0)
    # Post-stimulus processing
    return oegchi

#-------------------------------------------------------------------------------------------------
# Astrocyte class
#-------------------------------------------------------------------------------------------------
class Astrocyte(object):
    def __init__(self,model,**kwargs):
        self._model = model
        if model=='lra':
            self.pars = lra_parameters(**kwargs)
        elif model in ['lra2', 'lra_fit']:
            self.pars = lra_parameters(**kwargs)
        elif model=='chi':
            self.pars = chi_parameters(**kwargs)
        elif model=='gchi':
            self.pars = gchi_parameters(**kwargs)
        elif model=='pkc':
            self.pars = egchi_parameters(**kwargs)

        # Assign initial conditions
        self.ICs = self.pars['ICs']

    def popen(self,ca,ip3):
        # Return open probability
        m_inf = Hill(ip3,self.pars['d1'],1)
        n_inf = Hill(ca,self.pars['d5'],1)
        Q2 = self.pars['d2']*(ip3 + self.pars['d1'])/(ip3 + self.pars['d3'])
        # Q2 = self.pars['d2'] * (ip3 + self.pars['d1']) / (ip3 + self.pars['d1'])
        h_inf = Hill(Q2,ca,1)
        self.po  = (m_inf*n_inf*h_inf)**3

    def popen_2(self,ca,ip3):
        # Return open probability (assuming d1==d3)
        m_inf = Hill(ip3,self.pars['d1'],1)
        n_inf = Hill(ca,self.pars['d5'],1)
        Q2 = self.pars['d2']
        h_inf = Hill(Q2,ca,1)
        self.po  = (m_inf*n_inf*h_inf)**3

    def bias(self,twin=[0.,1.],dt=0.1):
        # Return [t,bias(t)] as used in the simulations
        t = np.arange(twin[0],twin[-1]+dt,dt)
        bias = np.zeros(np.size(t))
        if 'gchi' in self._model:
            bias[np.fmod(t,self.pars['T'])<self.pars['pw']] = self.pars['yb']
        return np.vstack((t,bias))

    def integrate(self,
                  algparams=None,
                  normalized=False,**kwargs):
        self._algparams = algparams
        self._normalized = normalized
        if (not self._algparams):
            self._algparams = su.solver_opts(method='gsl',**kwargs)
        if self._model=='lra':
            self.sol = lr_model(self.pars,self._algparams)
        elif self._model=='chi':
            self.sol = chi_model(self.pars,self._algparams)
        elif self._model=='gchi':
            self.sol = gchi_model(self.pars,self._algparams)
        elif self._model=='pkc':
            self.sol = egchi_model(self.pars, self._algparams)

        if normalized:
            # Create in parallel an attribute "mm" (min/max) whereto save min and max of data to reverse normalization
            self.mm = {'ca': np.asarray([np.amin(self.sol['ca']), np.amax(self.sol['ca'])])}
            self.sol['ca'] = normalize(self.sol['ca'])
            if self._model in ['chi','gchi','pkc']:
                self.mm['ip3'] = np.asarray([np.amin(self.sol['ip3']), np.amax(self.sol['ip3'])])
                self.sol['ip3'] = normalize(self.sol['ip3'])
            if 'pkc' in self._model:
                self.mm['dag'] = np.asarray([np.amin(self.sol['dag']), np.amax(self.sol['dag'])])
                self.mm['pkc'] = np.asarray([np.amin(self.sol['pkc']), np.amax(self.sol['pkc'])])
                self.sol['dag'] = normalize(self.sol['dag'])
                self.sol['pkc'] = normalize(self.sol['pkc'])

    def ip3_production(self, solution=None):
        if (not solution):
            if hasattr(self, 'sol'):
                solution = cp.deepcopy(self.sol)
                if hasattr(self, 'mm'):
                    solution['ca'] = unnormalize(self.sol['ca'], self.mm['ca'])
                    solution['ip3'] = unnormalize(self.sol['ip3'], self.mm['ip3'])
            else:
                raise ValueError('A valid solution dictionary with \'ca\' and \'ip3\' keys must be provided')
        prod = {}
        if 'chi' in self._model:
            prod['Jdelta'] = self.pars['vdelta']*Hill(solution['ca'],self.pars['Kdelta'],2)*(1-Hill(solution['ip3'],self.pars['kappad'],1))
            prod['Jbeta']  = np.zeros(np.size(solution['ca']))
        if 'gchi' in self._model:
            prod['Jbeta'] = self.pars['vbeta']*solution['rec']
        return prod

    def ip3_degradation(self, solution=None):
        if (not solution):
            if hasattr(self, 'sol'):
                solution = cp.deepcopy(self.sol)
                if hasattr(self, 'mm'):
                    solution['ca'] = unnormalize(self.sol['ca'], self.mm['ca'])
                    solution['ip3'] = unnormalize(self.sol['ip3'], self.mm['ip3'])
            else:
                raise ValueError('A valid solution dictionary with \'ca\' and \'ip3\' keys must be provided')
        degr = {}
        degr['J3k'] = self.pars['v3k']*Hill(solution['ca'],self.pars['Kd'],4)*Hill(solution['ip3'],self.pars['K3'],1)
        degr['J5p'] = self.pars['r5p']*solution['ip3']
        return degr

    def __del__(self):
        del self.pars
        del self.ICs

if __name__ == "__main__":
    # -------------------------------------------------------------------------------------------------
    # Testing Astrocyte class
    # -------------------------------------------------------------------------------------------------
    options = su.solver_opts(t0=0.0, tfin=30.0, dt=1e-4, atol=1e-8, rtol=1e-6, method="gsl_msadams")

    # ##LRA model
    astro = Astrocyte(model='lra',ip3=0.8)
    astro.integrate(algparams=options,normalized=True)
    plt.plot(astro.sol['ts'], astro.sol['ca'], 'y-', astro.sol['ts'], astro.sol['h'], 'm-')

    # ##CHI model
    # astro = Astrocyte(model='chi',vbeta=1)
    # astro.integrate(algparams=options)
    # plt.plot(astro.sol['ts'], astro.sol['ip3'], 'c-', astro.sol['ts'], astro.sol['ca'], 'y-', astro.sol['ts'], astro.sol['h'], 'm-')

    # # ##G-CHI model
    # astro = Astrocyte(model='gchi',vbeta=1.0,T=2,pw=1,yb=10)
    # astro.integrate(algparams=options)
    # plt.plot(astro.sol['ts'], astro.sol['ip3'], 'c-', astro.sol['ts'], astro.sol['ca'], 'y-', astro.sol['ts'], astro.sol['h'], 'm-')

    # ##eG-CHI model
    # astro = Astrocyte(model='pkc',vbias=0.05)
    # ca_data = svu.loaddata('../data/codazzi_data.pkl')[0]
    # po_data = svu.loaddata('../data/fit_po.pkl')[0]
    # lr_data = svu.loaddata('../data/fit_lra.pkl')[0]
    # chi_data = svu.loaddata('../data/chi_fit_0_bee.pkl')[0]
    #
    # c0 = 7.0
    # rc = lr_data[0]
    # ver = 0.7 * lr_data[1]
    # vbeta = 0.2
    # vdelta = 0.05 * chi_data[1]
    # v3k = 0.05 * chi_data[2]
    # r5p = 0.5 * chi_data[3]
    #
    # # Reconstruct fit
    # astro = Astrocyte(model='pkc',
    #                          d1=po_data[0], d2=po_data[1], d3=po_data[0], d5=po_data[2], a2=po_data[3],
    #                          c0=c0, c1=0.5, rl=0.01, Ker=0.1, rc=rc, ver=ver,
    #                          yrel=1.5, vbeta=vbeta, vdelta=vdelta, v3k=v3k, r5p=r5p,
    #                          ICs=np.asarray([0.9, 0.6, 0.1, 0.5, 0.1, 0.1]))
    #
    # astro.integrate(algparams=options)
    # plt.plot(astro.sol['ts'], astro.sol['ip3'], 'c-', astro.sol['ts'], astro.sol['ca'], 'y-',
    #          astro.sol['ts'], astro.sol['dag'], 'r-', astro.sol['ts'], astro.sol['pkc'], 'm-')
    #
    ## Sample use of the routine
    plt.show()