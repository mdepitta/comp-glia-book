"""
solvers_utils.py
Library that defines solver related data structures.

Methods:
- solver_opts : create dictionary for settings of "continuous solvers such as
                Euler or RK4 solvers (in related H and CPP files).

v1.0
Maurizio De Pitta', The University of Chicago, March 1st, 2015.
"""

import numpy as np
from scipy import *
import sys, os
sys.path.append(os.path.join(os.path.expanduser('~'),'Ch6.DePitta.GPCR/pycustommodules'))
import general_utils as gu


def solver_opts(method='euler',**kwargs):
    """
    Define a standard dictionary to use in C/C++ simulators.
    
    Use:
    options = solver_opts(...)
    
    Input:
    - method    : {'euler'} | 'rk4' | 'gsl' | 'gsl_rk8pd' | 'gsl_msadams'
    - ndims     : {1} | Integer   Number of equations (i.e. dimension) of the system to integrate

    **kwargs:
        - t0        : initial instant of integration [s]
        - tfin      : final instant of integration [s]
        - transient : transient to drop from output result [s]
        - tbin      : time bin of solution [s]
        - dt        : step of integration [s]
        - solver    : string {"euler"} | "rk4"
    
    Output:
    - options   : dictionary of solver settings (keys are inputs).

    v1.3
    Added options for GSL solvers.
    Append 'solver' key at the end for all methods.
    Maurizio De Pitta', INRIA Rhone-Alpes, November 1st, 2017.

    """


    ## User-defined parameters
    if method in ['euler','rk4']:
        opts = {'t0'        : 0.0,
                'tfin'      : 20.0,
                'transient' : 0.0,
                'tbin'      : 1e-2,
                'dt'        : 1e-3
                }
    elif method in ['gsl', 'gsl_rk8pd', 'gsl_msadams']:
        opts = {'t0'     : 0.0,
                'tfin'   : 1.0,
                'dt'     : 1e-3,
                'atol'   : 1e-8,
                'rtol'   : 1e-6
                }
    else:
        # This is the case of solver 'None','none','steady_state'
        opts = {'nmax' : 1000, # Max number of iterations for the solver
                'atol' : 1e-10,# Absolute tolerance on error
                'rtol' : 0.0   # Relative tolerance on error (default: 0.0: not considered)
                }

    opts = gu.varargin(opts, **kwargs)
    for k,item in opts.iteritems():
        if (not isinstance(item, (basestring, bool)))&(not hasattr(item, '__len__'))&(item != None):
            opts[k] = float(item)

    if method in ['gsl','gsl_rk8pd']:
        opts['nstep'] = (opts['tfin']-opts['t0'])//opts['dt'] + 1 # +1 if when we start counting from zero

    if (not method) or (method in ['none','steady_state']):
        for p in ['nmax']:
            opts[p] = np.intc(opts[p])

    # Include solver specifier in the final dictionary
    opts['solver'] = method
    return opts

#-----------------------------------------------------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    algparams = solver_opts(method='gsl_rk8pd',ndims=7)
    print algparams