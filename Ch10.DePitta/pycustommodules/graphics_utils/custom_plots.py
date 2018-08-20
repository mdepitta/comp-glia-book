import numpy as np
from matplotlib import pyplot as plt

from .. import general_utils as gu

def stp_indexes(u,r):
    '''
    Compute short-term facilitation and depression indexes.

    Inputs:
    - u : Array-like, 'u' synaptic variable
    - r : Array-like, synaptic release 'r'

    Returns:
    - fac_idx,dep_idx,rec_idx :  facilitation/depression/recovery indexes for 'r' array

    v1.0
    Maurizio De Pitta', The University of Chicago, March 1, 2016.
    '''

    real_idx = lambda log_idx : log_idx.nonzero()[0]+1

    du = np.diff(u)
    dr = np.diff(r)

    # Logical indexes
    fac_idx = (dr>0.)&(du>0.)
    dep_idx = dr<0.
    rec_idx = ~fac_idx&~dep_idx

    return real_idx(fac_idx),real_idx(dep_idx),real_idx(rec_idx)

def stp_raster(t_spk,u,r,syn_index,**kwargs):
    '''
    Synaptic raster plot for short-term plasticity (invokes stp_indexes).

    INPUT ARGUMENTS:
    - t_spk     : Array-like (1xN), spike timing
    - u         : Array-like (1xN), 'u' synaptic variable
    - r         : Array-like (1xN), synaptic release 'r'
    - syn_index : Array-like (1xN), neuron indexes in t_spk
    - **kwargs  :
      - color   : Dictionary of colors for rasterization, e.g. {'fac' : <color>, 'dep': <color>, 'rec': <color>}
      - order   : String (1x3) for plotting order, e.g. {'fdr'}
      - lw      : LineWidth
      - ax      : {None} | Axis Object wherein to plot raster

    OUTPUT:
    - ax        : Axis object for raster plot

    v1.0
    Maurizio De Pitta', The University of Chicago, March 1, 2016.
    '''

    opts = {'color' : {'fac' : 'm', 'dep': 'g', 'rec': 'k'},
            'order' : 'fdr',
            'lw'    : .6,
            'ax'    : None,
            'N'     : 0}

    # User-defined parameters
    opts = gu.varargin(opts,**kwargs)

    if opts['ax']==None:
        _,ax = plt.subplots(1,1)
    else:
        ax = opts['ax']

    if 'N' in kwargs:
        N = int(opts['N'])
    else:
        N = int(max(syn_index)+1)

    for i in xrange(N):
        fac_idx,dep_idx,rec_idx = stp_indexes(u[syn_index==i],r[syn_index==i])
        for _,v in enumerate(opts['order']):
            if v=='f':
                ax.vlines(t_spk[syn_index==i][fac_idx], i-.5, i+.5, color=opts['color']['fac'], linewidth=opts['lw'])
            elif v=='d':
                ax.vlines(t_spk[syn_index==i][dep_idx], i-.5, i+.5, color=opts['color']['dep'], linewidth=opts['lw'])
            else :
                ax.vlines(t_spk[syn_index==i][rec_idx], i-.5, i+.5, color=opts['color']['rec'], linewidth=opts['lw'])
    ax.set_ylim([-.5,N-0.5])

    return ax