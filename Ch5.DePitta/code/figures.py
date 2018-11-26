from __future__ import division
import numpy as np
import gc # Garbage collector to contain possible memory leaks due to non-optimized code
import scipy.signal as sps
import scipy.io as io
import sys, os
base_dir = '/Ch4.DePitta'
sys.path.append(os.path.join(os.path.expanduser('~'),base_dir+'/pycustommodules'))
import pycustommodules.general_utils as gu
import pycustommodules.solvers.solver_utils as su
import pycustommodules.save_utils as svu
from pycustommodules.graphics_utils import plot_utils as pu

import matplotlib.pylab as plt

import astrocyte_models as models
import fit_data as fd
import egchi_bif as bif

# Some useful Lambdas
# Percentage variation
percentdiff = lambda data : (data-data[0])/data[0]*100

# Data normalization
def normalize(data,min=None,max=None):
    if (not min) and (not max):
        return (data-np.amin(data))/(np.amax(data)-np.amin(data))
    elif min!=None and (not max):
        return (data - min) / (np.amax(data) - min)
    elif max!=None and (not min):
        return (data - np.amin(data)) / (max - np.amin(data))
    else:
        return (data - min) / (max - min)

#-----------------------------------------------------------------------------------------------------------------------
# Figure 1
#-----------------------------------------------------------------------------------------------------------------------
def chi_model():
    """
    Provide sample fits for sole ChI model.

    """

    # Load MPL default style
    plt.style.use('figures.mplstyle')

    #--------------------------------------------------------------------------
    # Plot open probability
    #--------------------------------------------------------------------------
    p_open = svu.loaddata('../data/p_open.pkl')[0]
    po_data = svu.loaddata('../data/fit_po.pkl')[0]

    # Load significant data
    astro = models.Astrocyte(model='lra',
                             d1=po_data[0], d2=po_data[1], d3=po_data[0], d5=po_data[2], a2=po_data[3])

    # # p_open vs. IP3 was obtained with 0.25 uM free Ca2+
    astro.popen(0.25, p_open['ip3'][0])
    po_ip3 = astro.po
    # print po_ip3

    # p_open vs. Ca2+ was obtained with 1 uM IP3
    astro.popen(p_open['ca'][0], 1.0)
    po_ca = astro.po

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 4), sharex=False,
                           gridspec_kw={'top': 0.98, 'bottom': 0.15,
                                        'left': 0.1, 'right': 0.98})

    ax[0].plot(p_open['ip3'][0],p_open['ip3'][1],'k^',ms=10,ls='none',label=r'Ca$^{2+}$=25 nM')
    ax[0].plot(p_open['ca'][0],p_open['ca'][1],'ko',ms=10,ls='none', label=r'IP$_3$=1 $\mu$M')
    ax[0].plot(p_open['ip3'][0],po_ip3,'k--')
    ax[0].plot(p_open['ca'][0], po_ca, 'k-')
    # Adjust axes
    pu.adjust_spines(ax[0], ['left', 'bottom'])
    xticks = 10**np.arange(-2,2.1,1)
    ax[0].set(xlim = (0.008,500),
              ylim = (0., 0.5)  , yticks = np.arange(0.,0.51,0.1),
              xscale='log', xlabel=r'Ca$^{2+}$, IP$_3$ ($\mu$M)', ylabel='IP$_3$R open probability')
    ax[0].set(xticks = xticks, xticklabels = ([str(xt) for xt in xticks]))
    # Add legend
    ax[0].legend(loc='lower right', facecolor='w', edgecolor='none')

    #--------------------------------------------------------------------------
    # Plot ChI model from data
    #--------------------------------------------------------------------------
    # Load relevant data
    ca_data = svu.loaddata('../data/herzog_data.pkl')[0]
    lr_data = svu.loaddata('../data/fit_lra.pkl')[0]
    c0 = 5.0
    iset = 2

    colors = {'data' : '0.65',
              'ca'   : 'k',
              'ip3'  : pu.hexc((44,160,44)),
              'h'    : pu.hexc((23,190,207))
              }

    # Colors for different data sets
    filename = 'chi_fit_0_bee'
    x_rescaled = svu.loaddata('../data/' + filename + '.pkl')[0]
    # Fitted data
    astro = models.Astrocyte(model='chi',
                             d1=po_data[0], d2=po_data[1], d3=po_data[0], d5=po_data[2], a2=po_data[3],
                             c0=c0, c1=0.5, rl=0.1, Ker=0.1, rc=lr_data[0], ver=lr_data[1],
                             vbeta=x_rescaled[0], vdelta=x_rescaled[1], v3k=x_rescaled[2], r5p=x_rescaled[3],
                             ICs=np.asarray([x_rescaled[6], x_rescaled[4], x_rescaled[5]]))
    options = su.solver_opts(t0=ca_data['time'][iset][0], tfin=ca_data['time'][iset][-1], dt=1e-4, atol=1e-8,
                             rtol=1e-6, method="gsl_msadams")
    astro.integrate(algparams=options, normalized=True)

    # Original data
    ax[1].plot(ca_data['time'][iset], ca_data['original'][iset], '.', ms=4, color=colors['data'], ls='none')
    # Fit data
    ax[1].plot(astro.sol['ts'], astro.sol['ca'], '-', color=colors['ca'], label='C')
    ax[1].plot(astro.sol['ts'], astro.sol['h'], '-', color=colors['h'], label='h')
    ax[1].plot(astro.sol['ts'], astro.sol['ip3'],'-', color=colors['ip3'], label='I')

    ax[1].set(xlim = (0,30.), xticks = np.arange(0,31,5),
              ylim = (0., 1.01)  , yticks = np.arange(0.,1.1,0.2),
              xlabel='Time (s)', ylabel='$ChI$ variables (n.u.)')
    pu.adjust_spines(ax[1], ['left', 'bottom'])
    pu.adjust_ylabels(ax, x_offset=-0.1)
    # Add legend
    ax[1].legend(loc='lower right', facecolor='w', edgecolor='none')

    del astro
    plt.savefig('../Figures/f1_chi_model.eps', dpi=600)
    plt.show()

#-----------------------------------------------------------------------------------------------------------------------
# Figure 2
#-----------------------------------------------------------------------------------------------------------------------
def gchi_model():
    """
    Simulation of the G-ChI model with some pulse function
    """

    # Load relevant data from fits
    po_data = svu.loaddata('../data/fit_po.pkl')[0]
    lr_data = svu.loaddata('../data/fit_lra.pkl')[0]
    chi_data = svu.loaddata('../data/chi_fit_0_bee.pkl')[0]

    # Sample parameters
    c0 = 10.0
    ver = 10.0
    vbeta = 0.8
    v3k = chi_data[3]
    # Simulation duration
    t0 = 0.0
    tfin = 10.0
    # Generate model
    astro = models.Astrocyte(model='gchi',
                             d1=po_data[0], d2=po_data[1], d3=po_data[0], d5=po_data[2], a2=po_data[3],
                             c0=c0, c1=0.5, rl=0.1, Ker=0.1, rc=lr_data[0], ver=ver,
                             vbeta=vbeta, vdelta=2*chi_data[1], v3k=v3k, r5p=chi_data[3],
                             ICs=np.asarray([0.0,0.05, 0.05, 0.9]))
    options = su.solver_opts(t0=t0, tfin=30., dt=1e-4, atol=1e-8,
                             rtol=1e-6, method="gsl_msadams")
    # First run to seek resting state
    astro.integrate(algparams=options, normalized=False)
    # Update model with new ICs and add pulse train
    options = su.solver_opts(t0=t0, tfin=tfin, dt=1e-4, atol=1e-8,
                             rtol=1e-6, method="gsl_msadams")
    astro.ICs = np.r_[astro.sol['rec'][-1],astro.sol['ip3'][-1],astro.sol['ca'][-1],astro.sol['h'][-1]]
    astro.pars['T'] = (tfin-t0)/(3.0*tfin)
    astro.pars['pw'] = 0.15*astro.pars['T']
    astro.pars['yb'] = 8.
    # Effective simulation
    astro.integrate(algparams=options, normalized=True)
    # Retrieve bias
    bias = astro.bias(twin=[t0,tfin],dt=0.01)
    # Retrieve production and degradation terms
    prod = astro.ip3_production()
    degr = astro.ip3_degradation()

    #---------------------------------------------------------------------
    # Plot figure
    #---------------------------------------------------------------------
    # Load MPL default style
    plt.style.use('figures.mplstyle')
    colors = {'data' : '0.65',
              'ca'   : 'k',
              'ip3'  : pu.hexc((44,160,44)),
              'Jbeta'  : pu.hexc((255,127, 14)), # Orange
              'Jdelta' : pu.hexc((31, 119,180)), # Tourquois
              'J3k'    : pu.hexc((214, 39, 40)), # DarkRed
              'J5p'    : pu.hexc((227,119,194))  # Magenta
              }
    fig, ax = plt.subplots(nrows=5, ncols=1, figsize=(8, 6.5),
                           sharex=False,
                           gridspec_kw={'height_ratios': [0.5, 3, 3, 3, 3],
                                        'top': 0.96, 'bottom': 0.12,
                                        'left': 0.14, 'right': 0.95})

    # Plot stimulation
    ax[0].plot(bias[0],bias[1],'k-')
    ax[0].set(xlim = (t0,tfin), xticks=[],
              ylim = (-0.2, 1.05*astro.pars['yb']), yticks = [],
              ylabel='Glu\n($\mu$M)')
    pu.adjust_spines(ax[0], ['left'])

    # Plot receptor activation
    ax[1].plot(astro.sol['ts'],astro.sol['rec'],'k-')
    pu.adjust_spines(ax[1], ['left'])
    ax[1].set(xlim = (t0,tfin), xticks=[],
              ylim = (0., 0.3), yticks = np.arange(0.,0.31,0.1),
              ylabel='Ast. Rec.\n$\Gamma_A$')

    # Plot Ca2+ / IP3
    ax[2].plot(astro.sol['ts'], astro.sol['ca'], '-', color=colors['ca'], label='C')
    ax[2].plot(astro.sol['ts'], astro.sol['ip3'],'-', color=colors['ip3'], label='I')
    pu.adjust_spines(ax[2], ['left'])
    ax[2].set(xlim = (t0,tfin), xticks=[],
              ylim = (-0.05, 1.01), yticks = np.arange(0,1.1,0.5),
              ylabel='Ca$^{2+}$, IP$_3$\n(n.u.)')
    ax[2].legend(loc='lower right', facecolor='w', edgecolor='none')

    # Plot IP3 production
    ax[3].plot(astro.sol['ts'], prod['Jbeta'], '-', color=colors['Jbeta'], label=r'J$_\beta$')
    ax[3].plot(astro.sol['ts'], prod['Jdelta'],'-', color=colors['Jdelta'], label=r'J$_\delta$')
    ax[3].plot(astro.sol['ts'], prod['Jbeta']+prod['Jdelta'], 'k--')
    pu.adjust_spines(ax[3], ['left'])
    ax[3].set(xlim = (t0,tfin), xticks=[],
              ylim = (-0.02, 0.2), yticks = np.arange(0,0.21,0.1),
              ylabel='IP$_3$ prod.\n($\mu$M/s)')
    ax[3].legend(loc='lower right', facecolor='w', edgecolor='none')

    # Plot IP3 degradation
    ax[4].plot(astro.sol['ts'], degr['J5p'], '-', color=colors['J5p'], label=r'J$_{5P}$')
    ax[4].plot(astro.sol['ts'], degr['J3k'], '-', color=colors['J3k'], label=r'J$_{3K}$')
    ax[4].plot(astro.sol['ts'], degr['J5p']+degr['J3k'], 'k--')
    pu.adjust_spines(ax[4], ['left','bottom'])
    ax[4].set(xlim = (t0,tfin), xticks=np.arange(t0,tfin+1,2.),
              ylim = (-0.05, 0.22), yticks = np.arange(0,0.21,0.1),
              xlabel='Time (s)', ylabel='IP$_3$ degr.\n($\mu$M/s)')
    ax[4].legend(loc='lower right', facecolor='w', edgecolor='none')

    # Final adjustment of y-labels
    pu.adjust_ylabels(ax, x_offset=-0.07)

    del astro
    plt.savefig('../Figures/f2_gchi_model.eps', dpi=600)
    plt.show()

def compute_thr(**kwargs):
    """
    Simulation routine to compute CICR threshold
    """

    # Load relevant data from fits
    po_data = svu.loaddata('../data/fit_po.pkl')[0]
    lr_data = svu.loaddata('../data/fit_lra.pkl')[0]
    chi_data = svu.loaddata('../data/chi_fit_0_bee.pkl')[0]

    # Control parameter set
    pars = models.gchi_parameters(d1 = po_data[0], d2 = po_data[1], d3 = po_data[0], d5 = po_data[2], a2 = po_data[3],
                                  c0 = 10.0, c1 = 0.5, rl = 0.1, Ker = 0.1, rc = lr_data[0], ver = 10.0,
                                  vbeta = 0.8, vdelta = 2 * chi_data[1], v3k = chi_data[3], r5p = chi_data[3])
    # Custom parameters
    pars = gu.varargin(pars,**kwargs)

    # Simulation duration
    t0 = 0.0
    tfin = 20.0

    # Generate model
    astro = models.Astrocyte(model='gchi',
                             d1=pars['d1'], d2=pars['d2'], d3=pars['d3'], d5=pars['d5'], a2=pars['a2'],
                             c0=pars['c0'], c1=pars['c1'], rl=pars['rl'], Ker=pars['Ker'], rc=pars['rc'], ver=pars['ver'],
                             vbeta=pars['vbeta'], vdelta=pars['vdelta'], v3k=pars['v3k'], r5p=pars['r5p'],
                             ICs=np.asarray([0.0,0.05, 0.05, 0.9]))
    options = su.solver_opts(t0=t0, tfin=30., dt=1e-4, atol=1e-8,rtol=1e-6, method="gsl_msadams")
    # First run to seek resting state
    astro.integrate(algparams=options, normalized=False)
    # Update model with new ICs and add pulse train
    options = su.solver_opts(t0=t0, tfin=tfin, dt=1e-4, atol=1e-8,
                             rtol=1e-6, method="gsl_msadams")
    astro.ICs = np.r_[astro.sol['rec'][-1],astro.sol['ip3'][-1],astro.sol['ca'][-1],astro.sol['h'][-1]]
    astro.pars['T'] = tfin
    astro.pars['pw'] = tfin
    yb = np.arange(0.5,5.0,0.05)
    # yb = np.arange(1.0, 3.0, 0.3)
    Nt = np.size(yb)
    # Allocate results
    ca_traces = [None] * Nt
    ip3_traces = [None] * Nt
    bias_traces = [None] * Nt
    for i in xrange(Nt):
        astro.pars['yb'] = yb[i]
        # Effective simulation
        astro.integrate(algparams=options, normalized=False)
        # Save results
        ca_traces[i] = astro.sol['ca']
        ip3_traces[i] = astro.sol['ip3']
        bias_traces[i] = astro.bias(twin=[t0,tfin],dt=0.01)

    traces = {'yb': yb, 'ts': astro.sol['ts'], 'ca': ca_traces, 'ip3': ip3_traces, 'bias': bias_traces}
    del astro
    gc.collect()
    return traces,pars


def compute_latency(traces, dthr=0.5):
    """
    Compute latency of CICR for constant puff application. Detect CICR on the first value of dC/dt > dthr.
    Traces are as those provided by compute_thr().

    Return 3x1 tuple with (bias values, latency, IP3 injected).
    """

    Nt = np.size(traces['yb'])
    thr = []
    for i in xrange(Nt):
        index = np.where(np.diff(traces['ca'][i]) / np.diff(traces['ts']) > dthr)[0]
        if np.size(index) > 0:
            thr.append((traces['yb'][i], traces['ts'][index[0]],
                        np.sum(np.multiply(traces['ip3'][i][:index[0]], traces['ts'][:index[0]]))))
    bv, latency, inj = zip(*thr)
    return bv, latency, inj

def gchi_threshold(sim=False):
    """
    Simulation of the G-ChI model with some pulse function
    """

    if sim:
        print "Control"
        tr_0, p = compute_thr()
        print "r5p"
        tr_5p,_ = compute_thr(r5p=1.5*p['r5p'])
        print "v3k"
        tr_3k, _ = compute_thr(v3k=1.5*p['v3k'])
        print "vbeta"
        tr_vb, _ = compute_thr(vbeta=1.5 * p['vbeta'])
        print "vdelta"
        tr_vd, _ = compute_thr(vdelta=1.5 * p['vbeta'])

        print "Saving data"
        svu.savedata([tr_0,tr_vb,tr_vd,tr_3k,tr_5p],'../data/cicr_thr.pkl')

        print "done!"

    #---------------------------------------------------------------------
    # Plot figure
    #---------------------------------------------------------------------
    # Load MPL default style
    plt.style.use('figures.mplstyle')

    # Load data
    tr_0, tr_vb, tr_vd, tr_3k, tr_5p = svu.loaddata('../data/cicr_thr.pkl')

    #---------------------------------------------------------------------------------------------------------
    # Showing threshold
    #---------------------------------------------------------------------------------------------------------
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(6, 5.5),
                           sharex=False,
                           gridspec_kw={'height_ratios': [0.75, 2, 2],
                                        'top': 0.96, 'bottom': 0.12,
                                        'left': 0.14, 'right': 0.95})
    Nt = np.size(tr_0['yb'])
    t0 = -0.2
    tfin = 5.0
    mask = tr_0['ts'] < tfin
    maskb = tr_0['bias'][0][0] < tfin
    for i in xrange(10,51,5):
        # Plot adding some points on the left to show instant of stimulation
        ax[0].plot(np.r_[t0,0,tr_0['bias'][i][0][maskb]], np.r_[0,0,tr_0['bias'][i][1][maskb]])
        ax[1].plot(np.r_[t0,tr_0['ts'][mask]],np.r_[tr_0['ip3'][i][0],tr_0['ip3'][i][mask]])
        ax[2].plot(np.r_[t0,tr_0['ts'][mask]],np.r_[tr_0['ca'][i][0],tr_0['ca'][i][mask]])

    # Annotate
    ax[0].plot(0., 3.6, 'v', color='k', ms=10, clip_on = False)
    ax[1].plot(0., 0.065, 'v', color='k', ms=10, clip_on = False)
    ax[2].plot(0., 0.25, 'v', color='k', ms=10, clip_on = False)

    # Adjust axes
    pu.adjust_spines(ax[0], ['left'])
    ax[0].set(xlim = (t0,tfin), xticks=[],
              ylim = (-0.05, 3.05), yticks = np.arange(0,3.1,1.0),
              ylabel='Glu\n($\mu$M)')
    pu.adjust_spines(ax[1], ['left'])
    ax[1].set(xlim = (t0,tfin), xticks=[],
              ylim = (0., 0.2), yticks = np.arange(0,0.21,0.1),
              ylabel='IP$_3$ ($\mu$M)')
    pu.adjust_spines(ax[2], ['left','bottom'])
    ax[2].set(xlim = (t0,tfin), xticks=np.arange(0,5.1),
              ylim = (0., 2.5), yticks = np.arange(0,2.51,1.0),
              xlabel='Time (s)', ylabel='Ca$^{2+}$ ($\mu$M)')
    # Final adjustment of y-labels
    pu.adjust_ylabels(ax, x_offset=-0.07)

    # Save figure
    plt.savefig('../Figures/f3_gchi_thr_0.eps', dpi=600)

    #---------------------------------------------------------------------------------------------------------
    # Compute latency and IP3 injected
    #---------------------------------------------------------------------------------------------------------
    bv, latency, inj = {}, {}, {}
    # Control
    key = 'ctrl'
    bv[key], latency[key], inj[key] = compute_latency(tr_0, dthr=0.5)
    key = 'r5p'
    bv[key], latency[key], inj[key] = compute_latency(tr_5p, dthr=0.5)
    key = 'v3k'
    bv[key], latency[key], inj[key] = compute_latency(tr_3k, dthr=0.5)
    key = 'vb'
    bv[key], latency[key], inj[key] = compute_latency(tr_vb, dthr=0.5)
    key = 'vd'
    bv[key], latency[key], inj[key] = compute_latency(tr_vd, dthr=0.5)

    colors = {'ctrl'   : 'k',
              'vb'     : pu.hexc((255,127, 14)), # Orange
              'vd'     : pu.hexc((31, 119,180)), # Tourquois
              'v3k'    : pu.hexc((214, 39, 40)), # DarkRed
              'r5p'    : pu.hexc((227,119,194))  # Magenta
              }

    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(5, 5.5),
                           sharex=False,
                           gridspec_kw={'height_ratios': [1, 1],
                                        'top': 0.96, 'bottom': 0.12,
                                        'left': 0.14, 'right': 0.95})
    plt.subplots_adjust(hspace=0.3)

    keys = ['ctrl','vb','vd','v3k','r5p']
    label = ['reference',
             '$O_{\\beta} \uparrow$ ($\\times 1.5$)',
             '$O_\delta \uparrow$ ($\\times 1.5$)',
             '$O_{3K} \uparrow$ ($\\times 1.5$)',
             '$\Omega_{5P} \uparrow$ ($\\times 1.5$)']
    for _,k in enumerate(keys):
        if k=='ctrl':
            ax[0].plot(bv[k], latency[k], color=colors[k], ls='--', label=label[_], zorder=10)
            ax[1].plot(latency[k], inj[k], color=colors[k], ls='--', zorder=10)
        else:
            ax[0].plot(bv[k],latency[k], color=colors[k], ls='-', label=label[_])
            ax[1].plot(latency[k], inj[k], color=colors[k], ls='-')

    # Adjust axes
    pu.adjust_spines(ax[0], ['left','bottom'])
    ax[0].set(xlim = (0.,5.0), xticks=np.arange(0,5.1),
          ylim = (0.,4.), yticks = np.arange(0,4.1,1.0),
          xlabel='Glu ($\mu$M)',ylabel='Latency (s)')
    ax[1].set(xlim = [0.,4.0], xticks=np.arange(0.,4.1,1.0),
              ylim = [100,10000], yscale='log',
              xlabel='Latency (s)',
              ylabel='IP$_3$ Time Integral (mM$\cdot$s)')
    pu.adjust_spines(ax[1], ['left','bottom'])
    ax[1].set(yticks = [100,1000,10000], yticklabels = [str(i) for i in [0.1,1,10]])
    # # Final adjustment of y-labels
    pu.adjust_ylabels(ax, x_offset=-0.09)
    # Add legend
    ax[0].legend(loc='upper right', facecolor='w', edgecolor='none')

    # # Save figure
    plt.savefig('../Figures/f3_gchi_thr_1.eps', dpi=600)

    plt.show()

#-----------------------------------------------------------------------------------------------------------------------
# Figure 4
#-----------------------------------------------------------------------------------------------------------------------
def pkc_effect(sim=False):
    """
    Shows ExGChI and compare with data from Codazzi et al (CON 2001)
    """

    # Load relevant data
    ca_data = svu.loaddata('../data/codazzi_data.pkl')[0]
    po_data = svu.loaddata('../data/fit_po.pkl')[0]
    lr_data = svu.loaddata('../data/fit_lra.pkl')[0]
    chi_data = svu.loaddata('../data/chi_fit_0_bee.pkl')[0]

    if sim:
        print "1. Simulate PKC dynamics according to Codazzi"
        # yrel = 1.508
        yrel = 1.48
        # yrel = 1.448
        c0 = 5.0
        rc = 0.8 * lr_data[0]
        ver = lr_data[1]
        vbeta = 1.0
        vdelta = chi_data[1]
        v3k = 1.4 * chi_data[2]
        r5p = chi_data[3]

        vk     = 1.0
        vkd    = 0.28
        OmegaKD= 0.33
        vd     = 0.45
        OmegaD = 0.26

        astro = models.Astrocyte(model='pkc',
                                 d1=po_data[0], d2=po_data[1], d3=po_data[0], d5=po_data[2], a2=po_data[3],
                                 c1=0.5, rl=0.01, Ker=0.1, rc=rc, ver=ver, Kkc=0.5, Kkd=0.1, Kdd=0.1,
                                 c0=c0, yrel=yrel, vbeta=vbeta, vdelta=vdelta, v3k=v3k, r5p=r5p,
                                 vk=vk, vkd=vkd, OmegaKD=OmegaKD, vd=vd, OmegaD=OmegaD,
                                 ICs=np.asarray([0.0, 0.05, 0.05, 0.9, 0.0, 0.0]))

        # Preliminary integration to start from resting conditions
        options = su.solver_opts(t0=0.0, tfin=60., dt=1e-4, atol=1e-4,rtol=1e-4, method="gsl_msadams")
        astro.integrate(algparams=options, normalized=False)
        # Update model with new ICs and add pulse train
        options = su.solver_opts(t0=0.0, tfin=ca_data['time'][-1], dt=1e-4, atol=1e-4,rtol=1e-4, method="gsl_msadams")
        astro.ICs = np.r_[astro.sol['rec'][-1],astro.sol['ip3'][-1],astro.sol['ca'][-1],astro.sol['h'][-1],astro.sol['dag'],astro.sol['pkc']]
        astro.integrate(algparams=options, normalized=False)
        sol_fit = astro.sol
        # _,sol_0,sol_1,sol_2 = svu.loaddata('../data/pkc_qual_fit.pkl')
        # svu.savedata([sol_fit,sol_0,sol_1,sol_2],'../data/pkc_qual_fit.pkl')

        # Second simulation
        print "2. Show effect of O_K"
        # First standard astrocyte without
        astro.pars['yrel'] = 1.55
        astro.pars['vkd'] = 0.0
        astro.pars['vk']  = 0.0
        options = su.solver_opts(t0=0.0, tfin=30, dt=1e-4, atol=1e-4,rtol=1e-4, method="gsl_msadams")
        astro.integrate(algparams=options, normalized=False)
        sol_0 = astro.sol

        astro.pars['vkd']  = 0.28
        astro.pars['vk']   = 1.0
        astro.integrate(algparams=options, normalized=False)
        sol_1 = astro.sol

        astro.pars['vk']   = 3.0
        astro.integrate(algparams=options, normalized=False)
        sol_2 = astro.sol

        # Save data
        svu.savedata([sol_fit,sol_0,sol_1,sol_2],'../data/pkc_qual_fit.pkl')

        # # Third simulation (requires PyDSTool)
        print "3. Computing oscillation period (requires PyDSTool)"
        # Effective computation
        pkc0 = bif.compute_period(vkd=0.001, vk=0.001)
        svu.savedata([pkc0], '../data/pkc_period_Ok.pkl')
        pkc1 = bif.compute_period(vkd=0.28, vk=1.0)
        svu.savedata([pkc0, pkc1], '../data/pkc_period_Ok.pkl')
        pkc2 = bif.compute_period(vkd=0.28, vk=3.0)
        svu.savedata([pkc0, pkc1, pkc2], '../data/pkc_period_Ok.pkl')
        pkc3 = bif.compute_period(vkd=0.84, vk=1.0)
        svu.savedata([pkc0, pkc1, pkc2, pkc3], '../data/pkc_period_Ok.pkl')
        pkc4 = bif.compute_period(vkd=0.84, vk=3.0)
        svu.savedata([pkc0, pkc1, pkc2, pkc3, pkc4], '../data/pkc_period_Ok.pkl')
    #---------------------------------------------------------------------
    # Plot figure
    #---------------------------------------------------------------------
    # Load MPL default style
    plt.style.use('figures.mplstyle')

    # Load data
    sol, pkc0, pkc1, pkc2 = svu.loaddata('../data/pkc_qual_fit.pkl')

    colors = {'ctrl': '0.7'}

    #---------------------------------------------------------------------------------------------------------
    # Fit with experimental data
    #---------------------------------------------------------------------------------------------------------
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(6, 5.5),
                           sharex=False,
                           gridspec_kw={'height_ratios': [1, 1],
                                        'top': 0.96, 'bottom': 0.12,
                                        'left': 0.15, 'right': 0.95})

    # # Plot experimental data
    ax[0].plot(ca_data['original']['ca'][:,0], normalize(ca_data['original']['ca'][:,1]), 'k-', label='Ca$^{2+}$')
    ax[0].plot(ca_data['original']['pkc'][:,0], normalize(ca_data['original']['pkc'][:,1]), 'r-', label='cPKC$^*$')

    # Plot solution
    tclip = [3,27]
    # # Mask
    cmn, _ = fd.relextrema(sol['ca'])
    pmn, _ = fd.relextrema(sol['pkc'])
    cbase = (sol['ca'][cmn].max()-sol['ca'][cmn].min())/2
    # pbase = (sol['pkc'][pmn].max() - sol['pkc'][pmn].min())/2
    cbase += sol['ca'][cmn].min()
    pbase = sol['pkc'][pmn].min()
    cmask = np.where((tclip[0]<=sol['ts']) & (sol['ts']<tclip[1]) & (sol['ca']>cbase))[0]
    pmask = np.where((tclip[0]<=sol['ts']) & (sol['ts']<tclip[1]) & (sol['pkc']>pbase))[0]
    # # Effective plotting
    ax[1].plot(sol['ts'][cmask]-sol['ts'][cmask[0]], normalize(sol['ca'][cmask],min=cbase), 'k-')
    ax[1].plot(sol['ts'][pmask]-sol['ts'][cmask[0]], normalize(sol['pkc'][pmask],min=pbase), 'r-')

    # Fit Check
    # ax[1].plot(sol['ts'], normalize(sol['ca']), 'k-')
    # ax[1].plot(sol['ts'], normalize(sol['pkc']), 'r-')
    # ax[0].plot(sol['ts'], sol['ca'], 'k-', sol['ts'][cmn], sol['ca'][cmn], 'ko', sol['ts'][[0,-1]], [cbase, cbase], 'b-', ms=6)
    # ax[0].plot(sol['ts'], sol['pkc'], 'r-', sol['ts'][pmn], sol['pkc'][pmn], 'ro', ms=6)

    # Adjust axes
    pu.adjust_spines(ax[0], ['left','bottom'])
    ax[0].set(xlim = (0.,140), xticks=np.arange(0,141,35),
          ylim = (-0.02,1.05), yticks = np.arange(0,1.1,0.5),
          ylabel='Experimental\nTraces (n.u.)')
    pu.adjust_spines(ax[1], ['left', 'bottom'])
    ax[1].set(xlim = (0.,20), xticks=np.arange(0.,30,5),
              ylim=(-0.02, 1.05), yticks=np.arange(0, 1.1, 0.5),
              xlabel='Time (s)',
              ylabel='Simulated\nTraces (n.u.)')
    # # Final adjustment of y-labels
    pu.adjust_ylabels(ax, x_offset=-0.08)
    # Add legend
    ax[0].legend(loc='upper right', facecolor='w', edgecolor='none')

    # Save figure
    plt.savefig('../Figures/f4_pkc_model.eps', dpi=600)

    # ---------------------------------------------------------------
    # Build effect of O_k on period
    # ---------------------------------------------------------------
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(6, 4.5),
                           sharex=False,
                           gridspec_kw={'height_ratios': [1, 1, 1],
                                        'top': 0.96, 'bottom': 0.12,
                                        'left': 0.17, 'right': 0.95})

    # C
    ax[0].plot(pkc0['ts'], pkc0['ca'], '-', color=colors['ctrl'])
    ax[0].plot(pkc1['ts'], pkc1['ca'], 'k-')
    ax[0].plot(pkc2['ts'], pkc2['ca'], 'k-.')
    # D
    # First line is a dummy one to add proper legend in this plot rather than on the top one
    # ax[1].plot(pkc0['ts'][[0,1]],pkc0['dag'][[0,1]],'-',color=colors['ctrl'],label='$O_{K}=0$')
    ax[1].plot(pkc0['ts'], pkc0['dag'], '-', color=colors['ctrl'], label='$O_{K}=0$')
    ax[1].plot(pkc1['ts'], pkc1['dag'], 'k-',label='$O_{K}=1.0$ $\mu$M$^{-1}$s$^{-1}$')
    ax[1].plot(pkc2['ts'], pkc2['dag'], 'k-.',label='$O_{K}=3.0$ $\mu$M$^{-1}$s$^{-1}$')
    # P
    ax[2].plot(pkc0['ts'], pkc0['pkc'], '-', color=colors['ctrl'], label='$O_{K}=0$')
    ax[2].plot(pkc1['ts'], pkc1['pkc'], 'k-', label='$O_{K}=1.0$ $\mu$M$^{-1}$s$^{-1}$')
    ax[2].plot(pkc2['ts'], pkc2['pkc'], 'k-.', label='$O_{K}=3.0$ $\mu$M$^{-1}$s$^{-1}$')

    # Adjust axes
    pu.adjust_spines(ax[0], ['left'])
    ax[0].set(xlim = (0.,30), xticks=[],
          ylim = (-0.02,0.6), yticks = np.arange(0,0.61,0.2),
          ylabel='Ca$^{2+}$\n($\mu$M)')
    pu.adjust_spines(ax[1], ['left'])
    ax[1].set(xlim = (0.,30), xticks=[],
          ylim = (-0.02,0.41), yticks = np.arange(0,0.4,0.1),
          ylabel='DAG\n($\mu$M)')
    pu.adjust_spines(ax[2], ['left', 'bottom'])
    ax[2].set(xlim = (0.,30), xticks=np.arange(0.,31,10),
              ylim=(-0.01, 0.075), yticks=np.arange(0, 0.071, 0.02), yticklabels=[str(int(i*1e3)) for i in np.arange(0, 0.071, 0.02)],
              xlabel='Time (s)',
              ylabel='cPKC$^{*}$\n(nM)')
    # # Final adjustment of y-labels
    pu.adjust_ylabels(ax, x_offset=-0.08)
    # Add legend
    ax[2].legend(loc='upper right', facecolor='w', edgecolor='none')

    # Save figure
    plt.savefig('../Figures/f4_pkc_osc_1.eps', dpi=600)

    #---------------------------------------------------------------
    # Plot period of oscillations
    #---------------------------------------------------------------
    p0,p1,p2,p3,p4 = svu.loaddata('../data/pkc_period_Ok.pkl')

    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(6, 4.5),
                           sharex=False,
                           gridspec_kw={'height_ratios': [1, 2],
                                        'top': 0.96, 'bottom': 0.12,
                                        'left': 0.17, 'right': 0.95})

    for i,p in enumerate([p0,p1,p2,p3,p4]):
        # idx = (p['T'] > 3.5)
        idx = np.argsort(p['T'][(p['T'] > 3.5)])[::-1]
        if i==3: idx = idx[0:-2]
        p['T'] = p['T'][idx]
        p['yrel'] = p['yrel'][idx]
    mask = lambda period: np.argsort(1.0 / period)
    ax[1].plot(p0['yrel'][mask(p0['T'])],p0['T'][mask(p0['T'])],'-',color=colors['ctrl'], label='O$_{KD}$=0')
    ax[1].plot(p1['yrel'][mask(p1['T'])], p1['T'][mask(p1['T'])], 'k-', zorder=5,  label='O$_{KD}$=0.28 $\mu$Ms$^{-1}$')
    ax[1].plot(p2['yrel'][mask(p2['T'])], p2['T'][mask(p2['T'])], 'b-', zorder=5,  label='O$_{KD}$=0.84 $\mu$Ms$^{-1}$')
    ax[1].plot(p3['yrel'][mask(p3['T'])], p3['T'][mask(p3['T'])], 'k-.')
    ax[1].plot(p4['yrel'][mask(p4['T'])], p4['T'][mask(p4['T'])], 'b-.')

    # Adjust axes
    ax[0].set_visible(False)
    pu.adjust_spines(ax[1], ['left','bottom'])
    ax[1].set(xlim = (1.4,1.9), xticks=np.arange(1.4,1.95,0.1),
          ylim = (3.5,15), yticks = np.arange(4,15.1,3),
          xlabel='Glu ($\mu$M)',
          ylabel='Oscillation Period (s)')

    # Add Legend
    ax[1].legend(loc='upper right', facecolor='w', edgecolor='none')

    # Save figure
    plt.savefig('../Figures/f4_pkc_osc_2.eps', dpi=600)

    plt.show()

#-----------------------------------------------------------------------------------------------------------------------
# Effective plotting
#-----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    gc.collect()
    chi_model()
    gchi_model()
    gchi_threshold(sim=True)
    # gchi_threshold(sim=False)
    pkc_effect(sim=True)
    # pkc_effect(sim=False)
    gc.collect()
