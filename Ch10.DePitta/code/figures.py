from __future__ import division
import numpy as np
import numpy.ma as ma  # Masked arrays
import copy as cp
import sys, os
# Change
base_dir = '/Ch10.DePitta/code'
sys.path.append(os.path.join(os.path.expanduser('~'),base_dir))
import pycustommodules.general_utils as gu
import pycustommodules.solvers.solver_utils as su
import pycustommodules.save_utils as svu
import pycustommodules.spk_generators as spg
from pycustommodules.graphics_utils import plot_utils as pu
from pycustommodules.graphics_utils import custom_plots as cpt  # contains routines for raster plots

# Models
import gliotransmission_models as gtm

# Plotting
import matplotlib as mpl
import matplotlib.pylab as plt
import matplotlib.patches as patches
import matplotlib.offsetbox as osb

#-----------------------------------------------------------------------------------------------------------------------
# Some useful color specifications (from Tableaux)
#-----------------------------------------------------------------------------------------------------------------------
DarkOrange = pu.hexc((255,127,14))
DarkGreen  = pu.hexc((44,160,44))
DarkBlue   = pu.hexc((31,119,180))
Red        = pu.hexc((214,39,40))
DarkRed    = pu.hexc((177,3,24))
Purple     = pu.hexc((148,103,189))
DarkGray   = pu.hexc((127,127,127))
DarkYellow = pu.hexc((219,161,58))
DarkMagenta= pu.hexc((227,119,194))
LightOrange= pu.hexc((255,198,140))
LightBlue  = pu.hexc((175,216,255))
LightMagenta = pu.hexc((247,182,210))
LightGreen = pu.hexc((155,223,138))
MediumGreen= pu.hexc((103,191,92))


#-----------------------------------------------------------------------------------------------------------------------
# Some useful Lambdas
#-----------------------------------------------------------------------------------------------------------------------
# u0
u_0 = lambda x,u0_star,alpha : u0_star + (alpha-u0_star)*x
# Plasticity threshold
u_thr = lambda x : x/(x+1)
# U0 at steady state
u0_inf = lambda nua,omg,oma,ua,js,alpha,u0 : (oma*omg*u0 + (omg*u0 + alpha*js*oma)*ua*nua)/(oma*omg + (js*oma + omg)*ua*nua)
# Asymptote for GRE rate
u0_asymptote = lambda x,js,oma,omg,alpha : u_thr(x)+js*oma/omg*(u_thr(x)-alpha)
# Threshold GRE rate
nu_thr = lambda omg,oma,ua,js,alpha,ut,u0 : omg*oma*(ut-u0)/(js*oma*ua*(alpha-ut)-omg*ua*(ut -u0))
# ppr = lambda u0,omd,omf,dt : (1+(1-u0)*np.exp(-omf*dt))*(1-np.exp(-omd*dt)+(1-u0)*np.exp(-omd*dt))

#-----------------------------------------------------------------------------------------------------------------------
# Figure 2 / Model compartments
#-----------------------------------------------------------------------------------------------------------------------
def model_compartments(sim=True,data_dir='../data/',fig_dir='../Figures/',format='pdf',dpi=600):
    """
    Show how each model compartment work.

    Input:
    - sim      : Run simulations
    - data_dir : Directory whereto save data to produce figures
    - fig_dir  : Directory whereto save figures
    - format   : Format string
    - dpi      : dpi of saved figure

    Return : void type

    Maurizio De Pitta', Basque Center for Applied Mathematics, June 23, 2018.
    """
    # Load MPL default style
    plt.style.use('figures.mplstyle')

    # General figure specifications
    figsize = (6.,5.5)
    TopMargin = 0.98
    BottomMargin = 0.12
    LeftMargin = 0.13
    RightMargin = 0.98

    # Color specifications
    colors = {'xs': Red,
              'us': DarkGray,
              'rs': Purple,
              'ca': DarkYellow,
              'h' : DarkGray,
              'xg': 'black',
              'ra': DarkBlue,
              'ga': DarkOrange,
              'Gs': 'black'}

    # Model parameters and ranges
    tRange = {'exo': [0, 3.],
              'gtr': [0, 60.],
              'ecs': [0, 5.]}
    tTicks = {'exo': np.arange(0,1.1*tRange['exo'][1]),
              'gtr': np.arange(0,1.1*tRange['gtr'][1],20.),
              'ecs': np.arange(0,1.1*tRange['ecs'][1])}
    yRange = {'aps': [0.,1],
              'ux' : [-0.02,1.02],
              'r'  : [0.,0.4],
              'ca' : [0.,2.5],
              'xa' : [0.,1.02],
              'ra' : [0.,0.6],
              'ga' : [-5.,200],
              'Gs' : [-0.02,0.7]}
    yTicks = {'aps': [],
              'ux' : [0.,0.5,1],
              'r'  : [0.,0.2,0.4],
              'ca' : np.arange(0.,1.1*yRange['ca'][1],0.5),
              'xa' : [0.,0.5,1],
              'ra' : [0.,0.3,0.6],
              'ga' : [0.,100.,200],
              'Gs' : np.arange(0.,yRange['Gs'][1],0.2)}

    #-------------------------------------------------------------------------------------------------------------------
    # Synaptic compartment
    #-------------------------------------------------------------------------------------------------------------------
    if sim:
        print "1. Synaptic compartment"
        exo = gtm.gtrs(model='exo_spk',
                       u0=0.2, taud=0.5, tauf=1.,
                       ICs=[1.0,0.,0.0])
        nus = [7,2.5]
        aps_pre = 0.2+spg.input_spikes(1,1.5,nus[0],0,stimulus='periodic')[0]
        aps_pre += np.random.normal(0,0.2/nus[0],len(aps_pre))
        aps_pre = np.r_[aps_pre,2.0+spg.input_spikes(1,1.0,nus[1],0,stimulus='periodic')[0]]
        exo.stimulation(twin=[0., tRange['exo'][-1]], Nsyn=1, Nvar=2,
                        stimulus='fixed', spikes_pre=aps_pre)
        exo.simulate(reconstruct=True,recompile=0)
        # Save data
        svu.savedata([exo.sol,exo.pars],data_dir+'compartment_syn.pkl')

        print "2. Calcium-dependent gliotransmitter exocytosis compartment"
        options = su.solver_opts(t0=0.0, tfin=tRange['gtr'][-1], dt=1e-3, atol=1e-8, rtol=1e-8, method="gsl_rk2imp")
        astro = gtm.gtrs(model='lra_exo',ip3=0.4,noise=True,ICs=[0.05,0.9,1.0])
        astro.simulate(algparams=options,reconstruct=True,recompile=0)
        # Save data
        svu.savedata([astro.sol,astro.pars],data_dir+'compartment_gtr.pkl')

        print "3. Extracellular gliotransmitter time course and receptor activation"
        options = su.solver_opts(t0=0.0, tfin=tRange['ecs'][-1], dt=1e-4, atol=1e-8, rtol=1e-8, method="gsl_rk2imp")
        ecs = gtm.gtrs(model='asn_spk', Ou=1e2, Ku=10., taue=1./30., Og=0.1)
        ecs.stimulation(rate=0, rate_gtr=1, twin_gtr=[0.,5.],stimulus_gtr='poisson')
        ecs.simulate(algparams=options,reconstruct=False,recompile=0)
        # Save data
        svu.savedata([ecs.sol,ecs.pars],data_dir+'compartment_ecs.pkl')

    # -------------------------------------------------------------------------------------------------------------------
    # Effective plotting
    # -------------------------------------------------------------------------------------------------------------------
    # SYN compartment
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=figsize, sharex=True,
                           gridspec_kw={'height_ratios': [1, 1, 1],
                                        'hspace': 0.1,
                                        'top': TopMargin, 'bottom': BottomMargin,
                                        'left': LeftMargin, 'right': RightMargin})
    # Load data
    sol,p = svu.loaddata(data_dir+'compartment_syn.pkl')

    ax[0].stem(sol['spk'],np.ones(len(sol['spk'])),linefmt='k-',markerfmt=" ",basefmt=" ")
    ax[1].plot(sol['ts'][0],sol['xs'][0],'-',color=colors['xs'])
    ax[1].plot(sol['ts'][0], sol['us'][0], '-', color=colors['us'])
    ml, sl, bl = ax[2].stem(sol['spk'],sol['r'],linefmt='-',markerfmt=" ",basefmt=" ")
    plt.setp(sl,color=colors['rs'])

    # Adjust spines
    pu.adjust_spines(ax[0], [])
    pu.adjust_spines(ax[1], ['left'])
    pu.adjust_spines(ax[2], ['left','bottom'])
    # Adjust axes
    ax[0].set(xlim=tRange['exo'],xticks=[],ylim=yRange['aps'],yticks=[],ylabel='APs')
    ax[1].set(xlim=tRange['exo'],xticks=[],
              ylim=yRange['ux'],yticks=yTicks['ux'])
    ax[2].set(xlim=tRange['exo'],xticks=tTicks['exo'],
              xticklabels=[str(t) for t in tTicks['exo']],xlabel='Time (s)',
              ylim=yRange['r'],yticks=yTicks['r'],ylabel='$r_S$')
    # Adjust y-labels
    pu.adjust_ylabels(ax, x_offset=-0.09)
    # Add a two-color y-label to the intermediate axis (...matplotlib is notoriously pretty straight forward on this sort of things...)
    x_box = osb.TextArea("$x_S$", textprops=dict(color=colors['xs'], size=16, rotation=90,ha='left',va='bottom'))
    c_box = osb.TextArea(",", textprops=dict(color='k', size=16, rotation=90,ha='left',va='bottom'))
    u_box = osb.TextArea("$u_S$", textprops=dict(color=colors['us'], size=16, rotation=90,ha='left',va='bottom'))
    ybox = osb.VPacker(children=[x_box, c_box, u_box], align='bottom', pad=0, sep=0)
    anchored_ybox = osb.AnchoredOffsetbox(loc=8, child=ybox, pad=0., frameon=False, bbox_to_anchor=(-0.1,0.2),
                                          bbox_transform=ax[1].transAxes, borderpad=0.)
    ax[1].add_artist(anchored_ybox)
    # Save figure
    plt.savefig(fig_dir+'compartment_exo.'+format, dpi=dpi,transparent=True)

    # GTR compartment
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=figsize, sharex=True,
                           gridspec_kw={'height_ratios': [1, 1, 1],
                                        'hspace': 0.1,
                                        'top': TopMargin, 'bottom': BottomMargin,
                                        'left': LeftMargin, 'right': RightMargin})
    # Load data
    sol, p = svu.loaddata(data_dir + 'compartment_gtr.pkl')
    # You must take the first element regardless
    # Calcium dynamics
    ax[0].plot(sol['t'], sol['ca'],'-',color=colors['ca'])
    ax[0].plot(sol['t'][[0,-1]], p['cthr']*np.ones(2),'--', color=DarkRed)
    # x_a
    ax[1].plot(sol['tg'][0], sol['xg'][0],'-', color=colors['xg'])
    # r_a
    ml, sl, bl = ax[2].stem(sol['gre'], sol['ra'], linefmt='-',markerfmt=" ",basefmt=" ")
    plt.setp(sl, color=colors['ra'])

    # Adjust spines
    for i in xrange(2): pu.adjust_spines(ax[i], ['left'])
    pu.adjust_spines(ax[2], ['left','bottom'])
    # Adjust axes
    ax[0].set(xlim=tRange['gtr'],xticks=[],
           ylim=yRange['ca'],yticks=yTicks['ca'],ylabel=r'Ca$^{2+}$ ($\mu$M)')
    ax[1].set(xlim=tRange['gtr'],xticks=[],
           ylim=yRange['xa'],yticks=yTicks['xa'],ylabel='$x_A$')
    ax[2].set(xlim=tRange['gtr'],xticks=tTicks['gtr'],
              xticklabels=[str(t) for t in tTicks['gtr']],xlabel='Time (s)',
              ylim=yRange['ra'],yticks=yTicks['ra'],ylabel='$r_A$')
    # Adjust y-labels
    pu.adjust_ylabels(ax, x_offset=-0.09)
    # Save figure
    plt.savefig(fig_dir+'compartment_gtr.'+format, dpi=dpi,transparent=True)

    # ECS compartment
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=figsize, sharex=True,
                           gridspec_kw={'height_ratios': [1, 1, 1],
                                        'hspace': 0.1,
                                        'top': TopMargin, 'bottom': BottomMargin,
                                        'left': LeftMargin, 'right': RightMargin})
    # Load data
    sol, p = svu.loaddata(data_dir + 'compartment_ecs.pkl')
    # You must take the first element regardless
    # r_a
    ml, sl, bl = ax[0].stem(sol['gre'], sol['ra'], linefmt='-', markerfmt=" ", basefmt=" ")
    plt.setp(sl, color=colors['ra'])
    # g_a
    ax[1].plot(sol['t'], sol['ga'], '-', color=colors['ga'])
    # Gamma_s
    ax[2].plot(sol['t'], sol['gammas'], '-', color=colors['Gs'])

    # Adjust spines
    for i in xrange(2): pu.adjust_spines(ax[i], ['left'])
    pu.adjust_spines(ax[2], ['left','bottom'])
    # Adjust axes
    ax[0].set(xlim=tRange['ecs'],xticks=[],
           ylim=yRange['ra'],yticks=yTicks['ra'],ylabel='$r_A$')
    ax[1].set(xlim=tRange['ecs'],xticks=[],
           ylim=yRange['ga'],yticks=yTicks['ga'],ylabel='$G_A$ ($\mu$M)')
    ax[2].set(xlim=tRange['ecs'],xticks=tTicks['ecs'],
              xticklabels=[str(t) for t in tTicks['ecs']],xlabel='Time (s)',
              ylim=yRange['Gs'],yticks=yTicks['Gs'],ylabel=r'$\Gamma_S$')
    # Adjust y-labels
    pu.adjust_ylabels(ax, x_offset=-0.09)
    # Save figure
    plt.savefig(fig_dir+'compartment_ecs.'+format, dpi=dpi,transparent=True)

    # Show plots
    plt.show()

#-----------------------------------------------------------------------------------------------------------------------
# Figure 3 / Full model
#-----------------------------------------------------------------------------------------------------------------------
def full_model(sim=True,data_dir='../data/',fig_dir='../Figures/',format='pdf',dpi=600):
    """
    Full model exemplification.

    Input:
    - sim      : Run simulations
    - data_dir : Directory whereto save data to produce figures
    - fig_dir  : Directory whereto save figures
    - format   : Format string
    - dpi      : dpi of saved figure

    Return : void type

    Maurizio De Pitta', Basque Center for Applied Mathematics, June 23, 2018.
    """
    # Load MPL default style
    plt.style.use('figures.mplstyle')

    # General figure specifications
    figsize = (6.5,6.0)
    TopMargin = 0.98
    BottomMargin = 0.11
    LeftMargin = 0.12
    RightMargin = 0.78

    # Color specifications
    colors = {'xs': Red,
              'us': DarkGray,
              'rs': Purple,
              'ca': DarkYellow,
              'h' : DarkGray,
              'xg': 'black',
              'ra': DarkBlue,
              'ga': DarkOrange,
              'Gs': 'black'}

    # Model parameters and ranges
    tRange = [0., 20.]
    tTicks = np.arange(tRange[0],1.1*tRange[1],5.)
    yRange = {'ca' : [0.,2.5],
              'xa' : [0.,1.02],
              'ra' : [0.,0.6],
              'ga' : [-5.,200],
              'Gs' : [-0.02,0.7],
              'u0' : [0.,1.]}
    yTicks = {'ca' : np.arange(0.,1.1*yRange['ca'][1],0.5),
              'xa' : [0.,0.5,1],
              'ra' : [0.,0.3,0.6],
              'ga' : [0.,100.,200],
              'Gs' : np.arange(0.,yRange['Gs'][1],0.2),
              'u0' : [0,0.5,1.]}

    #-------------------------------------------------------------------------------------------------------------------
    # Synaptic compartment
    #-------------------------------------------------------------------------------------------------------------------
    if sim:
        print "1. Full model simulation"
        options = su.solver_opts(t0=tRange[0], tfin=tRange[1], dt=1e-4, atol=1e-8, rtol=1e-8, method="gsl_rk2imp")
        asn = gtm.gtrs(model='asn_spk', noise=True, ip3=0.4,
                       Ou=1e2, Ku=10., taue=1./30., Og=0.1)
        asn.stimulation()
        asn.simulate(algparams=options,reconstruct=False,recompile=0)
        # Save data
        svu.savedata([asn.sol,asn.pars],data_dir+'full_model.pkl')

    # -------------------------------------------------------------------------------------------------------------------
    # Effective plotting
    # -------------------------------------------------------------------------------------------------------------------
    fig, ax = plt.subplots(nrows=4, ncols=1, figsize=figsize, sharex=True,
                           gridspec_kw={'height_ratios': [1, 1, 1, 2],
                                        'hspace': 0.15,
                                        'top': TopMargin, 'bottom': BottomMargin,
                                        'left': LeftMargin, 'right': RightMargin})
    # Load data
    sol,p = svu.loaddata(data_dir+'full_model.pkl')

    # # Calcium dynamics
    ax[0].plot(sol['t'], sol['ca'],'-',color=colors['ca'])
    ax[0].plot(sol['t'][[0,-1]], p['cthr']*np.ones(2),'--', color=DarkRed)
    # G_a
    ax[1].plot(sol['t'], sol['ga'],'-',color=colors['ga'])
    # Gamma_S
    ax[2].plot(sol['t'], sol['gammas'], '-', color=colors['Gs'])

    # Generate background
    npts = 100
    time = np.linspace(tRange[0],tRange[1],npts)
    u0 = np.linspace(yRange['u0'][0],yRange['u0'][1],npts)
    T,U0 = np.meshgrid(time,u0,indexing='xy')
    # Generate divergent colormap
    cmap = pu.CustomColormap([0,1], yTicks['u0'], [LightMagenta,'white',MediumGreen], cmap_name='custom')
    # Plot background
    pm = ax[3].pcolormesh(T, U0, U0, cmap=cmap, vmin=yRange['u0'][0], vmax=yRange['u0'][1], shading='gouraud')
    ax[3].plot(sol['t'],-0.02+u_0(sol['gammas'],0.5,0),'-',color=DarkMagenta)
    ax[3].plot(sol['t'],0.02+u_0(sol['gammas'],0.5,1), '-', color=DarkGreen)
    ax[3].plot(sol['t'],u_0(sol['gammas'], 0.5, 0.5), '--', color=DarkGray)

    ax_pos = ax[3].get_position().get_points()
    cax = fig.add_axes(np.r_[0.83, ax_pos[0][1], 0.05, 0.3175])
    cbar = fig.colorbar(pm, cax=cax)

    # Adjust spines
    for i in xrange(3): pu.adjust_spines(ax[i], ['left'])
    pu.adjust_spines(ax[3], ['left','bottom'])
    # Adjust axes
    ax[0].set(xlim=tRange,xticks=[],
           ylim=yRange['ca'],yticks=yTicks['ca'],ylabel=r'Ca$^{2+}$ ($\mu$M)')
    ax[1].set(xlim=tRange,xticks=[],
           ylim=yRange['ga'],yticks=yTicks['ga'],ylabel='$G_A$ ($\mu$M)')
    ax[2].set(xlim=tRange,xticks=[],
              ylim=yRange['Gs'],yticks=yTicks['Gs'],ylabel=r'$\Gamma_S$')
    ax[3].set(xlim=tRange,xticks=tTicks,
              xticklabels=[str(t) for t in tTicks],xlabel='Time (s)',
              ylim=yRange['u0'],yticks=yTicks['u0'],ylabel='$u_0$')
    # Adjust y-labels
    pu.adjust_ylabels(ax, x_offset=-0.1)
    # Adjust colorbar
    cbar.ax.set_ylabel(r'Gt. Type, $\alpha$')

    # Save figure
    plt.savefig(fig_dir+'full_model.'+format, dpi=dpi,transparent=True)

    plt.show()

#-----------------------------------------------------------------------------------------------------------------------
# Figure 4 / Synaptic Filtering
#-----------------------------------------------------------------------------------------------------------------------
def step_simulation(pars,rate,twin,N_syn,gtr=False,dt=1e-3,gres=[]):
    """
    Simulate the tripartite synapse model, for step variations of synaptic rate.

    Maurizio De Pitta', Basque Center of Applied Mathematics, June 21, 2018.
    """

    rate = np.asarray(rate,dtype=float)
    twin = np.asarray(twin,dtype=float)
    assert rate.size==twin.size, "rate and twin must be of same length"
    tfin = np.cumsum(np.r_[0,twin])
    spk_sol = []
    spk_idx = []
    spk_tsp = []
    sol_gms = []
    sol_tim = []
    time = np.arange(0.,tfin[-1],dt)
    if not gtr:
        pars['ICs'] = [1.0,0.,0.0]
        exo = gtm.gtrs(model='exo_spk',**pars)
        for i,tf in enumerate(twin):
            if i>0: exo.ICs = exo.sol['LCs'] # Start from last point of previous simulation
            exo.stimulation(rate=rate[i],twin=[0., tf], Nsyn=N_syn, Nvar=2,stimulus='poisson')
            exo.simulate(reconstruct=False,recompile=0)
            spk_tsp = np.r_[spk_tsp,tfin[i]+exo.sol['spk']]
            spk_sol = np.r_[spk_sol,exo.sol['y']]
            spk_idx = np.r_[spk_idx,exo.sol['is']]
    else:
        pars['ICr'] = [1.0,0.0,0.0,1.0]
        asn = gtm.gtrs(model='asn_spk',**pars)
        # Auxiliary parameters to handle stimulation of GTR
        gres = np.asarray(gres)
        for i,tf in enumerate(twin):
            if i>0:
                asn.ICs = np.r_[asn.sol['ga'][-1],asn.sol['gammas'][-1]]
                asn.ICr = asn.sol['LCr'] # Start from last point of previous simulation
            gres_ = gres-tfin[i]
            gres_ = gres_[(gres_>=0)&(gres_<=tf)]
            rate_gtr = 0.
            twin_gtr = [0., tf]
            if len(gres_)>0:
                stimulus_gtr = 'fixed'
            else:
                stimulus_gtr = 'poisson'
            asn.stimulation(rate=rate[i],twin=[0., tf], Nsyn=N_syn, Nvar=2,stimulus='poisson',
                            rate_gtr=rate_gtr,twin_gtr=twin_gtr,stimulus_gtr=stimulus_gtr,pre_gtr=gres_)
            options = su.solver_opts(t0=0.0, tfin=tf, dt=1e-4, atol=1e-8, rtol=1e-8, method="gsl_rk2imp")
            asn.simulate(algparams=options,reconstruct=False,recompile=0)
            print "Screening gamma_S: step: ",i,"-> value: ",max(asn.sol['gammas'])," --> u0: ",u_0(max(asn.sol['gammas']),asn.pars['u0'],asn.pars['alpha'])
            spk_tsp = np.r_[spk_tsp,tfin[i]+asn.sol['spk']]
            spk_sol = np.r_[spk_sol,asn.sol['y']]
            spk_idx = np.r_[spk_idx,asn.sol['is']]
            # Continuous solution
            sol_tim = np.r_[sol_tim,tfin[i]+asn.sol['t']]
            sol_gms = np.r_[sol_gms,asn.sol['gammas']]

    # Reconstruct y solution
    sol = []
    for i in xrange(int(max(spk_idx))+1):
        index = spk_idx==i
        tmp = gtm.reconstruct_solution(spk_tsp[index], spk_sol[index], [], [0.,tfin[-1]], 0.0, pars['taui'], 'y', dt=dt)
        keep = np.searchsorted(tmp[0],time)
        sol.append(tmp[1][keep])
    sol = np.mean(np.vstack(sol),axis=0)

    if not gtr:
        return time,sol
    else:
        keep = np.searchsorted(sol_tim, time)
        gammas = sol_gms[keep]
        return time,sol,gammas

def syn_filtering(sim=True,data_dir='../data/',fig_dir='../../Figures/',format='pdf',dpi=600):
    # Load MPL default style
    plt.style.use('figures.mplstyle')

    # Color specifications
    colors = {'y' : DarkGray,
              'ra': DarkBlue,
              'Gs': 'black'}

    # General figure specifications
    figsize = (6.0,5.5)
    TopMargin = 0.98
    BottomMargin = 0.12
    LeftMargin = 0.12
    RightMargin = 0.96

    # Stimulus specification
    rate = [1.5,10.,20.]
    twin = [10.,10.,10.]
    gres = np.arange(5.,26.,10.)

    if sim:
        time, sol, gammas, p = {}, {}, {}, {}
        Nsyn = 1000

        print "1. Simulation facilitating synapse"
        pars = gtm.exocytosis_parameters(u0=0.05,taud=0.3,tauf=1.0)
        p['fac'] = pars
        time['fac'],sol['fac'] = step_simulation(pars, rate, twin, Nsyn, gtr=False,dt=1e-3)
        # Transient saving of data
        svu.savedata([time,sol,gammas,p],'syn_filter.pkl')

        print "2. Simulation facilitating synapse + release-increasing gliotransmission"
        pars = gtm.asn_parameters(Ou=1e2, Ku=0.1, taue=1./30., Gtot=150,
                                  Og=0.1,taug=20.,alpha=1.0, **p['fac'])
        p['fac_ri'] = pars
        time['fac_ri'],sol['fac_ri'],gammas['fac_ri'] = step_simulation(pars, rate, twin, Nsyn, dt=1e-3,gtr=True, gres=gres)
        # Transient saving of data
        svu.savedata([time,sol,gammas,p],'syn_filter.pkl')

        print "3. Simulation depressing synapse"
        pars = gtm.exocytosis_parameters(u0=0.5,taud=0.5,tauf=0.3)
        p['dep'] = pars
        time['dep'],sol['dep'] = step_simulation(pars, rate, twin, Nsyn, gtr=False,dt=1e-3)
        # Transient saving of data
        svu.savedata([time,sol,gammas,p],'syn_filter.pkl')

        # time,sol,p = svu.loaddata('syn_filter.pkl')
        print "4. Simulation depressing synapse + release-decreasing gliotransmission"
        pars = gtm.asn_parameters(Ou=1e2, Ku=0.1, taue=1./30., Gtot=150,
                                  Og=0.7,taug=20.,alpha=0.0, **p['dep'])
        p['dep_rd'] = pars
        time['dep_rd'],sol['dep_rd'],gammas['dep_rd'] = step_simulation(pars, rate, twin, Nsyn, dt=1e-3,gtr=True, gres=gres)
        # Final saving of data
        svu.savedata([time,sol,gammas,p],'syn_filter.pkl')
    # -------------------------------------------------------------------------------------------------------------------
    # Effective plotting
    # -------------------------------------------------------------------------------------------------------------------
    # Load data
    time, sol, gammas, p = svu.loaddata('syn_filter.pkl')
    # Add zero point in all plots
    t0 = -1
    tstep = np.r_[-1,0.,np.cumsum(twin)]
    rate = np.r_[0.,rate,rate[-1]]

    # Ranges
    tRange = [-1.,np.sum(twin)]
    tTicks = np.arange(0.,31.,10.)
    yRange = {'nu' : [-1.,21.],
              'y_fac'  : [-0.1,2.7],
              'y_dep'  : [-0.1,2.5],
              'Gs' : [-0.02,1.01]}
    yTicks = {'nu' : np.arange(0.,21.,10.),
              'y_fac'  : np.arange(0.,2.8,.5),
              'y_dep'  : np.arange(0.,2.8,.5),
              'Gs' : [0.,0.5,1.]}

    #--------------------------------------------------------------------
    # Fac / Without gliotransmission
    #--------------------------------------------------------------------
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=figsize, sharex=True,
                           gridspec_kw={'height_ratios': [1, 2],
                                        'hspace': 0.15,
                                        'top': TopMargin, 'bottom': BottomMargin,
                                        'left': LeftMargin, 'right': RightMargin})
    # Rate
    ax[0].step(tstep,rate,'k-',where='post')
    # PSCs
    ax[1].plot(time['fac'],sol['fac'],'-',color='k')
    ax[1].plot(np.r_[t0,time['fac'][0]],sol['fac'][[0,0]],'-',color='k')
    # Adjust spines
    pu.adjust_spines(ax[0], ['left'])
    pu.adjust_spines(ax[1], ['left','bottom'])
    # Adjust axes
    ax[0].set(xlim=tRange,xticks=[],
              ylim=yRange['nu'],yticks=yTicks['nu'],ylabel=r'$\nu_S(t)$')
    ax[1].set(xlim=tRange,xticks=tTicks,
              xticklabels=[str(t) for t in tTicks],xlabel='Time (s)',
              ylim=yRange['y_fac'],yticks=yTicks['y_fac'],ylabel='$I(t)$ (pA)')
    # Adjust y-labels
    pu.adjust_ylabels(ax, x_offset=-0.08)
    # Save figure
    plt.savefig(fig_dir+'filtering_fac.'+format, dpi=dpi,transparent=True)

    #--------------------------------------------------------------------
    # Dep / Without gliotransmission
    #--------------------------------------------------------------------
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=figsize, sharex=True,
                           gridspec_kw={'height_ratios': [1, 2],
                                        'hspace': 0.15,
                                        'top': TopMargin, 'bottom': BottomMargin,
                                        'left': LeftMargin, 'right': RightMargin})
    # Rate
    ax[0].step(tstep,rate,'k-',where='post')
    # PSCs
    ax[1].plot(time['dep'],sol['dep'],'-',color='k')
    ax[1].plot(np.r_[t0,time['dep'][0]],sol['dep'][[0,0]],'-',color='k')
    # Adjust spines
    pu.adjust_spines(ax[0], ['left'])
    pu.adjust_spines(ax[1], ['left','bottom'])
    # Adjust axes
    ax[0].set(xlim=tRange,xticks=[],
              ylim=yRange['nu'],yticks=yTicks['nu'],ylabel=r'$\nu_S(t)$')
    ax[1].set(xlim=tRange,xticks=tTicks,
              xticklabels=[str(t) for t in tTicks],xlabel='Time (s)',
              ylim=yRange['y_dep'],yticks=yTicks['y_dep'],ylabel='$I(t)$ (pA)')
    # Adjust y-labels
    pu.adjust_ylabels(ax, x_offset=-0.08)
    # Save figure
    plt.savefig(fig_dir+'filtering_dep.'+format, dpi=dpi,transparent=True)

    #--------------------------------------------------------------------
    # Fac / with gliotransmission
    #--------------------------------------------------------------------
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=figsize, sharex=True,
                           gridspec_kw={'height_ratios': [.3, .7, 2],
                                        'hspace': 0.15,
                                        'top': TopMargin, 'bottom': BottomMargin,
                                        'left': LeftMargin, 'right': RightMargin})
    # GRE
    ml, sl, bl = ax[0].stem(gres, np.ones(len(gres)), linefmt='-', markerfmt=" ", basefmt=" ")
    plt.setp(sl, color=colors['ra'])
    # Gammas
    ax[1].plot(time['fac_ri'],gammas['fac_ri'],'-',color=colors['Gs'])
    ax[1].plot(np.r_[t0,time['fac_ri'][0]],gammas['fac_ri'][[0,1]],color=colors['Gs'])
    # PSCs
    ax[2].plot(time['fac_ri'],sol['fac_ri'],'-',color=colors['y'])
    ax[2].plot(np.r_[t0,time['fac_ri'][0]],sol['fac_ri'][[0,0]],'-',color=colors['y'])
    # Adjust spines
    pu.adjust_spines(ax[0], [])
    pu.adjust_spines(ax[1], ['left'])
    pu.adjust_spines(ax[2], ['left','bottom'])
    # Adjust axes
    ax[0].set(xlim=tRange,xticks=[],
              ylim=[0.,1.],yticks=[],ylabel='GREs')
    ax[1].set(xlim=tRange,xticks=[],
              ylim=yRange['Gs'],yticks=yTicks['Gs'],ylabel=r'$\Gamma_S$')
    ax[2].set(xlim=tRange,xticks=tTicks,
              xticklabels=[str(t) for t in tTicks],xlabel='Time (s)',
              ylim=yRange['y_fac'],yticks=yTicks['y_fac'],ylabel='$I(t)$ (pA)')
    # Adjust y-labels
    pu.adjust_ylabels(ax, x_offset=-0.08)
    # Save figure
    plt.savefig(fig_dir+'filtering_fac_ri.'+format, dpi=dpi,transparent=True)

    #--------------------------------------------------------------------
    # Dep / with gliotransmission
    #--------------------------------------------------------------------
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=figsize, sharex=True,
                           gridspec_kw={'height_ratios': [.3, .7, 2],
                                        'hspace': 0.15,
                                        'top': TopMargin, 'bottom': BottomMargin,
                                        'left': LeftMargin, 'right': RightMargin})
    # GRE
    ml, sl, bl = ax[0].stem(gres, np.ones(len(gres)), linefmt='-', markerfmt=" ", basefmt=" ")
    plt.setp(sl, color=colors['ra'])
    # Gammas
    ax[1].plot(time['dep_rd'],gammas['dep_rd'],'-',color=colors['Gs'])
    ax[1].plot(np.r_[t0,time['dep_rd'][0]],gammas['dep_rd'][[0,1]],color=colors['Gs'])
    # PSCs
    ax[2].plot(time['dep_rd'],sol['dep_rd'],'-',color=colors['y'])
    ax[2].plot(np.r_[t0,time['dep_rd'][0]],sol['dep_rd'][[0,0]],'-',color=colors['y'])
    # Adjust spines
    pu.adjust_spines(ax[0], [])
    pu.adjust_spines(ax[1], ['left'])
    pu.adjust_spines(ax[2], ['left','bottom'])
    # Adjust axes
    ax[0].set(xlim=tRange,xticks=[],
              ylim=[0.,1.],yticks=[],ylabel='GREs')
    ax[1].set(xlim=tRange,xticks=[],
              ylim=yRange['Gs'],yticks=yTicks['Gs'],ylabel=r'$\Gamma_S$')
    ax[2].set(xlim=tRange,xticks=tTicks,
              xticklabels=[str(t) for t in tTicks],xlabel='Time (s)',
              ylim=yRange['y_dep'],yticks=yTicks['y_dep'],ylabel='$I(t)$ (pA)')
    # Adjust y-labels
    pu.adjust_ylabels(ax, x_offset=-0.08)
    # Save figure
    plt.savefig(fig_dir+'filtering_dep_rd.'+format, dpi=dpi,transparent=True)

    plt.show()

#-----------------------------------------------------------------------------------------------------------------------
# Figure 5 / Paired-pulse plasticity
#-----------------------------------------------------------------------------------------------------------------------
def ppp_resolve(u,r,syn_index):
    """
    Provide a 3-valued vector for all t_spk: 1: PPF; -1: PPD; 0: RFD

    Maurizio De Pitta', Basque Center for Applied Mathematics, June 27, 2018.
    """
    N = int(max(syn_index)+1)
    ip = np.zeros(np.size(syn_index))
    for i in xrange(N):
        idx = np.where(syn_index==i)[0]
        fac,dep,rec = cpt.stp_indexes(u[idx],r[idx])
        ip[idx[fac]] = 1
        ip[idx[dep]] = -1
        # The unassigned ones are r
    return ip

def ppp_histogram(twin,t_spk,u,r,syn_index,binw=20e-3,average=True):
    """
    Provide rates of different PPP events during stimulation. Uses binning by binw.

    Maurizio De Pitta', Basque Center of Applied Mathematics, June 27, 2018.
    """

    ip = ppp_resolve(u,r,syn_index)
    t_bins = np.arange(twin[0],twin[1],binw)
    v_bins = np.arange(-1,3)
    # Find which parts of t_spk are within the bins: this works because t_spk are sorted and so are the bins
    idx,bins = np.histogram(t_spk,t_bins,range=twin)
    fac,dep,rfd = np.zeros(len(idx)),np.zeros(len(idx)),np.zeros(len(idx))
    i = 0
    for j,v in enumerate(idx):
        ppp,_ = np.histogram(ip[i:i+v],v_bins)
        i += v
        dep[j] = ppp[0]
        rfd[j] = ppp[1]
        fac[j] = ppp[2]
    if average:
        # Average over all synapses
        Nsyn = int(max(syn_index)+1)
        dep /= Nsyn
        fac /= Nsyn
        rfd /= Nsyn
    return bins[:-1],dep,fac,rfd

def pp_plasticity(sim=True,data_dir='../data',fig_dir='../../Figures/',format='pdf',dpi=600):
    """
    Generate the raster plots with plasticity color coding.

    Maurizio De Pitta', Basque Center of Applied Mathematics, June 27, 2018.
    """
    # Load MPL default style
    plt.style.use('figures.mplstyle')

    # General figure specifications
    figsize = (6.5,6.0)
    TopMargin = 0.98
    BottomMargin = 0.11
    LeftMargin = 0.12
    RightMargin = 0.96

    # Color specifications
    colors = {'dep': DarkOrange,
              'fac': DarkBlue,
              'rec': 'black',
              'ra' : 'black'}

    # Model parameters and ranges
    Nsyn = 1000
    tRange = [0., 60.]
    tTicks = np.arange(tRange[0],1.1*tRange[1],10.)
    yRange = {'si' : [-0.5,19.5],
              'pp' : [0.,3.0]}
    yTicks = {'si' : np.arange(0.,20,5),
              'pp' : np.arange(0.,3.1,1.)}

    sol,pars = {},{}
    gres = [2.,15.,17.]
    if sim:
        print "1. Plasticity of FAC synapse"
        pars['fac'] = gtm.exocytosis_parameters(u0=0.05,taud=0.3,tauf=1.0)
        exo = gtm.gtrs(model='exo_spk',**pars['fac'])
        exo.stimulation(rate=2.5,twin=tRange,Nsyn=Nsyn,Nvar=2)
        exo.simulate(reconstruct=False,recompile=0)
        # Save data
        sol['fac'] = cp.copy(exo.sol)
        svu.savedata([sol,pars],data_dir+'pp_plasticity.pkl')

        print "2. Plasticity of FAC synapse + release-increasing gliotransmission"
        options = su.solver_opts(t0=tRange[0], tfin=tRange[1], dt=1e-4, atol=1e-8, rtol=1e-8, method="gsl_rk2imp")
        pars['fac_ri'] = gtm.asn_parameters(Ou=1e2, Ku=0.1, taue=1./30., Gtot=150, Og=0.1, taug=20., alpha=1.0, **pars['fac'])
        tsn = gtm.gtrs(model='asn_spk',**pars['fac_ri'])
        tsn.stimulation(rate=0.,twin=tRange,Nsyn=Nsyn,Nvar=2,
                        stimulus='fixed', spikes_pre=np.vstack((sol['fac']['spk'], sol['fac']['is'])),
                        stimulus_gtr='fixed',pre_gtr=gres)
        tsn.simulate(algparams=options,reconstruct=False,recompile=0)
        # Save data
        sol['fac_ri'] = cp.copy(tsn.sol)
        svu.savedata([sol, pars], data_dir + 'pp_plasticity.pkl')

        print "3. Plasticity of DEP synapse"
        pars['dep'] = gtm.exocytosis_parameters(u0=0.5,taud=0.5,tauf=0.3)
        exo = gtm.gtrs(model='exo_spk',**pars['dep'])
        exo.stimulation(rate=2.5,twin=tRange,Nsyn=Nsyn,Nvar=2)
        exo.simulate(reconstruct=False,recompile=0)
        # Save data
        sol['dep'] = cp.copy(exo.sol)
        svu.savedata([sol,pars],data_dir+'pp_plasticity.pkl')

        print "4. Plasticity of DEP synapse + release-decreasing gliotransmission"
        options = su.solver_opts(t0=tRange[0], tfin=tRange[1], dt=1e-4, atol=1e-8, rtol=1e-8, method="gsl_rk2imp")
        pars['dep_rd'] = gtm.asn_parameters(Ou=1e2, Ku=0.1, taue=1./30., Gtot=150, Og=0.7,taug=20.,alpha=0.0, **pars['dep'])
        tsn = gtm.gtrs(model='asn_spk',**pars['dep_rd'])
        tsn.stimulation(rate=0.,twin=tRange,Nsyn=Nsyn,Nvar=2,
                        stimulus = 'fixed', spikes_pre=np.vstack((sol['dep']['spk'],sol['dep']['is'])),
                        stimulus_gtr='fixed',pre_gtr=gres)
        tsn.simulate(algparams=options,reconstruct=False,recompile=0)
        # Save data
        sol['dep_rd'] = cp.copy(tsn.sol)
        svu.savedata([sol, pars], data_dir + 'pp_plasticity.pkl')

    # -------------------------------------------------------------------------------------------------------------------
    # Effective plotting
    # -------------------------------------------------------------------------------------------------------------------
    # Load data
    sol, pars = svu.loaddata(data_dir+'pp_plasticity.pkl')
    # time bins for solutions
    binw = 200e-3
    #--------------------------------------------------------------------
    # Fac / Without gliotransmission
    #--------------------------------------------------------------------
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=figsize, sharex=True,
                           gridspec_kw={'height_ratios': [2, 1],
                                        'hspace': 0.15,
                                        'top': TopMargin, 'bottom': BottomMargin,
                                        'left': LeftMargin, 'right': RightMargin})

    cpt.stp_raster(sol['fac']['spk'], sol['fac']['u'], sol['fac']['r'], sol['fac']['is'],
                   N=20, ax=ax[0], color=colors)

    # Plot histogram
    drf = [1,0,2]
    tbins,dep,fac,rec = ppp_histogram(tRange,sol['fac']['spk'], sol['fac']['u'], sol['fac']['r'], sol['fac']['is'], binw=binw, average=True)
    dep /= binw
    fac /= binw
    rec /= binw
    ax[1].step(tbins,dep,where='pre',color=colors['dep'],zorder=drf[0])
    ax[1].step(tbins,rec,where='pre',color=colors['rec'],zorder=drf[1])
    ax[1].step(tbins,fac,where='pre',color=colors['fac'],zorder=drf[2])
    # Adjust spines
    pu.adjust_spines(ax[0], ['left'])
    pu.adjust_spines(ax[1], ['left','bottom'])
    # Adjust axes
    ax[0].set(xlim=tRange,xticks=[],
              ylim=yRange['si'],yticks=yTicks['si'],ylabel='Synapse Index')
    ax[1].set(xlim=tRange,xticks=tTicks,
              xticklabels=[str(t) for t in tTicks],xlabel='Time (s)',
              ylim=yRange['pp'],yticks=yTicks['pp'],ylabel='PPP (pairs/s)')
    # Adjust y-labels
    pu.adjust_ylabels(ax, x_offset=-0.08)
    # Save figure
    plt.savefig(fig_dir+'pp_fac.'+format, dpi=dpi,transparent=True)

    #--------------------------------------------------------------------
    # Dep / Without gliotransmission
    #--------------------------------------------------------------------
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=figsize, sharex=True,
                           gridspec_kw={'height_ratios': [2, 1],
                                        'hspace': 0.15,
                                        'top': TopMargin, 'bottom': BottomMargin,
                                        'left': LeftMargin, 'right': RightMargin})

    cpt.stp_raster(sol['dep']['spk'], sol['dep']['u'], sol['dep']['r'], sol['dep']['is'],
                   N=20, ax=ax[0], color=colors)

    # Plot histogram
    drf = [2,0,1]
    tbins,dep,fac,rec = ppp_histogram(tRange,sol['dep']['spk'], sol['dep']['u'], sol['dep']['r'], sol['dep']['is'], binw=binw, average=True)
    dep /= binw
    fac /= binw
    rec /= binw
    ax[1].step(tbins,dep,where='pre',color=colors['dep'],zorder=drf[0])
    ax[1].step(tbins,rec,where='pre',color=colors['rec'],zorder=drf[1])
    ax[1].step(tbins,fac,where='pre',color=colors['fac'],zorder=drf[2])
    # Adjust spines
    pu.adjust_spines(ax[0], ['left'])
    pu.adjust_spines(ax[1], ['left','bottom'])
    # Adjust axes
    ax[0].set(xlim=tRange,xticks=[],
              ylim=yRange['si'],yticks=yTicks['si'],ylabel='Synapse Index')
    ax[1].set(xlim=tRange,xticks=tTicks,
              xticklabels=[str(t) for t in tTicks],xlabel='Time (s)',
              ylim=yRange['pp'],yticks=yTicks['pp'],ylabel='PPP  (pairs/s)')
    # Adjust y-labels
    pu.adjust_ylabels(ax, x_offset=-0.08)

    # Add legend
    d = patches.Patch(color=colors['dep'], label='PPD')
    f = patches.Patch(color=colors['fac'], label='PPF')
    r = patches.Patch(color=colors['rec'], label='RFD')
    ax[1].legend(handles=[d,f,r],loc='upper right', facecolor='w', edgecolor='none')

    # # Save figure
    plt.savefig(fig_dir+'pp_dep.'+format, dpi=dpi,transparent=True)

    #--------------------------------------------------------------------
    # Fac / With RI gliotransmission
    #--------------------------------------------------------------------
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=figsize, sharex=True,
                           gridspec_kw={'height_ratios': [0.3, 2, 1],
                                        'hspace': 0.15,
                                        'top': TopMargin, 'bottom': BottomMargin,
                                        'left': LeftMargin, 'right': RightMargin})
    # Plot
    ml, sl, bl = ax[0].stem(gres, np.ones(len(gres)), linefmt='-', markerfmt=" ", basefmt=" ")
    plt.setp(sl, color=colors['ra'])
    # Plot Raster
    cpt.stp_raster(sol['fac_ri']['spk'], sol['fac_ri']['u'], sol['fac_ri']['r'], sol['fac_ri']['is'],
                   N=20, ax=ax[1], color=colors)
    # Plot histogram
    drf = [2,0,1]
    tbins,dep,fac,rec = ppp_histogram(tRange,sol['fac_ri']['spk'], sol['fac_ri']['u'], sol['fac_ri']['r'], sol['fac_ri']['is'], binw=binw, average=True)
    dep /= binw
    fac /= binw
    rec /= binw
    ax[2].step(tbins,dep,where='pre',color=colors['dep'],zorder=drf[0])
    ax[2].step(tbins,rec,where='pre',color=colors['rec'],zorder=drf[1])
    ax[2].step(tbins,fac,where='pre',color=colors['fac'],zorder=drf[2])
    # Adjust spines
    pu.adjust_spines(ax[0], [])
    pu.adjust_spines(ax[1], ['left'])
    pu.adjust_spines(ax[2], ['left','bottom'])
    # Adjust axes
    ax[0].set(xlim=tRange,xticks=[],
              ylim=[0.,1.],yticks=[],ylabel='GREs')
    ax[1].set(xlim=tRange,xticks=[],
              ylim=yRange['si'],yticks=yTicks['si'],ylabel='Synapse Index')
    ax[2].set(xlim=tRange,xticks=tTicks,
              xticklabels=[str(t) for t in tTicks],xlabel='Time (s)',
              ylim=yRange['pp'],yticks=yTicks['pp'],ylabel='PPP  (pairs/s)')
    # Adjust y-labels
    pu.adjust_ylabels(ax, x_offset=-0.08)
    # Save figure
    plt.savefig(fig_dir+'pp_fac_ri.'+format, dpi=dpi,transparent=True)

    #--------------------------------------------------------------------
    # Dep / With RD gliotransmission
    #--------------------------------------------------------------------
    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=figsize, sharex=True,
                           gridspec_kw={'height_ratios': [0.3, 2, 1],
                                        'hspace': 0.15,
                                        'top': TopMargin, 'bottom': BottomMargin,
                                        'left': LeftMargin, 'right': RightMargin})
    # Plot
    ml, sl, bl = ax[0].stem(gres, np.ones(len(gres)), linefmt='-', markerfmt=" ", basefmt=" ")
    plt.setp(sl, color=colors['ra'])
    # Plot Raster
    cpt.stp_raster(sol['dep_rd']['spk'], sol['dep_rd']['u'], sol['dep_rd']['r'], sol['dep_rd']['is'],
                   N=20, ax=ax[1], color=colors)
    # Plot histogram
    drf = [1,0,2]
    tbins,dep,fac,rec = ppp_histogram(tRange,sol['dep_rd']['spk'], sol['dep_rd']['u'], sol['dep_rd']['r'], sol['dep_rd']['is'], binw=binw, average=True)
    dep /= binw
    fac /= binw
    rec /= binw
    ax[2].step(tbins,dep,where='pre',color=colors['dep'],zorder=drf[0])
    ax[2].step(tbins,rec,where='pre',color=colors['rec'],zorder=drf[1])
    ax[2].step(tbins,fac,where='pre',color=colors['fac'],zorder=drf[2])
    # Adjust spines
    pu.adjust_spines(ax[0], [])
    pu.adjust_spines(ax[1], ['left'])
    pu.adjust_spines(ax[2], ['left','bottom'])
    # Adjust axes
    ax[0].set(xlim=tRange,xticks=[],
              ylim=[0.,1.],yticks=[],ylabel='GREs')
    ax[1].set(xlim=tRange,xticks=[],
              ylim=yRange['si'],yticks=yTicks['si'],ylabel='Synapse Index')
    ax[2].set(xlim=tRange,xticks=tTicks,
              xticklabels=[str(t) for t in tTicks],xlabel='Time (s)',
              ylim=yRange['pp'],yticks=yTicks['pp'],ylabel='PPP  (pairs/s)')
    # Adjust y-labels
    pu.adjust_ylabels(ax, x_offset=-0.08)
    # Save figure
    plt.savefig(fig_dir+'pp_dep_rd.'+format, dpi=dpi,transparent=True)

    plt.show()

#-----------------------------------------------------------------------------------------------------------------------
# Figure 6 / Mean-field model against spiking model
#-----------------------------------------------------------------------------------------------------------------------
def step_gre_stimulation(pars,rate_syn,rate_gtr,twin,N_syn,N_gtr,dt=1e-3,gtr=True):
    """
    Simulate ASN model for step variations of GRE rate, at constant synaptic rate.

    Maurizio De Pitta', Basque Center of Applied Mathematics, June 28, 2018.
    """
    rate_gtr = np.asarray(rate_gtr,dtype=float)
    twin = np.asarray(twin,dtype=float)
    assert rate_gtr.size==twin.size, "rate of GTR and twin must be of same length"
    tfin = np.cumsum(np.r_[0,twin])
    time = np.arange(0.,tfin[-1],dt)
    if not gtr:
        spk_sol = []
        spk_idx = []
        spk_tsp = []
        # Added for debugging purposes
        pars['ICs'] = [1.0,0.,0.0]
        exo = gtm.gtrs(model='exo_spk',**pars)
        for i,tf in enumerate(twin):
            if i>0: exo.ICs = exo.sol['LCs'] # Start from last point of previous simulation
            exo.stimulation(rate=rate_syn,twin=[0., tf], Nsyn=N_syn, Nvar=2,stimulus='poisson')
            exo.simulate(reconstruct=False,recompile=0)
            spk_tsp = np.r_[spk_tsp,tfin[i]+exo.sol['spk']]
            spk_sol = np.r_[spk_sol,exo.sol['y']]
            spk_idx = np.r_[spk_idx,exo.sol['is']]

        # Reconstruct y solution
        sol = []
        for i in xrange(int(max(spk_idx))+1):
            index = spk_idx==i
            tmp = gtm.reconstruct_solution(spk_tsp[index], spk_sol[index], [], [0.,tfin[-1]], 0.0, pars['taui'], 'y', dt=dt)
            keep = np.searchsorted(tmp[0],time)
            sol.append(tmp[1][keep])
        y = np.mean(np.vstack(sol),axis=0)
        sol = {'t': time, 'y': y}
    else:
        y  = []
        ga = []
        gs = []
        ra = []
        gre = []
        for j in xrange(N_gtr):
            print "Astro :",str(j),"/",str(N_gtr)
            spk_sol = []  # y
            spk_idx = []  # is
            spk_tsp = []  # spk
            sol_ga  = []  # G_A
            sol_gms = []  # Gamma_S
            sol_tim = []  # Time of continuous solution
            sol_gre = []  # gre
            sol_ra  = []  # ra
            pars['ICr'] = [1.0,0.0,0.0,1.0]
            asn = gtm.gtrs(model='asn_spk',**pars)
            for i,tf in enumerate(twin):
                if i>0:
                    asn.ICs = np.r_[asn.sol['ga'][-1],asn.sol['gammas'][-1]]
                    asn.ICr = asn.sol['LCr'] # Start from last point of previous simulation
                asn.stimulation(rate=rate_syn,twin=[0., tf], Nsyn=N_syn, Nvar=2,stimulus='poisson',
                                rate_gtr=rate_gtr[i],twin_gtr=[0., tf],stimulus_gtr='poisson')
                options = su.solver_opts(t0=0.0, tfin=tf, dt=1e-4, atol=1e-8, rtol=1e-8, method="gsl_rk2imp")
                asn.simulate(algparams=options,reconstruct=False,recompile=0)
                print "Screening gamma_S: step: ",i,"-> value: ",max(asn.sol['gammas'])," --> u0: ",u_0(max(asn.sol['gammas']),asn.pars['u0'],asn.pars['alpha'])
                spk_tsp = np.r_[spk_tsp,tfin[i]+asn.sol['spk']]
                spk_sol = np.r_[spk_sol,asn.sol['y']]
                spk_idx = np.r_[spk_idx,asn.sol['is']]
                # Continuous solution
                sol_tim = np.r_[sol_tim,tfin[i]+asn.sol['t']]
                sol_gms = np.r_[sol_gms,asn.sol['gammas']]
                sol_ga  = np.r_[sol_ga,asn.sol['ga']]
                sol_gre = np.r_[sol_gre,tfin[i]+asn.sol['gre']]
                sol_ra = np.r_[sol_ra,asn.sol['ra']]

            # Reconstruct y solution
            sol = []
            for i in xrange(int(max(spk_idx))+1):
                index = spk_idx==i
                tmp = gtm.reconstruct_solution(spk_tsp[index], spk_sol[index], [], [0.,tfin[-1]], 0.0, pars['taui'], 'y', dt=dt)
                keep = np.searchsorted(tmp[0],time)
                sol.append(tmp[1][keep])
            sol = np.mean(np.vstack(sol),axis=0)
            y.append(sol)

            # Append gs and ga
            keep = np.searchsorted(sol_tim, time)
            ga.append(sol_ga[keep])
            gs.append(sol_gms[keep])

            # Append gre and ra
            gre.append(sol_gre)
            ra.append(sol_ra)

        # Do final averaging
        sol = {}
        sol['t']  = time
        sol['y']  = np.mean(np.vstack(y),axis=0)
        sol['ga'] = np.mean(np.vstack(ga),axis=0)
        sol['gammas'] = np.mean(np.vstack(gs),axis=0)
        gre = np.concatenate(gre)
        ra  = np.concatenate(ra)
        idx = np.argsort(gre)
        sol['gre'] = gre[idx]
        sol['ra']  = ra[idx]

    return sol

def step_gre_stimulation_mf(pars,rate_syn,rate_gtr,twin,dt=1e-3,gtr=True):
    """
    Simulate the mean field ASN model for step variations of GRE rate, at constant synaptic rate.

    Maurizio De Pitta', Basque Center of Applied Mathematics, June 28, 2018.
    """
    rate_gtr = np.asarray(rate_gtr,dtype=float)
    twin = np.asarray(twin,dtype=float)
    assert rate_gtr.size==twin.size, "rate of GTR and twin must be of same length"
    tfin = np.cumsum(np.r_[0,twin])
    if not gtr:
        t = []
        x = []
        u = []
        y = []
        # Added for debugging purposes
        pars['ICs'] = [1.0,pars['u0']*pars['psc0'],pars['u0']]
        exo = gtm.gtrs(model='exo_ave',**pars)
        for i,tf in enumerate(twin):
            options = options = su.solver_opts(t0=0.0, tfin=tf, dt=dt, atol=1e-8, rtol=1e-8, method="gsl_rk2imp")
            if i>0: exo.ICs = np.r_[exo.sol['x'][-1],exo.sol['y'][-1],exo.sol['u'][-1]], # Start from last point of previous simulation
            exo.stimulation(rate=rate_syn,twin=[0., tf], Nvar=2)
            exo.simulate(algparams=options,recompile=0)
            t = np.r_[t,tfin[i]+exo.sol['t']]
            x = np.r_[x,exo.sol['x']]
            y = np.r_[y,exo.sol['y']]
            u = np.r_[u,exo.sol['u']]

        # Reconstruct full solution
        sol = {'t' : t, 'x' : x, 'y' : y, 'u' : u}
    else:
        t = []
        x = []
        u = []
        y = []
        xa = []
        ga = []
        gs = []
        # pars['ICs'] = np.asarray([0.,0.])
        pars['ICr'] = np.asarray([1.0, 0., pars['u0'],1.0])
        asn = gtm.gtrs(model='asn_ave',**pars)
        for i,tf in enumerate(twin):
            if i > 0:
                asn.ICs = np.r_[asn.sol['ga'][-1], asn.sol['gammas'][-1]]
                asn.ICr = np.r_[asn.sol['x'][-1],asn.sol['y'][-1],asn.sol['u'][-1],asn.sol['xa'][-1]]  # Start from last point of previous simulation
            options = su.solver_opts(t0=0.0, tfin=tf, dt=1e-4, atol=1e-8, rtol=1e-8, method="gsl_bsimp")
            asn.stimulation(rate=rate_syn, twin=[0., tf], Nvar=2,rate_gtr=rate_gtr[i], twin_gtr=[0., tf])
            asn.simulate(algparams=options, recompile=0)
            t = np.r_[t,tfin[i]+asn.sol['t']]
            x = np.r_[x,asn.sol['x']]
            y = np.r_[y,asn.sol['y']]
            u = np.r_[u,asn.sol['u']]
            xa = np.r_[xa,asn.sol['xa']]
            ga = np.r_[ga,asn.sol['ga']]
            gs = np.r_[gs,asn.sol['gammas']]

        # Reconstruct full solution
        sol = {'t': t, 'x': x, 'y': y, 'u': u,
               'xa': xa, 'ga': ga, 'gammas': gs}

    return sol

def mf_model(sim=True,data_dir='../data/',fig_dir='../../Figures/',format='pdf',dpi=600):
    """
    Fit spiking dynamics with mean-field dynamics.

    Maurizio De Pitta', Basque Center of Applied Mathematics, June 29, 2018.
    """

    # Load MPL default style
    plt.style.use('figures.mplstyle')

    # General figure specifications
    figsize = (6.5,6.5)
    TopMargin = 0.98
    BottomMargin = 0.11
    LeftMargin = 0.12
    RightMargin = 0.96

    # Color specifications
    colors = {'nu_s': Red,
              'ra_s': LightBlue,
              'ga_s': LightOrange,
              'Gs_s': DarkGray,
              'rd_s': LightMagenta,
              'ri_s': LightGreen,
              'nu_m': DarkRed,
              'ra_m': DarkBlue,
              'ga_m': DarkOrange,
              'Gs_m': 'black',
              'rd_m': DarkMagenta,
              'ri_m': DarkGreen}

    # Model parameters and ranges
    tRange = [0., 60.]
    tTicks = np.arange(tRange[0],1.1*tRange[1],10.)
    yRange = {'nua' : [-0.02,0.52],
              'ra'  : [0.,0.7],
              'ga'  : [-0.1,5.1],
              'gs'  : [-0.02,1.0],
              'y'   : [0.,1.]}
    yTicks = {'nua' : np.arange(0.,0.55,.1),
              'ra'  : np.arange(0.,0.8,.2),
              'ga'  : np.arange(0.,5.1),
              'gs'  : np.arange(0.,1.01,0.2),
              'y'   : np.arange(0,1.01,0.5)}

    sol,pars = {},{}
    N_syn = 10
    N_gtr = 100
    rate_syn = 2.0
    rate_gtr = [0.,0.1,0.5]
    twin = [5.,30.,25.]
    if sim:
        print "1. Release-increasing gliotransmission / Spiking model"
        pars['ri'] = gtm.asn_parameters(u0=0.3,taud=0.3,tauf=0.5, Ou=0.0,
                                        taue=1./40., Gtot=150, Og=0.3, taug=20., alpha=1.0)
        sol['ri'] = step_gre_stimulation(pars['ri'], rate_syn, rate_gtr, twin, N_syn, N_gtr, dt=1e-3, gtr=True)
        svu.savedata([sol, pars], data_dir + 'mf_model.pkl')

        print "2. Release-increasing gliotransmission / Mean-field model"
        sol['ri_mf'] = step_gre_stimulation_mf(pars['ri'],rate_syn, rate_gtr, twin, dt=1e-3, gtr=True)
        svu.savedata([sol, pars], data_dir + 'mf_model.pkl')

        print "3. Release-decreasing gliotransmission / Spiking model"
        pars['rd'] = gtm.asn_parameters(u0=0.3,taud=0.3,tauf=0.5,Ou=0.0,
                                        taue=1./40., Gtot=150, Og=0.3, taug=20., alpha=0.0)
        sol['rd'] = step_gre_stimulation(pars['rd'], rate_syn, rate_gtr, twin, N_syn, N_gtr, dt=1e-3, gtr=True)
        svu.savedata([sol, pars], data_dir + 'mf_model.pkl')

        print "4. Release-decreasing gliotransmission / Mean-field model"
        sol['rd_mf'] = step_gre_stimulation_mf(pars['rd'],rate_syn, rate_gtr, twin, dt=1e-3, gtr=True)
        svu.savedata([sol, pars], data_dir + 'mf_model.pkl')

    # -------------------------------------------------------------------------------------------------------------------
    # Effective plotting
    # -------------------------------------------------------------------------------------------------------------------
    # Load data
    sol,p = svu.loaddata(data_dir+'mf_model.pkl')

    #--------------------------------------------------------------------
    # Fac / With RI gliotransmission
    #--------------------------------------------------------------------
    fig, ax = plt.subplots(nrows=5, ncols=1, figsize=figsize, sharex=True,
                           gridspec_kw={'height_ratios': [1, 1, 1, 1, 2],
                                        'hspace': 0.22,
                                        'top': TopMargin, 'bottom': BottomMargin,
                                        'left': LeftMargin, 'right': RightMargin})

    # # Plot spiking model
    ax[0].step(np.r_[0.,np.cumsum(twin)],np.r_[rate_gtr,rate_gtr[-1]],'-',color=colors['nu_m'],where='post')
    ax[1].plot(sol['ri']['gre'],sol['ri']['ra'],'o',color=colors['ra_s'])
    ax[2].plot(sol['ri']['t'], sol['ri']['ga'], '-', color=colors['ga_s'])
    ax[3].plot(sol['ri']['t'], sol['ri']['gammas'],'-', color=colors['Gs_s'])
    ax[4].plot(sol['ri']['t'], sol['ri']['y'],'-', color=colors['ri_s'])
    ax[4].plot(sol['rd']['t'], sol['rd']['y'], '-', color=colors['rd_s'])
    # Plot MF model
    ax[1].plot(sol['ri_mf']['t'], p['ri']['ua']*sol['ri_mf']['xa'], '-', color=colors['ra_m'])
    ax[2].plot(sol['ri_mf']['t'], sol['ri_mf']['ga'], '-', color=colors['ga_m'])
    ax[3].plot(sol['ri_mf']['t'], sol['ri_mf']['gammas'],'-', color=colors['Gs_m'])
    ax[4].plot(sol['ri_mf']['t'], sol['ri_mf']['y'],'-', color=colors['ri_m'], label=r'RIG ($\alpha=1$)')
    ax[4].plot(sol['rd_mf']['t'], sol['rd_mf']['y'],'-', color=colors['rd_m'], label=r'RDG ($\alpha=0$)')
    # Adjust spines
    for i in xrange(4): pu.adjust_spines(ax[i], ['left'])
    pu.adjust_spines(ax[4], ['left','bottom'])
    # Adjust axes
    ax[0].set(xlim=tRange,xticks=[],
              ylim=yRange['nua'],yticks=yTicks['nua'],ylabel=r'$\nu_A(t)$ (Hz)')
    ax[1].set(xlim=tRange,xticks=[],
              ylim=yRange['ra'],yticks=yTicks['ra'],ylabel='$r_A(t)$')
    ax[2].set(xlim=tRange,xticks=[],
              ylim=yRange['ga'],yticks=yTicks['ga'],ylabel=r'$G_A(t)$ ($\mu$M)')
    ax[3].set(xlim=tRange,xticks=[],
              ylim=yRange['gs'],yticks=yTicks['gs'],ylabel=r'$\Gamma_S(t)$')
    ax[4].set(xlim=tRange,xticks=tTicks,
              xticklabels=[str(t) for t in tTicks],xlabel='Time (s)',
              ylim=yRange['y'],yticks=yTicks['y'],ylabel='$I(t)$ (pA)')
    # Adjust y-labels
    pu.adjust_ylabels(ax, x_offset=-0.08)
    # Add legend
    ax[4].legend(loc='center right', facecolor='w', edgecolor='none')

    # Save figure
    plt.savefig(fig_dir+'mf_model.'+format, dpi=dpi,transparent=True)

    plt.show()

#-----------------------------------------------------------------------------------------------------------------------
# Figure 7 / Mean-field U0
#-----------------------------------------------------------------------------------------------------------------------
def spk_u0(pars,nua_vals,**kwargs):

    p = {'nspk': 30,
         'dt'  : 1e-3,}
    ## User-defined parameters
    pars = gu.varargin(pars,**kwargs)

    # Set initial conditions
    pars['ICs'] = [0.,0.]
    pars['ICr'] = [1.0, 0.0, 0.0, 1.0]
    asn = gtm.gtrs(model='asn_spk', **pars)
    u0_vals = np.zeros((np.size(nua_vals),2))
    for i, nua in enumerate(nua_vals):
        print "nua val: ",str(i),"/",str(np.size(nua_vals))
        tfin = np.round(3.0*p['nspk']/nua,0)
        if tfin > 300:
            p['dt'] = 1e-2
        elif (tfin>90)&(tfin<=300):
            p['dt'] = 1e-3
        else:
            p['dt'] = 1e-4
        options = su.solver_opts(t0=0.0, tfin=tfin, dt=p['dt'], atol=1e-8, rtol=1e-8, method="gsl_rk2imp")
        asn.stimulation(rate=0., twin=[0., 1.], Nsyn=1, Nvar=1, stimulus='poisson',
                        rate_gtr=nua, twin_gtr=[0., tfin], stimulus_gtr='poisson')
        asn.simulate(algparams=options, reconstruct=False, recompile=0)
        keep = np.searchsorted(asn.sol['t'],asn.sol['gre'][-p['nspk']:])
        gs = asn.sol['gammas'][keep]
        u0 = u_0(gs,pars['u0'],pars['alpha'])
        u0_vals[i] = np.r_[np.mean(u0),np.std(u0)]
    return u0_vals.T

def mf_brp(sim=True,data_dir='../data/',fig_dir='../../Figures/',format='pdf',dpi=600):
    """
    Mean-field U_0 vs. spiking values of U_0.

    Maurizio De Pitta', Basque Center of Applied Mathematics, June 29, 2018.
    """
    # Load MPL default style
    plt.style.use('figures.mplstyle')

    # General figure specifications
    figsize = (5.,4.5)
    TopMargin = 0.98
    BottomMargin = 0.15
    LeftMargin = 0.15
    RightMargin = 0.95

    # Color specifications
    colors = {'u0'  : DarkRed,
              'uthr': DarkBlue,
              'dep' : LightOrange,
              'fac' : LightBlue,
              'nut' : DarkMagenta}

    # Model parameters and ranges
    nuaRange = np.asarray([0.0095, 1.04])
    nuaTicks = np.asarray([0.01,0.1,1.])
    uRange = np.asarray([0.,1])
    uTicks = np.arange(0.,1.1,0.2)

    u0sol,pars = {},{}
    npts = 10
    nua_vals = np.logspace(-2.,0.,npts)
    # npts = 10
    # nua_vals = np.asarray([0.1,1.])
    if sim:
        print "1. Release-increasing gliotransmission / Spiking model"
        pars['ri'] = gtm.asn_parameters(u0=0.05,taud=0.3,tauf=1.0,
                                        Ou=0.0,taue=1./40., Gtot=150, Og=0.3, taug=20., alpha=1.0)
        u0sol['ri'] = spk_u0(pars['ri'], nua_vals)
        svu.savedata([u0sol, pars], data_dir + 'mf_brp.pkl')

        print "2. Release-decreasing gliotransmission / Spiking model"
        pars['rd'] = gtm.asn_parameters(u0=0.7,taud=0.5,tauf=0.3,
                                        Ou=0.0,taue=1./40., Gtot=150, Og=0.3, taug=20., alpha=0.0)
        u0sol['rd'] = spk_u0(pars['rd'], nua_vals)
        svu.savedata([u0sol, pars], data_dir + 'mf_brp.pkl')

    # -------------------------------------------------------------------------------------------------------------------
    # Effective plotting
    # -------------------------------------------------------------------------------------------------------------------
    # Load data
    upts,p = svu.loaddata(data_dir+'mf_brp.pkl')
    # Generate substrate to plot MF curve of U0_inf
    npts = 100
    nua = np.logspace(np.log10(nua_vals[0]),np.log10(nua_vals[-1]),npts)
    #--------------------------------------------------------------------
    # Fac / With RI gliotransmission
    #--------------------------------------------------------------------
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize, sharex=True,
                           gridspec_kw={'top': TopMargin, 'bottom': BottomMargin,
                                        'left': LeftMargin, 'right': RightMargin})

    U0 = u0_inf(nua,1/p['ri']['taug'],1/p['ri']['taua'],p['ri']['ua'],p['ri']['js'],p['ri']['alpha'],p['ri']['u0'])
    ut = u_thr(p['ri']['tauf']/p['ri']['taud'])
    nut = nu_thr(1/p['ri']['taug'],1/p['ri']['taua'],p['ri']['ua'],p['ri']['js'],p['ri']['alpha'],ut,p['ri']['u0'])
    fac = patches.Rectangle([nuaRange[0],uRange[0]],nut-nuaRange[0],ut,color=colors['fac'],ec='none',fc=colors['fac'],zorder=0)
    dep = patches.Rectangle([nut, ut], nuaRange[1]-nut, uRange[1]-ut, color=colors['dep'], ec='none',fc=colors['dep'],zorder=0)
    ax.add_patch(fac)
    ax.add_patch(dep)
    # Plot thresholds
    ax.plot(nuaRange[[0,-1]],ut*np.ones(2),'-',color=colors['uthr'])
    ax.plot(nut*np.ones(2),uRange,'--',color=colors['nut'])
    # Plot U0
    ax.plot(nua,U0,'-',color=colors['u0'],zorder=2)
    ax.errorbar(nua_vals,upts['ri'][0],yerr=upts['ri'][1],
                fmt='ko',ms=6,capsize=3,capthick=1.8,zorder=1)
    # Adjust spines
    pu.adjust_spines(ax, ['left','bottom'])
    # Adjust axes
    ax.set_xscale('log', nonposx='clip')
    ax.set(xlim=nuaRange,xticks=nuaTicks,xlabel=r'$\nu_A(t\rightarrow\infty)$ (Hz)',
              ylim=uRange,yticks=uTicks,ylabel=r'$U_0(t\rightarrow\infty)$')
    ax.set(xticklabels = [str(nu) for nu in nuaTicks])

    # # Save figure
    plt.savefig(fig_dir+'mf_brp_fac.'+format, dpi=dpi,transparent=True)

    #--------------------------------------------------------------------
    # Dep / With RD gliotransmission
    #--------------------------------------------------------------------
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize, sharex=True,
                           gridspec_kw={'top': TopMargin, 'bottom': BottomMargin,
                                        'left': LeftMargin, 'right': RightMargin})

    U0 = u0_inf(nua,1/p['rd']['taug'],1/p['rd']['taua'],p['rd']['ua'],p['rd']['js'],p['rd']['alpha'],p['rd']['u0'])
    ut = u_thr(p['rd']['tauf'] / p['rd']['taud'])
    nut = nu_thr(1/p['rd']['taug'],1/p['rd']['taua'],p['rd']['ua'],p['rd']['js'],p['rd']['alpha'],ut,p['rd']['u0'])

    fac = patches.Rectangle([nut,uRange[0]],nuaRange[1]-nut,ut,color=colors['fac'],ec='none',fc=colors['fac'],zorder=0)
    dep = patches.Rectangle([nuaRange[0], ut], nut-nuaRange[0], uRange[1]-ut, color=colors['dep'], ec='none',fc=colors['dep'],zorder=0)
    ax.add_patch(fac)
    ax.add_patch(dep)
    # Plot thresholds
    ax.plot(nuaRange[[0, -1]], ut*np.ones(2), '-', color=colors['uthr'])
    ax.plot(nut*np.ones(2),uRange,'--',color=colors['nut'])
    # Plot U_0
    ax.plot(nua,U0,'-',color=colors['u0'],zorder=2)
    ax.errorbar(nua_vals,upts['rd'][0],yerr=upts['rd'][1],
                fmt='ko', ms=6, capsize=3, capthick=1.8,zorder=1)
    # Adjust spines
    pu.adjust_spines(ax, ['left','bottom'])
    # Adjust axes
    ax.set_xscale('log', nonposx='clip')
    ax.set(xlim=nuaRange,xticks=nuaTicks,xlabel=r'$\nu_A(t\rightarrow\infty)$ (Hz)',
              ylim=uRange,yticks=uTicks,ylabel=r'$U_0(t\rightarrow\infty)$')
    ax.set(xticklabels = [str(nu) for nu in nuaTicks])
    # # Save figure
    plt.savefig(fig_dir+'mf_brp_dep.'+format, dpi=dpi,transparent=True)

    plt.show()

#-----------------------------------------------------------------------------------------------------------------------
# Figure 8 / Heat maps of plasticity
#-----------------------------------------------------------------------------------------------------------------------
def generate_domains(npts,ratioRange=[0.1,10.],u0Range=[0.,1.],grid=False):
    ratio = np.logspace(np.log10(ratioRange[0]), np.log10(ratioRange[-1]), npts)
    u0 = np.linspace(0, 1., npts)
    u_theta = u_thr(ratio)
    if grid:
        R,U0 = np.meshgrid(ratio,u0,indexing='xy') # Generate Cartesian coordinates
        return ratio,u0,u_theta,R,U0
    else:
        return ratio,u0,u_theta

def synapse_space(npts,ratioRange=[0.1,10.],u0Range=[0.,1.],colors={},domains=False):
    assert type(colors)==dict, "Colors must be a dictionary with 'dep' and 'fac' keys"
    ratioRange = np.asarray(ratioRange)
    u0Range = np.asarray(u0Range)
    ratio, u0, u_theta = generate_domains(npts=npts,ratioRange=ratioRange,u0Range=u0Range,grid=False)
    if 'dep' not in colors:
        colors['dep'] = DarkOrange
    if 'fac' not in colors:
        colors['fac'] = DarkBlue
    dep = patches.Polygon(np.vstack((np.r_[ratioRange[0],ratio,ratioRange[[1,1,0]]],np.r_[u_theta[0],u_theta,u_theta[-1],u0Range[[-1,-1]]])).T,
                          closed=True,facecolor=colors['dep'])
    fac = patches.Polygon(np.vstack((np.r_[ratio,ratioRange[[1,1,0,0]]],np.r_[u_theta,u_theta[-1],u0Range[[0,0]],u_theta[0]])).T,
                          closed=True,facecolor=colors['fac'])
    if domains:
        return dep,fac,ratio,u0,u_theta
    else:
        return dep,fac

def nu_threhsold(pars, npts, ratioRange=[0.1, 10.], u0Range=[0., 1.], grid=True):
    """
    Generate 2D data to map nu_threshold.

    Maurizio De Pitta', Basque Center of Applied Mathematics, June 17, 2018.
    """

    u0a, NUT_, NUT = {}, {}, {}
    alpha = [0, 1]
    # Generate domains
    ratio, u0, u_theta = generate_domains(npts=npts, ratioRange=ratioRange, u0Range=u0Range, grid=False)
    # Generate grid
    R, U0 = np.meshgrid(ratio, u0)
    for i, a in enumerate(alpha):
        if a == 0:
            tag = 'dep'
        else:
            tag = 'fac'
        u0a[tag] = u0_asymptote(ratio, pars['js'], 1 / pars['taua'], 1 / pars['taug'], a)
        NUT_[tag] = nu_thr(1 / pars['taug'], 1 / pars['taua'], pars['ua'], pars['js'], a, u_thr(R), U0)
        # Manipulate NUT and mask values  <0 by NaNs
        NUT_[tag][NUT_[tag] < 0] = np.nan
        NUT[tag] = ma.array(NUT_[tag], mask=np.isnan(NUT_[tag]))
    if grid:
        return NUT, u0a, R, U0
    else:
        return NUT, u0a

def adjust_synaptic_map(ax,ratioRange,u0Range):
    # Adjust axes in a fixed way for all figures
    pu.adjust_spines(ax, ['left', 'bottom'])
    xticks = 10 ** np.arange(-1, 1.1, 1)
    ax.set(xlim=ratioRange,ylim=u0Range,
           xscale='log',
           xlabel=r'$\Omega_d/\Omega_f$ ratio', ylabel='Synaptic BRP ($U_0$)')
    ax.set(xticks=xticks, xticklabels=([str(xt) for xt in xticks]),
              yticks=np.arange(0., 1.1, 0.5))

def threhsold_heatmap(pars, NPTS, ratioRange, u0Range, nuRange, colors, ax=None):
    dep, fac, ratio, u0, u_theta = synapse_space(NPTS, ratioRange=ratioRange, u0Range=u0Range, domains=True,
                                                 colors=colors)
    # Compute threshold frequency
    NUT, u0a, R, U0 = nu_threhsold(pars, NPTS, ratioRange=ratioRange, u0Range=u0Range, grid=True)

    if not ax: fig, ax = plt.subplots(figsize=(5, 5))
    # Add graphic objects
    ax.add_patch(dep)
    ax.add_patch(fac)
    cmap = mpl.cm.viridis
    # cmap.set_bad('w', 1.)
    ax.pcolormesh(R, U0, NUT['dep'], cmap=cmap, vmin=nuRange[0], vmax=nuRange[-1], shading='gouraud',
                  zorder=10)
    ax.pcolormesh(R, U0, NUT['fac'], cmap=cmap, vmin=nuRange[0], vmax=nuRange[-1], shading='gouraud',
                  zorder=10)
    # Add lines
    ax.plot(ratio, u0a['dep'], '--', ratio, u0a['fac'], '--', color=colors['asymptote'], zorder=20)
    ax.plot(ratio, u_theta, '-', color=colors['thr'], zorder=20)

    # Adjust axes
    adjust_synaptic_map(ax, ratioRange, u0Range)

    return ax

def pswitch_domains(dir='../Figures/',format='pdf',dpi=600):
    """
    Generate heat maps of different domains of short-term plasticity

    Maurizio De Pitta', Basque Center for Applied Mathematics, June 20, 2018
    """
    # Load MPL default style
    plt.style.use('figures.mplstyle')

    # General figure specifications
    figsize = (5,4)
    TopMargin = 0.98
    BottomMargin = 0.15
    LeftMargin = 0.14
    RightMargin = 0.95

    # Color specifications
    colors = {'dep': DarkOrange,
              'fac': DarkBlue,
              'thr': 'white',
              'asymptote': DarkRed,
              'cmap': 'viridis'}

    # Model parameters and ranges
    ratioRange = [0.1,10.]
    u0Range = [0., 1.]
    nuRange = [0., 5.]
    pars = gtm.asn_parameters(taua=1.,taug=30.,tau_e=1./40,js=0.2)

    #-------------------------------------------------------------------------------------------------------------------
    # General synaptic space
    #-------------------------------------------------------------------------------------------------------------------
    print "1. Synaptic space"
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize, sharex=False,
                           gridspec_kw={'top': TopMargin, 'bottom': BottomMargin,
                                        'left': LeftMargin, 'right': RightMargin})
    # Generate patches
    NPTS = 100
    dep,fac,ratio,u0,u_theta = synapse_space(NPTS,ratioRange=ratioRange,u0Range=u0Range,domains=True)
    # Add patches
    ax.add_patch(dep)
    ax.add_patch(fac)
    ax.plot(ratio,u_theta,'-',color=colors['thr'])

    # Adjust axes
    adjust_synaptic_map(ax, ratioRange, u0Range)
    plt.savefig(dir + 'heatmap_synspace.' + format, dpi=dpi)

    # -------------------------------------------------------------------------------------------------------------------
    # Control scenario with Gt. release
    # -------------------------------------------------------------------------------------------------------------------
    print "2. Control space"

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize, sharex=False,
                           gridspec_kw={'top': TopMargin, 'bottom': BottomMargin,
                                        'left': LeftMargin, 'right': RightMargin})
    # General settings
    NPTS = 200
    colors['dep'] = LightOrange
    colors['fac'] = LightBlue
    # Generate domains
    threhsold_heatmap(pars, NPTS, ratioRange, u0Range, nuRange, colors, ax=ax)
    ax.set_rasterized(True)  # This option is needed to save transparency
    plt.savefig(dir + 'heatmap_ctrl.' + format, dpi=dpi,transparent=True)

    # -------------------------------------------------------------------------------------------------------------------
    # 2-fold increase of taua
    # -------------------------------------------------------------------------------------------------------------------
    print "3. Two-fold increase of taua"

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize, sharex=False,
                           gridspec_kw={'top': TopMargin, 'bottom': BottomMargin,
                                        'left': LeftMargin, 'right': RightMargin})
    # Generate parameters
    p = cp.copy(pars)
    p['taua'] *= 2
    # Generate heatmap
    threhsold_heatmap(p, NPTS, ratioRange, u0Range, nuRange, colors, ax=ax)
    ax.set_rasterized(True)
    plt.savefig(dir + 'heatmap_taua.' + format, dpi=dpi,transparent=True)

    # -------------------------------------------------------------------------------------------------------------------
    # 3-fold decrease of taug
    # -------------------------------------------------------------------------------------------------------------------
    print "4. Three-fold decrease of taug"

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize, sharex=False,
                           gridspec_kw={'top': TopMargin, 'bottom': BottomMargin,
                                        'left': LeftMargin, 'right': RightMargin})
    # Generate parameters
    p = cp.copy(pars)
    p['taug'] /= 3
    # Generate heatmap
    threhsold_heatmap(p, NPTS, ratioRange, u0Range, nuRange, colors, ax=ax)
    ax.set_rasterized(True)
    plt.savefig(dir + 'heatmap_taug.' + format, dpi=dpi,transparent=True)

    # -------------------------------------------------------------------------------------------------------------------
    # 5-fold decrease of js
    # -------------------------------------------------------------------------------------------------------------------
    print "5. Five-fold decrease of js"

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize, sharex=False,
                           gridspec_kw={'top': TopMargin, 'bottom': BottomMargin,
                                        'left': LeftMargin, 'right': RightMargin})
    # Generate parameters
    p = cp.copy(pars)
    p['js'] /= 5
    # Produce heatmap
    threhsold_heatmap(p, NPTS, ratioRange, u0Range, nuRange, colors, ax=ax)
    ax.set_rasterized(True)
    plt.savefig(dir+'heatmap_js.'+format, dpi=dpi,transparent=True)

    # -------------------------------------------------------------------------------------------------------------------
    # Colorbar
    # -------------------------------------------------------------------------------------------------------------------
    print "6. Colorbar"

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(figsize[0],1), sharex=False,
                           gridspec_kw={'top': 0.98, 'bottom': 0.7,
                                        'left': 0.1, 'right': 0.98})
    # Create a standalone colormap (taken from https://matplotlib.org/examples/api/colorbar_only.html)
    cmap = mpl.cm.viridis
    cmap.set_over(cmap(255))
    norm = mpl.colors.Normalize(vmin=nuRange[0], vmax=nuRange[-1])
    # ColorbarBase derives from ScalarMappable and puts a colorbar
    # in a specified axes, so it has everything needed for a
    # standalone colorbar.  There are many more kwargs, but the
    # following gives a basic continuous colorbar with ticks
    # and labels.
    cb = mpl.colorbar.ColorbarBase(ax,cmap=cmap,
                                   norm=norm,ticks=np.arange(0,5.1),
                                   extend='max',
                                   orientation='horizontal')
    cb.set_label(r'Threshold GRE rate, $\nu_\theta$ (Hz)')

    # Save figure
    plt.savefig(dir+'heatmap_cbar.'+format, dpi=dpi,transparent=True)

    # -------------------------------------------------------------------------------------------------------------------
    # Show figures
    # -------------------------------------------------------------------------------------------------------------------
    plt.show()

#-----------------------------------------------------------------------------------------------------------------------
# Figure A1 / Error on Gamma_S approximation
#-----------------------------------------------------------------------------------------------------------------------
def gammas_error(fig_dir='../../Figures/',format='pdf',dpi=600):
    # Load MPL default style
    plt.style.use('figures.mplstyle')

    Gamma_inf = lambda nua,oma,omg,js,ua : js*oma*ua*nua/(oma*omg + (js*oma+omg)*ua*nua)
    Gamma_approx = lambda nua, omg, js, ua: (np.exp(js*ua)-1)*nua / ((omg + nua)*np.exp(js*ua) - nua)

    # Color specifications
    colors = {'gs_inf'  : DarkRed,
              'gs_apx'  : DarkGray}

    # General figure specifications
    figsize = (5.,4.5)
    TopMargin = 0.98
    BottomMargin = 0.15
    LeftMargin = 0.15
    RightMargin = 0.95

    # Axes
    nuaRange = [0.01,1.]
    nuaTicks = [0.01,0.1,1.0]
    yRange = {'gs' : [0.,1.],
              'err' : [-18,8.]}
    yTicks = {'gs' : [0.,0.5,1.],
              'err' : np.arange(-20, 11.,5)}

    # Generate substrate
    npts = 100
    nua = np.logspace(np.log10(nuaRange[0]),np.log10(nuaRange[1]),npts)
    # Model parameters
    pars = gtm.asn_parameters(taue=1./40., Gtot=150, Og=0.3, taug=20.)

    # Generate Curves
    Gs_i = Gamma_inf(nua,1/pars['taua'],1/pars['taug'],pars['js'],pars['ua'])
    Gs_a = Gamma_approx(nua,1/pars['taug'],pars['js'],pars['ua'])
    Err = (Gs_a - Gs_i)/Gs_i*100

    #--------------------------------------------------------------------
    # Effective plotting
    #--------------------------------------------------------------------
    # Plot Theory vs. approximation
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize, sharex=True,
                           gridspec_kw={'top': TopMargin, 'bottom': BottomMargin,
                                        'left': LeftMargin, 'right': RightMargin})

    ax.plot(nua,Gs_i,'-',color=colors['gs_inf'],
            label=r'$\Gamma_S(t\rightarrow \infty)$ (exact)')
    ax.plot(nua,Gs_a,'-',color=colors['gs_apx'],
            label=r'$\Gamma_S(t\rightarrow \infty)$ (approx)')
    # Adjust spines
    pu.adjust_spines(ax, ['left','bottom'])
    # Adjust axes
    ax.set_xscale('log', nonposx='clip')
    ax.set(xlim=nuaRange,xticks=nuaTicks,xlabel=r'$\nu_A(t\rightarrow\infty)$ (Hz)',
           ylim=yRange['gs'],yticks=yTicks['gs'],ylabel=r'$\Gamma_S(t\rightarrow\infty)$')
    ax.set(xticklabels = [str(nu) for nu in nuaTicks])
    pu.adjust_ylabels(ax, x_offset=-0.1)
    # Add legend
    ax.legend(loc='lower right', facecolor='w', edgecolor='none')
    # # Save figure
    plt.savefig(fig_dir+'gs_approx.'+format, dpi=dpi,transparent=True)



    # Plot error
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=figsize, sharex=True,
                           gridspec_kw={'top': TopMargin, 'bottom': BottomMargin,
                                        'left': LeftMargin, 'right': RightMargin})
    ax.plot(nua,Err,'-',color='k')
    # Adjust spines
    pu.adjust_spines(ax, ['left','bottom'])
    # Adjust axes
    ax.set_xscale('log', nonposx='clip')
    ax.set(xlim=nuaRange,xticks=nuaTicks,xlabel=r'$\nu_A(t\rightarrow\infty)$ (Hz)',
           ylim=yRange['err'],yticks=yTicks['err'],ylabel='% Error')
    ax.set(xticklabels = [str(nu) for nu in nuaTicks])
    pu.adjust_ylabels(ax, x_offset=-0.1)
    # # Save figure
    plt.savefig(fig_dir+'gs_error.'+format, dpi=dpi,transparent=True)

    plt.show()

#-----------------------------------------------------------------------------------------------------------------------
# Check on validity of MF approximation
#-----------------------------------------------------------------------------------------------------------------------
def cv_ux(nus_, u0_, omd_, omf_):
    u_inf = lambda nus, u0, omd, omf : u0*(omf+nus)/(omf+u0*nus)
    xi_u  = lambda nus, u0, omd, omf : (nus*(omf*(1-u0))**2)/((2*omf+u0*(2-u0)*nus)*(omf+nus)**2)
    xi_x  = lambda nus, u0, omd, omf : nus*(u_inf(nus, u0, omd, omf)**2)/(2*omd + nus*(u_inf(nus, u0, omd, omf)*(2-u_inf(nus, u0, omd, omf))))

    return np.sqrt(xi_u(nus_,u0_,omd_,omf_)*xi_x(nus_,u0_,omd_,omf_))

def cv_xg(nua_, ua_, oma_, js_, omg_):
    beta = np.exp(-js_*ua_)
    xi_xa = lambda nua, ua, oma : nua*(ua**2)/(2*oma + nua*(ua*(2-ua)))
    xi_gs = lambda nua, omg : (omg**2)/(omg+(1-beta)*nua)/(omg+(1+beta)*nua)

    return np.sqrt(xi_xa(nua_, ua_, oma_)*xi_gs(nua_,omg_))

def cv_check():
    # Generte substrates
    npts = 50
    nus = np.logspace(np.log10(0.1),np.log10(30),npts)
    nua = np.logspace(np.log10(0.01), np.log10(1), npts)

    maxcv = 0

    print "Checking Cauchy-Schwarz bound on synaptic dynamics: "
    ntest = 20
    u0 = np.logspace(-2,0,npts)
    omd = np.linspace(0.5,10,npts)
    omf = np.linspace(0.5,2.,npts)
    for i in xrange(ntest):
        for j in xrange(ntest):
            for k in xrange(ntest):
                ux_cv = cv_ux(nus,u0[i],omd[j],omf[k])
                # print np.amax(ux_cv)
                maxcv = np.amax(np.r_[maxcv,ux_cv])
    print "<",str(maxcv)

    print "Checking Cauchy-Schwarz bound on gtr dynamics: "
    maxcv = 0.0
    ntest = 30
    ua  = np.linspace(0.05,0.9,npts)
    oma = np.logspace(np.log10(0.6),np.log10(5),npts)
    js  = np.logspace(np.log10(0.0024),np.log10(630),npts)
    omg = np.logspace(np.log10(0.005),np.log10(0.17),npts)
    for i in xrange(ntest):
        for j in xrange(ntest):
            for k in xrange(ntest):
                for w in xrange(ntest):
                    xg_cv = cv_xg(nua,ua[i],oma[j],js[k],omg[w])
                    # print np.amax(xg_cv)
                    maxcv = np.amax(np.r_[maxcv,xg_cv])
    print "<",str(maxcv)

#-----------------------------------------------------------------------------------------------------------------------
# Generate figures
#-----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    data_dir = '../data/'
    fig_dir = '../Figures/'
    
    ## Figure 2
    model_compartments(sim=True, data_dir=data_dir, fig_dir=fig_dir, format='svg', dpi=600)
    # model_compartments(sim=False, data_dir=data_dir, fig_dir=fig_dir, format='svg', dpi=600)

    ## Figure 3
    full_model(sim=True, data_dir=data_dir, fig_dir=fig_dir, format='pdf', dpi=600)
    # full_model(sim=False, data_dir=data_dir, fig_dir=fig_dir, format='pdf', dpi=600)

    ## Figure 4
    syn_filtering(sim=True, data_dir=data_dir, fig_dir=fig_dir, format='svg', dpi=600)
    # syn_filtering(sim=False, data_dir=data_dir, fig_dir=fig_dir, format='svg', dpi=600)

    ## Figure 5
    pp_plasticity(sim=True, data_dir=data_dir, fig_dir=fig_dir, format='svg', dpi=600)
    # pp_plasticity(sim=False, data_dir=data_dir, fig_dir=fig_dir, format='svg', dpi=600)

    ## Figure 6
    mf_model(sim=True, data_dir=data_dir, fig_dir=fig_dir, format='pdf', dpi=600)
    # mf_model(sim=False, data_dir=data_dir, fig_dir=fig_dir, format='pdf', dpi=600)

    ## Figure 7
    mf_brp(sim=True, data_dir=data_dir, fig_dir=fig_dir, format='pdf', dpi=600)
    # mf_brp(sim=False, data_dir=data_dir, fig_dir=fig_dir, format='pdf', dpi=600)

    ## Figure 8
    pswitch_domains(dir=fig_dir,format='svg',dpi=1200)

    ## Figure 9
    gammas_error(fig_dir=fig_dir, format='eps', dpi=600)

    ## Check validity of MF
    cv_check()