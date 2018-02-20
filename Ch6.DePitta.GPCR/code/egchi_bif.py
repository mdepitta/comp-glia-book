import PyDSTool as dst
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

# Extended ChI model
import sys, os

sys.path.append(os.path.join(os.path.expanduser('~'), 'Dropbox/Ongoing.Projects'))
import pycustommodules.general_utils as gu
import pycustommodules.save_utils as svu

def gchidp_pars(**kwargs):
    pars = {'yrel': 0, 'Op': 0.3, 'OmegaP': 1.8,
            'vbias': 0, 'vbeta': 3, 'vdelta': 0.5, 'kappad': 1, 'Kdelta': 0.5,
            'v3k': 2, 'Kd': 0.5, 'K3': 1, 'r5p': 0.1,
            'vkd': 0.28, 'Kkc': 0.5, 'OmegaKD': 0.33, 'vk': 1.0,
            'vd' : 1.0, 'Kdc': 0.3, 'Kdd': 0.005, 'OmegaD': 0.1,
            'd1' : 0.1, 'd2': 2.1, 'd3': 0.9967, 'd5': 0.2, 'a2': 0.4, 'c1': 0.4, 'c0': 4, 'rc': 7, 'rl': 0.05,
            'ver': 0.9, 'KER': 0.1}
    # Custom parameters
    pars = gu.varargin(pars, **kwargs)
    return pars


def gchidp_model(pars):
    """Creates PyDSTool DSargs object for the 'ChI' model vector field
    """
    auxfuncs = {'hill'  : (['x', 'k', 'n'], 'pow(x,n)/(pow(x,n)+pow(k,n))'),
                'minf'  : (['ip3'], 'hill(ip3,d1,1)'),
                'ninf'  : (['ca'], 'hill(ca,d5,1)'),
                'Jchan' : (['ca', 'h', 'ip3'], 'rc*pow(minf(ip3)*ninf(ca)*h,3)*(c0-(1+c1)*ca)'),
                'Jleak' : (['ca'], 'rl*(c0-(1+c1)*ca)'),
                'Jpump' : (['ca'], 'ver*hill(ca,KER,2)'),
                'vglu'  : (['gammaa'], 'vbeta*gammaa'),
                'vplcd' : (['ca', 'ip3'], 'vdelta/(1+ip3/kappad)*hill(ca,Kdelta,2)'),
                'v3keff': (['ca', 'ip3'], 'v3k*hill(ca,Kd,4)*hill(ip3,K3,1)'),
                'v5peff': (['ip3'], 'r5p*ip3'),
                'vpkc'  : (['ca', 'dag'], 'vkd*dag*hill(ca,Kkc,1)'),
                'vdagk' : (['ca', 'dag'], 'vd*hill(ca,Kdc,2)*hill(dag,Kdd,2)')}

    rhs = {'gammaa': 'Op*(1-gammaa)*yrel-OmegaP*(1+vk*pkc/OmegaP)*gammaa',
           'ip3'   : 'vbias+vglu(gammaa)+vplcd(ca,ip3)-v3keff(ca,ip3)-v5peff(ip3)',
           'ca'    : 'Jchan(ca,h,ip3)-Jpump(ca)+Jleak(ca)',
           'h'     : '(a2*d2*(ip3+d1)/(ip3+d3))*(1-h)-a2*ca*h',
           'dag'   : 'vbias+vglu(gammaa)+vplcd(ca,ip3)-vpkc(ca,dag)-vdagk(ca,dag)-OmegaD*dag',
           'pkc'   : 'vpkc(ca,dag)-OmegaKD*pkc'}
    ICs = {'gammaa': 0, 'ip3': 0, 'ca': 0, 'h': 0.9, 'dag': 0, 'pkc': 0}

    DSargs = dst.args(name='GChIDP')
    DSargs.pars = pars
    DSargs.fnspecs = auxfuncs
    DSargs.varspecs = rhs
    DSargs.ics = ICs

    return DSargs

def compute_period(vkd=0.28,vk=1.0):
    """
    Compute period of oscillations as a function of yrel for given value of O_K

    Return dictionary: {'yrel': ndarray, 'T': ndarray}
    """
    # First check / Default set
    po_data = svu.loaddata('../data/fit_po.pkl')[0]
    lr_data = svu.loaddata('../data/fit_lra.pkl')[0]
    chi_data = svu.loaddata('../data/chi_fit_0_bee.pkl')[0]

    # yrel = 1.45
    # yrel = 1.4715  # vk=0.5
    # yrel = 1.448     # vk=0.2
    # yrel   = 2.0
    yrel   = 1.3
    c0     = 5.0
    rc     = 0.8*lr_data[0]
    ver    = lr_data[1]
    vbeta  = 1.0
    vdelta = chi_data[1]
    v3k    = 1.4*chi_data[2]
    r5p    = chi_data[3]

    # DAG metabolism parameters
    OmegaKD = 0.33
    vd      = 0.45
    OmegaD  = 0.26
    pars = gchidp_pars(d1=po_data[0], d2=po_data[1], d3=po_data[0], d5=po_data[2], a2=po_data[3],
                    c1=0.5, rl=0.01, KER=0.1, rc=rc, ver=ver, Kkc=0.5, Kkd=0.1, Kdd=0.1,
                    c0=c0,yrel=yrel,vbeta=vbeta,vdelta=vdelta,v3k=v3k,r5p=r5p,
                    vkd=vkd,vk=vk,OmegaKD=OmegaKD,vd=vd,OmegaD=OmegaD)
    # Generate ODE system
    DSargs = gchidp_model(pars)
    DSargs.tdomain = [0,20]                         # set the range of integration
    ode  = dst.Generator.Vode_ODEsystem(DSargs)     # an instance of the 'Generator' class
    traj = ode.compute('steady')                    # Integrate ODE
    sol = traj.sample()                             # Retrieve solution
    # Intermediate plotting
    # plt.plot(sol['t'],sol['gammaa'],'k-')
    # plt.plot(sol['t'],sol['ca'],'k-')
    # plt.plot(sol['t'],sol['dag'],'r-')
    # plt.plot(sol['t'],sol['pkc'],'b-')

    # Prepare the system to start close to a steady state
    p0 = traj(DSargs.tdomain[1])
    ode.set(pars = {'yrel': yrel} )       # Lower bound of the control parameter 'i'
    ode.set(ics =  p0)                 # Start from steady-state

    # Generate PyCont object
    PC = dst.ContClass(ode)            # Set up continuation class
    # EP-C
    PCargs = dst.args(name='EQ1', type='EP-C', force=True)       # 'EP-C' stands for Equilibrium Point Curve. The branch will be labeled 'EQ1'.
    PCargs.freepars     = ['yrel']                   # control parameter(s) (it should be among those specified in DSargs.pars)
    PCargs.MaxNumPoints = 300                        # The following 3 parameters are set after trial-and-error
    PCargs.MaxStepSize  = 0.01
    PCargs.MinStepSize  = 1e-4
    PCargs.StepSize     = 1e-3
    PCargs.VarTol       = 1e-5
    PCargs.FuncTol      = 1e-5
    PCargs.TestTol      = 1e-5
    PCargs.LocBifPoints = 'all'                      # detect limit points / saddle-node bifurcations
    PCargs.SaveEigen    = True                       # to tell unstable from stable branches
    PC.newCurve(PCargs)
    # Compute equilibria
    print "Computing EQ-C"
    PC['EQ1'].forward()
    # PC['EQ1'].backward()
    print "done!"

    svu.savedata([PC],'tmp_cont.pkl')
    PC = svu.loaddata('tmp_cont.pkl')[0]

    # LC-C
    PCargs = dst.args(name='LC1', type='LC-C', force=True)
    PCargs.initpoint    = 'EQ1:H1'
    PCargs.freepars = ['yrel']
    PCargs.MaxNumPoints = 200
    PCargs.MaxStepSize  = 0.05
    PCargs.MinStepSize  = 0.01
    PCargs.StepSize     = 1e-3
    # PCargs.VarTol       = 1e-6
    # PCargs.FuncTol      = 1e-4
    # PCargs.TestTol      = 1e-5
    PCargs.LocBifPoints = 'all'

    #PCargs.NumSPOut = 50;
    #PCargs.SolutionMeasures = ['max','min']
    PCargs.SolutionMeasures = 'all'
    PCargs.SaveEigen = True
    PC.newCurve(PCargs)
    # Compute Orbits
    print "Computing LC-C"
    PC['LC1'].forward()
    # PC['LC1'].backward()
    print "done!"

    # Retrieve solution
    T = PC['LC1'].sol['_T']
    yrel = PC['LC1'].sol['yrel']
    # Retrieve stabled and unstable branches
    idx = np.asarray([PC['LC1'].sol.labels['LC'][i]['stab'] for i in PC['LC1'].sol.labels['LC']], 'S')
    plt.plot(yrel[idx == 'S'],T[idx == 'S'],'ko')
    plt.show()
    return {'yrel': yrel[idx == 'S'], 'T': T[idx == 'S']}

if __name__ == "__main__":
    # Effective computation
    pkc0 = compute_period(vkd=0.001,vk=0.001)
    svu.savedata([pkc0], '../data/pkc_period_Ok.pkl')
    pkc1 = compute_period(vkd=0.28,vk=1.0)
    svu.savedata([pkc0,pkc1], '../data/pkc_period_Ok.pkl')
    pkc2 = compute_period(vkd=0.28,vk=3.0)
    svu.savedata([pkc0,pkc1,pkc2],'../data/pkc_period_Ok.pkl')

    pkc0, pkc1, pkc2 = svu.loaddata('../data/pkc_period_Ok.pkl')
    pkc3 = compute_period(vkd=0.84,vk=1.0)
    svu.savedata([pkc0, pkc1, pkc2, pkc3], '../data/pkc_period_Ok.pkl')
    pkc4 = compute_period(vkd=0.84,vk=3.0)
    svu.savedata([pkc0, pkc1, pkc2, pkc3, pkc4], '../data/pkc_period_Ok.pkl')