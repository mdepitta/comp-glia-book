"""
spk_generators.py

Module for the generation of spike trains of given statistics.
Each generator produce a spike train 2xN [t1,1,...t1,N1,t2,1,...,t2,N2,...,tN,1,...,tN,NN]

Contents:
- spk_poisson  : Poisson spike train generator (also with refractoriness)
- spk_periodic : Periodic spike train generator (also with refractoriness)

Other routines:
- compute_frequency : compute frequency of a spike train
- compute_averr     : compute average released resources
- isi_stat          : compute mean, std and cv of ISIs

Maurizio De Pitta', The University of Chicago, September 1st, 2014.
"""

from __future__ import division
from numpy import *
import numpy.matlib as matlib
from scipy import *
import os, sys
base_dir = '/Ch10.DePitta/code'
sys.path.append(os.path.join(os.path.expanduser('~'),base_dir))
import pycustommodules.general_utils as gu

import matplotlib.pylab as plt

def cut_isi(isi,tmax):
    """
    Drop spike times above tmax
    
    Maurizio De Pitta', The University of Chicago, September 10th, 2014.
    """
    tspk = cumsum(isi)
    return tspk[tspk<=tmax]

def make_current(val,tspk):
    """
    Build a 3-row array: 1st row: spike instants; 2nd row: values; 3rd row: source neuron

    Maurizio De Pitta', The University of Chicago, September 10th, 2014.    
    """
    if (isscalar(val) or (val.size==1)):
        vals = val*ones(tspk[0].size)
    else:
        vals = val
    return vstack((tspk[0],vals,tspk[1]))
    
def compute_frequency(tspk,T):
    """
    Compute firing frequency within time period T, for spike train instants tspk

    Maurizio De Pitta', The University of Chicago, September 10th, 2014.    
    """
    isi = diff(tspk)
    isi = isi[isi>0]
    if not size(isi):
        return 0
    else:
        nspk = size(isi)+1
        return nspk/double(T)

def compute_averr(rr,N,T):
    """
    Compute average released resources out of a rr train of released resources
    for N synapses in the time window T

    Maurizio De Pitta', The University of Chicago, September 29th, 2014.    
    """
    
    return sum(rr)/(N*double(T))
    
def isi_stat(tspk,correction=False):
    """
    Compute mean, std and CV of the ISI of a spike train specified in tspk.

    v1.1 - Added correction option
    Maurizio De Pitta', The University of Chicago, August 24, 2016.

    v1.0
    Maurizio De Pitta', The University of Chicago, October 2nd, 2014.    
    """

    isi = diff(tspk)
    isi = isi[isi>0]
    if not size(isi):
        return 0,0,0
    else:
        m = mean(isi)
        s = std(isi)
        if correction and (s==0.):
            s=m # when there is only one spike you may legitimately assume that std(ISI)>=ISI
        cv = s/m
        return m,s,cv    

def isi_dist(tspk,bins=10,range=None):
    """
    Compute ISI distribution of spike train specified in tspk.
    It is essentially a call to np.histogram.

    Input arguments:
    tspk   : array with spike times as produced by the simulator
    bins   : number of bins or sequence of bin edges (as in np.histogram)
    range  : {None} | [min_isi,max_isi] interval for consider min/max to consider

    Return:
    pd        : p.d.f. values
    bin_edges : bin_edges for pdf

    Maurizio De Pitta', The University of Chicago, April 22, 2016.
    """

    isi = diff(tspk)
    pd,bin_edges = histogram(isi,bins=bins,range=range,density=True)
    return pd,bin_edges

def thin_isi(isi_homog,tau_r,N_spikes,rate,T_rate):
    """
    Thinning of homogeneous Poisson process to generate inhomogeneous Poisson process.

    Input arguments:
    - isi_homog  : Array-like of ISIs of Poisson process at max(rate)
    - tau_r      : Float. Refractory period.
    - N_spikes   : Int. Number of spikes used to compute ISIs.
    - rate       : Array-like of floats for rates (time ordered)
    - T_rate     : Array-like of floats for durations intervals of each rate value w.r.t. preceding one.

    NOTE: This method can be used also to generate inhomogenous Poisson processes for ANY time-varying rate.
    As far as rate and T_rate are properly given.

    v3.0
    Corrected handling of absolute refractory (meant as case of inhomogenous Poisson process where the rate is zero for
    tau_r right after each spike. Improved array-based coding.
    Maurizio De Pitta', Basque Center for Applied Mathematics, August 19, 2018.

    v2.0
    Added handling of multiple rates and non-homogeneous Poisson processes.
    Maurizio De Pitta', The University of Chicago, August 27, 2016.

    v1.0
    Maurizio De Pitta', The University of Chicago, September 9th, 2014.
    """

    tspk_homog = cumsum(isi_homog)  # Spike sequence at maximum rate

    x = random.rand(N_spikes)
    rate_ = rate[0]*ones(N_spikes)
    if tau_r>0:
        index = logical_and(tspk_homog<T_rate[0],isi_homog<tau_r)
        rate_[index] = 0.0
    if len(rate)>1:
        for i,t in enumerate(T_rate[:-1]):
            rate_[tspk_homog>=t] = rate[i+1]
            if tau_r > 0:
                index = logical_and(logical_and(tspk_homog>=t,tspk_homog<T_rate[i+1]),isi_homog<tau_r)
                rate_[index] = 0.0
    z = rate_/amax(rate_)
    tspk_ = tspk_homog[less(x,z)]
    isi = diff(tspk_)                      # Difference between successive spikes
    return isi

def spk_poisson(T=10,N=1,**kwargs):
    """
    Generate Poisson-distributed spike train.

    Use:
    tspk = spk_poisson(T=10,N=1,**kwargs)

    Input arguments:
    - T        : period [0,T] for the spike train
    - N        : number of neurons
    - **kwargs :
        - rate : average rate of Poisson spike [Hz]
        - trp  : refractory period [s]

    Output:
    - tspk     : spike train instants (sorted) with corresponding neuron

    v1.1
    Corrected minor bug in handling tau_r and rate array.
    Maurizio De Pitta', Basque Center of Applied Mathematics, August 19, 2018.

    v1.0
    Maurizio De Pitta', The University of Chicago, September 9th, 2014.
    """
    pars = {'rate': 100, 'trp': 2e-3}
    pars = gu.varargin(pars,**kwargs)
    # Check that T and rate are of the same size
    assert size(pars['rate']) == size(T),"Size of stimulus ending times (T) (" + str(size(T)) + ") does not match size of stimulus rates (" + str(size(pars['rate'])) + ")"
    # The following assures that rate and T are Numpy arrays
    if hasattr(pars['rate'],'__len__'):
        pars['rate'] = asarray(pars['rate'])
        T = asarray(T)
    else:
        pars['rate'] = asarray([pars['rate']])
        T = asarray([T])

    # Compute N_spikes : allows for some fluctuations in the number of spikes in order to avoid border effects
    N_spikes = int(max(sum(T)*(1.+random.random())*amax(pars['rate']),2))  # Average number of spikes per trial >=2 due to thinning
    spk_train,indexes = empty(0),empty(0)
    for i in xrange(N):
        # For an homogeneous Poisson process the following works. For non-homogeneous processes then it will prepare it
        # for thinning
        isi_homog = -log(random.rand(N_spikes))/amax(pars['rate'])  # Random ISIs at maximum rate
        if (pars['trp']>0)or(size(pars['rate'])>1):
            tspk = cut_isi(thin_isi(isi_homog,pars['trp'],N_spikes,pars['rate'],cumsum(T)),sum(T))
        else:
            tspk = cut_isi(isi_homog,sum(T))
        spk_train = concatenate((spk_train,tspk))
        indexes = concatenate((indexes,i*ones(tspk.size)))
    # Sort spikes
    si = spk_train.argsort()
    # Stack and sort
    spk_train = vstack((spk_train[si],indexes[si]))
    return spk_train

def spk_poisson_tvrate(rate_fun,T=10.0,N=1,mode='spike',**kwargs):
    """
    Generate Poisson spike train for arbitrary time-varying rate.

    Inputs:
    rate_fun   : Non-negative rate function in the form of f(x)
    - T        : Period [0,T] for the spike train
    - N        : Number of neurons
    mode       : {'spike'} | 'interval' If 'spike', recompute every spike; if 'interval' rate is compute on the whole
                 [0,T] interval
    kwargs     :

    Output:
    - spk_train: spike train instants (sorted) with corresponding neuron

    Maurizio De Pitta', Basque Center of Applied Mathematics, June 11th, 2018.
    """

    pars = {'npts': 100, 'dt': 1e-3, 'eps': 5e-3} # number of points to evaluate rate
    pars = gu.varargin(pars, **kwargs)
    if 'dt' in kwargs.keys():
        pars['npts'] = T/kwargs['dt']
    if mode=='interval':
        time = linspace(0.,T,int(pars['npts']))
        rate = rate_fun(time)
        dt = time[1]-time[0]
        sp_index = less(random.rand(N,pars['npts']),tile((rate*dt),(N,1)))
        spk_train = tile(time,(N,1))[sp_index]
        indexes = (tile(arange(0, N),(pars['npts'], 1))).T[sp_index]
        si = spk_train.argsort()
        # Stack and sort
        spk_train = vstack((spk_train[si], indexes[si]))
    else:
        # mode=='spike'
        spk_train, indexes = zeros(1), zeros(1)
        for i in xrange(N):
            t = 0
            isi = []
            xscale = float(T - t)
            yscale = amax(rate_fun(linspace(t, T, pars['npts'])))
            while ((T-t)>pars['eps']):
                tval = xscale*random.rand(1)
                yval = yscale*random.rand(1)
                if yval<rate_fun(tval):
                    isi.append(tval)
                    t = t + tval
                    xscale = float(T - t)
                    yscale = amax(rate_fun(linspace(t, T, pars['npts'])))
            # Generate final spike train
            tspk = cut_isi(isi, T)
            spk_train = concatenate((spk_train, r_[0,tspk]))
            indexes = concatenate((indexes, i * ones(len(isi)+1)))
        spk_train = vstack((spk_train,indexes))

    return spk_train

def spk_periodic(T=10,N=1,**kwargs):
    """
    Generate periodic spike train


    Use:
    tspk = spk_periodic(T=10,N=1,**kwargs)
    
    Input arguments:
    - T        : period [0,T] for the spike train
    - N        : number of neurons
    - **kwargs : 
        - rate : average rate of Poisson spike [Hz]
        - trp  : refractory period [s]
    
    Output:
    - tspk     : spike train instants (sorted) with corresponding neuron
 
    Maurizio De Pitta', The University of Chicago, September 23rd, 2014.        
    """
    pars = {'rate': 100, 'trp': 2e-3}
    pars = gu.varargin(pars,**kwargs)
    N_spikes = int(T * pars['rate'])  # Average number of spikes per trial
    spk_train,indexes = empty(0),empty(0)
    for i in xrange(N):
        isi = cut_isi(ones(N_spikes)*(1/pars['rate']+pars['trp']),T)  # Periodic ISI at effective rate
        spk_train = concatenate((spk_train,isi))
        indexes = concatenate((indexes,i*ones(isi.size)))
    # Sort spikes
    si = spk_train.argsort()
    # Stack and sort
    spk_train = vstack((spk_train[si],indexes[si]))
    return spk_train

def input_spikes(N, tfin, rate, trp, stimulus='poisson', spikes_pre=None):
    update_time = lambda spk_vect,t0 : vstack((spk_vect[0,:] + t0,spk_vect[1,:]))
    if stimulus=='poisson':
        spk_train = spk_poisson(tfin,N,rate=rate,trp=trp)
    elif stimulus=='poisson_tv': # It is added here for completeness, but actually better to use the standalone routine
        # In this case rate is either a function
        spk_train = spk_poisson_tvrate(rate,T=tfin,N=N,mode='interval')
    elif stimulus=='periodic':
        # Check whether rate is given as scalar or as an array. In the latter case concatenate solutions
        if hasattr(rate, '__len__') and (not isinstance(rate, str)):
            if hasattr(tfin,'__len__'):
                # In this case tfin is a vector whose elements represent duration of each rate stimulus
                assert size(rate)==size(tfin),"Size of stimulus intervals (tfin) ("+str(size(tfin))+") does not match size of stimulus rates ("+str(size(rate))+")"
                t0 = cumsum(r_[0.0,tfin[:-1]]) # Initial instants to add to timings
                spk_train = concatenate([update_time(spk_periodic(tfin[i],N,rate=rate[i],trp=trp),t0[i]) for i in xrange(rate.size)],axis=1)
            else:
                spk_train = concatenate([update_time(spk_periodic(tfin/rate.size,N,rate=rate[i],trp=trp),i*tfin/rate.size) for i in xrange(rate.size)],axis=1)
        else:
            spk_train = spk_periodic(tfin,N,rate=rate,trp=trp)
    elif stimulus=='fixed':
        if size(shape(spikes_pre))<2:
            # Only a series of spikes is given and it is assumed to be for all the inputs
            spk_train = vstack((tile(sort(spikes_pre),(1,N))[0],tile(arange(N),(size(spikes_pre),1)).reshape((1,N*size(spikes_pre)),order='F')[0]))
        else:
            # The input as produced elsewhere is just reused (it must be a (2xM) array wth the second row having synaptic indexes
            # Issue a warning --> This will suggest where to look at in case of SEG FAULTS
            if max(spikes_pre[1])>N : print "WARNING: expected ",N," inputs but found ",max(spikes_pre[1])
            spk_train = spikes_pre
    return spk_train

#-----------------------------------------------------------------------------------------------------------------------
# Testing routines
#-----------------------------------------------------------------------------------------------------------------------
def sin_rate(x):
    return 10.*(sin(0.5*x))**2

def exp_rate(x):
    return 5*(1-exp(-2.0*x))

if __name__ == "__main__":
    # Verify Poisson distribution with time varying rate
    spk_train = spk_poisson_tvrate(sin_rate,T=10,N=10,mode='spike')
    plt.plot(spk_train[0],ones(shape(spk_train)[1]),'o')

    # # Verify Poisson distribution
    # N_trials = 1000
    # N_bins = 100
    # bins = linspace(50.,150.,N_bins)
    # h = zeros((N_trials,N_bins))
    #
    # rate = 100.
    # tfin = 1.
    # spk = input_spikes(N_trials,tfin,rate,0.,stimulus='poisson')
    # N_spk = [sum(spk[1]==i) for i in xrange(N_trials)]
    # f,b = histogram(N_spk,bins)
    #
    # plt.plot(b[:-1],f/sum(f,dtype=float),'ko')

    plt.show()