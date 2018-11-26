"""
general_utils.py

Module with general purpose routines.
- checkscalar(val)         : Convert val to array if scalar otherwise leave it unchanged (equivalent to numpy.asarray)
- containsAny(str,set)     : Check whether 'str' contains ANY of the chars in 'set'
- containAll(str,set)      : Check whether 'str' contains ALL of the chars in 'set'
- extendval(val,maxl)      : Reshape array to maxl size, taking scalar value or the first element of the array.
- find(condx,*args)        : Python equivalent to matlab find()
- find_nearest(value,array,**kwargs) : Find array elements nearest to value
- permutation_indices(data): Return index of permutation in a sorted array (alternative to argsort())
- normalize_elements(array): Normalize elements of an array so that all elements are in [0,1]
- rescale(val,va_range)    : Rescale value within val_range (already vectorized)
- toarray(val)             : Modified np.asarray(val) that also converts scalars to array consistently
- unique_rows(a,col)       : Unique rows of a matrix 'a' sorted according to elements in 'col' column  
- varargin(pars,**kwargs)  : Handling of varargin for pars dictionary
"""

import numpy as np

def checkscalar(val):
    """
    Convert val into an array if it is a scalar, otherwise leave it as array.

    Maurizio De Pitta', The University of Chicago, September 6th, 2014.
    """
    if isscalar(val):
        return np.asarray(val)
    return val

def extendval(val,maxl):
    """
    Conveniently extend a value val to an ml-element array if val is scalar or
    an array of size <maxl. In the case that val is an array of size maxl, it considers only the first value and extends
    it to maxl.

    v2.0
    Recoded to also accept scalars.
    Maurizio De Pitta', The University of Chicago, August 11th, 2016.

    v1.0
    Maurizio De Pitta', The University of Chicago, September 6th, 2014.
    """
    if not hasattr(val,'__len__'):
        val_aux = val*np.ones(maxl)
    else:
        if val.size!=maxl:
            val_aux = val[0]*np.ones(maxl)
        else:
            val_aux = val
    return val_aux

def containsAny(str, set):
    """Check whether 'str' contains ANY of the chars in 'set'"""
    return 1 in [c in str for c in set]

def containsAll(str, set):
    """Check whether 'str' contains ALL of the chars in 'set'"""
    return 0 not in [c in str for c in set]

def permutation_indices(data):
    """Return index of permutation in a sorted array"""
    return sorted(range(len(data)), key = data.__getitem__)
     
def find(condx,*args):
    """Find indices of an array that satistfy some condition"""
    """find(x>10), find(x>10,1), find(x>10,-1)"""
    indices = np.array([i for i,elem in enumerate(condx) if elem])
    if len(args):
        if args[0]<0:
                indices = indices[args[0]:]
        if args[0]>0:
                indices = indices[:args[0]]
        if np.size(args)==0:
            print "zero indices requested"
            raise IndexError
                #print "zero indices requested"
    return indices

def find_nearest(value,array,which_nearest='min',index=False):
    """Return nearest array elements to a given value w/in min(array) and max(array)
    Note that if the array is nonmonotonic, considers only the first values
    from the beginning of the array.
    Opts: 
    - which_nearest : {'min'}|'max'|'minmax'
    - index         : {True}|False return index in the array not value
    """
    # Preliminary check that 
    if value<np.amin(array) or value>np.amax(array):
        print "WARNING: Value out of array range"
        return np.array([])
    # A couple of functions to retrieve the correct index in case value is one of the elements in the array    
    imin = lambda val : find(array>=val,0)[0] if array[find(array>=val,0)[0]]==val else find(array>=val,0)[0]-1    
    imax = lambda val : find(array<=val,0)[-1] if array[find(array<=val,0)[-1]]==val else find(array<=val,0)[-1]+1
#    idx_min = find(array>value,0)[0]-1  # first index
#    idx_max = find(array<value,0)[-1]+1 # last index
    idx_min = imin(value)
    idx_max = imax(value)
    if which_nearest=='min':
        if index:
            return idx_min
        else:
            return array[idx_min]
    if which_nearest=='max':
        if index:
            return idx_max
        else:
            return array[idx_min]
    if which_nearest=="minmax":
        if index:
            return idx_min,idx_max
        else:
            return array[idx_min],array[idx_max]

def merge_dicts(dict_0,*nargs):
    """
    Merge arbitrary number of dictionaries with dict_0

    Input arguments:
    - dict_0   : First dictionary to merge
    - *nargs   : list of further dictionaries to merge

    Return merged dictionary (with unique keys).
    Each key contains the last presentation.

    Maurizio De Pitta', The University of Chicago, Apr 27, 2016.
    """
    if np.size(nargs)==0 : return dict_0
    if np.size(nargs)==1 : return dict(dict_0, **nargs[0])
    return merge_dicts(dict_0,merge_dicts(nargs[0],*nargs[1:]))

def normalize_elements(array_):
    """
    Normalize elements of an array

    Input arguments:
    - array_ : Array

    Return:
    - normalized_array : Array
    """
    return (array_-(array_).min())/np.ptp(array_)

def transpose_dict(dict_):
    """
    Transpose (in the sense of an array operation) all elements of a dictionary (check that they are numpy array
    otherwise issue an error)

    Input arguments:
    - dict_  : dictionary

    Return:
    - dict_T : dictionary with transposed elements
    """
    dict_T = {}
    for key,val in dict_.iteritems():
        assert hasattr(val, "__len__"), "Dictionary elements must be numpy arrays"
        assert val.ndim==2, "Dictionary elements must be numpy 2D arrays"
        dict_T[key] = val.T
    return dict_T

def rescale(value,val_range):
    """Rescale value in a given range.
       If range is a point returns the point scaled to 0. The output is a np.array.
       This routine is already suitable to be used for array inputs.
    """
    # First verifies that val_range is an array
    if type(value)!=type(np.zeros(1)):
        value = np.array([value])         # Put scalar as array      
    if type(val_range)!=type(np.zeros(1)):
        val_range = np.asarray(val_range) # Make val_range array
    if np.size(np.shape(val_range))<2:
        val_range = val_range[None,:]     # If val_range is just one entry make it as an array 2x1
    rx = lambda x,r : (x-r[0])/np.diff(r) if np.diff(r)!=0 else x-r[0]
    return np.asarray([float(rx(value[i],val_range[i,:])) for i in xrange(np.shape(val_range)[0])])      

def toarray(val):
    """
    Modified np.asarray that accept also scalars and convert them into finite shape arrays.
    """
    if np.isscalar(val):
        return np.asarray([val])
    else:
        return np.asarray(val)

def unique_rows(a,col=0):
    """ Find unique rows of a matrix a and sort according to elements in column col
    """
    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    _, idx = np.unique(b, return_index=True)
    # a = a[idx]
    # Retrieve complementary index to be dumped 
    idxc = tuple(set(range(np.size(b)))-set(idx))
    a = np.delete(a,idxc,axis=0)
    if col==None:
        return a
    else:
        return a[a[:,col].argsort(),]

    
def varargin(pars,**kwargs):
    """
    varargin-like option for user-defined parameters in any function/module
    Use:
    pars = varargin(pars,**kwargs)
    
    Input:
    - pars     : the dictionary of parameters of the calling function
    - **kwargs : a dictionary of user-defined parameters

    Output:
    - pars     : modified dictionary of parameters to be used inside the calling
                 (parent) function
    
    Maurizio De Pitta', The University of Chicago, August 27th, 2014.
    """
    for key,val in kwargs.iteritems():
        if key in pars:
            pars[key] = val
    return pars
    
if __name__ == "__main__":
    
#    ###########################################################################
#    # Test of find_nearest()    
#    ###########################################################################    
#    x = arange(10)
#    val = [5,7,8.9]
#    f = np.vectorize(find_nearest,excluded=[1,2,3])
#    v = f(val,x,which_nearest='minmax',index=False)
##    v=  np.array(zip(*v),dtype=tuple)
##    print type(v),shape(v)

    ###########################################################################
    # Test of rescale()    
    ###########################################################################
    # print rescale(5,[0.1,10.0])

    ###########################################################################
    # Test of merge_dicts()
    ###########################################################################
    a = {'0':0, '1':1}
    b = {'2':2, '3':3}
    c = {'4':4, '5':5}
    d = {'6':6, '7':7}

    print merge_dicts(a)
    print merge_dicts(a,b)
    print merge_dicts(a,b,c)
    print merge_dicts(a,b,c,d)