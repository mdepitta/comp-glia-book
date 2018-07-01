"""
save_utils.py

Module which contains several saving utilities to be used in python.
Routines:
- savedata (pickle) : save objects on file (by pickle)
- loaddata (pickle) : load objects from file (by pickle)

Maurizio De Pitta', The University of Chicago, August 18th, 2014.
"""

# import cPickle as pickle
import dill as pickle

def savedata(data,filename='./data.pkl'):
    """
    Save list of variables in 'data' in file 'filename
    Use:
    savedata(data, filename='./data.pkl')
    
    Input:
    - data     : [var1,var2,...] list of variables names (NO strings! but real names)
    - filename : string for file name (usually with extension 'pkl') 
    
    Maurizio De Pitta', The University of Chicago, August 18th, 2014.
    """
    with open(filename, "wb") as f:
        pickle.dump(len(data), f)
        for value in data:
            pickle.dump(value, f)

def loaddata(filename,verbose=False):
    """
    load list of variables in 'data' in file 'filename
    Use:
    data = loaddata(filename)
    
    Input:
    - filename : string for file name (usually with extension 'pkl')
    
    Output:
    - data : list of objects
    
    Maurizio De Pitta', The University of Chicago, August 18th, 2014.
    """

    data = []
    with open(filename, "rb") as f:
        for _ in range(pickle.load(f)):
            data.append(pickle.load(f))
    if verbose : print "loaded: ", data
    return data