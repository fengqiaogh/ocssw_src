'''
Created on Feb 9, 2016

@author: rhealy
retrieved from https://github.com/marcelhaas/python/blob/master/value_locate.py

'''
import numpy as np


def value_locate(refx, x):
    """
    VALUE_LOCATE locates the positions of given values within a
    reference array.  The reference array need not be regularly
    spaced.  This is useful for various searching, sorting and
    interpolation algorithms.
    The reference array should be a monotonically increasing or
    decreasing list of values which partition the real numbers.  A
    reference array of NBINS numbers partitions the real number line
    into NBINS+1 regions, like so:
    REF:           X[0]         X[1]   X[2] X[3]     X[NBINS-1]
        <----------|-------------|------|---|----...---|--------------->
    INDICES:  -1           0          1    2       3        NBINS-1
    VALUE_LOCATE returns which partition each of the VALUES falls
    into, according to the figure above.  For example, a value between
    X[1] and X[2] would return a value of 1.  Values below X[0] return
    -1, and above X[NBINS-1] return NBINS-1.  Thus, besides the value
    of -1, the returned INDICES refer to the nearest reference value
    to the left of the requested value.
    
    Example:
    >>> refx = [2, 4, 6, 8, 10]
    >>> x = [-1, 1, 2, 3, 5, 5, 5, 8, 12, 30]
    >>> print value_locate(refx, x)
    array([-1, -1,  0,  0,  1,  1,  1,  3,  4,  4])
    
    
    This implementation is likely no the most efficient one, as there is
    a loop over all x, which will in practice be long. As long as x is
    shorter than 1e6 or so elements, it should still be fast (~sec).
    
    
    """
    print ("TODO: check if refx is monotonically increasing.")
    
    refx = np.array(refx)
    x = np.array(x)
    loc = np.zeros(len(x), dtype='int')
    
    for i in range(len(x)):
        ix = x[i]
        ind = ((refx - ix) <= 0).nonzero()[0]
        if len(ind) == 0:
            loc[i] = -1
        else: loc[i] = ind[-1]
    
    return loc
    


if __name__ =="__main__":
    # test case
    refx = [2, 4, 6, 8, 10]
    x = [-1, 1, 2, 3, 5, 5, 5, 8, 12, 30]
    
    res = value_locate(refx, x)
    assert list(res) == [-1, -1,  0,  0,  1,  1,  1,  3,  4,  4]
    print ("Test(s) passed!")
    
    x= np.random.random(1e6)*20
    res = value_locate(refx, x)
    print ("Done!")