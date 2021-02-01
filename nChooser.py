#Returns the value of n choose r, relatively efficiently ( better than n!/((n-r)!*r!) )

import operator as op

def nChooser(n, r):
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, xrange(n, n-r, -1))
    denom = reduce(op.mul, xrange(1, r+1))
    return numer//denom