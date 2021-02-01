#Define a function to return a list of numbers spaced evenly on a logarithmic scale between the start and stop.
#The method is to define a regular array of powers, and then exponentiate.
#So if you take the log of the list, you get regular spacing.
#Note, start and stop must both be positive.  

import numpy as np
import math 

def logList(start,stop,n_elems):
    ratio=float(start)/float(stop)
    base=ratio**(1.0/(n_elems-1))
    powers = np.arange(0,n_elems)
    return np.sort(stop * base**powers)
