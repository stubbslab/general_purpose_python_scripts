#Returns a list of specified length that goes
# from given start to (and including) given end
import numpy as np

def getListOfLength(start, end, length):
    if length == 1:
        print 'Cannot generate a list that covers the range of length 1.  Returning start.'
        return np.array(start)
    else:
        #print 'end = ' + str(end)
        #print 'start = ' + str(start)
        #print 'length = ' + str(length)
        step_size = (end - start) / (length - 1)
        #print np.arange(start, end + step_size / 2.0 , step_size)
        return np.arange(start, end + step_size / 2.0 , step_size)

