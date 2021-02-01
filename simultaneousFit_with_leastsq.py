#This function simultaneously fits multiple data sets with various
# functions.  The functions are allowed to share some parameters

#The technique is a bit odd, but it seems to work.  Here's the description:
#The user passes in a set of 5 arrays of the same length and one array of any length.
# The length is the number of data sets to which we are fitting.
#The first argument is an array of arrays.
# Each ith array is the independent data (the 'x values') for data set i .
#The second argument is an array of arrays.
# Each ith array is the dependent data (the 'y values') for data set i .
#The third argument is an array of functions.
# Each ith function is the expected relation between the ith x values and ith y values.
# The other parameters that go into each function are what are fit over.
#The fourth argument is an array of arrays.
# Each ith array is the set of indeces that tell the user which fitted parameters go into
# function i.
# The order of parameters in the ith array describes the order that the parameters will be
# entered into the ith function.  So one should be careful.
#The fifth argument is an array that provides the best guess of the parameters.
# The length of this array is the number of free fit parameters.
#The sixth argument is an array of arrays.
# Each ith array is the dependent data uncertainties (the 'y error values') for data set i.

#For example, say we want to fit two lines with different slopes but the same intercept.
# We have two independent data sets: x_1 = [an array of x values],  x_2 = [an array of x values]
# We have two dependent data sets: y_1 = [an array of y values],  y_2 = [an array of y values]
# We suspect: y_1 = a_1 * x_1 + b_1
#             y_2 = a_2 * x_2 + b_2
#             b_1 = b_2 = b
# Let's assume we have the hypothesis that a_1 = 1.0, a_2 = 2.0, b = 0.0
# So we would run:
# getSimultaneousFitParams([x_1,x_2], [y_1,y_2], [lambda x, a, b: a * x + b, lambda x, a, b: a * x + b], [[1,0], [2,0]],  [0.0, 1.0, 2.0])
# Note that by not specifying y_errs, the function just weights everything equally
# The order of the parameters are [b, a_1, a_2].  Note that the index of 0 (which specifies the b parameter) is shared between the parameters.
# Also note that since the inputs into the function for data set 1 are [a_1, b], which have indeces [1, 0],
# that is the order in which the parameters must be given.  Ditto for parameter set 2. 

import math
import numpy as np
import scipy.optimize

def err(function, x_data, y_data, y_err):
    #y_exp = function(np.array(x_data))
    print 'entering err function' 
    #for x_elem in x_data:
    #    print 'x_elem = ' + str(x_elem)
    #    print 'function(x_elem) = ' + str(function(x_elem))
                
    y_exp = np.array([function(x_elem) for x_elem in x_data])
    #print 'managed to enter y_exp' 
    diff = y_data - y_exp
    normalized_diff = diff / y_err
    return normalized_diff 


def err_global(function_parameters, functions, param_indeces, x_data_arrays, y_data_arrays, y_err_arrays):

    #print 'len(x_data_arrays) = ' + str(len(x_data_arrays))
    print 'function_parameters = ' + str(function_parameters) 

    #print 'param_indeces = '
    #print param_indeces 
    param_loaded_functions = [lambda x: functions[i](x, *[ function_parameters[index] for index in param_indeces[i] ]) for i in range(len(x_data_arrays))]
    #print 'param_loaded_functions = '
    #print param_loaded_functions 
    errs = [ err(param_loaded_functions[i], x_data_arrays[i], y_data_arrays[i], y_err_arrays[i]) for i in range(len(x_data_arrays)) ]
    #print 'errs = '
    #print errs
    return np.concatenate(errs)

def getSimultaneousFitParams(x_data_arrays, y_data_arrays, functions, param_indeces,  guess_parameters, y_err_arrays = None ):

    if y_err_arrays is None or y_err_arrays is 'none':
        y_err_arrays = [ np.zeros(len(y_data_array)) + 1.0 for y_data_array in y_data_arrays ]
    p_best, ier = scipy.optimize.leastsq(err_global, guess_parameters, args = (functions, param_indeces, x_data_arrays, y_data_arrays, y_err_arrays))

    return (p_best, ier) 

    
    
