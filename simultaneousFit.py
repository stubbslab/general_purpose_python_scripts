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
import itertools

def err(function, x_data, y_data, y_err):
    #y_exp = function(np.array(x_data))
    #print 'entering err function' 
    #for x_elem in x_data:
    #    print 'x_elem = ' + str(x_elem)
    #    print 'function(x_elem) = ' + str(function(x_elem))
    
    #y_exp = np.array([function(x_elem) for x_elem in x_data])
    y_exp = np.array(function(x_data))

    diff = y_data - y_exp
    #print 'diff = '
    #print diff
    normalized_diff = diff / y_err
    return normalized_diff 


def err_global(function_parameters, functions, param_indeces, x_data_arrays, y_data_arrays, y_err_arrays):

    #print 'len(x_data_arrays) = ' + str(len(x_data_arrays))
    #print 'function_parameters = ' + str(function_parameters) 

    #print 'param_indeces = '
    #print param_indeces
    #print functions
    #print param_indeces 
    param_loaded_functions = [lambda x: functions[i](x, *[ function_parameters[index] for index in param_indeces[i] ]) for i in range(len(x_data_arrays))]
    #print 'param_loaded_functions = '
    #print param_loaded_functions
    errs = [ err(param_loaded_functions[i], x_data_arrays[i], y_data_arrays[i], y_err_arrays[i]) for i in range(len(x_data_arrays)) ]
    #print 'sum of squares of err_global for these function parameters are = ' + str(sum([elem ** 2.0 for elem in np.concatenate(errs)])) 
    return np.concatenate(errs)

#bounds must not be infinite
def simpleLeastSquareFinder(error_function, bounds, n_sample_points, args = (), space_sampling_method = 'random', ): 
    #if type(n_sample_points) is float or type(n_sample_points) is int:
    #    n_sample_points = [n_sample_points for bound in bounds[0]]
    #print 'n_sample_points = ' + str(n_sample_points)
    #for param_range in param_ranges:
    #    print 'len(param_range) = ' + str(len(param_range))
    #    print 'param_range = ' + str(param_range)
    #print 'bounds = ' + str(bounds)
    #print 'n_sample_points = ' + str(n_sample_points)
    space_sampling_method = space_sampling_method.lower()
    all_combos = []
    if space_sampling_method == 'grid':
        param_ranges = [np.arange(bounds[0][i], bounds[1][i] + (bounds[1][i] - bounds[0][i]) / (2.0 * (n_sample_points[i] - 1.0)), (bounds[1][i] - bounds[0][i]) / (n_sample_points[i] - 1.0) )
                        if n_sample_points[i] > 1
                        else [(bounds[1][i] - bounds[0][i]) / 2.0] 
                        for i in range(len(n_sample_points))]
        all_combos = itertools.product(*param_ranges)
        all_combos = list(all_combos)
        all_combos = [list(combo) for combo in all_combos]
        #print 'all_combos = ' + str(all_combos)
        #print 'len(all_combos) = ' + str(len(all_combos))
    elif space_sampling_method == 'random':
        all_combos = np.random.rand(n_sample_points, len(bounds[0]) )
        for i in range(len(bounds[0])):
            bound = [bounds[0][i], bounds[1][i]]
            all_combos[:,i] = all_combos[:,i] * (bound[1] - bound[0]) + bound[0]
        all_combos = all_combos.tolist() 
    else:
        all_combos = np.random.rand(n_sample_points, len(bounds[0]) )
        for i in range(len(bounds[0])):
            bound = [bounds[0][i], bounds[1][i]]
            all_combos[:,i] = all_combos[:,i] * (bound[1] - bound[0]) + bound[0]
        all_combos = all_combos.tolist() 
            
    calculated_least_squares = {}
    least_square = np.inf
    least_square_combo = -1
    current_index = -1
    
    for combo in list(all_combos):
        current_index = current_index + 1
        if current_index % 200 == 0:
            print 'On iteration ' + str(current_index) 
            #print 'Measuring least square value for combo ' + str(combo)
        #print 'combo = ' + str(combo)
        residuals = err_global(list(combo), *args)
        new_least_square = math.sqrt(sum( [ residual ** 2.0 for residual in residuals ] ))
        #print 'new_least_square = ' + str(new_least_square) 
        if new_least_square < least_square:
            least_square = new_least_square
            #print 'Assigning least_square combo to = ' + str(combo) 
            least_square_combo = combo
        calculated_least_squares[tuple(combo)] = new_least_square 
    #print 'param_ranges = ' + str(param_ranges) 
    #print 'least_square value is ' + str(least_square)
    #print 'at parameter set ' + str(least_square_combo) 

    return least_square_combo

def getSimultaneousFitParams(x_data_arrays, y_data_arrays, functions, param_indeces, guess_params = [], y_err_arrays = None, param_bounds = (-np.inf, np.inf), n_steps = 10, space_sampling_method = 'random'  ):
    
    if y_err_arrays is None or y_err_arrays is 'none':
        y_err_arrays = [ np.zeros(len(y_data_array)) + 1.0 for y_data_array in y_data_arrays ]
    #print 'scipy.optimize.least_squares(err_global, guess_parameters, args = (functions, param_indeces, x_data_arrays, y_data_arrays, y_err_arrays)) = '
    #print scipy.optimize.least_squares(err_global, guess_parameters, args = (functions, param_indeces, x_data_arrays, y_data_arrays, y_err_arrays))
    #print "scipy.optimize.least_squares(err_global, guess_parameters, args = (functions, param_indeces, x_data_arrays, y_data_arrays, y_err_arrays))['x']= "
    #print scipy.optimize.least_squares(err_global, guess_parameters, args = (functions, param_indeces, x_data_arrays, y_data_arrays, y_err_arrays))['x']
    #print 'guess_params = ' + str(guess_params) 
    if len(guess_params) == 0:
        #print 'here 1' 
        guess_params = simpleLeastSquareFinder(err_global, param_bounds, n_steps, args = (functions, param_indeces, x_data_arrays, y_data_arrays, y_err_arrays), space_sampling_method = space_sampling_method)
        #print 'Now guess_params = ' + str(guess_params) 
    for i in range(len(guess_params)):
        param = guess_params[i]
        if type(param_bounds[0]) is float or type(param_bounds[0]) is int:
            if param <= param_bounds[0]: guess_params[i] = param_bounds[0]
        else:
            if param <= param_bounds[0][i]: guess_params[i] = param_bounds[0][i]
        if type(param_bounds[1]) is float or type(param_bounds[1]) is int:
            if param >= param_bounds[1]: guess_params[i] = param_bounds[1]
        else:
            if param >= param_bounds[1][i]: guess_params[i] = param_bounds[1][i]
    #print 'guess_params = ' + str(guess_params)
    #print 'param_bounds = ' + str(param_bounds ) 
    simul_fit_results = scipy.optimize.least_squares(err_global, guess_params, bounds = param_bounds, args = (functions, param_indeces, x_data_arrays, y_data_arrays, y_err_arrays))
    p_best = simul_fit_results['x']
    residuals = simul_fit_results['fun']
    jacobian = simul_fit_results['jac']
    #print 'jacobian = '
    #print jacobian
    #print 'residuals = '
    #print residuals
    #print 'np.dot(jacobian.transpose(), jacobian) = supposed Hessian = '
    #print np.dot(jacobian.transpose(), jacobian)
    #print 'np.linalg.inv(np.dot(jacobian.transpose(), jacobian)) = supposed covariance matrix = '
    #print np.linalg.inv(np.dot(jacobian.transpose(), jacobian))

    hessian = np.dot(jacobian.transpose(), jacobian)
    try:
        cov_matrix = np.linalg.inv(np.dot(jacobian.transpose(), jacobian))
    except np.linalg.LinAlgError as err:
        print 'Could not invert hessian.  Got error: '
        print err
    #    print 'jacobian matrix is: '
    #    print jacobian
    #    print 'jacobian.transpose(): '
    #    print jacobian.transpose()
    #    print 'Hessian matrix is: '
    #    print hessian 
        cov_matrix = np.zeros(np.shape(hessian)) 
    return p_best, [cov_matrix[i][i] for i in range(len(p_best))]
    #return simul_fit_results 

    
    
