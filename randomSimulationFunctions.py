#Take in two lists (x-data and y-data for some points) and return
# the same lists, with a random rearrangement of which x-y pairings.

import random
from cantrips import removeListElement

def simulateNormallyDistributedData(x_data, function, y_errs, n_draws = -1, replace = 0):
    
    if n_draws < 0: n_draws = len(x_data) 

    new_x_data = []
    new_y_data = []
    new_y_errs = []

    for i in range(n_draws):
        x_index = random.randint(0, len(x_data) - 1)
        #print 'x_index = ' + str(x_index)
        #print 'y_index = ' + str(y_index) 
        new_x_data = new_x_data + [x_data[x_index]]
        new_y_data = new_y_data + [function(x_data[x_index]) + random.gauss(0.0, y_errs[x_index])]
        new_y_errs = new_y_errs + [y_errs[x_index]]
        
        if not replace:
            x_data = removeListElement(x_data[:], x_index)
            y_errs =removeListElement(y_errs[:], x_index) 

    return new_x_data, new_y_data, new_y_errs

def randomSortNSigma(x_data, y_data, y_errs, exp_values, n_draws = -1, replace = 0):
    if n_draws < 0: n_draws = len(x_data)
    
    n_sigma = [(y_data[i] - exp_values[i]) / y_errs[i] for i in range(len(x_data)) ]

    new_x_data = []
    new_y_data = []
    new_y_errs = []
    for i in range(n_draws):
        x_index = random.randint(0, len(x_data) - 1)
        n_sigma_index = random.randint(0, len(n_sigma) - 1)

        new_x_data = new_x_data + [x_data[x_index]]
        new_n_sigma = n_sigma[n_sigma_index]
        new_y_data = new_y_data + [new_n_sigma * y_errs[x_index] + exp_values[x_index]]
        new_y_errs = new_y_errs + [y_errs[x_index]]

        if not replace:
            x_data = removeListElement(x_data[:], x_index)
            n_sigma =removeListElement(n_sigma[:], n_sigma_index)
            y_errs =removeListElement(y_errs[:], x_index)

    return new_x_data, new_y_data, new_y_errs  


def randomSortData(x_data, y_data, n_draws = -1, y_errs = None, replace = 0):

    if n_draws < 0: n_draws = len(x_data) 

    new_x_data = []
    new_y_data = []
    if not y_errs is None: new_y_errs = []

    for i in range(n_draws):
        #print 'i = ' + str(i) 
        x_index = random.randint(0, len(x_data) - 1)
        y_index = random.randint(0, len(x_data) - 1)
        #print 'x_index = ' + str(x_index)
        #print 'y_index = ' + str(y_index) 
        new_x_data = new_x_data + [x_data[x_index]]
        new_y_data = new_y_data + [y_data[y_index]]
        if not y_errs is None: new_y_errs = new_y_errs + [y_errs[y_index]]
        
        if not replace:
            x_data = removeListElement(x_data[:], x_index)
            y_data =removeListElement(y_data[:], y_index)
            if not y_errs is None: y_errs =removeListElement(y_errs[:], y_index) 

    if not y_errs is None: return new_x_data, new_y_data, new_y_errs
    else: return new_x_data, new_y_data 
