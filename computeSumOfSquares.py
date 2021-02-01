#Computes and returns the 'sum of squares' residual for som data set, given some funciton
import math

def computeSumOfSquares(x_data, y_data, function, y_errs = [] ):
    if len(y_errs) == 0:
        y_errs = [1.0 for x in x_data]
    resids = [function(x_data[i]) - y_data[i] for i in range(len(x_data))]
    squared_normalized_dists = [(resids[i] / y_errs[i]) ** 2.0 for i in range(len(x_data)) ]

    sum_of_squares =  sum(squared_normalized_dists)
    return math.sqrt(sum_of_squares) 
