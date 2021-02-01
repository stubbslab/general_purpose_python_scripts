import numpy as np
import scipy.interpolate as interpolate

#Note: this function must take in a list, no
def computePositionsRecursivelyFromGroups(initial_positions, min_sep):
    #We have two return conditions: either we have only 1 or two points remaining
    initial_positions = np.array(initial_positions).tolist() 
    n_positions = len(initial_positions)
    #print 'initial_positions = ' + str(initial_positions) 
    if n_positions == 1:
        return initial_positions
    elif n_positions == 2:
        if initial_positions[1] - initial_positions[0] > min_sep:
            return initial_positions
        else:
            mean = (initial_positions[0] + initial_positions[1]) / 2.0
            return [mean - min_sep / 2.0, mean + min_sep / 2.0]
    #Otherwise, we have not yet reached the end of our recursive series
    else:
        left_positions = computePositionsRecursivelyFromGroups(initial_positions[0:n_positions / 2], min_sep = min_sep)
        right_positions = computePositionsRecursivelyFromGroups(initial_positions[n_positions / 2:], min_sep = min_sep)
        if right_positions[0] - left_positions[-1] > min_sep:
            return left_positions + right_positions
        else:
            center = (left_positions[-1] + right_positions[0]) / 2.0
            left_positions[-1] = center - min_sep
            right_positions[0] = center + min_sep
            reversed_left_positions = left_positions[:]
            reversed_left_positions.reverse()
            for pos_index in range(1, min(len(left_positions), len(right_positions))):
                if left_positions[-1-pos_index] > left_positions[-1-pos_index+1] - min_sep:
                    left_positions[-1-pos_index] = left_positions[-1-pos_index+1] - min_sep
                if right_positions[pos_index] < right_positions[pos_index] + min_sep:
                    right_positions[pos_index] = right_positions[pos_index] + min_sep
            #deal with possibility that odd-length arrays could have there last element not checked
            if len(right_positions) > 1:
                if right_positions[-1] < right_positions[-2] + min_sep: right_positions[-1] = right_positions[-2] + min_sep
            if len(left_positions) > 1:
                if left_positions[0] > left_positions[1]-min_sep: left_positions[0] = left_positions[1] - min_sep
            return left_positions + right_positions 
                

def computePositionsRecursivelyFromFractions(start, positionSpecifications):
    #if in put is just a single number, than we have reached the end of the list.
    # thus, the value is just the actual end point
    if len(positionSpecifications) is 1:
        return positionSpecifications
    else:
        new_fraction = positionSpecifications[0]
        calculated_positions = (computePositionsRecursivelyFromFractions(start, positionSpecifications[1:]))
        return [start + new_fraction * (calculated_positions[0] - start)] + calculated_positions

def getKnotXPositions(knot_args, spline_technique, min_sep):
    n_knots = len(knot_args) / 2
    knot_xs = []
    knot_ys = knot_args[n_knots:]
    
    if spline_technique in ['fixed_separations', 'fixed_seps','fixedseparations', 'fixedseps']:
        knot_start_x = knot_args[0]
        knot_seps = knot_args[1:n_knots]
        knot_xs = [knot_start_x]
        current_knot_value = knot_start_x
        for sep in knot_seps:
            current_knot_value = current_knot_value + sep
            knot_xs = knot_xs + [current_knot_value]
    elif spline_technique in ['recursive_positions', 'rec_positions', 'recursive_pos', 'rec_pos']:
        knot_start_x = knot_args[0]
        knot_fracs = knot_args[1:n_knots]
        knot_xs = [knot_start_x] + computePositionsRecursivelyFromFractions(knot_start_x, knot_fracs)
    elif spline_technique in ['recursive_groups', 'rec_groups', 'rec_group','recursivegroups', 'recgroups','recursion_groups']:
        knot_initial_xs = knot_args[0:n_knots] 
        knot_xs = computePositionsRecursivelyFromGroups(knot_initial_xs, min_sep)
    else:
        #If we simply allow the values to vary freely, we need sort them, as the spline has trouble when the xs are out of order
        knot_xs = knot_args[0:n_knots]
        knot_xs, knot_ys = (list(knot_vals) for knot_vals in zip(*sorted(zip(knot_xs, knot_ys))))

    return [knot_xs, knot_ys]


def calculateSplineWithFixedNumberOfKnots(xs, shift, min_x, max_x, min_sep, spline_technique, *knot_args):
    xs = np.array(xs) 
    knot_args = list(knot_args)
    
    spline_technique = spline_technique.lower()
    knot_xs = []
    knot_xs, knot_ys = getKnotXPositions(knot_args, spline_technique, min_sep)
    
    x_range = max_x - min_x
    end_anchor_scales = [0.01, 0.1]
    #I tried to keep it under control at the ends of the spline, but I don't think that really helps 
    #my_spline = interpolate.InterpolatedUnivariateSpline([min_x - x_range * (end_anchor_scales[1]),
    #                                                      min_x - x_range * (end_anchor_scales[0])]
    #                                                      + knot_xs
    #                                                      + [max_x + x_range * (end_anchor_scales[0]),
    #                                                         max_x + x_range * (end_anchor_scales[1])],
    #                                                      [0.0,0.0] + knot_ys + [0.0,0.0])

    
    my_spline = interpolate.InterpolatedUnivariateSpline(knot_xs, knot_ys)
    
    ys = my_spline(xs) + shift 
    
    return ys 

def getIterableSplineFunction(min_x, max_x, min_sep, spline_technique):
    return lambda xs, shift, *knot_args: calculateSplineWithFixedNumberOfKnots(xs, shift, min_x, max_x, min_sep, spline_technique, *knot_args)


