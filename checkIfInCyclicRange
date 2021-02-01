

#Checks if a parameter is within some range, where the range is cyclic
# (eg, could be an angle with cyclicity 2 pi rad or an RA on the sky
# with 360 deg)
#It is on the user to make sure that all entered values are of the same units.

def getCyclicRanges(orig_range, cyclicity, cyclic_start = 0):
    num_cyclities = [int((elem - cyclic_start) / cyclicity) for elem in orig_range]
    n_ranges = max(num_cyclicities) - min(num_cyclicities)
    cyclic_range_bounds = [(orig_range[0] - cyclic_start) % cyclicity]
    for i in range(n_ranges):
        cyclic_range_bounds = cyclic_range+bounds + [cyclicity, 0.0]

    cyclic_range_bounds = cyclic_range_bounds + [(orig_range[1] - cyclic_start) % cyclicity]

    cyclic_range_bounds = [bound + cyclic_start for bound in cyclic_range_bounds]

    cyclic_ranges = [[bound[i], bound[i+1]] for i in range(len(cyclic_range_bounds) - 1) ]
    return cyclic_ranges

def checkIfInCyclicRange(val, orig_range, cyclicity, cyclic_start = 0):
    
    cyclic_ranges = getCyclicRanges(orig_range, cyclicity, cyclic_start = cyclic_start)
    for range in ranges:
        if val in range:
            return 1
    return 0 

    
