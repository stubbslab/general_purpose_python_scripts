from ApproximateFourierTransform import ApproximateFourierTransform
from computeDistributionOfFourierModesOfFunction import measureLikelihoodOfFunctionGivenData
from readInArchiveSpread import readInArchiveSpread 
import numpy as np
import math
import scipy.optimize as optimize

#You're guesses MUST sit above and below the contour you're hoping to hit.  Otherwise, the algorithm will not converge. 

def measureFourierBounds (xs, ys, y_errs, single_var_function, init_low_guess, init_high_guess, 
                          archive_spread_file = 'standard_rand_fit_results_file_n100.npy', windowing_funct = 'rect', max_n_draw = 100,
                          n_chi_sqr_sigma_tolerance = 3.0, frequencies_to_use = 'all', n_refinements = 3, n_needed_to_reject = 2):
    archive_data = readInArchiveSpread(file_name = archive_spread_file)

    chi_sqr_index = 1
    freq_index = 2
    params_index = 3
    random_draws_index = 4
    

    frequencies = archive_data[freq_index]
    if frequencies_to_use is 'all':
        frequencies_to_use = frequencies
    original_random_draws = archive_data[random_draws_index]

    n_draws_used_in_archive = len(original_random_draws[0]) 
    archive_chi_sqr = archive_data[chi_sqr_index]
    #asymptotic_target_chi_sqr = archive_chi_sqr + math.sqrt(2 * archive_chi_sqr) * target_n_chi_sqr_sigmas_from_archive
    #print 'frequencies = ' + str(frequencies)
    #print 'len(frequencies) = ' + str(len(frequencies))
    print 'archive_chi_sqr = ' + str(archive_chi_sqr)
    print 'archive_chi_sqr = ' + str(archive_chi_sqr)

    chi_sqr_at_0_for_1_trial = 3474.9 #Measured from repeated computation 

    extra_factor_for_guessing_peaks_of_dist = 2.0 
    #min_chi_sqr_separation = lambda n_draws_determining_chi_sqr: ((math.sqrt(2 * len(frequencies)) * target_n_chi_sqr_sigmas_from_archive)
    #                                                               * ( extra_factor_for_guessing_peaks_of_dist / np.sqrt(n_draws_determining_chi_sqr) )
    #                                                             )
    #My formula (without theoretical motivation at present) is chiSqr(n) = chiSqr(infinity) * (1.0 + 1.0 / n * CONSTANT)
    const_coef = (chi_sqr_at_0_for_1_trial - archive_chi_sqr) * n_draws_used_in_archive / (archive_chi_sqr *  n_draws_used_in_archive - chi_sqr_at_0_for_1_trial)
    asymptotic_target_chi_sqr = chi_sqr_at_0_for_1_trial / (1.0 + const_coef)
    print 'const_coef = ' + str(const_coef)
    print 'asymptotic_target_chi_sqr = ' + str(asymptotic_target_chi_sqr) 
    
    targ_chi_sqr_bounds_from_number_of_draws = lambda n_draws_determining_chi_sqr: [ (asymptotic_target_chi_sqr * (1.0 + const_coef / n_draws_determining_chi_sqr) 
                                                                                     - math.sqrt(2 * asymptotic_target_chi_sqr * (1.0 + const_coef / n_draws_determining_chi_sqr) )
                                                                                      * n_chi_sqr_sigma_tolerance),
                                                                                     (asymptotic_target_chi_sqr * (1.0 + const_coef / n_draws_determining_chi_sqr) 
                                                                                     + math.sqrt(2 * asymptotic_target_chi_sqr * (1.0 + const_coef / n_draws_determining_chi_sqr) )
                                                                                      * n_chi_sqr_sigma_tolerance)
                                                                                     ]
    

    true_fourier_trans = ApproximateFourierTransform(xs, ys, y_errs, frequencies = frequencies_to_use, windowing_funct = windowing_funct)
    true_mags = true_fourier_trans.normalized_coef_mags 
    
    widths = [params[1] for params in archive_data[params_index]]

    low_param_bound = init_low_guess
    high_param_bound = init_high_guess

    has_converged = 0
    boundary_found = 0
    n_refinements_done = 0

    while not (has_converged):
        print 'Here 1. '
        if boundary_found:
            n_refinements_done = n_refinements_done + 1
            print 'On refinement ' + str(n_refinements_done)
        #Now regenerated random data until a peak can be determined with sufficient confidence to say which side of our desired contour it is on.
        new_param = (high_param_bound + low_param_bound) / 2.0
        new_funct_ys = single_var_function(xs, new_param)
        peak_determined_to_be_beyond_bounds = 0
        randomly_determined_modes = []
        n_rand_samples = 0
        reached_max_n_draw = 0
        n_above = 0
        n_below = 0
        while not (peak_determined_to_be_beyond_bounds or reached_max_n_draw):
            n_rand_samples = n_rand_samples + 1
            noised_funct_ys = np.random.normal(new_funct_ys, y_errs) 
            randomly_determined_modes = randomly_determined_modes + [ApproximateFourierTransform(xs, noised_funct_ys, y_errs, frequencies = frequencies, windowing_funct = windowing_funct).normalized_coef_mags]
            mean_mags = np.mean(np.array(randomly_determined_modes), axis = 0)
            new_n_sigma_differences = [ ((true_mags[i] - mean_mags[i]) / widths[i]) for i in range(len(true_mags)) ]
            new_chi_sqr = sum([elem ** 2.0 for elem in new_n_sigma_differences])
            
            target_chi_sqr_bounds = targ_chi_sqr_bounds_from_number_of_draws(n_rand_samples)
            print 'For n_rand_samples = ' + str(n_rand_samples) +', target_chi_sqr_bounds = ' + str(target_chi_sqr_bounds) + ' and new_chi_sqr = ' + str(new_chi_sqr) 
            max_tolerable_chi_sqr = target_chi_sqr_bounds[1]
            min_tolerable_chi_sqr = target_chi_sqr_bounds[0]
            
            if new_chi_sqr > max_tolerable_chi_sqr:
                n_above = n_above + 1
                print 'new_chi_sqr larger than max target_chi_sqr boundary'
                if n_above >= n_needed_to_reject:
                    high_param_bound = new_param
                    peak_determined_to_be_beyond_bounds = 1
            elif new_chi_sqr < min_tolerable_chi_sqr:
                n_below = n_below + 1
                print 'new_chi_sqr smaller than min target_chi_sqr boundary'
                if n_below >= n_needed_to_reject:
                    low_param_bound = new_param
                    peak_determined_to_be_beyond_bounds = 1
            elif n_rand_samples > max_n_draw:
                print 'With parameter ' + str(new_param) + ', we have needed to iterate over ' + str(max_n_draw) + ' times to determine if peak is above or below target. '
                reached_max_n_draw = 1
                low_param_bound = new_param #We are interested in the UPPER BOUND of the region.  So this encountered wall defines our new lower value.
                if not(boundary_found): 
                    boundary_found = 1
                    wall_hit_param = new_param
            else:
                print 'Parameter produces chi square in tolerable range, and we can sample new parameters.  Continuing...' 
        if n_refinements_done >= n_refinements:
            has_converged = 1
    print 'The boundary wall was reached with param ' + str(wall_hit_param) 
    print 'With refinments, we determined that the upper boundary lies somewhere between ' + str(low_param_bound) + ' and ' + str(high_param_bound)
    
    return [low_param_bound, high_param_bound]
        
    
