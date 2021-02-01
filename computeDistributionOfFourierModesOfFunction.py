import numpy as np
import math
from DirectoryArchive import DirectoryArchive 
from ApproximateFourierTransform import ApproximateFourierTransform
from StandardPeriodigram import StandardPeriodigram 
from cantrips import getCPBFunctFromArray
import matplotlib.pyplot as plt
from showParallel1DTempPlots import showParallel1DTempPlots
from nChooser import nChooser
import scipy.integrate as integrate
import scipy.special as special
import scipy.optimize as optimize
import scipy.interpolate as interpolate
import time 

def calcSingleModes(xs, ys, y_errs, frequencies, windowing_funct, apply_normalization, compute_convolution = 0):
    fourier_funct = StandardPeriodigram(xs, ys, y_errs,
                                        frequencies = frequencies, windowing_funct = windowing_funct,
                                        apply_normalization = apply_normalization, compute_convolution = compute_convolution)
    sin_coefs = fourier_funct.sin_coefs
    cos_coefs = fourier_funct.cos_coefs
    total_mags = fourier_funct.normalized_coef_mags
    frequencies = fourier_funct.frequencies

    return sin_coefs, cos_coefs, total_mags 

def computeDistributionOfFourierModesOfFunction(xs, ys, y_errs, function, calc_true_modes = 1, n_terms = 100, frequencies = None, return_full = 0, windowing_funct = 'rect', pre_calculated_function = 0, sort_by_frequency = 0, apply_normalization = 0, resample_from_data = 0, calc_single_set_function_modes = 0, compute_convolution = 0 ):

    funct_sin_coefs = []
    funct_cos_coefs = []
    funct_total_coef_mags = []
        
    if pre_calculated_function:
        true_function_vals = function
    else: 
        true_function_vals = [function(x) for x in xs]
    
    for i in range(n_terms):
        start = time.time()
        if resample_from_data:
            randomized_vals_for_resampling = np.random.normal(ys, y_errs)
        else: 
            randomized_vals_for_resampling = np.random.normal(true_function_vals, y_errs)
        resampled_wMean = sum([randomized_vals_for_resampling[j] * (1.0 / y_errs[j]) ** 2.0 for j in range(len(ys)) ]) / sum((1.0 / y_err)**2.0 for y_err in y_errs)
        #print 'resampled wMean = ' + str(resampled_wMean )
        randomized_vals_for_resampling = [val for val in randomized_vals_for_resampling - resampled_wMean ]
        
        #windowed_function_vals = [ hamming_vals[i] * calculated_function_vals[i] for i in range(len(hamming_vals)) ]
        fourier_funct = StandardPeriodigram(xs, randomized_vals_for_resampling, y_errs, frequencies = frequencies, apply_normalization = apply_normalization, windowing_funct = windowing_funct, compute_convolution = compute_convolution )
        funct_sin_coefs = funct_sin_coefs + [fourier_funct.sin_coefs]
        funct_cos_coefs = funct_cos_coefs + [fourier_funct.cos_coefs]
        funct_total_coef_mags = funct_total_coef_mags + [fourier_funct.normalized_coef_mags]
        end = time.time()
        if i % 400 == 399: print 'Random draw ' + str(i) + ' took: ' + str(end - start) + 's.' 

    #windowed_ys = [ hamming_vals[i] * ys[i] for i in range(len(hamming_vals)) ]
    if sort_by_frequency:
        sin_coefs_set_by_frequency = [[funct_sin_coef_set[i] for funct_sin_coef_set in funct_sin_coefs] for i in range(len(frequencies)) ]
        cos_coefs_set_by_frequency = [[funct_cos_coef_set[i] for funct_cos_coef_set in funct_cos_coefs] for i in range(len(frequencies)) ]
        total_coef_mags_set_by_frequency = [[funct_total_coef_set[i] for funct_total_coef_set in funct_total_coef_mags] for i in range(len(frequencies)) ]
        return_set = [frequencies, sin_coefs_set_by_frequency, cos_coefs_set_by_frequency, total_coef_mags_set_by_frequency]
    else:
        return_set = [frequencies, funct_sin_coefs, funct_cos_coefs, funct_total_coef_mags]
    
    if calc_true_modes:
        true_sin_coefs, true_cos_coefs, true_total_coef_mags = calcSingleModes(xs, ys, y_errs, frequencies, windowing_funct, apply_normalization, compute_convolution = compute_convolution)
        return_set = return_set + [true_sin_coefs, true_cos_coefs, true_total_coef_mags]

    #for i in range(len(frequencies)):
    #    frequency = frequencies[i]
    #    funct_sin_coefs_at_f = [single_fit[i] for single_fit in funct_sin_coefs]
    #    funct_cos_coefs_at_f = [single_fit[i] for single_fit in funct_cos_coefs]
    #    true_sin_coef_at_f = true_sin_coefs[i]
    #    true_cos_coef_at_f = true_cos_coefs[i]
    #    frac_sin_coefs_above_true = frac_sin_coefs_above_true + [float(len([elem for elem in funct_sin_coefs_at_f if elem > true_sin_coef_at_f])) / float(len(funct_sin_coefs_at_f))]
    #    frac_cos_coefs_above_true = frac_cos_coefs_above_true + [float(len([elem for elem in funct_cos_coefs_at_f if elem > true_cos_coef_at_f])) / float(len(funct_cos_coefs_at_f))]
    if calc_single_set_function_modes:
        single_randomized_funct_vals = np.random.normal(true_function_vals, y_errs)
        single_resampled_wMean = sum([randomized_vals_for_resampling[j] * (1.0 / y_errs[j]) ** 2.0 for j in range(len(ys)) ]) / sum((1.0 / y_err)**2.0 for y_err in y_errs)
        single_randomized_funct_vals = [val for val in single_randomized_funct_vals - single_resampled_wMean ]
        funct_single_sin_coefs, funct_single_cos_coefs, funct_single_total_mags = calcSingleModes(xs, single_randomized_funct_vals, y_errs, frequencies, windowing_funct, apply_normalization)
        return_set = return_set + [funct_single_sin_coefs, funct_single_cos_coefs, funct_single_total_mags]
            
        
    #if return_full
    print 'len(return_set) = ' + str(len(return_set)) 
    return return_set
    #else:
    #    return [frequencies, frac_sin_coefs_above_true, frac_cos_coefs_above_true]


def determineBestFitCBPFunctionsByFrequency(measured_mags_by_frequency, model_CPB_funct, extra_funct, best_guess_minimization_params = 'normal', max_iterations = 100000):
    n_redraw = len(measured_mags_by_frequency[0])
    #exp_frac_to_have_n_below_funct = lambda n_below, n_total: integrate.quad(lambda x: (n_total) * nChooser(n_total - 1, n_below) * x ** (n_below+1) * (1.0 - x) ** (n_total - 1 - n_below), 0, 1 )[0]
    #as it turns out, this can be given as just the actual fraction below.  Who would have thought
    exp_frac_to_have_n_below_funct = lambda n_below, n_total: float(n_below+1) / float(n_total+1)
    #This quantity we just approximate by interpolation from some number, as there is a limit where the computer cannot handle the numbers
    exp_square_frac_to_have_n_below_funct = lambda n_below, n_total: integrate.quad(lambda x: (n_total) * nChooser(n_total - 1, n_below) * x ** (n_below+2) * (1.0 - x) ** (n_total - 1 - n_below), 0, 1 )[0]
    max_number_of_points_to_compute_uncertainty = 100
    if n_redraw < max_number_of_points_to_compute_uncertainty:
        number_of_points_to_compute_uncertainty = n_redraw
    else: 
        number_of_points_to_compute_uncertainty = max_number_of_points_to_compute_uncertainty
    uncertainty_of_exp_frac_interp = interpolate.interp1d(np.linspace(0, n_redraw - 0.5, number_of_points_to_compute_uncertainty) / n_redraw,
                                                                      [ exp_square_frac_to_have_n_below_funct(n_below, number_of_points_to_compute_uncertainty)
                                                                      - exp_frac_to_have_n_below_funct(n_below, number_of_points_to_compute_uncertainty) ** 2.0
                                                                      for n_below in range(number_of_points_to_compute_uncertainty) ],
                                                                      kind = 'linear', )
    exp_frac_to_have_n_below_array = [exp_frac_to_have_n_below_funct(i, n_redraw) for i in range(n_redraw)]
    
    uncertainty_of_exp_frac_array = [uncertainty_of_exp_frac_interp(i / float(n_redraw)) for i in range(n_redraw)]
    mag_vs_exp_frac_by_frequency = [ (sorted(mags), exp_frac_to_have_n_below_array, uncertainty_of_exp_frac_array) for mags in measured_mags_by_frequency] 
    
    funct_to_minimize = lambda  xs, ys, y_errs, funct, params: sum( ((funct(xs, params) - ys) / y_errs ) ** 2.0)
    print 'Determining best fit CPB parameters...'
    minimization_results = []
    best_fit_param_vals_by_frequency = []
    for mag_vs_exp_frac in mag_vs_exp_frac_by_frequency:
        start = time.time()
        bounds = ((min(mag_vs_exp_frac[0]), None), (0.0001, None))
        
        if best_guess_minimization_params in ['normal','gaussian','norm','gauss']:
            best_guess_minimization_params_instance = [np.mean(mag_vs_exp_frac[0]), np.std(mag_vs_exp_frac[0])]
        elif best_guess_minimization_params in ['alt_normal', 'alt_norm', 'alt_guass','alt_gaussian','alternative']:
            best_guess_minimization_params_instance = [np.std(mag_vs_exp_frac[0]), 2.0 ]
        else:
            best_guess_minimization_params_instance = best_guess_minimization_params
        #print 'best_guess_minimization_params_instance = ' + str(best_guess_minimization_params_instance)
        new_min_results = optimize.minimize( lambda params: funct_to_minimize(mag_vs_exp_frac[0], mag_vs_exp_frac[1], mag_vs_exp_frac[2], model_CPB_funct, params),best_guess_minimization_params_instance, options=dict({'maxiter':max_iterations}), bounds = bounds )
        #print 'best_fit_results = ' + str(new_min_results.x) 
        minimization_results = minimization_results + [new_min_results]
        
        end = time.time()
        #print 'Took ' + str(end- start) + ' s to do one minization. '
    best_fit_params_by_frequency = [res.x for res in minimization_results]
    #print best_fit_params_by_frequency 
    #best_fit_param_vals_by_frequency = [optimize.minimize( lambda params: funct_to_minimize(mag_vs_exp_frac[0], mag_vs_exp_frac[1], mag_vs_exp_frac[2], model_CPB_funct, params),
    #                                                        best_guess_minimization_params ).x
    #                                     for mag_vs_exp_frac in mag_vs_exp_frac_by_frequency ]
    print 'Done determining best CPB parameters. '

    best_fit_cpb_functs = [(lambda xs, best_fit_cpb_vals = best_fit_cpb_vals: model_CPB_funct (xs, best_fit_cpb_vals) )
                               for best_fit_cpb_vals in best_fit_params_by_frequency]
    best_fit_extra_functs = [(lambda xs, best_fit_cpb_vals = best_fit_cpb_vals: extra_funct (xs, best_fit_cpb_vals) )
                                for best_fit_cpb_vals in best_fit_params_by_frequency]
    return [mag_vs_exp_frac_by_frequency, best_fit_params_by_frequency, best_fit_extra_functs, best_fit_cpb_functs] 
    
    

def measureLikelihoodOfFunctionGivenData(xs, ys, y_errs, function_to_decompose,
                                         best_guess_minimization_params = None, model_CPB_funct = 'err_funct', model_prob_density_funct = 'normal', resample_from_data = 1,
                                         n_redraw = 100, frequencies = None, windowing_funct = 'rect', pre_calculated_function = 0, measure_prob_from_total_mag = 1, apply_normalization = 1 ):
    frequencies, sin_coefs_set_by_frequency, cos_coefs_set_by_frequency, total_coef_mags_set_by_frequency, true_sin_coefs, true_cos_coefs, true_total_coef_mags, funct_single_draw_sin_coefs, funct_single_draw_cos_coefs, funct_single_draw_total_coef_mags  = computeDistributionOfFourierModesOfFunction(xs, ys, y_errs, function_to_decompose, n_terms = n_redraw, frequencies = frequencies, return_full = 1, windowing_funct = windowing_funct, pre_calculated_function = pre_calculated_function, sort_by_frequency = 1, apply_normalization = apply_normalization, resample_from_data = resample_from_data, calc_single_set_function_modes = 1)

    if model_CPB_funct is 'err_funct':     
        err_funct = lambda x, params: 0.5 * (1.0 + special.erf ((x - params[0])/ (params[1] * math.sqrt(2.0))))
        truncated_err_funct = lambda x, params: (err_funct(x, params) - err_funct(0.0, params)) / (1.0 - err_funct(0.0, params))
        model_CPB_funct = truncated_err_funct
    #determine a bunch of normal fits to the data and then calculate the probability 
    if model_prob_density_funct is 'normal':
        err_funct = lambda x, params: 0.5 * (1.0 + special.erf ((x - params[0])/ (params[1] * math.sqrt(2.0))))
        truncated_gaussian = lambda x, params: 1.0 / (math.sqrt(2.0 * math.pi) * params[1]) * np.exp(-(x-params[0]) ** 2.0 / (2.0 * params[1] ** 2.0)) / (1.0 - err_funct(0.0, params)) 
        model_prob_density_funct = truncated_gaussian
    elif model_prob_density_funct is 'unity_normal': 
        err_funct = lambda x, params: 0.5 * (1.0 + special.erf ((x - params[0])/ (params[1] * math.sqrt(2.0))))
        unity_gaussian = lambda x, params: np.exp(-(x-params[0]) ** 2.0 / (2.0 * params[1] ** 2.0)) / (1.0 - err_funct(0.0, params)) 
        model_prob_density_funct = unity_gaussian
        

    funct_to_minimize = lambda  xs, ys, y_errs, funct, params: sum( ((funct(xs, params) - ys) / y_errs ) ** 2.0)
    
    if measure_prob_from_total_mag:
        best_fit_cpb_results = determineBestFitCBPFunctionsByFrequency(total_coef_mags_set_by_frequency, model_CPB_funct, model_prob_density_funct, best_guess_minimization_params)
        if resample_from_data: 
            total_log_prob = sum([np.log(best_fit_cpb_results[2][i](funct_single_draw_total_coef_mags[i])) for i in range(len(funct_single_draw_total_coef_mags))] )
            raw_chi_sqr = sum( (funct_single_draw_total_coef_mags[i] - best_fit_cpb_results[1][i][0]) ** 2.0 /(best_fit_cpb_results[1][i][1]) ** 2.0 for i in range(len(funct_single_draw_total_coef_mags)) )
        else: 
            total_log_prob = sum([np.log(best_fit_cpb_results[2][i](true_total_coef_mags[i])) for i in range(len(true_total_coef_mags))] )
            raw_chi_sqr = sum( (true_total_coef_mags[i] - best_fit_cpb_results[1][i][0]) ** 2.0 /(best_fit_cpb_results[1][i][1]) ** 2.0 for i in range(len(true_total_coef_mags)) )
        return [total_log_prob, raw_chi_sqr] + [frequencies, true_total_coef_mags] + best_fit_cpb_results 
    else:
        best_fit_cpb_results_sin = determineBestFitCBPFunctionsByFrequency(sin_coefs_set_by_frequency, model_CPB_funct, model_prob_density_funct, best_guess_minimization_params)
        best_fit_cpb_results_cos = determineBestFitCBPFunctionsByFrequency(cos_coefs_set_by_frequency, model_CPB_funct, model_prob_density_funct, best_guess_minimization_params) 
        total_log_prob = ( sum([np.log(best_fit_cpb_results_sin[2][i](true_total_coef_mags[i])) for i in range(len(true_sin_coefs))] )
                           + sum([np.log(best_fit_cpb_results_cos[2][i](true_total_coef_mags[i])) for i in range(len(true_cos_coefs))] )
                         )
        raw_chi_sqr = ( sum( (true_sin_coefs[i] - best_fit_cpb_results_sin[1][i][0]) ** 2.0 /(best_fit_cpb_results_sin[1][i][1]) ** 2.0 for i in range(len(true_sin_coefs)) )
                        + sum( (true_cos_coefs[i] - best_fit_cpb_results_cos[1][i][0]) ** 2.0 /(best_fit_cpb_results_cos[1][i][1]) ** 2.0 for i in range(len(true_cos_coefs)) )
                      )
        return [total_log_prob, raw_chi_sqr] + best_fit_cpb_results_sin + best_fit_cpb_results_cos

def batchComputeDistributionOfFourierModesGivenFunction(xs, ys, y_errs, function_to_decompose, frequencies,
                                                        best_guess_minimization_params = 'normal', model_CPB_funct = 'err_funct', model_prob_density_funct = 'normal',
                                                        freq_batch_size = 100, n_redraw = 100, windowing_funct = 'rect', compute_convolution = 0, 
                                                        pre_calculated_function = 0, measure_prob_from_total_mag = 1, apply_normalization = 1,  show_single_batch_fit = 0):

    #Note that we cannot compute the normalization, as it is likely dependent on the overall frequency region, of which we are computing only a small section 

    if model_CPB_funct is 'err_funct':     
        err_funct = lambda x, params: 0.5 * (1.0 + special.erf ((x - params[0])/ (params[1] * math.sqrt(2.0))))
        truncated_err_funct = lambda x, params: (err_funct(x, params) - err_funct(0.0, params)) / (1.0 - err_funct(0.0, params))
        model_CPB_funct = truncated_err_funct
    elif model_CPB_funct is 'zero_alt_err_funct':
        alt_err_funct = lambda x, params: special.gammainc(1.0 / params[1], (x/params[0]) ** params[1]) / (special.gamma(1.0 / params[1]))
        model_CPB_funct = alt_err_funct 
                                                                                                                 
    #determine a bunch of normal fits to the data and then calculate the probability 
    if model_prob_density_funct is 'normal':
        err_funct = lambda x, params: 0.5 * (1.0 + special.erf ((x - params[0])/ (params[1] * math.sqrt(2.0))))
        truncated_gaussian = lambda x, params: 1.0 / (math.sqrt(2.0 * math.pi) * params[1]) * np.exp(-(x-params[0]) ** 2.0 / (2.0 * params[1] ** 2.0)) / (1.0 - err_funct(0.0, params)) 
        model_prob_density_funct = truncated_gaussian
    elif model_prob_density_funct is 'unity_normal': 
        err_funct = lambda x, params: 0.5 * (1.0 + special.erf ((x - params[0])/ (params[1] * math.sqrt(2.0))))
        unity_gaussian = lambda x, params: np.exp(-(x-params[0]) ** 2.0 / (2.0 * params[1] ** 2.0)) / (1.0 - err_funct(0.0, params)) 
        model_prob_density_funct = unity_gaussian

    
    true_total_coef_mags = []
    best_fit_cpb_results_params = []
    true_sin_coefs = []
    true_cos_coefs = []
    best_fit_cpb_results_params_sin = []
    best_fit_cpb_results_params_cos = []

    still_batching = 1
    batch_start_index = 0
    
    while still_batching:
        print 'Working on batch with start frequency index ' + str(batch_start_index) + ' of ' + str(len(frequencies)) 
        if batch_start_index + freq_batch_size < len(frequencies):
            frequency_batch = frequencies[batch_start_index:batch_start_index + freq_batch_size]
            batch_start_index = batch_start_index + freq_batch_size
        else:
            frequency_batch = frequencies[batch_start_index:]
            still_batching = 0
        frequency_batch, sin_coefs_set_by_frequency_batch, cos_coefs_set_by_frequency_batch, total_coef_mags_set_by_frequency_batch, true_sin_coefs_batch, true_cos_coefs_batch, true_total_coef_mags_batch, funct_single_draw_sin_coefs_batch, funct_single_draw_cos_coefs_batch, funct_single_draw_total_coef_mags_batch  = computeDistributionOfFourierModesOfFunction(xs, ys, y_errs, function_to_decompose, n_terms = n_redraw, frequencies = frequency_batch, return_full = 1, windowing_funct = windowing_funct, pre_calculated_function = pre_calculated_function, sort_by_frequency = 1, apply_normalization = apply_normalization, resample_from_data = 0, calc_single_set_function_modes = 1, compute_convolution = compute_convolution)

        if measure_prob_from_total_mag:
            best_fit_cpb_results_batch = determineBestFitCBPFunctionsByFrequency(total_coef_mags_set_by_frequency_batch, model_CPB_funct, model_prob_density_funct, best_guess_minimization_params)
            best_fit_cpb_results_params_batch = best_fit_cpb_results_batch[1] 
            true_total_coef_mags = true_total_coef_mags + true_total_coef_mags_batch
            best_fit_cpb_results_params = best_fit_cpb_results_params + best_fit_cpb_results_params_batch
            if show_single_batch_fit:
                displayed_frequency = frequency_batch[len(total_coef_mags_set_by_frequency_batch) / 2]
                fitted_points_to_show = best_fit_cpb_results_batch[0][len(total_coef_mags_set_by_frequency_batch) / 2]
                fitted_params_to_show = best_fit_cpb_results_params_batch[len(total_coef_mags_set_by_frequency_batch) / 2]
        else:
            best_fit_cpb_results_sin_batch = determineBestFitCBPFunctionsByFrequency(sin_coefs_set_by_frequency_batch, model_CPB_funct, model_prob_density_funct, best_guess_minimization_params)
            best_fit_cpb_results_cos_batch = determineBestFitCBPFunctionsByFrequency(cos_coefs_set_by_frequency_batch, model_CPB_funct, model_prob_density_funct, best_guess_minimization_params)
            best_fit_cpb_results_sin_params_batch = best_fit_cpb_results_sin_batch[1]
            best_fit_cpb_results_cos_params_batch = best_fit_cpb_results_cos_batch[1] 
            true_sin_coefs = true_sin_coefs + true_sin_coefs_batch
            true_cos_coefs = true_cos_coefs + true_cos_coefs_batch
            best_fit_cpb_results_params_sin = best_fit_cpb_results_params_sin + best_fit_cpb_results_sin_params_batch
            best_fit_cpb_results_params_cos = best_fit_cpb_results_params_cos + best_fit_cpb_results_cos_params_batch

    if show_single_batch_fit:
        plt.scatter(fitted_points_to_show[0], fitted_points_to_show[1])
        plt.errorbar(fitted_points_to_show[0], fitted_points_to_show[1], yerr= fitted_points_to_show[2], fmt = None)
        plt.plot(fitted_points_to_show[0], model_CPB_funct(fitted_points_to_show[0], fitted_params_to_show ))
        plt.xlabel('Measured Fourier Power')
        plt.ylabel('Expected Portion of probability Distribution Below Point')
        plt.title('Sample determination of Distribution of Zero for Frequency ' + str(displayed_frequency))
        plt.xlim([0.0, max (fitted_points_to_show[0]) * 1.01])
        plt.ylim([0.0, max (fitted_points_to_show[1]) * 1.01])
        plt.text(max (fitted_points_to_show[0]) * 1.01 * 0.5, max (fitted_points_to_show[1]) * 1.01 * 0.75, str(fitted_params_to_show))
        plt.show() 
    if measure_prob_from_total_mag:
        return_set = [frequencies, true_total_coef_mags, best_fit_cpb_results_params]
    else:
        return_set = [frequencies, true_sin_coefs, true_cos_coefs, best_fit_cpb_results_sin_params, best_fit_cpb_results_cos_params]
    return return_set 

    
def showDistributionOfFourierModesOfFunction(xs, ys, y_errs, function_to_decompose, disp_y_range = None, n_redraw = 100, frequencies = None, pre_calculated_function = 0, 
                                             white_space = 0.0, n_ys_in_range = 100, windowing_funct = 'rect', display_direct_cpbs = 1, show_n_sigma = -1.0, resample_from_data = 0, 
                                             model_CPB_funct = 'err_funct', funct_to_display = 'unity_normal', best_guess_minimization_params = None, compute_convolution = 0,
                                             disp_y_limits = [0.0, np.inf], levels = None, colormap = None, n_xticks = None, compute_total_coef_at_freq = 1, figsize = (10,10),
                                             apply_normalization = 0,  xlabels = ['frequency', 'frequency'], ylabels = ['magnitude of coefficient', 'cos coefficient'], titles = ['', ''],
                                             suptitle = '',show = 1, save = 0, return_mags = 0, funct_name = 'UnNamed', ind_var_name = 'UnNamed', dep_var_name = 'UnNamed', other_plot_name = ''):
    dir_archive = DirectoryArchive()

    frequencies, sin_coefs_set_by_frequency, cos_coefs_set_by_frequency, total_coef_mags_set_by_frequency, true_sin_coefs, true_cos_coefs, true_total_coef_mags, funct_single_draw_sin_coefs, funct_single_draw_cos_coefs, funct_single_draw_total_coef_mags  = computeDistributionOfFourierModesOfFunction(xs, ys, y_errs, function_to_decompose, n_terms = n_redraw, frequencies = frequencies, return_full = 1, windowing_funct = windowing_funct, pre_calculated_function = pre_calculated_function, sort_by_frequency = 1, apply_normalization = apply_normalization, resample_from_data = resample_from_data,calc_single_set_function_modes = 1, compute_convolution = compute_convolution)

    if display_direct_cpbs:
        #Display the directly measured CPBs
        if compute_total_coef_at_freq:
            total_mag_funct_by_frequency = [ lambda sampling_total_coef_mag, i=i: (getCPBFunctFromArray(sorted(total_coef_mags_set_by_frequency[i]))(sampling_total_coef_mag))
                                             for i in range(len(frequencies)) ]
        else: 
            sin_mag_funct_by_frequency = [ lambda sampling_sin_coefs, i=i: (getCPBFunctFromArray(sorted(sin_coefs_set_by_frequency[i]))(sampling_sin_coefs))
                                           for i in range(len(frequencies)) ]
            cos_mag_funct_by_frequency = [ lambda sampling_cos_coefs, i=i: (getCPBFunctFromArray(sorted(cos_coefs_set_by_frequency[i]))(sampling_cos_coefs))
                                           for i in range(len(frequencies)) ]
    else:
        if model_CPB_funct is 'err_funct':     
            err_funct = lambda x, params: 0.5 * (1.0 + special.erf ((x - params[0])/ (params[1] * math.sqrt(2.0))))
            truncated_err_funct = lambda x, params: (err_funct(x, params) - err_funct(0.0, params)) / (1.0 - err_funct(0.0, params))
            model_CPB_funct = truncated_err_funct
        #determine a bunch of normal fits to the data and then calculate the probability 
        if funct_to_display is 'normal':
            err_funct = lambda x, params: 0.5 * (1.0 + special.erf ((x - params[0])/ (params[1] * math.sqrt(2.0))))
            truncated_gaussian = lambda x, params: 1.0 / (math.sqrt(2.0 * math.pi) * params[1]) * np.exp(-(x-params[0]) ** 2.0 / (2.0 * params[1] ** 2.0)) / (1.0 - err_funct(0.0, params)) 
            funct_to_display = truncated_gaussian
        elif funct_to_display is 'unity_normal': 
            err_funct = lambda x, params: 0.5 * (1.0 + special.erf ((x - params[0])/ (params[1] * math.sqrt(2.0))))
            unity_gaussian = lambda x, params: np.exp(-(x-params[0]) ** 2.0 / (2.0 * params[1] ** 2.0)) / (1.0 - err_funct(0.0, params)) 
            funct_to_display = unity_gaussian
        if compute_total_coef_at_freq:
            bestFitCPBFunctionResults = determineBestFitCBPFunctionsByFrequency(total_coef_mags_set_by_frequency, model_CPB_funct, funct_to_display, best_guess_minimization_params)
            total_mag_best_fit_params_by_frequency = bestFitCPBFunctionResults[1]
            total_mag_funct_by_frequency = bestFitCPBFunctionResults[2]
            #total_mag_funct_by_frequency = determineBestFitCBPFunctionsByFrequency(total_coef_mags_set_by_frequency, model_CPB_funct, funct_to_display, best_guess_minimization_params)[2]
        else: 
            sin_mag_funct_by_frequency = determineBestFitCBPFunctionsByFrequency(sin_coefs_set_by_frequency, model_CPB_funct, funct_to_display, best_guess_minimization_params)[2]
            cos_mag_funct_by_frequency = determineBestFitCBPFunctionsByFrequency(cos_coefs_set_by_frequency, model_CPB_funct, funct_to_display, best_guess_minimization_params)[2]
        

    if compute_total_coef_at_freq:
        f, axarr = plt.subplots(1,1, squeeze = False, figsize = figsize)
    else:
        f, axarr = plt.subplots(2,1, squeeze = False, figsize = figsize)

    if disp_y_range is None:
        min_sin_coef = np.inf
        max_sin_coef = -np.inf
        min_cos_coef = np.inf
        max_cos_coef = -np.inf
        min_total_mag_coef = np.inf
        max_total_mag_coef = -np.inf
        fourier_freqs = []
        fourier_mags = []
        for i in range(len(frequencies)):
            sin_coefs_at_f = sin_coefs_set_by_frequency[i]
            cos_coefs_at_f = cos_coefs_set_by_frequency[i]
            total_mag_coef_at_f = total_coef_mags_set_by_frequency[i]
            for sin_coef in sin_coefs_at_f:
                if sin_coef > max_sin_coef: max_sin_coef = sin_coef
                if sin_coef < min_sin_coef: min_sin_coef = sin_coef
            for cos_coef in cos_coefs_at_f:
                if cos_coef > max_cos_coef: max_cos_coef = cos_coef
                if cos_coef < min_cos_coef: min_cos_coef = cos_coef
            for total_coef_mag in total_mag_coef_at_f:
                if total_coef_mag > max_total_mag_coef: max_total_mag_coef = total_coef_mag
                if total_coef_mag < min_total_mag_coef: min_total_mag_coef = total_coef_mag

            if return_mags:
                mean_mag_at_f = np.mean(total_mag_coef_at_f)
                fourier_freqs = fourier_freqs + [frequencies[i]]
                fourier_mags = fourier_mags + [mean_mag_at_f]
                
        sin_coef_range = [min_sin_coef, max_sin_coef]
        cos_coef_range = [min_cos_coef, max_cos_coef]
        total_mag_coef_range = [min_total_mag_coef, max_total_mag_coef]
        if  display_direct_cpbs:
            total_mag_coef_range[0] = min([total_mag_coef_range[0]] + true_total_coef_mags)
            total_mag_coef_range[1] = max([total_mag_coef_range[1]] + true_total_coef_mags) 
        else:
            total_mag_coef_range[0] = min([total_mag_coef_range[0]]
                                      + true_total_coef_mags
                                      + [params_at_freq[1] * show_n_sigma + params_at_freq[0] for params_at_freq in total_mag_best_fit_params_by_frequency])
            total_mag_coef_range[1] = max([total_mag_coef_range[1]]
                                      + true_total_coef_mags
                                      + [params_at_freq[1] * show_n_sigma + params_at_freq[0] for params_at_freq in total_mag_best_fit_params_by_frequency])
        if total_mag_coef_range[0] < disp_y_limits[0]: total_mag_coef_range[0] = disp_y_limits[0]
        if total_mag_coef_range[1] > disp_y_limits[1]: total_mag_coef_range[1] = disp_y_limits[1]
        print 'total_mag_coef_range = ' + str(total_mag_coef_range) 
    else:
        sin_coef_range = disp_y_range
        cos_coef_range = disp_y_range
        total_mag_coef_range = disp_y_range 
    
    if n_xticks is None:
        xticks = None
    else:
        xticks = np.linspace(min(frequencies), max(frequencies), n_xticks)
    if show or save: 
        if compute_total_coef_at_freq:
            im_for_colorbar = showParallel1DTempPlots(frequencies, total_mag_coef_range, total_mag_funct_by_frequency, axarr = axarr, axarr_indeces = [0,0],
                                                  white_space = white_space, n_ys_in_range = n_ys_in_range, show = 0, save = 0,
                                                  xlabel = xlabels[0], ylabel = ylabels[0], title = titles[0], 
                                                  levels = levels, colormap = colormap, xticks = xticks, show_colorbar = 0, return_single_plot = 1)
            if (show_n_sigma > 0.0 and not(display_direct_cpbs)):
                axarr[0,0].scatter(frequencies, [params_at_freq[1] * show_n_sigma + params_at_freq[0] for params_at_freq in total_mag_best_fit_params_by_frequency], c = 'r')
            axarr[0,0].scatter(frequencies, true_total_coef_mags, c = 'white')
            if resample_from_data:
                axarr[0,0].scatter(frequencies, funct_single_draw_total_coef_mags, c = 'black')
            
        else: 
            showParallel1DTempPlots(frequencies, sin_coef_range, sin_mag_funct_by_frequency, axarr = axarr, axarr_indeces = [0,0],
                                white_space = white_space, n_ys_in_range = n_ys_in_range, show = 0, save = 0,
                                xlabel = xlabels[0], ylabel = ylabels[0], title = titles[0], 
                                levels = levels, colormap = colormap, xticks = xticks, show_colorbar = 0)
            im_for_colorbar = showParallel1DTempPlots(frequencies, cos_coef_range, cos_mag_funct_by_frequency, axarr = axarr, axarr_indeces = [1,0],
                                                  white_space = white_space, n_ys_in_range = n_ys_in_range, show = 0, save = 0,
                                                  xlabel = xlabels[1], ylabel = ylabels[1], title = titles[1], 
                                                  levels = levels, colormap = colormap, xticks = xticks, show_colorbar = 0, return_single_plot = 1)
            axarr[0,0].scatter(frequencies, true_sin_coefs, c = 'white')
            axarr[1,0].scatter(frequencies, true_cos_coefs, c = 'white')

            if resample_from_data:
                axarr[0,0].scatter(frequencies, funct_single_draw_sin_coefs, c = 'black')
                axarr[1,0].scatter(frequencies, funct_single_draw_cos_coefs, c = 'black')

        f.subplots_adjust(right=0.9)
        cbar_ax = f.add_axes([0.92, 0.08, 0.05, 0.7])
        f.colorbar(im_for_colorbar, cax = cbar_ax)
        f.suptitle(suptitle)
        
    if show:
        plt.show()
        
    if save:
        plot_dir = dir_archive.getPlotDirectory() 
        plt.savefig(plot_dir + 'fittedFunctionPlots/fourierFittedPlots/' + 'fourierSpecOfData_' + ind_var_name + '_vs_' + dep_var_name + '_over_' + str(n_redraw) + '_rDrawsFromFunct_' + funct_name + '_' + other_plot_name + '.png') 

    if return_mags:
        return [fourier_freqs, fourier_mags]
    else:
        return 1
