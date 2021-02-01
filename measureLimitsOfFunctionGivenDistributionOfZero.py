import numpy as np
import math
from StandardPeriodigram import StandardPeriodigram
import scipy.optimize as optimize
import scipy.special as special 
import matplotlib.pyplot as plt
from showParallel1DTempPlots import showParallel1DTempPlots
plt.rc('font', family='serif')
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-small')

def plotFunctionVsDistributionOfZero(xs, y_errs, function_to_limit, specific_params_for_funct, mag_limit_params, limits_file, 
                                     axarr = None, axarr_indeces = (0,0), return_fig = 0, show_periodogram_of_data = 1, 
                                     white_space = 0.0, n_ys_in_range = 100, n_sigma_limits_colors = ['g','b','m','r','orange'], freq_range_to_show = 'all',
                                     normalization_coef = 1.0 / (276.509398007 * 2.0), n_xticks = 10, xtick_scaling = 1.0, xtick_labels = [], tick_label_size = 14.0, 
                                     funct_vals_color = 'goldenrod', levels = np.linspace(0.0, 1.0, 20), colormap = None, vert_lines_freq_indeces = [], 
                                     apply_normalization_to_periodigram = 1, apply_normalization_to_params = 0, fitting_funct_type = 'alt_normal', show_fourier_amplitude = 0, 
                                     xlabels = ['Periodogram frequency scaled by ' + r'$h$, ' + r'$f_n$' + ' (1/Gyr)'], ylabels = ['Periodogram Value Approximating Fourier Power'], titles = [''], label_size = 18.0, 
                                     show_temp_plots = 0, compute_minimization_value = 0, figsize = [16,8], smoothing_integer = 1, 
                                     show = 1, save = 0, save_name = '/Users/sasha/Desktop/PeriodigramVsPeriodigramFrom0.png', 
                                     transparency = 0.5):
    results_for_minimization = np.genfromtxt(limits_file, delimiter = ',').transpose()
    frequencies = results_for_minimization[0]

    smoothing_integer = int(smoothing_integer)

    if freq_range_to_show in ['all','ALL','All','full','FULL','full',None]:
        freq_range_to_show = [0, len(frequencies)]
    
    if fitting_funct_type in ['normal', 'norm','guass','gaussian']: 
        mus = results_for_minimization[2]
        sigmas = results_for_minimization[3]
        if apply_normalization_to_params:
            mus = [mu * normalization_coef for mu in mus]
            sigmas = [sig * normalization_coef for sig in sigmas]
        n_sigma = mag_limit_param
        mag_limits_by_limit_param = [mus[i] + n_sigma * sigmas[i] for i in range(len(mus))]
        
    elif fitting_funct_type in ['alt_normal', 'alt_norm','alt_guass','alt_gaussian']:
        alphas = results_for_minimization[2]
        betas = results_for_minimization[3]
        if apply_normalization_to_params:
            alphas = [alpha * normalization_coef for alpha in alphas]
        mag_limits_by_limit_params = []
        smoothed_mag_limits_by_limit_params = []
        smoothing_indeces = [[max(0, i - smoothing_integer // 2), min(len(betas), i + smoothing_integer // 2 + smoothing_integer % 2)]
                                 for i in range(freq_range_to_show[0], freq_range_to_show[1])]
        
        for mag_limit_param in mag_limit_params:
            
            mag_limits_by_limit_param = []
            for i in range(int(np.min(smoothing_indeces)), int(np.max(smoothing_indeces)) ):
                limit_param = optimize.brenth(lambda lim: special.gammainc(1.0 / betas[i], (lim/alphas[i]) ** betas[i]) - mag_limit_param, 0.0, 1.0)
                mag_limits_by_limit_param = mag_limits_by_limit_param + [limit_param]
            #print [(min(len(mag_limits_by_limit_param), int(i + smoothing_integer / 2)) - max(0, int(i - smoothing_integer / 2))) for i in range(len(mag_limits_by_limit_param))]
            
            effective_zero_index = np.min(smoothing_indeces)
            smoothed_mag_limits = [( sum(mag_limits_by_limit_param[smoothing_index[0]- effective_zero_index:smoothing_index[1] - effective_zero_index])
                                      / (smoothing_index[1] - smoothing_index[0]) )
                                   for smoothing_index in smoothing_indeces]
            mag_limits_by_limit_params = mag_limits_by_limit_params + [mag_limits_by_limit_param]
            smoothed_mag_limits_by_limit_params = smoothed_mag_limits_by_limit_params + [smoothed_mag_limits]

    

    funct_vals = function_to_limit (xs, specific_params_for_funct)


    if fitting_funct_type in ['normal', 'norm','guass','gaussian']: 
        total_mag_funct_by_frequency = [lambda xs, i = i: np.exp(- (np.array(xs) - mus[i]) ** 2.0 / (2.0 * sigmas[i] ** 2.0))
                                          for i in range(len(mus)) ]
    elif fitting_funct_type in ['alt_normal', 'alt_norm','alt_guass','alt_gaussian']:
        #total_mag_funct_by_frequency = [lambda xs, i = i: [special.gammainc(1.0 / betas[i], (x/alphas[i]) ** betas[i]) for x in xs]
        #                                  for i in range(len(mus)) ]
        #total_mag_funct_by_frequency = [lambda xs, i=i: [betas[i]  / (4.0 * alphas[i] * special.gamma(1.0 / betas[i])) * np.exp(-(np.array(xs)/alphas[i])**betas[i] ) for x in xs]
        #                                  for i in range(len(mus)) ]
        total_mag_funct_by_frequency = [lambda xs, i=i: [np.exp(-(np.array(x)/alphas[i])**betas[i] ) for x in xs]
                                          for i in range(len(alphas)) ]

    if freq_range_to_show is 'all':
        select_freqs = frequencies[:]
        #select_mus = mus[:]
        #select_sigmas = sigmas[:]
        select_mag_limits_by_limit_params = [mag_limits_by_limit_param [:] for mag_limits_by_limit_param in mag_limits_by_limit_params]
        select_smoothed_mag_limits_by_limit_params  = [smoothed_mag_limits[:] for smoothed_mag_limits in smoothed_mag_limits_by_limit_params]
        select_total_mag_funct_by_frequency  = total_mag_funct_by_frequency[:]

        #select_normalized_mags = funct_periodigram.normalized_coef_mags[:] 
    else:
        select_freqs = frequencies[freq_range_to_show[0]:freq_range_to_show[1]]
        #select_mus = mus[freq_range_to_show[0]:freq_range_to_show[1]]
        #select_sigmas = sigmas[freq_range_to_show[0]:freq_range_to_show[1]]
        select_total_mag_funct_by_frequency  = total_mag_funct_by_frequency[freq_range_to_show[0]:freq_range_to_show[1]]

        #select_mag_limits_by_limit_params = [mag_limits_by_limit_param [freq_range_to_show[0]:freq_range_to_show[1]] for mag_limits_by_limit_param in mag_limits_by_limit_params]
        select_mag_limits_by_limit_params = mag_limits_by_limit_params
        #select_smoothed_mag_limits_by_limit_params = [smoothed_mag_limits[freq_range_to_show[0]:freq_range_to_show[1]] for smoothed_mag_limits in smoothed_mag_limits_by_limit_params]
        select_smoothed_mag_limits_by_limit_params =  smoothed_mag_limits_by_limit_params

        #select_normalized_mags = funct_periodigram.normalized_coef_mags[freq_range_to_show[0]:freq_range_to_show[1]]
    print ('len(select_total_mag_funct_by_frequency) = ' + str(len(select_total_mag_funct_by_frequency)))

    select_smoothed_fourier_limits_by_limit_params = [[(limit * 2.0) ** 0.5 for limit in smoothed_mag_limits] for smoothed_mag_limits in select_smoothed_mag_limits_by_limit_params]

    partial_periodigram = StandardPeriodigram(xs, funct_vals, y_errs,
                                              frequencies = select_freqs, apply_normalization = apply_normalization_to_periodigram, compute_convolution = 0)
    select_normalized_mags = partial_periodigram.normalized_coef_mags
    select_normalized_fourier_modes = [(mag * 2.0) ** 0.5 for mag in select_normalized_mags]

    print ('len(select_normalized_mags) = ' + str(len(select_normalized_mags)))
    print ('len(select_mag_limits_by_limit_params) = ' + str(len(select_mag_limits_by_limit_params)))
    print ('np.shape(mag_limits_by_limit_params) = ' + str(np.shape(mag_limits_by_limit_params)))
    print ('np.shape(select_mag_limits_by_limit_params) = ' + str(np.shape(select_mag_limits_by_limit_params)))
    total_mag_coef_range = [0.0, max(max(select_normalized_mags), np.max(select_mag_limits_by_limit_params + []))]
    total_fourier_coef_range = [0.0, max(max(select_normalized_fourier_modes), np.max(select_smoothed_fourier_limits_by_limit_params + []), 0.045)]
    if (show or save or return_fig):
        if axarr is None: 
            f, axarr = plt.subplots(1, 1, squeeze = False, figsize = figsize)
            axarr_indeces = (0,0)
        
        freqs_to_show = [freq / (2.0 * np.pi) for freq in select_freqs]
        #freqs_to_show = [freq for freq in select_freqs]
        
        if show_temp_plots:
            artificial_functs = [lambda xs: np.array(xs) ** 0.0,  lambda xs: np.array(xs) ** 1.0, lambda xs: np.array(xs) ** 2.0, lambda xs: np.array(xs) ** 3.0, lambda xs: np.array(xs) ** 4.0]
            #showParallel1DTempPlots(freqs_to_show, total_mag_coef_range, artificial_functs, levels = levels, show_colorbar = 0, n_ys_in_range = 100)
            #im_for_colorbar = showParallel1DTempPlots(freqs_to_show, total_mag_coef_range, artificial_functs, axarr = axarr, axarr_indeces = [0,0],
            #                                          n_ys_in_range = n_ys_in_range, show = 0, save = 0,
            #                                          levels = levels, show_colorbar = 0, colormap = colormap)
            #im_for_colorbar = showParallel1DTempPlots(freqs_to_show, total_mag_coef_range, artificial_functs, axarr = axarr, axarr_indeces = [0,0],
            #                                          white_space = white_space, n_ys_in_range = n_ys_in_range, show = 1, save = 0,
            #                                          xlabel = xlabels[0], ylabel = ylabels[0], title = titles[0], 
            #                                          levels = levels, colormap = colormap, xticks = np.linspace(freqs_to_show[0], freqs_to_show[-1], 10), show_colorbar = 0, return_single_plot = 1)
            #plt.show() 
            im_for_colorbar = showParallel1DTempPlots(freqs_to_show, total_mag_coef_range, select_total_mag_funct_by_frequency, axarr = axarr, axarr_indeces = axarr_indeces,
                                                      white_space = white_space, n_ys_in_range = n_ys_in_range, show = 0, save = 0,
                                                      xlabel = xlabels[0], ylabel = ylabels[0], title = titles[0], 
                                                      levels = levels, colormap = colormap, xticks = np.linspace(freqs_to_show[0], freqs_to_show[-1], 10), show_colorbar = 1, return_single_plot = 1)

        base_for_fill = [0.0 for freq in freqs_to_show]
        for i in range(len(select_smoothed_mag_limits_by_limit_params )):
            if show_fourier_amplitude:
                select_smoothed_fourier_limits = select_smoothed_fourier_limits_by_limit_params[i]
                axarr[axarr_indeces[0],axarr_indeces[1]].fill_between(freqs_to_show, base_for_fill, select_smoothed_fourier_limits , edgecolor = 'black', facecolor = n_sigma_limits_colors[i], alpha = transparency)
                base_for_fill = select_smoothed_fourier_limits
            else:
                select_smoothed_mag_limits = select_smoothed_mag_limits_by_limit_params[i]
                print ('select_smoothed_mag_limits = ' + str(select_smoothed_mag_limits))
                axarr[axarr_indeces[0],axarr_indeces[1]].fill_between(freqs_to_show, base_for_fill, select_smoothed_mag_limits , facecolor = n_sigma_limits_colors[i], alpha = transparency)
                base_for_fill = select_smoothed_mag_limits

        if show_periodogram_of_data: 
            if show_fourier_amplitude:
                axarr[axarr_indeces[0],axarr_indeces[1]].scatter(freqs_to_show, select_normalized_fourier_modes, c = funct_vals_color, marker = '.')
            else: 
                axarr[axarr_indeces[0],axarr_indeces[1]].scatter(freqs_to_show, select_normalized_mags, c = funct_vals_color, marker = '.')
        
 
        axarr[axarr_indeces[0],axarr_indeces[1]].set_xlabel(xlabels[0], labelpad = 15.0, fontsize = label_size)
        axarr[axarr_indeces[0], axarr_indeces[1]].set_ylabel(ylabels[0], labelpad = 10.0, fontsize = label_size)

        if len(xtick_labels) > 0:
            axarr[axarr_indeces[0],axarr_indeces[1]].set_xticks([float(tick) * xtick_scaling for tick in xtick_labels])
            axarr[axarr_indeces[0],axarr_indeces[1]].set_xticklabels(xtick_labels)
        axarr[axarr_indeces[0],axarr_indeces[1]].tick_params(axis='both', which='major', labelsize=tick_label_size)
            
        axarr[axarr_indeces[0],axarr_indeces[1]].set_title(titles[0])
        if show_temp_plots:
            axarr[axarr_indeces[0],axarr_indeces[1]].set_xlim(min(freqs_to_show), max(freqs_to_show))
        else:
            #axarr[0,0].set_xlim(min(freqs_to_show) - 100.0 * (freqs_to_show[1] - freqs_to_show[0]), max(freqs_to_show) + 100.0 * (freqs_to_show[1] - freqs_to_show[0]))
            axarr[axarr_indeces[0],axarr_indeces[1]].set_xlim(min(freqs_to_show), max(freqs_to_show))
        

        if show_fourier_amplitude:
            axarr[axarr_indeces[0],axarr_indeces[1]].set_ylim(total_fourier_coef_range)
        else:
            axarr[axarr_indeces[0],axarr_indeces[1]].set_ylim(total_mag_coef_range)
        

        if len(vert_lines_freq_indeces) > 0:
            for vert_line_pos in vert_lines_freq_indeces:
                freq_for_vert_line = freqs_to_show[vert_line_pos]
                axarr[axarr_indeces[0],axarr_indeces[1]].axvline(x = freq_for_vert_line , color = 'grey', linewidth = 2)

        if show:
            plt.show()
        if save:
            my_dpi = 96
            plt.savefig(save_name, figsize = (figsize[0], figsize[1]), dpi = my_dpi) 
    if compute_minimization_value:
        print ('The value for the limiting function is: ')
        functToMinimize = functionToMinimize(xs, y_errs, function_to_limit, select_freqs, select_mag_limits_by_limit_param, apply_normalization_to_periodigram)
        print (functToMinimize(specific_params_for_funct))

    if return_fig:
        return axarr 

    
    

def difference_between_Periodigram_and_mag_limits(funct_params, xs, y_errs, function_to_limit, frequencies, mag_limits_by_n_sigma, apply_normalization):
    return (np.array(StandardPeriodigram(xs, function_to_limit(xs, funct_params), y_errs, 
                                                      frequencies = frequencies, apply_normalization = apply_normalization, compute_convolution = 0).normalized_coef_mags)
                         - np.array(mag_limits_by_n_sigma))
        
def positive_max_difference_between_Periodigram_and_mag_limits (funct_params, xs, y_errs, function_to_limit, frequencies, mag_limits_by_n_sigma, apply_normalization):
    #print 'Inside function to minimize. '
    #print 'funct_params = ' + str(funct_params)
    return_val = abs(max(difference_between_Periodigram_and_mag_limits(funct_params, xs, y_errs, function_to_limit, frequencies, mag_limits_by_n_sigma, apply_normalization)))
    #print 'value at funct_params = ' + str(return_val)
    return return_val 

def functionToMinimize(xs, y_errs, function_to_limit, frequencies, mag_limits_by_n_sigma, apply_normalization ):

    #square_difference_between_Periodigram_and_mag_limits = lambda funct_params: abs(max(np.array(StandardPeriodigram(xs, function_to_limit(xs, funct_params), y_errs, 
    #                                                                                           frequencies = frequencies, apply_normalization = apply_normalization, compute_convolution = 0).normalized_coef_mags)
    #                                                                                - np.array(mags_limit_by_n_sigma)))
    returnable_max_square_difference = lambda funct_params: positive_max_difference_between_Periodigram_and_mag_limits (funct_params, xs, y_errs, function_to_limit, frequencies, mag_limits_by_n_sigma, apply_normalization)
    return returnable_max_square_difference


def measurePeakCPBValOfFunctionGivenDistributionOfZero(xs, y_errs, function_to_limit, val_to_calc, mag_limit_param, limits_file, file_to_save, 
                                                       apply_normalization_to_periodigram = 1, apply_normalization_to_params = 0,
                                                       normalization_coef = 1.0 / (276.509398007 * 2.0), fitting_funct_type = 'alt_normal', 
                                                       freq_bounds_over_which_to_find_peak = None):
    results_for_minimization = np.genfromtxt(limits_file, delimiter = ',').transpose()
    frequencies = results_for_minimization[0]
    if fitting_funct_type in ['normal', 'norm','guass','gaussian']: 
        mus = results_for_minimization[2]
        sigmas = results_for_minimization[3]
        if apply_normalization_to_params:
            mus = [mu * normalization_coef for mu in mus]
            sigmas = [sig * normalization_coef for sig in sigmas]
        n_sigma = mag_limit_param
        mag_limits_by_limit_param = [mus[i] + n_sigma * sigmas[i] for i in range(len(mus))]
        
    elif fitting_funct_type in ['alt_normal', 'alt_norm','alt_guass','alt_gaussian']:
        alphas = results_for_minimization[2]
        betas = results_for_minimization[3]
        if apply_normalization_to_params:
            alphas = [alpha * normalization_coef for alpha in alphas]
        print ('mag_limit_param = ' + str(mag_limit_param) )
        mag_limits_by_limit_param = []
        for i in range(len(betas)):
            limit_param = optimize.brenth(lambda lim: special.gammainc(1.0 / betas[i], (lim/alphas[i]) ** betas[i]) - mag_limit_param, 0.0, 1.0)

            mag_limits_by_limit_param = mag_limits_by_limit_param + [limit_param]

    if freq_bounds_over_which_to_find_peak is None:
        freq_bounds_over_which_to_find_peak = [0, len(frequencies)]

    periodigram_vs_limits_difference = difference_between_Periodigram_and_mag_limits(val_to_calc, xs, y_errs, function_to_limit,
                                                                                     frequencies[freq_bounds_over_which_to_find_peak[0]:freq_bounds_over_which_to_find_peak[1]],
                                                                                     mag_limits_by_limit_param[freq_bounds_over_which_to_find_peak[0]:freq_bounds_over_which_to_find_peak[1]],
                                                                                     apply_normalization_to_periodigram) 
    freq_index_of_max_n_sigma_difference = np.argmax(periodigram_vs_limits_difference)
    max_n_sigma_difference = max(periodigram_vs_limits_difference)
    

    CMB_val_of_max_n_sigma_difference = 0.0

    return CMB_val_of_max_n_sigma_difference 
    

def measureLimitsOfFunctionGivenDistributionOfZero(xs, y_errs, function_to_limit, initialization_val, mag_limit_param, limits_file, file_to_save, bounds, 
                                                   args = (), apply_normalization_to_periodigram = 1, apply_normalization_to_params = 0,
                                                   normalization_coef = 1.0 / (276.509398007 * 2.0), fitting_funct_type = 'normal', 
                                                   tol = 1e-5, freq_width_to_calc = 100, freq_bounds_over_which_to_find_peak = None, freq_bounds_over_which_to_check_peak = None):
    #constraints look like ([lower_bound_for_var1, lower_bound_for_var2, ...],[upper_bound_for_var1, upper_bound_for_var2, ...])
    results_for_minimization = np.genfromtxt(limits_file, delimiter = ',').transpose()
    frequencies = results_for_minimization[0]
    if fitting_funct_type in ['normal', 'norm','guass','gaussian']: 
        mus = results_for_minimization[2]
        sigmas = results_for_minimization[3]
        if apply_normalization_to_params:
            mus = [mu * normalization_coef for mu in mus]
            sigmas = [sig * normalization_coef for sig in sigmas]
        n_sigma = mag_limit_param
        mag_limits_by_limit_param = [mus[i] + n_sigma * sigmas[i] for i in range(len(mus))]
        
    elif fitting_funct_type in ['alt_normal', 'alt_norm','alt_guass','alt_gaussian']:
        alphas = results_for_minimization[2]
        betas = results_for_minimization[3]
        #print 'alphas = ' + str(alphas)
        #print 'betas = ' + str(betas)
        if apply_normalization_to_params:
            alphas = [alpha * normalization_coef for alpha in alphas]
        print ('mag_limit_param = ' + str(mag_limit_param) )
        mag_limits_by_limit_param = []
        for i in range(len(betas)):
            #print 'i = ' + str(i) 
            #print 'alphas[i] = ' + str(alphas[i])
            #print 'betas[i] = ' + str(betas[i])
            limit_param = optimize.brenth(lambda lim: special.gammainc(1.0 / betas[i], (lim/alphas[i]) ** betas[i]) - mag_limit_param, 0.0, 1.0)

            #print 'limit_param = ' + str(limit_param) 
            mag_limits_by_limit_param = mag_limits_by_limit_param + [limit_param]
    #print 'frequencies = '
    #print frequencies.tolist() 
    #print 'mus = ' + str(mus)
    #print 'sigmas = ' + str(sigmas)

    if freq_bounds_over_which_to_find_peak is None:
        freq_bounds_over_which_to_find_peak = [0, len(frequencies)]
    if freq_bounds_over_which_to_check_peak is None:
        freq_bounds_over_which_to_check_peak = [0, len(frequencies)]

    freq_index_of_initial_max_n_sigma_difference = np.argmax(difference_between_Periodigram_and_mag_limits(initialization_val, xs, y_errs, function_to_limit,
                                                                                                           frequencies[freq_bounds_over_which_to_find_peak[0]:freq_bounds_over_which_to_find_peak[1]],
                                                                                                           mag_limits_by_limit_param[freq_bounds_over_which_to_find_peak[0]:freq_bounds_over_which_to_find_peak[1]]
                                                                                                           , apply_normalization_to_periodigram)) + freq_bounds_over_which_to_find_peak[0]
    print ('freq_index_of_initial_max_n_sigma_difference = ' + str(freq_index_of_initial_max_n_sigma_difference) )

    freq_range_on_which_to_minimize = [max(freq_index_of_initial_max_n_sigma_difference  - freq_width_to_calc / 2, 0),
                                       min(freq_index_of_initial_max_n_sigma_difference  + freq_width_to_calc / 2, len(frequencies))
                                       ]
     
    
    select_freqs = frequencies[freq_range_on_which_to_minimize[0]:freq_range_on_which_to_minimize[1]]
    
    #select_mus = mus[freq_range_on_which_to_minimize[0]:freq_range_on_which_to_minimize[1]]
    #select_funct_vals = funct_vals[freq_range_to_show[0]:freq_range_to_show[1]]
    #select_sigmas = sigmas[freq_range_on_which_to_minimize[0]:freq_range_on_which_to_minimize[1]]
    select_mag_limits_by_limit_param = mag_limits_by_limit_param[freq_range_on_which_to_minimize[0]:freq_range_on_which_to_minimize[1]]
    #select_total_mag_funct_by_frequency  = total_mag_funct_by_frequency[freq_range_to_show[0]:freq_range_to_show[1]]
    

    #print 'mag_limits_by_limit_param = ' + str(mag_limits_by_limit_param) 
    functToMinimizeOverSmallFreqRange = functionToMinimize(xs, y_errs, function_to_limit, select_freqs, select_mag_limits_by_limit_param, apply_normalization_to_periodigram)
    #print 'functToMinimize(init_guess) = ' + str(functToMinimize(init_guess) )
    #min_res = optimize.minimize(functToMinimize, init_guess, args = args, tol = tol, constraints = constraints)
    min_res = optimize.minimize_scalar(functToMinimizeOverSmallFreqRange, method = 'bounded', bounds = bounds, args = args)

    min_param_val = min_res['x']

    freq_index_of_final_max_n_sigma_difference = np.argmax(difference_between_Periodigram_and_mag_limits(min_param_val, xs, y_errs, function_to_limit,
                                                                                                         frequencies[freq_bounds_over_which_to_check_peak[0]:freq_bounds_over_which_to_check_peak[1]],
                                                                                                           mag_limits_by_limit_param[freq_bounds_over_which_to_check_peak[0]:freq_bounds_over_which_to_check_peak[1]]
                                                                                                           , apply_normalization_to_periodigram)) + freq_bounds_over_which_to_find_peak[0]
    print ('freq_index_of_final_max_n_sigma_difference = ' + str(freq_index_of_initial_max_n_sigma_difference)    )

    
    return min_res 
    


class ComparatorBetweenFunctionAndDistributionOfZero:

    def computeLikelihoodsOfSignalFromZeros(self, xs, ys, y_errs,
                                            freq_bounds_over_which_to_find_peak = None, apply_normalization = 1, normalization_coef = 1.0 / (276.509398007 * 2.0)):
        if freq_bounds_over_which_to_find_peak in ['all','ALL','All','full','FULL','full',None]:
            freq_bounds_over_which_to_find_peak = [0, len(self.frequencies)]
        freqs_to_check = self.frequencies[freq_bounds_over_which_to_find_peak[0]: freq_bounds_over_which_to_find_peak[1]]
 
        SP = StandardPeriodigram(xs, ys, y_errs,
                                 frequencies = freqs_to_check, apply_normalization = apply_normalization)
        normalized_coefs = SP.normalized_coef_mags
        relevant_functs = self.CPB_functs[freq_bounds_over_which_to_find_peak[0]: freq_bounds_over_which_to_find_peak[1]]
        zero_CPB_vals_from_funct_SP = [relevant_functs[i](normalized_coefs[i]) for i in range(len(freqs_to_check))]
        return (normalized_coefs, zero_CPB_vals_from_funct_SP)


    def __init__(self, limits_file,
                 apply_normalization_to_periodigram = 1, apply_normalization_to_params = 0,
                 normalization_coef = 1.0 / (276.509398007 * 2.0), fitting_funct_type = 'alt_normal'):
        self.limits_file = limits_file

        self.results_for_minimization = np.genfromtxt(self.limits_file, delimiter = ',').transpose()
        self.frequencies = self.results_for_minimization[0]
    
        if fitting_funct_type in ['normal', 'norm','guass','gaussian']: 
            self.mus = results_for_minimization[2]
            self.sigmas = results_for_minimization[3]
            if apply_normalization_to_params:
                self.mus = [mu * normalization_coef for mu in mus]
                self.sigmas = [sig * normalization_coef for sig in sigmas]
            
        
        elif fitting_funct_type in ['alt_normal', 'alt_norm','alt_guass','alt_gaussian']:
            self.alphas = self.results_for_minimization[2]
            self.betas = self.results_for_minimization[3]
            if apply_normalization_to_params:
                self.alphas = [alpha * normalization_coef for alpha in self.alphas]
            
            self.CPB_functs = [lambda mag, i=i: (special.gammainc(1.0 / self.betas[i], (mag/self.alphas[i]) ** self.betas[i])) for i in range(len(self.betas))]

