import matplotlib.pyplot as plt
import numpy as np

#bins can be either a set of bin edges OR an integer number of bins
# if instructed to match bounds and an integer is given, then will
# generate a matching set of bins according to full range of input arrays. 
def plotGroupedHistograms(data_sets, axarr = None, axarr_hist_index = None, bins = 20, bounds = None, labels = [], match_bounds = 1, layout_dim = None, show = 0, save = 0, fig_size = None,
                          xlabels = '', right_ylabels = '', left_ylabels = '', titles = '', add_arrow = 1):

    if axarr is None:
        f, axarr = plt.subplots(1,1, squeeze = False)
    if axarr_hist_index is None:
        axarr_hist_index = 0
    
    default_labeling_specs = {'val_y_pos': 0.0, 'label_x_offset':0.0, 'label_y_pos':40.0, 'label_color':'r', 'arrow_color':'red', 'arrow_shrink':0.05 }
    for i in range(len(data_sets)):
        data_set = data_sets[i]
        if type(data_set) is type(np.array((1))):
            data_sets[i] = data_set.tolist() 
    if not(type(data_sets[0]) is list):
        print 'data_sets must be a list of lists.  Attempting to correct by putting given list into a list.'
        data_sets = [data_sets]
    if match_bounds and type(bins) is int:
        if bounds is None:
            max_val = max(data_sets[0]) 
            min_val = min(data_sets[0]) 
            for data_set in data_sets[1:]:
                new_min_val = min(data_set)
                new_max_val = max(data_set)
                if new_min_val < min_val: min_val = new_min_val
                if new_max_val > max_val: max_val = new_max_val
            bounds = [min_val, max_val]
        bins = np.linspace(bounds[0], bounds[1], bins+1) 
        
    if layout_dim is None: layout_dim = [len(data_sets),2]

    #f, axarr = plt.subplots(layout_dim[0], layout_dim[1], squeeze = False)

    hist_peaks = []
    for i in range(len(data_sets)):
        data_set = data_sets[i]
        hist_counts, hist_bins, hist_patches = axarr[i, axarr_hist_index].hist(data_set, bins = bins)
        hist_peaks = hist_peaks + [max(hist_counts) ]

        if type(titles) is str: title = titles
        else: title = titles[i]
        if type(xlabels) is str: xlabel = xlabels
        else: xlabel = xlabels[i]
        if type(left_ylabels) is str: ylabel = left_ylabels
        else: ylabel = left_ylabels[i]
        if type(right_ylabels) is str: right_ylabel = right_ylabels
        else: right_ylabel = right_ylabels[i] 
        if i ==0: 
            axarr[i, axarr_hist_index].set_title(title)
        if i == len(data_sets) - 1:
            axarr[i, axarr_hist_index].set_xlabel(xlabel)
        axarr[i, axarr_hist_index].set_ylabel(ylabel) 
        if not(type(bins) is int):
            axarr[i, axarr_hist_index].set_xlim([min(bins), max(bins)])
        extra_ax = axarr[i, axarr_hist_index].twinx()
        extra_ax.get_yaxis().set_ticks([])
        extra_ax.set_ylabel(right_ylabel, rotation = 270, labelpad = 15 )

        if type(bins) is list or type(bins) is np.ndarray:
            bin_right_wall = max(bins)
            bin_left_wall = min(bins)
            n_above = sum([1 for elem in data_set if elem > bin_right_wall])
            n_below = sum([1 for elem in data_set if elem < bin_left_wall]) 
            if n_above > 0:
                above_specs = {'val_x_pos':bin_right_wall,'val_y_pos': hist_peaks[i] * 0.9, 'label_x_offset': -(bin_right_wall - bin_left_wall) / 8.0,
                               'label_y_pos': hist_peaks[i] * 0.9, 'label_color':'g', 'arrow_color':'g', 'arrow_shrink':0.05}
                
                axarr[i, axarr_hist_index].annotate('(' + str(n_above) + ' larger vals)', xy = (above_specs['val_x_pos'], above_specs['val_y_pos']),
                                                    xytext = (above_specs['val_x_pos'] + above_specs['label_x_offset'], above_specs['label_y_pos'] ),
                                                    color = above_specs['label_color'], arrowprops = dict(facecolor = above_specs['arrow_color'], shrink = above_specs['arrow_shrink']),
                                                        horizontalalignment='right', verticalalignment = 'center' )
            if n_below > 0:
                below_specs = {'val_x_pos':bin_left_wall,'val_y_pos': hist_peaks[i] * 0.9, 'label_x_offset': (bin_right_wall - bin_left_wall) / 8.0,
                               'label_y_pos': hist_peaks[i] * 0.9, 'label_color':'g', 'arrow_color':'g', 'arrow_shrink':0.05}

                axarr[i, axarr_hist_index].annotate('(' + str(n_below) + ' smaller vals)', xy = (below_specs['val_x_pos'], below_specs['val_y_pos']),
                                                    xytext = (below_specs['val_x_pos'] + below_specs['label_x_offset'], below_specs['label_y_pos']),
                                                    color = below_specs['label_color'], arrowprops = dict(facecolor = below_specs['arrow_color'], shrink = below_specs['arrow_shrink']),
                                                        horizontalalignment='left', verticalalignment = 'center' )
    for label in labels:
        set_index = label['set_to_label']
        val = label['val']
        label_text = label['text']
        specs = label['specs']
        if specs is None or specs in ['default', 'use_default', 'def', 'use_def', 'usedefault','useDefault']:
            specs = default_labeling_specs
            specs['label_y_pos'] = hist_peaks[set_index] * 0.2
        if add_arrow: 
            axarr[set_index, axarr_hist_index].annotate(label_text, xy = (val, specs['val_y_pos']), xytext = (val + specs['label_x_offset'], specs['label_y_pos']),
                                          color = specs['label_color'], arrowprops = dict(facecolor = specs['arrow_color'], shrink = specs['arrow_shrink'])  )
    if save:
        #Save it
        print 'save it' 
    if show:
        plt.show()
    if not show:
        return axarr
                    
            
    
