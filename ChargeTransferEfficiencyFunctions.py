from astropy.io import fits
import numpy as np
import cantrips as c
import time
import scipy.optimize as optimize 
import PISCO_Basic_Reduction_V1_0 as piscobr
import matplotlib.pyplot as plt
from scipy import signal


def makeCTETailLookUpTable(line_number, max_downstream_contributions):
    lookUpTable ={j:c.nChooser(line_number+j, line_number) for j in range(0, max_downstream_contributions)}
    return lookUpTable
    
def computeContributionOfSingleLine(line, ctes, line_number, n_lines, max_downstream_contributions, round_to_int = 1):
    #top_contribution = [[0.0 for elem in line] for j in range(line_number)]
    lookUpTable = makeCTETailLookUpTable(line_number, max_downstream_contributions)
    #print ('lookUpTable = ' + str(lookUpTable)) 
    middle_contribution = [line * (1.0 - ctes) ** j * ctes ** (line_number + 1) * 10.0 ** ( c.safeLog10BigN(lookUpTable[j]))  for j in range(0, min(n_lines - line_number, max_downstream_contributions))]
    #print ('middle_contribution = ' + str(middle_contribution))
    if round_to_int: 
        middle_contribution = np.round(middle_contribution) 
        total_contribution = np.sum(middle_contribution, axis = 0)
        electrons_off = line - total_contribution
        #print ('electrons_off = ' + str(electrons_off))
        electrons_off = 0.0 * electrons_off
        #print ('middle_contribution[0] = ' + str(middle_contribution[0] )) 
        middle_contribution[0] = (middle_contribution[0] + electrons_off).tolist()
        middle_contribution = middle_contribution.tolist()
        
    #bottom_contribution = [[0.0 for elem in line] for j in range(min(n_lines - line_number, max_downstream_contributions), n_lines - line_number)]
    
    #contribution = top_contribution + middle_contribution + bottom_contribution
    contribution = middle_contribution 
    return contribution 


def trueSimulateReadoutWithImperfectCTE(original_image, cte_funct, n_transfers_done = None, readout_axis = 1, output_image = None, readout_image = None, round_to_int = 1, max_downstream_contributions = 10):
    print ('Starting cte...')
    if readout_axis is 0:
        original_image = original_image.transpose()
    shape = np.shape(original_image) 
    readout_axis = 1
    if readout_image is None:
        readout_image = np.zeros(shape)

    new_image = np.copy(original_image)
    start = time.time() 
    for transfer_n in range(shape[0]):
        if transfer_n % 20 ==0:
            end = time.time()
            print ('Took ' + str(end - start) + 's since previous readout.' )
            print ('Simulating transfer ' + str(transfer_n))
            start = time.time()
            print ('np.shape(new_image) = ' + str(np.shape(new_image))) 
        shifted_image = cte_funct(new_image) * new_image 
        if round_to_int: shifted_image = np.round(shifted_image) 
        left_image = new_image - shifted_image
        new_image = left_image[0:-1, :] + shifted_image[1:, :]
        readout_image[transfer_n, :] = shifted_image[0, :]
        #print ('new_image = ' + str(new_image))
        #print ('readout_image = ' + str(readout_image)) 

    return readout_image 
        
    if n_transfers_done is None:
        n_transfers_done = 0 
    readout_image[n_transfers_done, :] = shifted_image[0, :]
        

    #if n_transfers_done == (np.shape(readout_image)[0] - 1):
    #    return readout_image
    #else:
    #    return simulateReadoutWithImperfectCTE(new_image, cte_funct, n_transfers_done = None, readout_axis = 1, output_image = None, readout_image = readout_image, n_transfers_done = 0, round_to_int = 1, max_downstream_contributions = 10)

    return shifted_image 

#Valid only for CONSTANT cte.  Full processing for variable cte should be done using above function 
def fastSimulateReadoutWithImperfectCTE(original_image, cte_funct, n_transfers = None, readout_axis = 1, output_image = None, n_transfers_done = 0, round_to_int = 1, max_downstream_contributions = 10):
    print ('Starting cte...') 
    if readout_axis is 0:
        original_image = original_image.transpose()
    readout_axis = 1 
    original_image = np.array(original_image)
    output_image = original_image * 0.0
    original_image = original_image.tolist()
    output_img = output_image.tolist() 
    if n_transfers == None:
        n_transfers = len(original_image)
    #n_transfers_done = len(original_image) - n_transfers
    
    n_lines = np.shape(original_image)[0]
    for i in range(n_lines):
        if i % 500 == 0: print ('Computing line ' + str(i) + ' of ' + str(n_lines)) 
        original_line = original_image[i]
        ctes = cte_funct(np.array(original_line)) 
        #if i >= 400: print ('ctes.tolist() = ' + str(ctes.tolist()) )
        total_contribution = computeContributionOfSingleLine(original_line, ctes, i, n_lines, max_downstream_contributions, round_to_int = round_to_int)
        #total_contribution = [[0.0 for elem in original_line] for j in range(i)] + [[elem * (1.0 - cte) ** j * (cte) ** (i+1) * c.nChooser(i+j, i) for elem in original_line] for j in range(0, min(n_lines - i, max_downstream_contributions))] + [[0.0 for elem in original_line] for j in range(min(n_lines - i, max_downstream_contributions), n_lines-i)]
        #print ('total_contribution = ' + str(total_contribution)) 
        output_image[i:min(n_lines, i+max_downstream_contributions)] = output_image[i:min(n_lines, i+max_downstream_contributions)] + total_contribution
    #print ('Took ' + str(end - start) + 's') 
    #
    #if output_image is None:
    #    output_image = [[0.0 * elem for elem in original_image[0]] for i in range(n_transfers)]
    #old recursive method 
    #if n_transfers == 0:
    #    return np.array(output_image) 
    #else:
    #    #print ('original_image = ' + str(original_image) )
    #    #print ('output_image = ' + str(output_image)) 
    #    output_image[n_transfers_done] = [np.ceil(elem * cte) if round_to_int else elem * cte for elem in original_image[0]]
    #    original_image = stepImageForward(original_image, cte, round_to_int = round_to_int)
    #    output_image = simulateReadoutWithImperfectCTE(original_image, cte, n_transfers = n_transfers - 1, readout_axis = readout_axis, output_image = output_image, n_transfers_done = n_transfers_done + 1, round_to_int = round_to_int)
    #    return np.array(output_image )
    return output_image 
    

def stepImageForward(original_image, cte, round_to_int = 1):
    #print ('round_to_int = ' + str(round_to_int)) 
    original_image = np.array(original_image) 
    original_image = original_image.tolist()
    if len(np.shape(original_image[0])) == 0:
        original_image = [[elem] for elem in original_image]
    #original_image = []
    #original_image = np.array(original_image) 
    #original_image = original_image.tolist()
    if len(original_image) == 0:
        return original_image
    if len(np.shape(original_image[0])) > 0: 
        line_len = len(original_image[0])
    else:
        line_len = 1 
    residual_image = original_image.copy() 
    next_contribution = [0.0 for i in range(line_len)]
    for i in c.niceReverse(list(range(-len(residual_image), 0))):
        #print ('i = ' + str(i)) 
        #print ('next_contribution = ' + str(next_contribution)) 
        new_line = (np.floor(np.array(residual_image[i])  * (1-cte)) + np.round(np.array(next_contribution )) if round_to_int else np.array(residual_image[i])  * (1-cte) + np.array(next_contribution )) 
        #print ('new_line = ' + str(new_line) ) 
        next_contribution = (np.ceil(np.array(residual_image[i]) * cte) if round_to_int else (np.array(residual_image[i]) * cte))
        residual_image[i] = new_line.tolist() 
        #print ('residual_image = ' + str(residual_image)) 
    return residual_image

def sumResidualOfImage(residual_image, cte):
    residual_image = np.array(residual_image)
    residual_image = residual_image.tolist()

    #print ('residual_image to be summed = ' + str(residual_image)) 
    total_contribution = [[elem * cte ** (i+1) for elem in residual_image[i]] for i in range(len(residual_image))]
    total_contribution = np.sum(total_contribution, axis = 0) 
    return np.array(total_contribution)

#This is different from stepping an image forward, because here we are interested in the contributions
# of each line to the residual image (the image of pixel values right before the line in question moved through them).  Thus, we only want to compute contribution to subsequent lines 
#def stepResidualImageForward(residual_image, cte):
    
    

#For each row, I need to know its contect before it is hit by the next row.  So I build an image of how the image appeared before it was hit by whatever row is presently being read out.  That I can build iteratively.
#  Currently, this procedure fails because of, I think, digitization noise being exponentially enlarged.  Not sure how to get around this, unfortunately... 
def precicelyCorrectReadoutInImageArray(image_to_correct, cte, max_downstream_contributions = 10, readout_axis = 1, round_to_int = 1):
    start = time.time() 
    if readout_axis is 0:
        image_to_correct = image_to_correct.transpose() 
    real_image = np.zeros(np.shape(image_to_correct)) 
    current_step = 0
    n_lines = np.shape(image_to_correct)[0]
    residual_image = np.zeros(np.shape(image_to_correct))
    for i in range(n_lines):
        if i % 10 == 0: print ('CORRECTING LINE ' + str(i)) 
        #print ('image_to_correct = ' + str(image_to_correct))
        #print ('real_image = ' + str(real_image)) 
        #print ('residual_image = ' + str(residual_image)) 
        original_line = image_to_correct[i]
        #if i > 0: 
        #    residual_contribution = simulateReadoutWithImperfectCTE(real_image, cte, n_transfers = i+1, readout_axis = readout_axis, round_to_int = round_to_int)[i]
        #    #addition_from_transfer = sumResidualOfImage(residual_contribution, cte)
        #else:
        #    residual_contribution = 0.0 * original_line
        residual_line = residual_image[i]
        real_line_scaled_by_cte = original_line - residual_line
        real_line_scaled_by_cte[real_line_scaled_by_cte < 0] = 0 
        real_line = (np.round(real_line_scaled_by_cte * (1.0 / cte) ** (i+1)) if round_to_int else real_line_scaled_by_cte * (1.0 / cte) ** (i+1))
        #for j in range(i+1):
        #    real_line_scaled_by_cte = (np.round(real_line_scaled_by_cte * (1.0 / cte)) if round_to_int else real_line_scaled_by_cte * (1.0 / cte)) 
        #real_line = real_line_scaled_by_cte
        #print ('original_line = ' + str(original_line))
        #print ('residual_image = ' + str(residual_image)) 
        #print ('residual_line = ' + str(residual_line))
        #print ('real_line_scaled_by_cte = ' + str(real_line_scaled_by_cte)) 
        #print ('real_line = ' + str(real_line) )
        real_image[i] = real_line
        residual_contribution = computeContributionOfSingleLine(real_line, cte, i, n_lines, max_downstream_contributions = max_downstream_contributions, round_to_int = round_to_int)
        #if i > 0: 
        #    residual_contribution_above = [[0.0 for j in range(len(original_line))] for earlier_line in range(i)]
        #    residual_contribution_below = simulateReadoutWithImperfectCTE(np.array([real_line] + [[0.0 for j in range(len(original_line))] for later_lines in range(i+1, n_lines)]), cte).tolist()
        #    #print ('[residual_contribution_above, residual_contribution_below] = ' + str([residual_contribution_above, residual_contribution_below])  ) 
        #    residual_contribution = np.array(residual_contribution_above + residual_contribution_below) 
        #else:
        #    residual_contribution = simulateReadoutWithImperfectCTE(np.array([real_line] + [[0.0 for j in range(len(original_line))] for later_lines in range(i+1, n_lines)]), cte)
        #print ('[residual_contributionA, residual_contributionB] = ' + str([residual_contributionA, residual_contributionB]))
        #print ('residual_contribution = ' + str(residual_contribution)) 
        residual_image[i:min(n_lines, i+max_downstream_contributions)] = residual_image[i:min(n_lines, i+max_downstream_contributions)] + residual_contribution
        #current_step = current_step + 1
        #residual_image = stepImageForward(residual_image, cte)
        #for j in range(i+1): 
        #    residual_image[j] = residual_image[j] + real_line * (cte) ** (i-j) * (1.0 - cte)
        #print ('residual_image before walk = ' + str(residual_image)) 

    if readout_axis is 0:
        image_to_correct = np.transpose(image_to_correct)
        real_image = np.transpose(real_image)

    end = time.time()
    print ('took ' + str(end - start) + 's to correct image of shape ' + str(np.shape(image_to_correct))) 
    return real_image 


#Based on the iterative method used for the HST outlined in https://arxiv.org/abs/0909.0507

def perturbativeCorrectCTI(image_to_correct, cte_funct, n_iterations = 2, readout_axis = 1, round_to_int = 1, max_downstream_contributions = 10):

    current_image = image_to_correct.copy() 
    for iter in range(n_iterations):
        more_overscanned_image = trueSimulateReadoutWithImperfectCTE(current_image, cte_funct, readout_axis = readout_axis, output_image = None, round_to_int = round_to_int, max_downstream_contributions = max_downstream_contributions)
        current_image = image_to_correct + current_image - more_overscanned_image

    return current_image 

def measurePISCOCTE(imgs_list, img_strs = None, 
                    target_dir ='', bias_file = 'BIAS.fits', already_cropped = 0, PISCO_mosaic_extensions = 8, binning = 1,
                    lowx_bin_2x2 = 40, highx_bin_2x2 = 1510, lowy_bin_2x2 = [365,429,300,375], highy_bin_2x2 = [2505,2569,2440,2515],
                    correlation_length = 1, n_max_for_correllation = 200, figsize = [15, 5], bands = ['g','r','i','z'], read_in_from_file = 1, 
                    save_file_name = 'cteTransferCurves.png', save_arrays_to_fits = 0, save_cte_fig = 0, show_cte_fig = 0, direction_flips = [[0,0],[1,0],[0,1],[1,1]] ):
    if img_strs is None:
        img_strs = ['' for img in imgs_list]
    if binning == 1:
        lowx = lowx_bin_2x2 * 2
        highx = highx_bin_2x2 * 2 
        lowy = [low * 2 for low in lowy_bin_2x2]
        highy = [high * 2 for high in highy_bin_2x2]
    else:
        lowx = lowx_bin_2x2
        highx = highx_bin_2x2
        lowy = lowy_bin_2x2
        highy = highy_bin_2x2
 
    if already_cropped: n_mosaic_extensions = PISCO_mosaic_extensions
    else: n_mosaic_extensions = 0
    #Need to work on images one at a time, as we otherwise will run out of memory
    img_medians = [[0,0,0,0, 0.0, 0.0] for img in imgs_list]
    #horizontal_shifted_imgs = [[[] for band in bands] for img_str in imgs_list]
    #vertical_shifted_imgs = [[[] for band in bands] for img_str in imgs_list]
    connected_deviations_by_row = [[[] for band in bands] for img in imgs_list]
    connected_deviations_by_col = [[[] for band in bands] for img in imgs_list]
    for i in range(len(imgs_list)):
        img_str = img_strs[i] 
        if read_in_from_file:
            img_file = imgs_list[i] 
            if bias_file is None:
                imgs, header = c.readInDataFromFitsFile(img_file, target_dir, n_mosaic_image_extensions = PISCO_mosaic_extensions)
                imgs = piscobr.PISCOStitch(imgs) 
            else:
                imgs, header = piscobr.OverscanCorrect(img_file, target_dir, binning, 
                                                      oscan_prefix = 'os_', return_data = 1,
                                                      save_data = 0, overscan_fit_order = 1,
                                                      overscan_buffer = 10)
                imgs = piscobr.PISCOStitch(imgs) 
                bias = piscobr.loadimageFB(target_dir + bias_file)
                imgs = imgs - bias
            cropped_imgs = [imgs[j][lowy[j]:highy[j], lowx:highx] for j in range(len(imgs))]
        else:
            cropped_imgs = [imgs_list[i][j][lowy[j]:highy[j], lowx:highx] for j in range(len(imgs_list[i]))]
            imgs = imgs_list 
        cropped_imgs = [np.transpose(img) for img in cropped_imgs]
        #for img in cropped_imgs:
        #    plt.imshow(img)
        #    plt.show() 
        cropped_medians = [np.median(img) for img in cropped_imgs]
        img_medians[i] = cropped_medians
        cropped_cte_measurements = [c.measureImageCorrelations(img, correlation_pixel_length = correlation_length,  n_max_for_correllation = n_max_for_correllation) for img in cropped_imgs]
        if save_arrays_to_fits: 
            [c.saveDataToFitsFile(cropped_imgs[j], 'TEST_CTE_CROPPED_' + bands[j] + img_str, target_dir) for j in range(len(cropped_imgs))]
        #print ('[cropped_cte_measurements[2][0][0,0],cropped_cte_measurements[0][1][0,0]] = ' + str([cropped_cte_measurements[0][0][0,0],cropped_cte_measurements[0][1][0,0]])) 
            [[c.saveDataToFitsFile(cropped_cte_measurements[j][0], 'TEST_CTE_HORIZ_' + bands[j] + img_str, target_dir),c.saveDataToFitsFile(cropped_cte_measurements[j][1], 'TEST_CTE_VERTICAL_' + ['g','r','i','z'][j] + img_str, target_dir)]  for j in range(len(cropped_imgs))]
        #print('[cropped_cte_measurements[0][0][3555 - 731, 267 - 21], cropped_cte_measurements[0][1][3555 - 731, 267 - 21]] = ' + str([cropped_cte_measurements[0][0][3555 - 731, 267 - 21], cropped_cte_measurements[0][1][3555 - 731, 267 - 21]])) 
        print ('Computing means of cropped images...')
        #horizontal_shifted_imgs[i] = [cropped_cte_measurement[0] for cropped_cte_measurement in cropped_cte_measurements] 
        #vertical_shifted_imgs[i] = [cropped_cte_measurement[1] for cropped_cte_measurement in cropped_cte_measurements] 
        connected_deviations_by_row[i] = [cropped_cte_measurement[2] for cropped_cte_measurement in cropped_cte_measurements] 
        connected_deviations_by_col[i] = [cropped_cte_measurement[3] for cropped_cte_measurement in cropped_cte_measurements] 

    img_medians, connected_deviations_by_row, connected_deviations_by_col, sorted_img_strs = c.safeSortOneListByAnother([med[-1] for med in img_medians], [img_medians, connected_deviations_by_row, connected_deviations_by_col, img_strs])
    f, axarr = plt.subplots(len(img_strs), len(bands), figsize = figsize, sharex = True, sharey = True, squeeze = False)
    plt.subplots_adjust(hspace = 0.0, wspace = 0.0)
    [axarr[0][band_num].set_xlabel('Pixel') for band_num in range(len(bands))]
    [axarr[-1][band_num].set_xlabel('Pixel') for band_num in range(len(bands))] 
    fit_funct = lambda pix, a, p, b: np.array(pix / a + b) ** p
    all_fits = [[[] for band in bands] for i in range(len(imgs_list))]
    for i in range(len(sorted_img_strs)):
        img_str = sorted_img_strs[i]
        axarr[i][0].set_ylabel(img_str)
        for band_num in range(len(bands)):
            n_rows = len(connected_deviations_by_row[i][band_num])
            n_cols = len(connected_deviations_by_col[i][band_num])
            print ('[n_cols, n_rows] = ' + str([n_cols, n_rows]))
            direction_flip = direction_flips[band_num]
            by_row_to_plot = connected_deviations_by_row[i][band_num]
            if direction_flip[1]: by_row_to_plot.reverse() 
            by_col_to_plot = connected_deviations_by_col[i][band_num]
            if direction_flip[0]: by_col_to_plot.reverse() 
            axarr[i][band_num].plot(range(n_rows), c.smoothList(by_row_to_plot, params = [20], averaging = 'mean'), c = 'b')
            #print ('n_rows = ' + str(n_rows))
            #print ('connected_deviations_by_row[i][band_num] = ' + str(connected_deviations_by_row[i][band_num])) 
            #by_row_fit = optimize.curve_fit(fit_funct, np.array(list(range(n_rows))), connected_deviations_by_row[i][band_num], p0 = [n_rows, 0.5, 0.1], maxfev = 2000)
            by_row_fit = np.polyfit(np.array(list(range(n_rows))), by_row_to_plot, 2)
            all_fits[i][band_num] = by_row_fit 
            print ('[a0,a1,a2] = ' + str(by_row_fit))
            axarr[i][band_num].plot(range(n_rows), np.poly1d(by_row_fit)(np.array(list(range(n_rows)))), c = 'r')
            #axarr[i][band_num].plot(range(n_rows), fit_funct(np.array(list(range(n_rows))), *by_row_fit[0]), c = 'r')
            #axarr[i][band_num].plot(range(n_rows), fit_funct(np.array(list(range(n_rows))), *[n_rows, 0.5, 0.1]), c = 'orange') 
            #by_col_plot = axarr[i][band_num].plot(range(len(connected_deviations_by_col[i][band_num])), connected_deviations_by_col[i][band_num], c = 'r')
            axarr[i][band_num].plot(range(n_cols), c.smoothList(by_col_to_plot, params = [20], averaging = 'mean'), c = 'g')
            axarr[i][band_num].text(0, 1.0, r'Median ADU$\simeq$' + str(c.round_to_n(img_medians[i][band_num], 3)), fontsize = 8.0)
            axarr[i][band_num].text(0, 0.9, r'[a2,a1,a0]$\simeq$' + str([c.round_to_n(term, 3) for term in by_row_fit]), fontsize = 8.0)
            axarr[i][band_num].set_ylim(-0.05, 1.15)
            if i == 0: axarr[i][band_num].set_title('PISCO ' + bands[band_num])
    plt.tight_layout() 
    if save_cte_fig:  
        f.savefig(save_file_name)
    if show_cte_fig:
        plt.show()
    plt.close('all') 
    return all_fits 
    
