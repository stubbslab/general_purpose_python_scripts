import numpy as np
from cantrips import readInFileLineByLine
from cantrips import recursiveReplaceMultipleSpaces
import scipy.integrate as integrate
from scipy.interpolate import interp1d
from AstronomicalParameterArchive import AstronomicalParameterArchive 

def getStarFixedParams(target_stars = 'all'):
    total_flux_file_common_dir = '/Users/sasha/Documents/Harvard/physics/stubbs/narayanWDs/WDdata/out/'
    all_cal_star_info = {'SDSS-J010322':{'RA':'01:03:22.19', 'Dec':'−00:20:47.73', 'uSDSS_synth':[18.631, 0.005], 'g_synth':[19.067, 0.006],'r_synth':[19.560, 0.005] ,'i_synth':[19.921, 0.005],'total_spectrum_file':total_flux_file_common_dir + 'sdssj010322-20131129-total/sdssj010322-20131129-total_spec_model.dat'},
                         'SDSS-J022817':{'RA':'02:28:17.17', 'Dec':'−08:27:16.41', 'uSDSS_synth':[19.775, 0.011], 'g_synth':[19.806, 0.011] ,'r_synth':[20.163, 0.006],'i_synth':[20.467, 0.007],'total_spectrum_file':total_flux_file_common_dir + 'sdssj022817-20131013-total/sdssj022817-20131013-total_spec_model.dat'},
                         'SDSS-J024854':{'RA':'02:48:54.96', 'Dec':'+33:45:48.30', 'uSDSS_synth':[18.115, 0.007], 'g_synth':[18.357, 0.008] ,'r_synth': [18.738, 0.005],'i_synth': [19.042, 0.003],'total_spectrum_file':total_flux_file_common_dir + 'sdssj024854-20151011.5-total/sdssj024854-20151011.5-total_spec_model.dat'},
                         #'SDSS-J041053':{'RA':'04:10:53.63', 'Dec':'+06:30:27.75', 'uSDSS_synth':[], 'g_synth':[], 'r_synth': [], 'i_synth':[],'total_spectrum_file':},
                         #'WD0554':      {'RA':'05:57:01.30', 'Dec':'−16:35:12.00', 'uSDSS_synth':[], 'g_synth':[], 'r_synth': [], 'i_synth':[],'total_spectrum_file':},
                         'SDSS-J072752':{'RA':'07:27:52.76', 'Dec':'+32:14:16.10', 'uSDSS_synth':[17.570, 0.002], 'g_synth':[17.976, 0.003] ,'r_synth': [18.447, 0.002],'i_synth':[18.797, 0.002],'total_spectrum_file':total_flux_file_common_dir + 'sdssj072752-20150124-total/sdssj072752-20150124-total_spec_model.dat'},
                         'SDSS-J081508':{'RA':'08:15:08.78', 'Dec':'+07:31:45.80', 'uSDSS_synth':[19.358, 0.008] , 'g_synth':[19.700, 0.005],'r_synth':[20.177, 0.006],'i_synth':[20.535, 0.007],'total_spectrum_file':total_flux_file_common_dir + 'sdssj081508-20130214-total/sdssj081508-20130214-total_spec_model.dat'},
                         'SDSS-J102430':{'RA':'10:24:30.93', 'Dec':'−00:32:07.03', 'u_synth':[18.588, 0.009] , 'g_synth':[18.896, 0.010],'r_synth':[19.309, 0.007],'i_synth':[19.631, 0.012],'total_spectrum_file':total_flux_file_common_dir + 'sdssj102430-20170223-total/sdssj102430-20170223-total_spec_model.dat'},
                         'SDSS-J111059':{'RA':'11:10:59.43', 'Dec':'−17:09:54.10', 'u_synth':[17.447, 0.004] , 'g_synth':[17.841, 0.003],'r_synth':[18.305, 0.002] ,'i_synth':[18.653, 0.002],'total_spectrum_file':total_flux_file_common_dir + 'sdssj111059-20150124-total/sdssj111059-20150124-total_spec_model.dat'},
                         'SDSS-J111127':{'RA':'11:11:27.30', 'Dec':'+39:56:28.00', 'u_synth':[17.930, 0.004], 'g_synth':[18.398, 0.003],'r_synth':[18.926, 0.003],'i_synth': [19.307, 0.003],'total_spectrum_file':total_flux_file_common_dir + 'sdssj111127-20150124-total/sdssj111127-20150124-total_spec_model.dat'},
                         'SDSS-J120650':{'RA':'12:06:50.41', 'Dec':'+02:01:42.46', 'u_synth':[18.553, 0.004], 'g_synth':[18.663, 0.004],'r_synth':[19.055, 0.005],'i_synth':[19.377, 0.006],'total_spectrum_file':total_flux_file_common_dir + 'sdssj120650-20130310-total/sdssj120650-20130310-total_spec_model.dat'},
                         'SDSS-J121405':{'RA':'12:14:05.11', 'Dec':'+45:38:18.50', 'u_synth':[17.378, 0.003] , 'g_synth':[17.740, 0.004],'r_synth':[18.227, 0.003],'i_synth':[18.593, 0.002],'total_spectrum_file':total_flux_file_common_dir + 'sdssj121405-20150218-total/sdssj121405-20150218-total_spec_model.dat'},
                         'SDSS-J130234':{'RA':'13:02:34.44', 'Dec':'+10:12:39.01', 'u_synth':[16.619, 0.002], 'g_synth':[17.016, 0.003],'r_synth': [17.503, 0.002],'i_synth':[17.865, 0.002],'total_spectrum_file':total_flux_file_common_dir + 'sdssj130234-20130215-total/sdssj130234-20130215-total_spec_model.dat'},
                         'SDSS-J131445':{'RA':'13:14:45.05', 'Dec':'−03:14:15.64', 'u_synth':[18.683, 0.005] , 'g_synth':[19.080, 0.007],'r_synth':[19.560, 0.006],'i_synth':[19.917, 0.004],'total_spectrum_file':total_flux_file_common_dir + 'sdssj131445-20130309-total/sdssj131445-20130309-total_spec_model.dat'},
                         'SDSS-J151421':{'RA':'15:14:21.27', 'Dec':'+00:47:52.79', 'u_synth':[15.464, 0.002], 'g_synth':[15.694, 0.004],'r_synth': [16.108, 0.002],'i_synth':[16.438, 0.002],'total_spectrum_file':total_flux_file_common_dir + 'sdssj151421-20130317-total/sdssj151421-20130317-total_spec_model.dat'},
                         'SDSS-J155745':{'RA':'15:57:45.40', 'Dec':'+55:46:09.70', 'u_synth':[16.983, 0.002] , 'g_synth':[17.447, 0.002] ,'r_synth': [17.975, 0.002],'i_synth':[18.355, 0.002],'total_spectrum_file':total_flux_file_common_dir + 'sdssj155745-20150124-total/sdssj155745-20150124-total_spec_model.dat'},
                         'SDSS-J163800':{'RA':'16:38:00.36', 'Dec':'+00:47:17.81', 'u_synth':[18.410, 0.005] , 'g_synth':[18.815, 0.006],'r_synth': [19.274, 0.004],'i_synth':[19.616, 0.004],'total_spectrum_file':total_flux_file_common_dir + 'sdssj163800-20150518-total/sdssj163800-20150518-total_spec_model.dat'},
                         #'SDSS-J172135':{'RA':'17:21:35.98', 'Dec':'+29:40:16.00', 'u_synth':[], 'g_synth':[], 'r_synth': [], 'i_synth':[],'total_spectrum_file':},
                         'SDSS-J181424':{'RA':'18:14:24.13', 'Dec':'+78:54:02.90', 'u_synth':[16.212, 0.002] , 'g_synth':[16.524, 0.002] ,'r_synth':[16.997, 0.002],'i_synth':[17.357, 0.002],'total_spectrum_file':total_flux_file_common_dir + 'sdssj181424-20150427-total/sdssj181424-20150427-total_spec_model.dat'},
                         'SDSS-J203722':{'RA':'20:37:22.17', 'Dec':'−05:13:03.03', 'u_synth':[18.412, 0.009] , 'g_synth':[18.643, 0.006],'r_synth': [19.055, 0.003],'i_synth':[19.382, 0.004],'total_spectrum_file':total_flux_file_common_dir + 'sdssj210150-20150518-total/sdssj210150-20150518-total_spec_model.dat'},
                         #'SDSS-J210150':{'RA':'21:01:50.66', 'Dec':'−05:45:50.97', 'u_synth':[], 'g_synth':[], 'r_synth': [], 'i_synth':[],'total_spectrum_file':},
                         'SDSS-J232941':{'RA':'23:29:41.33', 'Dec':'+00:11:07.80', 'u_synth':[18.156, 0.003] , 'g_synth':[18.147, 0.007],'r_synth':[18.468, 0.005],'i_synth':[18.752, 0.004],'total_spectrum_file':total_flux_file_common_dir + 'sdssj232941-20150917-20151106-total/sdssj232941-20150917-20151106-total_spec_model.dat'},
                         'SDSS-J235144':{'RA':'23:51:44.29', 'Dec':'+37:55:42.60', 'u_synth':[17.747, 0.004] , 'g_synth':[18.061, 0.006] ,'r_synth':[18.448, 0.004],'i_synth':[18.754, 0.003],'total_spectrum_file':total_flux_file_common_dir + 'sdssj235144-20150915-total/sdssj235144-20150915-total_spec_model.dat'},
                        }
    if target_stars in ['all','ALL']:
        return all_cal_star_info
    else:
        target_cal_star_info = {} 
        for target_star in target_stars:
            target_cal_star_info[target_star] = all_cal_star_info[target_star]

        return target_cal_star_info

#Reads in the spectrum of the Narayan stars 
def getStarSpectra(target_stars = 'all',
                   rows_to_ignore = 1, wavelength_col = 0, flux_col = 4, flux_err_col = 5, sep_char = ' ', spectrum_file_keyword = 'total_spectrum_file'):
    
    target_star_info = getStarFixedParams(target_stars = target_stars)
    print ('target_star_info = ' + str(target_star_info)) 
    target_stars = list(target_star_info.keys())
    spectra_of_stars = {target_star: [ recursiveReplaceMultipleSpaces(new_row).strip().split(' ') for new_row in readInFileLineByLine(target_star_info[target_star][spectrum_file_keyword])[rows_to_ignore:] ] for target_star in target_stars}
    #wavelengths should be expressed in nm, so we divide stored wavelength (expressed in Angstroms) by 10 
    spectra = {target_star:[[float(spectrum_row[wavelength_col])/10.0, float(spectrum_row[flux_col]), float(spectrum_row[flux_err_col])] for spectrum_row in spectra_of_stars[target_star]] for target_star in target_stars}
    spectra = {target_star:[[spectral_point[0] for spectral_point in spectra[target_star]], [spectral_point[1] for spectral_point in spectra[target_star]], [spectral_point[2] for spectral_point in spectra[target_star]]] for target_star in target_stars}

    return spectra

#bandpass curves should be in count transmission efficiency as a function of wavelength
# wavelengths should be in nm
# frequencies should be in Ghz
# lam = c / nu = => lam_nm = (c_km_s * km / s) / (nm nu_Ghz 10.0 ** 9.0 1/s) = c_km_s / nu_Ghz * (10.0 ** 12 nm / s) / (1.0 ** 9.0 nm / s) = 10.0 ** 3.0 c_km_s / nu_Ghz 
def getStarABMagnitudesGivenFilters(bandpass_curves, target_stars = 'all'):
    n_bandpass_curves = len(bandpass_curves)
    
    astro_arch = AstronomicalParameterArchive()
    c_km_s = astro_arch.getc() 
    target_spectra = getStarSpectra(target_stars = target_stars)
    target_stars = list(target_spectra.keys())
    bandpass_wavelength_bounds = [[min(bandpass_curve[0]), max(bandpass_curve[0])] for bandpass_curve in bandpass_curves]
    bandpass_interps = [interp1d(bandpass_curve[0], bandpass_curve[1]) for bandpass_curve in bandpass_curves]
    wavelength_from_freq = lambda freq: 1000.0 * c_km_s/ freq
    freq_from_wavelength = lambda wavelength: (1000.0 * c_km_s) / wavelength
    bandpass_freq_bounds = [[freq_from_wavelength(bound[1]), freq_from_wavelength(bound[0])] for bound in bandpass_wavelength_bounds]
    print ('bandpass_freq_bounds = ' + str(bandpass_freq_bounds)) 
    
    mags = {target_star: [] for target_star in target_stars}
    #numerators = {target_star:integrate.quad(lambda wavelength: , ) for target_star in target_stars}
    #denominators = {target_star:integrate.quad()  for target_star in target_stars} 
                
    for target_star in target_stars:
        target_spectrum = target_spectra[target_star]
        spectrum_wavelength_bounds = [min(target_spectrum[0]), max(target_spectrum[0])]
        print ('spectrum_wavelength_bounds = ' + str(spectrum_wavelength_bounds))  
        spectrum_freq_bounds = [freq_from_wavelength(spectrum_wavelength_bounds[1]), freq_from_wavelength(spectrum_wavelength_bounds[0])]
        print ('spectrum_freq_bounds = ' + str(spectrum_freq_bounds))
        numerators = [np.nan for i in range(n_bandpass_curves)]
        denominators = [np.nan for i in range(n_bandpass_curves)]
        for i in range(n_bandpass_curves):
            bandpass_freq_bound = bandpass_freq_bounds[i]
            bandpass_interp = bandpass_interps[i]
            if bandpass_freq_bound[1] < spectrum_freq_bounds[0] or bandpass_freq_bound[0] > spectrum_freq_bounds[1]:
                print ('for ' + str(i) + 'th filter, star spectrum unknown over full filter bandpass')
            else:
                spectrum_interp = interp1d(target_spectrum[0], target_spectrum[1])
                numerator = integrate.quad(lambda freq: 1.0 / freq * bandpass_interp(wavelength_from_freq(freq)) * spectrum_interp(wavelength_from_freq(freq)),
                                         max(bandpass_freq_bound[0], spectrum_freq_bounds[0]), min(bandpass_freq_bound[1], spectrum_freq_bounds[1]) )[0]
                numerators[i] = numerator 
                denominator = integrate.quad(lambda freq: 1.0 / freq * bandpass_interp(wavelength_from_freq(freq)), bandpass_freq_bound[0], bandpass_freq_bound[1])[0]
                denominators[i] = denominator 
                       
        mean_spectral_flux_densities = [numerators[i] / denominators[i] for i in range(len(numerators))] 
        mags[target_star] = [-2.5 * np.log10(mean_spectral_flux_density) - 48.60 for mean_spectral_flux_density in mean_spectral_flux_densities] 

    return mags 
    
    

def getStarColors(target_stars = 'all', color_pairs = [['g_synth','r_synth']]):
    target_star_info = getStarFixedParams(target_stars = target_stars)
    target_stars = list(target_star_info.keys()) 
    star_colors = {target_star:{(color_pair[0] + '-' + color_pair[1]):[target_star_info[target_star][color_pair[0]][0] - target_star_info[target_star][color_pair[1]][0],
                                                                       np.sqrt(target_star_info[target_star][color_pair[0]][1] ** 2.0 + target_star_info[target_star][color_pair[1]][1] ** 2.0)] 
                                 for color_pair in color_pairs}
                    for target_star in target_stars}

    return star_colors 
