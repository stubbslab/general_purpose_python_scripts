import numpy as np
import scipy.interpolate as interpolate
import scipy.integrate as integrate 
import math
import time 
import scipy.optimize as optimization 

class ApproximateFourierTransform:

    def computeFrequencies(self, interval, start_n, end_n, extra_f_factor = 1.0 ):
        return [2.0 * math.pi * n / interval for n in np.arange(start_n, end_n, 1.0 / extra_f_factor)]

    def computeSinCoefs(self):
        sin_coefs = []
        sin_coef_errs = []
        for frequency in self.frequencies:
            funct_to_minimize = lambda x, A: A * np.sin(frequency * (x - self.min_x - self.x_interval / 2.0) )
            A_coef, A_coef_err = optimization.curve_fit(funct_to_minimize, self.xs, self.windowed_ys, 0.0, self.yerrs)
            
            A_coef = A_coef[0]
            A_coef_err = A_coef_err[0][0]
            sin_coefs = sin_coefs + [A_coef]
            sin_coef_errs = sin_coef_errs + [A_coef_err]
        return [sin_coefs, sin_coef_errs] 
        #return [integrate.quad(lambda x: self.centered_funct(x) * np.sin(x * freq), - self.x_interval / 2.0, self.x_interval / 2.0 )[0] for freq in self.frequencies]

    def computeCosCoefs(self):
        cos_coefs = []
        cos_coef_errs = []
        for frequency in self.frequencies:
            funct_to_minimize = lambda x, A: A * np.cos(frequency * (x - self.min_x - self.x_interval / 2.0))
            A_coef, A_coef_err = optimization.curve_fit(funct_to_minimize, self.xs, self.windowed_ys, 0.0, self.yerrs)
            A_coef = A_coef[0]
            A_coef_err = A_coef_err[0][0]
            cos_coefs = cos_coefs + [A_coef]
            cos_coef_errs = cos_coef_errs + [A_coef_err]
        return [cos_coefs, cos_coef_errs] 

    def computeA0(self):
        return integrate.quad(lambda x: self.centered_funct(x) * np.cos(x * 0.0), - self.x_interval / 2.0, self.x_interval / 2.0 )[0]

    def getFourierFunction(self):
        return self.fourier_funct

    def getCalculatedFourierFunction(self, x_vals):
        return [self.fourier_funct(x_elem) for x_elem in x_vals]
    
    def __init__(self, xs, ys, yerrs, frequencies = None, apply_normalization = 0, windowing_funct = 'rect', extra_f_factor = 1.0, remove_zero = 1):
        self.max_x = max(xs)
        self.min_x = min(xs)
        self.x_interval = (self.max_x - self.min_x)
        self.n_xs = len(xs)
        self.weights = [1.0 / (err ** 2.0) for err in yerrs]
        self.weighted_mean = sum([ ys[i] * self.weights[i] for i in range(len(self.weights)) ]) / ( sum(self.weights) )
        self.xs = xs
        self.remove_zero = remove_zero 
        if self.remove_zero:
            self.ys = [y - self.weighted_mean for y in ys]
        else:
            self.ys = ys 
        self.yerrs = yerrs

        if windowing_funct.lower() in ['rect', 'rectangle','flat']:
            hamming_vals = [1.0 for x in xs]
        elif windowing_funct.lower() in ['hamming',' hamm']:
            hamm_const_1 = 25.0/46.0
            hamm_const_2 = 21.0/46.0
            min_x = min(xs)
            max_x = max(xs)
            n_xs = len(xs)
            hamming_vals = [hamm_const_1 - hamm_const_2 * np.cos( (2.0 * np.pi * ((np.array(x) - min_x) / (max_x - min_x)) * n_xs )/float(n_xs- 1.0) ) for x in xs]
        else:
            hamming_vals = [1.0 for x in xs]
            print 'Windowing function ' + windowing_funct + ' not recognized.  Setting to flat.'

        if frequencies is None:
            self.frequencies = self.computeFrequencies(self.x_interval, 1, self.n_xs / 2, extra_f_factor = extra_f_factor)
        elif type(frequencies) is int:
            self.frequencies = self.computeFrequencies(self.x_interval, 1, frequencies, extra_f_factor = extra_f_factor)
        else:
            self.frequencies = frequencies 

        self.windowed_ys = [hamming_vals[i] * self.ys[i] for i in range(len(hamming_vals))]
        
        #self.centered_funct = lambda shift_x: self.funct(shift_x + self.x_interval / 2.0 + self.min_x)

        #print 'Computing sin coefs...'
        start = time.time()
        self.sin_coefs, self.sin_coef_errs = self.computeSinCoefs()
        end = time.time()
        #print 'Took ' + str(end - start) + 's for that.  '  
        #print 'Computing cos coefs...'
        start = time.time()
        self.cos_coefs, self.cos_coef_errs = self.computeCosCoefs()
        end = time.time()
        #print 'Took ' + str(end - start) + 's for that.  '
        
        #print 'Computing coef mags...'
        start = time.time() 
        self.coef_mags = np.sqrt(np.array(self.sin_coefs) ** 2.0 + np.array(self.cos_coefs) ** 2.0).tolist()
        self.coef_mag_errs = np.sqrt( (np.array(self.sin_coefs) * np.array(self.sin_coef_errs)) ** 2.0
                                          + (np.array(self.cos_coefs) * np.array(self.cos_coef_errs)) ** 2.0 ).tolist()
        self.coef_mag_errs = [ self.coef_mag_errs[i] / self.coef_mags[i] if self.coef_mags[i] > 0.0 else 0.01 for i in range(len(self.coef_mags)) ]
        end = time.time()
        #print 'Took ' + str(end - start) + 's for that.  '

        #Chris does a normalization thing here that I am not sure is necessary.  But I will put it here, just so that I don't need to do it later:
        coef_variance = sum([coef_mag ** 2.0 for coef_mag in self.coef_mags])
        mean_y = np.mean(self.ys) 
        mu_resid_variance = (sum([(y - mean_y) ** 2.0 for y in self.ys])) / (len(self.ys) - 1.0)

        self.variance_ratio = mu_resid_variance / coef_variance

        if apply_normalization:
            self.normalized_coef_mags = [coef * math.sqrt(self.variance_ratio) for coef in self.coef_mags]
        else:
            self.normalized_coef_mags = [coef * 1.0 for coef in self.coef_mags]

        #self.A0 = self.computeA0()

        self.A0 = 0
        
        #self.fourier_mags = [math.sqrt(self.sin_coefs[i] ** 2.0 + self.cos_coefs[i] ** 2.0) for i in range(len(self.sin_coefs)) ]

        fourier_funct_single_x = lambda x: self.A0 / 2.0 + sum ([self.sin_coefs[i]
                                                                    * np.sin(self.frequencies[i] * (x - self.min_x - self.x_interval / 2.0) )
                                                                    + self.cos_coefs[i]
                                                                    * np.cos(self.frequencies[i] * (x - self.min_x - self.x_interval / 2.0) )
                                                                    for i in range(len(self.frequencies))])
        self.fourier_funct = lambda x: fourier_funct_single_x(x) if type(x) in [float, int] else [fourier_funct_single_x (x_elem) for x_elem in x]
