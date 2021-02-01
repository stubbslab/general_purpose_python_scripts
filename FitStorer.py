import numpy as np
import math
from scipy.optimize import curve_fit
from numpy import polyfit 


class FitStorer:

    def getDimension(self):
        return self.fit_dimension

    def generateFit(self, xs, ys, y_errs=None, start_params = None, fixed_vars = None):
        if self.fit_dimension < 1:
            print ('No fit loaded.  Will not try to generate a fit. ')
            return 0
        if self.fit_funct is 'polyfit':
            if y_errs is None:
                weights = None
            else:
                weights = [1 / y_err for y_err in y_errs]
            
            (self.fit_params, self.covariance) = polyfit (xs, ys, self.fit_order, w = weights, cov = True)

            print ('self.covariance = ')
            print (self.covariance)
            self.fit_string = 'poly_' + '+'.join([str(np.around(self.fit_params[i],5)) + '(' + str(np.around(self.covariance[i][i],5)) + ')'+ 'x^' + str(len(self.fit_params) - 1 - i) for i in range(len(self.fit_params)) ])
        else:
            #print self.fit_funct
            #print xs
            #print ys
            #print y_errs
            if self.bounds is None:
                self.fit_params, self.covariance = curve_fit (self.fit_funct, xs, ys, sigma = y_errs, p0 = self.p0, maxfev = 10000)
            else:
                self.fit_params, self.covariance = curve_fit (self.fit_funct, xs, ys, sigma = y_errs, p0 = self.p0, max_nfev = 10000, bounds = self.bounds )
            print (self.fit_funct )
            if type(self.in_funct) is str:
                print ('a')
                self.fit_string = self.in_funct + '_' + '_'.join([str(np.around(elem,5)) for elem in self.fit_params])
            else:
                print ('b')
                self.in_string = 'userFunction_' + '_'.join([str(elem) for elem in self.fit_params]) 

    def getFitValues(self, xs):
        if self.fit_params is '':
            print ('Fit params not yet assigned.  Returning 0.')
            return 0
        if self.fit_funct is 'polyfit':
            return [sum([ x ** (len(self.fit_params) - n - 1) * self.fit_params[n] for n in range(len(self.fit_params)) ]) for x in xs]
        else:
            return [self.fit_funct(x, *self.fit_params) for x in xs]

    #All input information of fit is contained in fit_information.  Should be:
    # {'funct':['poly','sine','shift_sine'], 'order': [order of polynomial fit], 'guess': [best guest parameters], 'dimension': [number of free parameters]}
    def __init__(self, fit_information):
        self.maxfev = 10000
        self.fit_information = fit_information 
        self.in_funct = self.fit_information['funct']
        
        #Default to no fit, which should stop any attempts to make a fit when no fit is here

        #Now try to see what sort of fit we have

        #Predefined (by string) fits
        if type(self.in_funct) == str:
            self.in_funct = self.in_funct.lower() 
            if self.in_funct == 'none':
                self.fit_funct = 'none'
                self.fit_dimension = 0
                self.fit_string = 'none'
            elif self.in_funct == 'poly' or self.in_funct == 'polyfit' or self.in_funct == 'polynomial' or self.in_funct == 'poly1d' or self.in_funct == 'polyfit1d' or self.in_funct == 'polynomial1d':
                self.fit_funct = 'polyfit'
                default_order = 0
                self.fit_order = 0
                if 'order' in self.fit_information.keys():
                    self.fit_order = self.fit_information['order']
                self.fit_dimension = 1
                self.fit_string = 'polynomial'
            elif self.in_funct == 'sine' or self.in_funct == 'sinusoidal' or self.in_funct == 'sine1d' or self.in_funct == 'sinusoidal1d':  
                self.fit_funct = lambda x, A, omega, phi: A * np.sin( omega * x + phi )
                self.fit_dimension = 1
                self.fit_string = 'sine'
            elif self.in_funct == 'shift_sine' or self.in_funct == 'shift_sinusoidal' or self.in_funct == 'shift_sine1d' or self.in_funct == 'shift_sinusoidal1d':
                self.fit_funct =  lambda x, shift, A, omega, phi: A * np.sin( omega * x + phi ) + shift
                self.fit_dimension = 1
                self.fit_string = 'shift_sine'
            elif self.in_funct == 'gauss' or self.in_funct == 'guassian' or self.in_funct == 'gauss1d' or self.in_funct == 'gaussian1d':
                self.fit_funct =  lambda x, A, mu, width: A * np.exp(-(x - mu) ** 2.0 / (2 * width ** 2.0))
                self.fit_dimension = 1
                self.fit_string = 'shift_gauss'
            elif self.in_funct == 'shift_gauss' or self.in_funct == 'shift_guassian' or self.in_funct == 'shift_gauss1d' or self.in_funct == 'shift_gaussian1d':
                self.fit_funct =  lambda x, shift, A, mu, width: A * np.exp(-(x - mu) ** 2.0 / (2 * width ** 2.0)) + shift
                self.fit_dimension = 1
                self.fit_string = 'shift_gauss'
            else:
                print ("Fit type '" + self.in_funct + "' not recognized. ")
        #No fit loaded
        elif self.in_funct is None:
            self.fit_funct = 'none'
            self.fit_dimension = 0
        #Also allow user to define their own function 
        else:
            self.fit_funct = self.in_funct
            self.fit_dimension = self.fit_information['dimension']

        if 'guess' in self.fit_information.keys():
            self.p0 = fit_information['guess']
        else:
            self.p0 = None

        if 'bounds' in self.fit_information.keys():
            self.bounds = fit_information['bounds']
        else:
            self.bounds = None
            
            
