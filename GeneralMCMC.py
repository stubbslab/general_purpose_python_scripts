#A general MCMC that works on an input fit function

import numpy as np
import VisitedMCMCPoint as vmp
import random
import scipy.optimize as optimize
from binData import binData
import matplotlib.pyplot as plt
import scipy.special as special

class MCMCFitter:

    def resampleParams(self, current_params, sample_step_sizes ):

        new_params = current_params[:]
        for param_index in range(len(current_params)):
            current_param = current_params[param_index]
            step_size = sample_step_sizes[param_index]
            single_param_bounds = self.bounds[param_index]
            if self.sample_funct[param_index].lower() == 'gauss':
                new_param = np.random.normal(current_param, step_size)
                if new_param < single_param_bounds[0]:
                    new_param = min(single_param_bounds[0] + (single_param_bounds[0] - new_param), single_param_bounds[1])
                if new_param > single_param_bounds[1]:
                    new_param = max(single_param_bounds[1] - (new_param - single_param_bounds[1]), single_param_bounds[1])
                new_params[param_index] = new_param

        return new_params

    def doSingleMCMCChain(self, start_params):
        start_likelihood = self.MCMC_fit_funct(start_params)

        current_params = start_params[:]
        current_likelihood = start_likelihood
        current_n_steps = 0
        visited_points = [vmp.VisitedMCMCPoint(current_params, current_likelihood, 1)]

        while current_n_steps < self.n_steps:
            new_params = self.resampleParams(current_params, [step * (self.large_step_multiplier if current_n_steps % self.large_step_frequency == self.large_step_frequency - 1 else 1) for step in self.sample_step_sizes])
            #new_fit_resids, new_fit_chisqr, new_fit_weighted_mean = self.printingFitFunct(fit_funct, zs, resids, errs, np.array(new_params))
            new_likelihood = self.MCMC_fit_funct(new_params)
            if np.isnan(new_likelihood):
                print ('!!!WARNING!!! new_likelihood is nan for params ' + str(new_params))
            if self.likelihood_from_chisqr:
                likelihood_rat = np.exp(-(new_likelihood - current_likelihood))
            else:
                likelihood_rat = new_likelihood / current_likelihood
            #print ('[new_likelihood, current_likelihood, likelihood_rat] = ' + str([new_likelihood, current_likelihood, likelihood_rat] ))
            likelihood_rat_comp = random.random()
            if current_n_steps % self.printing_freq == self.printing_freq - 1:
                print ('Working on MCMC chain step ' + str(current_n_steps))
            if likelihood_rat > likelihood_rat_comp :
                #Move to new point
                if current_n_steps % self.printing_freq == self.printing_freq - 1:
                    print ('Moving from params: ' + str(current_params) + ' to params: ' + str(new_params) + ' with respective likelihoods ' + str([current_likelihood, new_likelihood]))
                visited_points = visited_points + [vmp.VisitedMCMCPoint(new_params, new_likelihood, 1)]
                current_params, current_likelihood = [new_params, new_likelihood]
            else:
                #Stay at current point
                if current_n_steps % self.printing_freq == self.printing_freq - 1:
                    print ('Not moving from params: ' + str(current_params) + ' to params: ' + str(new_params) + ' with respective likelihoods ' + str([current_likelihood, new_likelihood]))
                visited_points[-1].n_visits = visited_points[-1].n_visits + 1
            current_n_steps = current_n_steps + 1

        return visited_points

    def fitMCMCPosteriors(self, ):

        fit_res = []
        n_params = len(self.fullMCMCOutputs[0][0].parameters)

        for i in range(n_params):
            param_vals, n_visits = [[], []]
            for single_chain in self.fullMCMCOutputs :
                chain_len = len(single_chain)
                new_param_vals = [chain_elem.parameters[i] for chain_elem in single_chain]
                new_n_visits = [chain_elem.n_visits for chain_elem in single_chain]
                param_vals = param_vals + new_param_vals
                n_visits = n_visits + new_n_visits
            plt.scatter(param_vals, n_visits)
            plt.show()
            full_x_vals, full_y_vals = binData(param_vals, n_visits, n_bins = self.n_fit_bins)
            full_y_vals = [full_y_vals[0][i] * full_y_vals[2][i] for i in range(len(full_y_vals[0]))]
            x_vals = full_x_vals[self.fit_buffer:len(full_y_vals) - self.fit_buffer]
            y_vals = full_y_vals[self.fit_buffer:len(full_y_vals) - self.fit_buffer]
            print ('[x_vals, y_vals] = ' + str([x_vals, y_vals]))
            print ('[full_x_vals, full_y_vals] = ' + str([full_x_vals, full_y_vals]))
            print ("self.posterior_fit_funct_str.lower() = " + str(self.posterior_fit_funct_str.lower()))

            if self.posterior_fit_funct_str.lower() == 'single_hump_shift':
                print ('Here1')
                fitting_function = lambda xs, A, mu, sigma, power, shift: A * np.exp(-np.abs(xs - mu)**power /((sigma) ** power)) + shift
                #fitting_function = lambda xs, A, mu, sigma, shift: A * np.exp(-np.abs(xs - mu)**2.0 /((np.sqrt(2.0) * sigma) ** 2.0)) + shift
                single_bounds = ([0.0, x_vals[0],  0.0, 1.0, -np.inf], [np.inf, x_vals[-1], (x_vals[-1] - x_vals[0]) / 1.0, 10.0, np.inf])
                recompute_param_function = lambda params: [params[0], params[1], np.sqrt(params[2] ** 2.0 * special.gamma(3.0 / params[3]) / special.gamma(1.0 / params[3])), params[3], params[4]]
                A = np.max(y_vals) * 1.0
                mu = x_vals[np.argmax(y_vals)]
                #It is true that int_0^infinity (dx x e ^ (-x^2 / (2.0 * sigma ** 2.0))) = sigma ** 2.0.  So assuming the distribution is gaussian and I know the peak, this is a good first guess of the width
                sigma = np.sqrt(np.sum([(x_vals[1] - x_vals[0]) * abs(mu - x_vals[point_index]) * y_vals[point_index] / A for point_index in range(len(x_vals))]))
                power = 2.0
                shift = 0.0
                guess_params = [A, mu, sigma, power, shift]

            plt.plot(x_vals, y_vals, c = 'k')
            plt.plot(x_vals, fitting_function(x_vals, *guess_params), c = 'r')
            plt.show()
            single_param_fit_res = optimize.curve_fit(fitting_function, x_vals, y_vals, p0 = guess_params, bounds = single_bounds )
            fit_res = fit_res + [recompute_param_function(single_param_fit_res[0])[1:3] ]

        return fit_res

    def doFullMCMC(self):

        MCMC_chain_outputs = [[] for start_params in self.start_params_set]
        for chain_num in range(len(self.start_params_set)):
            start_params = self.start_params_set[chain_num]
            new_MCMC_chain = self.doSingleMCMCChain(start_params)
            MCMC_chain_outputs[chain_num] = new_MCMC_chain

        return MCMC_chain_outputs

    def __init__(self, fit_funct, start_params_set, sample_step_sizes, n_steps,
                      large_step_frequency = 100, large_step_multiplier = 10, bounds = [[-np.inf, np.inf]], sample_funct = ['gauss'], likelihood_from_chisqr = 0,
                      posterior_fit_funct_str = 'single_hump_shift', n_fits_bins = 101, fit_buffer = 2, printing_freq = 100):
        self.printing_freq = printing_freq
        self.MCMC_fit_funct = fit_funct
        self.start_params_set = start_params_set
        self.sample_step_sizes = sample_step_sizes
        print ('self.sample_step_sizes = ' + str(self.sample_step_sizes))
        print ('bounds = ' + str(bounds))
        self.n_steps = n_steps
        self.large_step_frequency = large_step_frequency
        self.large_step_multiplier = large_step_multiplier
        if len(bounds) == 1:
            bounds = [bounds[0] for param in start_params_set[0]]
        self.bounds = bounds
        if len(sample_funct) == 1:
            sample_funct = [sample_funct[0] for param in start_params_set[0]]
        self.sample_funct = sample_funct
        self.likelihood_from_chisqr = likelihood_from_chisqr
        self.n_fit_bins = n_fits_bins
        self.fit_buffer = fit_buffer
        self.posterior_fit_funct_str = posterior_fit_funct_str

        self.fullMCMCOutputs = self.doFullMCMC()

        self.best_fit_res = self.fitMCMCPosteriors()
