import numpy as np
import math 
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interp1d
import scipy.integrate as integrate 

class ThreeDimensionalFunction:

    def integrate_los(self, xs = None, zs = None, integration_bin_centers = None, los_bin_step = 1.0): #,integration_rect_centers):
        if xs is None: xs = self.xs
        if zs is None: zs = self.zs
        if integration_bin_centers is None: integration_bin_centers = self.bin_centers
        n_centers = len(integration_bin_centers) 
        integration_bin_borders = ( [integration_bin_centers[0] - (integration_bin_centers[1] - integration_bin_centers[0]) / 2.0]
                                        + [(integration_bin_centers[i+1] + integration_bin_centers[i]) / 2.0 for i in range(n_centers - 1)]
                                        + [integration_bin_centers[n_centers-1] + (integration_bin_centers[n_centers-1] - integration_bin_centers[n_centers - 2]) /2.0 ] )
        integration_bin_borders.sort()
        integrated_array=np.zeros((len(xs),len(zs)))

        integrated_array = sum([ self.getRotatedCrossSection(self.funct, xs, integration_bin_centers[i], zs)
                                  * (integration_bin_borders[i + 1] - integration_bin_borders[i]) for i in range(len(integration_bin_centers))])
        onSkyInterpolator = RegularGridInterpolator((xs, zs), integrated_array, method = 'linear')

        return onSkyInterpolator 


    def getRotatedCrossSection(self, funct, xs, y_val, zs):
        xmesh, zmesh = np.meshgrid(xs, zs)
        return np.transpose(funct(xmesh, y_val, zmesh))

    def getRadialAverageMassFunction(self, Rs, phis, xs = None, zs = None, funct = None, integration_bin_centers = None):
        if xs is None: xs = self.xs
        if zs is None: zs = self.zs
        if funct is None: funct = self.funct 
        if integration_bin_centers is None: integration_bin_centers = self.bin_centers
        x_of_phi = lambda R, phi: R * math.cos(phi)
        z_of_phi = lambda R, phi: R * math.sin(phi)
        surface_mass_interpolator = self.integrate_los(xs = xs, zs = zs, integration_bin_centers = integration_bin_centers)
        R_seps = [Rs[1] - Rs[0]] + [(Rs[i+1] - Rs[i-1]) / 2.0 for i in range(1, len(Rs)-1) ] + [ Rs[-1] - Rs[-2] ]
        phi_seps = [phis[1] - phis[0]] + [(phis[i+1] - phis[i-1]) / 2.0 for i in range(1, len(phis)-1) ] + [ phis[-1] - phis[-2] ]
        ring_integrated_masses = [R_seps[i] * Rs[i] * sum([surface_mass_interpolator( (x_of_phi(Rs[i], phis[j]), z_of_phi(Rs[i], phis[j])) )
                                       * phi_seps[j] for j in range(len(phis)) ]) for i in range(len(Rs))]
        ring_integrated_average_masses = [sum([surface_mass_interpolator( (x_of_phi(Rs[i], phis[j]), z_of_phi(Rs[i], phis[j])) )
                                          * phi_seps[j] for j in range(len(phis)) ]) for i in range(len(Rs))]
        
        radial_interpolator = interp1d(Rs, ring_integrated_average_masses, kind = 'cubic')
        return radial_interpolator 

    def __init__(self, funct, xs, bin_centers, zs):
        self.funct = funct
        self.xs = xs
        self.bin_centers = bin_centers 
        self.zs = zs
        self.surfaceMassInterpolator = self.integrate_los(xs = self.xs, zs = self.zs, integration_bin_centers = self.bin_centers)
