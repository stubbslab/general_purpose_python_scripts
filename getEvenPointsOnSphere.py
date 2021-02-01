#Distributes roughly even number of points over a sphere, according to the
# method described in this paper:
# https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
#Note that it opts to distribute not exactly the requested number of
# points in favor of ensuring that the points are evenly distributed. 

import math
import numpy as np
from cantrips import normalizeVector
from cantrips import dotProduct

def getEvenPointsOnSphere(n_points, r = 1.0, phi_shift = 0.0, theta_shift = 0.0, return_format = 'vec', ref_z_axis = [0.0,0.0,1.0], ref_y_axis = [0.0,1.0,0.0]):
    r = float(r) 
    n_placed = 0
    point_area = 4 * math.pi / n_points
    point_side = math.sqrt(point_area)
    n_theta = int(round(math.pi / point_side))
    theta_sep = math.pi / float(n_theta)
    phi_sep = math.pi / float(n_theta)

    phi_theta_pairs = []

    for i in range(n_theta):
        theta_cent = (math.pi * (i + 0.5) / n_theta + theta_shift) % (math.pi) 
        n_phi = int(round(2 * math.pi * math.sin(theta_cent) / float(phi_sep)))
        for j in range(n_phi):
            phi_cent = (2 * math.pi * j / n_phi + phi_shift) % (2.0 * math.pi) 
            phi_theta_pairs = phi_theta_pairs + [[phi_cent, theta_cent]]
            n_placed = n_placed + 1
    print 'placed ' + str(n_placed) + ' points on sphere. '

    vectors =  [[r * math.cos(phi_theta[0]) * math.sin(phi_theta[1]),
                 r * math.sin(phi_theta[0]) * math.sin(phi_theta[1]),
                 r * math.cos(phi_theta[1]) ]
                 for phi_theta in phi_theta_pairs]
    #Now we want to reexpress these randomly chosen point on a sphere in terms of some use-specified base
    ref_y_axis = normalizeVector(ref_y_axis)
    ref_z_axis = normalizeVector(ref_z_axis)
    ref_x_axis = np.cross(ref_y_axis, ref_z_axis)

    new_coord_vectors = [[dotProduct(vector, ref_x_axis), dotProduct(vector, ref_y_axis), dotProduct(vector, ref_z_axis)] for vector in vectors]
    new_coord_phi_theta_pairs = [[math.acos(vector[0] / math.sqrt(vector[0] ** 2.0 + vector[1] ** 2.0 )) if vector[1] >= 0.0 else 2 * math.pi - math.acos(vector[0] / math.sqrt(vector[0] ** 2.0 + vector[1] ** 2.0 )), math.acos(vector[2] / r)] for vector in new_coord_vectors]
    

    if return_format == 'angles':
        return new_coord_phi_theta_pairs
    elif return_format == 'vec':
        return new_coord_vectors
    else:
        return new_coord_vectors
            
            
                      
