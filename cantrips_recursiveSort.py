#small functions that do basic computing things. 

import math
import numpy as np 

#Sorts list2 by list1, and will sort list2 in a consistent way based on duplicats in list1
# I just use quicksort (directly implemented), having the second array tailing the first 
def safeSortOneListByAnother(list1, list2):
    if len(list1) != len(list2):
        print 'Arrays to simultaneously sort must be of same size! Returning list2 unsorted. '
        return list2
    elif len(list1) <= 1:
        return list2
    else:
        pivot_index = len(list1)-1
        pivot_value = list1[pivot_index]
        pivotless_list1 = list1[0:pivot_index] + list1[pivot_index+1:]
        pivotless_list2 = list2[0:pivot_index] + list2[pivot_index+1:]
        left_list1 = []
        left_list2 = []
        right_list1 = []
        right_list2 = []
        for i in range(len(pivotless_list1)):
            list1_elem = pivotless_list1[i]
            list2_elem = pivotless_list2[i]
            if list1_elem < pivot_value:
                left_list1 = left_list1 + [list1_elem]
                left_list2 = left_list2 + [list2_elem]
            else:
                right_list1 = right_list1 + [list1_elem]
                right_list2 = right_list2 + [list2_elem]
        return safeSortOneListByAnother(left_list1, left_list2) + [list2[pivot_index]] + safeSortOneListByAnother(right_list1, right_list2)

def dotProduct(vec1, vec2):
    return sum([vec1[i] * vec2[i] for i in range(len(vec1))])

def normalizeVector(vec):
    orig_mag = math.sqrt(sum([elem ** 2.0 for elem in vec]))
    return [elem / orig_mag for elem in vec]

def removeListElement(orig_list, index):
     return orig_list[0:index] + orig_list[index + 1:len(orig_list)] if index < len(orig_list) - 1 else orig_list[0:len(orig_list)-1]

def insertListElement(orig_list, val, index):
    if isinstance(orig_list, np.ndarray):
        orig_list = orig_list.tolist() 
    return orig_list[:index] + [val] + orig_list[index:]

#Inserts vals at the specified indeces,
# where the indeces are the indeces in the FINAL LIST
def insertListElements(orig_list, index_to_val_dict):
    if not isinstance(index_to_val_dict, dict):
        print 'Second argument of insertListElements must be a dictionary, which it was not. '
        print 'Returning empty array. '
        return []
    if isinstance(orig_list, np.ndarray):
        orig_list = orig_list.tolist() 
    new_list = orig_list
    for index in sorted(index_to_val_dict.keys()):
        new_list = insertListElement(new_list, index_to_val_dict[index], index)

    return new_list

def combineTwoDicts(dict1, dict2):
    return_dict = dict1.copy()
    return_dict.update(dict2)
    return return_dict 

#Sort n lists simultaneously base on first list
def simulSortLists(*lists):
    s_list_of_lists = [list(tup) for tup in zip(*sorted(zip(*lists)))]
    return s_list_of_lists 

#get element of a list that is closest to some value:
def getClosestElement(list_of_elems, elem):
    closest_elem = list_of_elems[ (getIndexOfClosestElement(list_of_elems, elem)) ]
    return closest_elem

#get index of list element that is closest to some value
def getIndexOfClosestElement(list_of_elems, elem):
    mag_diff_array = (np.array(list_of_elems) - elem ) ** 2.0
    min_index = np.argmin(mag_diff_array)
    return min_index

#index multidimensional numpy array based on series of indeces (must be in correct order)
def indexMultiDimArrayWithUnknownNumberOfIndeces(full_array, indeces):
    array_to_index = np.copy(full_array)
    for index in indeces:
        array_to_index = array_to_index[index]

    return array_to_index

#Perform operation that is normally performed on single number
# to every element in numpy array.  Then return array of same
# size.

def operateOnNumpyArray(numpy_array, operation):
    flattened_array = numpy_array.flatten()
    operated_flat_array = np.array([operation(elem) for elem in  flattened_array.tolist()])
    operated_array = operated_flat_array.reshape(np.shape(numpy_array))
    return operated_array 

def getCPBValsFromArrayAtSpecifiedPoints(measured_values, values_at_which_to_sample_CPB):
    values_at_which_to_sample_CPB = np.array(values_at_which_to_sample_CPB)
    measured_values = np.array(measured_values)
    n_measured_vals = len(measured_values)
    sample_values_mesh, measured_values_mesh = np.meshgrid(values_at_which_to_sample_CPB, measured_values)

    measured_below_sample_point_mesh = sample_values_mesh > measured_values_mesh

    n_measured_vals_below_samples = np.sum(measured_below_sample_point_mesh, axis = 0)

    #n_below = sum(val_below_sampling_mesh, axis = 0)

    #n_below = np.zeros(np.shape(sampling_points))
    #for elem in list_to_determine_CPB:
    #    new_below = sampling_points > elem
    #    n_below = n_below + new_below
    return n_measured_vals_below_samples / float(n_measured_vals)

#Return the function that allows you to measured the
# Cumulative Probability Distribution (CPB) of the list passed into this function
# at some given values passed into the returned function.
#Note that list must bye SORTED 
def getCPBFunctFromArray(measured_values):
    return lambda values_at_which_to_sample_CPB: getCPBValsFromArrayAtSpecifiedPoints(measured_values, values_at_which_to_sample_CPB)

#For two coordinate meshes, return a grid of the same size that characterizes the areas defined by the meshes. 
def getAreaMesh(meshes):
    ext_meshes = []
    border_meshes = []
    for mesh in meshes:
        mesh_shape = np.shape(mesh) 
        ext_mesh = np.zeros(np.array(mesh_shape) + 2)
        ext_mesh[1:mesh_shape[0]+1, 1:mesh_shape[1]+1] = mesh
        ext_mesh[0,:] = ext_mesh[1,:] - (ext_mesh[2,:] - ext_mesh[1,:]) 
        ext_mesh[mesh_shape[0]+1,:] = ext_mesh[mesh_shape[0],:] + (ext_mesh[mesh_shape[0],:] - ext_mesh[mesh_shape[0]-1,:])
        ext_mesh[:,0] = ext_mesh[:,1] - (ext_mesh[:,2] - ext_mesh[:,1]) 
        ext_mesh[:,mesh_shape[1]+1] = ext_mesh[:,mesh_shape[1]] + (ext_mesh[:,mesh_shape[1]] - ext_mesh[:,mesh_shape[1]-1])
        ext_meshes = ext_meshes + [ext_mesh]
    edge_loc_mesh1 = (ext_meshes[1][1:,:] + ext_meshes[1][0:-1,:])[:,1:] / 2.0
    edge_loc_mesh2 = (ext_meshes[0][:,1:] + ext_meshes[0][:,0:-1])[1:,:] / 2.0
    edge_loc_meshes = [edge_loc_mesh1, edge_loc_mesh2]

    edge_size_mesh1 = (edge_loc_mesh1[1:,:] - edge_loc_mesh1[0:-1,:])[:,1:]
    edge_size_mesh2 = (edge_loc_mesh2[:,1:] - edge_loc_mesh2[:,0:-1])[1:,:]
    edge_size_meshes=[edge_size_mesh1, edge_size_mesh2]
    
    area_mesh = edge_size_mesh1 * edge_size_mesh2

    return area_mesh 
