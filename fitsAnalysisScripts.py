import cantrips as c
import numpy as np
import photutils
from photutils import aperture_photometry
from photutils import CircularAperture, CircularAnnulus
from astropy import wcs
from astropy.io import fits
import sys

#These processes require you pass numpy arrays, presumably read in from fits files
#You can use cantrips to load the files as:
# data, header = c.readInDataFromFitsFile(file_name, file_dir)
# Note: center_pix should be given in the form [y_pix, x_pix] because of how numpy arrays are indexed

def measureCircularAperture(data, center_x, center_y, rad, method='subpixel', subpixels = 5):
    aperture = CircularAperture([center_x,center_y], r=rad)
    phot_table = aperture_photometry(data, aperture, method = method, subpixels = subpixels)
    sum = np.array(phot_table['aperture_sum'].data) [0]
    area = np.pi * rad ** 2.0
    surface_brightness = sum / area

    return [sum , surface_brightness]

#Measure brightness by measuring brightness in an inner circle
# then measure the background via an aperture of larger radius.
def measureAnnulusAperture(data, center, rad_center, rad_background_inner, rad_background_outer,
                           method='subpixel', subpixels = 5):
    center_x, center_y = center
    inner_aperture = CircularAperture(center, r=rad_center)
    outer_aperture = CircularAnnulus(center, r_in=rad_background_inner, r_out=rad_background_outer)
    apers = [inner_aperture, outer_aperture]
    phot_table = aperture_photometry(data, apers)
    inner_sum = np.array(phot_table['aperture_sum_0'].data) [0]
    inner_area = np.pi * rad_center ** 2.0
    outer_sum = np.array(phot_table['aperture_sum_1'].data) [0]
    outer_area = np.pi * (rad_background_outer ** 2.0 - rad_background_inner ** 2.0)
    background = outer_sum / outer_area
    bg_sub_sum = inner_sum - inner_area * background
    return [inner_sum, inner_area, outer_sum, outer_area]

# Load the WCS information from a fits header, and use it
# to convert pixel coordinates to world coordinates.

def loadWCSFromHeader(header):
    # Load the FITS hdulist using astropy.io.fits
    w = wcs.WCS(header)
    return w
