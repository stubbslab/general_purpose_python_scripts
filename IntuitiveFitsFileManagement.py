import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import cantrips as can

class FitsObject:

    """
    A python object that reads in a fits data file, either an image or a table.  Written
        to make handling fits data more intuitive - the hdu/hdul stuff is all masked.
<<<<<<< HEAD
        The user can just query the object for the data or the header.  
=======
        The user can just query the object for the data or the header.
>>>>>>> 511b43e1351da7a86045ee0245caf8133e781f64
    """

    def showImage(self, n_mosaic_rows = 2, logscale = 1):
        if logscale:
            disp_data = np.log10(self.data)
        else:
            disp_data = self.data
        if self.fits_data_type != 'image':
            print ('This FitsObject contains fits data of type: "' + fits_data_type + '" - I cannot display that.')
        else:
            if self.n_mosaic_extensions <=1:
                plt.imshow(disp_data)
            else:
<<<<<<< HEAD
                f, axarr = plt.subplots(n_mosaic_rows, np.ceil(self.n_mosaic_extensions / 2))
                for i in range(self.n_mosaic_extensions):
                    axarr[i //2, i % 2].imshow(disp_data[i])
=======
                f, axarr = plt.subplots(n_mosaic_rows, int(np.ceil(self.n_mosaic_extensions / 2)))
                for i in range(self.n_mosaic_extensions):
                    axarr[i % 2, i // 2].imshow(disp_data[i])
>>>>>>> 511b43e1351da7a86045ee0245caf8133e781f64

        plt.show()

        return 1

<<<<<<< HEAD
    def saveFitsDataToFile(self, save_file_name, save_dir = None, overwrite = True):
        if save_dir == None:
            save_dir = self.load_dir
=======
    def saveFitsDataToFile(self, save_file_name, save_dir = None, overwrite = True, separate_mosaic = False, verbose = 1, headers = None):
        """
        Save the fits data contained in this class to a fits data file.  If this object contains fits image data,
            then this will make a fits image file.  Same for a fits table.
        If this FitsObject contains mosaic image data, the user has an option, specified by the separate_mosaic flag:
            they can either save the data as a single mosaic image (separate_mosaic = False) or as n_mosaic_extensions
            separate fits images, one for each mosaic extension (separate_mosaic = True).  Note that if the user
            wants to save a mosaic fits image to multiple files, they must give a list of file names in the
            save_file_name variable.
        """
        if save_dir == None:
            save_dir = self.target_dir
>>>>>>> 511b43e1351da7a86045ee0245caf8133e781f64
        if self.n_mosaic_extensions <= 1:
            if self.fits_data_type == 'image':
                col_names = []
            else:
                col_names = self.data.columns
<<<<<<< HEAD
            can.saveDataToFitsFile(self.data, save_file_name, save_dir, header = self.header, overwrite = overwrite, n_mosaic_extensions = self.n_mosaic_extensions, data_type = self.fits_data_type, col_names = col_names)
=======
            print ('self.n_mosaic_extensions = ' + str(self.n_mosaic_extensions))
            can.saveDataToFitsFile(self.data, save_file_name, save_dir, header = self.header, overwrite = overwrite, n_mosaic_extensions = self.n_mosaic_extensions, data_type = self.fits_data_type, col_names = col_names)
        else:
            if separate_mosaic:
                for i in range(self.n_mosaic_extensions):
                    if verbose: print ('Saving image from image extension ' + str(i+1) + ' of ' + str(self.n_mosaic_extensions))
                    ext_header = self.header[i+1]
                    ext_header['EXTNUMBR'] = i
                    can.saveDataToFitsFile(self.data[i], save_file_name[i], save_dir, header = ext_header, overwrite= overwrite, n_mosaic_extensions = 0)
            else:
                print('self.n_mosaic_extensions = ' + str(self.n_mosaic_extensions))
                can.saveDataToFitsFile(self.data, save_file_name, save_dir, header = self.header, overwrite= overwrite, n_mosaic_extensions = self.n_mosaic_extensions)
>>>>>>> 511b43e1351da7a86045ee0245caf8133e781f64

        return 1


    def __init__(self, fits_file, load_dir = '', fits_data_type = 'image', n_mosaic_extensions = 0):
        if fits_data_type in ['i','I','img','IMG','Img','image','IMAGE','Image']:
            fits_data_type = 'image'
        self.fits_file = fits_file
        self.target_dir = load_dir
        self.fits_data_type = fits_data_type
        self.n_mosaic_extensions = n_mosaic_extensions

        self.data, self.header = can.readInDataFromFitsFile(self.fits_file, self.target_dir, n_mosaic_image_extensions = self.n_mosaic_extensions, data_type = self.fits_data_type )

        if n_mosaic_extensions <=1:
            self.data = np.transpose(self.data)
        else:
<<<<<<< HEAD
            self.data = [np.transpose(self.data[i]) for i in range(len(n_mosaic_extensions))]
=======
            self.data = [np.transpose(self.data[i]) for i in range(self.n_mosaic_extensions)]
>>>>>>> 511b43e1351da7a86045ee0245caf8133e781f64
