#The astrometry.net software uses a healpix methodology to track the various index files.
# However, their convention is different from the default convention of the python healpy library.
# I have matched the astrometry.net healpix files (from this link: https://github.com/dstndstn/astrometry.net/blob/master/util/hp2.png)
#  to the healpix pixels of healpy
#  by running this code:
#  >>> import numpy as np
#  >>> import healpy as hp
#  >>> NSIDE = 2
#  >>> for i in range(48):
#  >>>    print ('pix ' + str(i) + ' Dec, RA:  ' + str(   np.array(hp.pix2ang(2, i)) * 180.0 / np.pi * (-1.0, 1.0) + (90.0, 0.0) ) )

#Downloaded index files should be written to /usr/local/Cellar/astrometry-net/0.78_6/data/

import numpy as np
import healpy as hp
import os
import AstronomicalParameterArchive
import subprocess

class HealpixLookupTable:

    # Get astrometry healpix from RA and Dec in degrees
    def getAstrometryHealpix(self, ra_deg, dec_deg):
        astro_arch = AstronomicalParameterArchive.AstronomicalParameterArchive()
        deg_to_rad = astro_arch.getDegToRad()
        ra_rad, dec_rad = [ra_deg * deg_to_rad, dec_deg * deg_to_rad]
        theta, phi = [-dec_rad + np.pi / 2.0, ra_rad]
        print ('[theta, phi] = ' + str([theta, phi]))
        healpy_pix = hp.ang2pix(self.nsides, theta, phi )
        astrometry_pix = self.healpixToAstrometryTable[healpy_pix]

        index_file = 'index-5000-' + str(astrometry_pix) + '.fits'
        index_website = 'http://data.astrometry.net/5000/'
        index_file_download_commands = 'wget'
        download_command_bash_dir_option = '-P'
        index_dir = '/usr/local/Cellar/astrometry-net/0.78_6/data/'
        print ('index_dir + index_file = ' + index_dir + index_file)
        if not(os.path.exists(index_dir + index_file)):
            copy_output = subprocess.run([index_file_download_commands,  index_website + index_file, download_command_bash_dir_option, index_dir])
            print ('copy_output = ' + str(copy_output))
        #print ('copy_output = ' + str(copy_output))

        return index_dir


    def __init__(self):
        self.healpixToAstrometryTable = {0:3, 1:7, 2:11, 3:15,
                                    4:1, 5:2, 6:5, 7:6, 8:7, 9:10, 10:13, 11:14,
                                    12:19, 13:0, 14:23, 15:4, 16:27, 17:8, 18:31, 19:12,
                                    20:18, 21:21, 22:22, 23:25, 24:26, 25:29, 26:30, 27:13,
                                    28:16, 29:35, 30:20, 31:39, 32:24, 33:43, 34:28, 35:47,
                                    36:33, 37:34, 38:37, 39:38, 40:41, 41:42, 42:45, 43:46,
                                    44:32, 45:36, 46:40, 47:44  }
        self.nsides = 2
