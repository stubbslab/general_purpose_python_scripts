import cantrips as cant
import matplotlib.pyplot as plt
import numpy as np

class rmapper:

    def showDataOnSky(self, z_range = [0.0, np.inf], color = 'k', marker = '.'):
        RAs_to_plot = [self.rm_RAs[i] for i in range(len(self.rm_RAs)) if (self.rm_zs[i] >= z_range[0] and self.rm_zs[i] <= z_range[1]) ]
        Decs_to_plot = [self.rm_Decs[i] for i in range(len(self.rm_Decs)) if (self.rm_zs[i] >= z_range[0] and self.rm_zs[i] <= z_range[1]) ]

        cant.plotStarsOnSky(RAs_to_plot, Decs_to_plot, color = color, marker = marker)

        return 1


    def __init__(self, data_file = 'redmapper_dr8_public_v6.3_catalog.fits', data_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/ClusterCatalogues/',
                 use_spec_zs = 0):

        data_table, col_names, header = cant.readInDataFromFitsFile(data_file, data_dir, data_type = 'table')
        self.spec_z_keyword = 'z_spec'
        self.phot_z_keyword = 'z_LAMBDA'
        self.RA_keyword = 'RA'
        self.Dec_keyword = 'Dec'

        self.rm_RAs = data_table[self.RA_keyword]
        self.rm_Decs = data_table[self.Dec_keyword]
        rm_phot_zs = data_table[self.phot_z_keyword]
        rm_spec_zs = data_table[self.spec_z_keyword]

        if use_spec_zs:
            self.rm_zs = rm_spec_zs
        else:
            self.rm_zs = rm_phot_zs
