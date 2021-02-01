import numpy as np
import scipy.interpolate as interpolate
from cantrips import readInFileLineByLine 

class FilterCurve:

    def getThroughput(self, wavelengths):
        wavelengths = np.array(wavelengths) 
        if len(np.shape(wavelengths)) > 0:  
            throughputs = [self.throughput_interp(wavelength) if (wavelength > self.wavelength_range[0] and wavelength < self.wavelength_range[1]) else 0.0
                           for wavelength in wavelengths
                           ]
        else:
            throughputs = self.throughput_interp(wavelengths) if (wavelengths > self.wavelength_range[0] and wavelengths < self.wavelength_range[1]) else 0.0
        return throughputs 

    def __init__(self, filter_file, filter_file_dir, split_char = ' '):
        rows = readInFileLineByLine(filter_file, filter_file_dir)
        rows = [row.split(split_char) for row in rows]
        self.wavelengths = [float(row[0]) for row in rows] 
        self.wavelength_range = [min(self.wavelengths), max(self.wavelengths)]
        self.throughputs = [float(row[1]) for row in rows]
        self.throughput_interp = interpolate.interp1d(self.wavelengths, self.throughputs)
        
                       
        
