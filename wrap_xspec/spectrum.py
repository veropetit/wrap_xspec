# Module with helper functions for xspec spectra
import os
import numpy as np
import xspec

class SpectrumLoader:
    '''This is a helper class to store the file information for a single spectrum'''
    def __init__(self, path, pi, back, rmf, arf) -> None:
        '''
        Class initialization for the spectrum class. 
        
        :param path: the path where the spectrum files are (absolute or relative)
        :param pi: the filename of the spectrum
        :param back: the filename of the backgroud spectrum
        :param rmf: the filename of the rmf file
        :param arf: the filename of the arf file
        '''
        self.path = path
        self.pi = pi
        self.back = back
        self.rmf = rmf
        self.arf = arf

    def xspec_load(self):
        '''
        A helper function to load a spectrum and responded in xspec. 
        It will change the directory to the self.path (xspec cannot handle the paths...)
        and return to the current directory afterwards, using the `os` package. 

        :param self: a spectrum object
        :rtype s: an Xspec spectrum object
        '''
        # changing directory, because pyXspec doesn't seem to be able to handle paths...
        print('Current path: {}'.format(os.getcwd()))
        path = os.getcwd()
        os.chdir(self.path)
        print('Moving to {}'.format(os.getcwd()))
        s = xspec.Spectrum(self.pi, 
                    backFile=self.back, 
                    respFile=self.rmf, 
                    arfFile=self.arf)
        # going back to the main directory
        os.chdir(path)
        print('Moving back {}'.format(os.getcwd()))
        return(s)
    



