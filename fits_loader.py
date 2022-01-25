from astropy.io import fits
import os


class FitsFileLoader():
    '''Class for handling fits files'''

    def __init__(self, spectra_directory=None):
        self.current_directory = os.getcwd()

        if spectra_directory is not None:
            self.spectra_directory = spectra_directory
        else:
            self.spectra_directory = os.path.join(self.current_directory, 'spectra/')
            if not os.path.exists(self.spectra_directory): self.spectra_directory = self.current_directory

        self.filenames = self.list_fits_filenames(self.spectra_directory)
        self.filepaths = self.list_fits_filepaths(self.spectra_directory)

    @staticmethod
    def list_fits_filepaths(spectra_directory=None):
        '''Lists all absolute paths of .fits files' in spectra_directory'''
        return [
            os.path.join(spectra_directory, filename)
            for filename in FitsFileLoader.list_fits_filenames(spectra_directory)
        ]
        
    @staticmethod
    def list_fits_filenames(spectra_directory=None):
        '''Lists all .fits filenames in spectra_directory'''
        return [
            filename
            for filename in os.listdir(spectra_directory)
            if filename.lower().endswith(".fits")
        ]

    @staticmethod
    def data_from_single_fits(filepath, hdu_index):
        '''Returns data from hdu[1] from fits file'''
        hdul = fits.open(filepath)
        return hdul[hdu_index].data

    
