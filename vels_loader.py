from astropy.io import fits
import os
import fits_loader

class VelsLoader():
    '''Class for handling fits files'''

    def __init__(self, vels_directory=None):
        self.current_directory = os.getcwd()

        if vels_directory is not None:
            self.vels_directory = vels_directory
        else:
            self.vels_directory = os.path.join(self.current_directory, 'vels/')
            if not os.path.exists(self.vels_directory): self.vels_directory = self.current_directory

        self.filenames = self.list_fits_filenames(self.vels_directory)
        self.filepaths = self.list_fits_filepaths(self.vels_directory)

    @staticmethod
    def list_fits_filepaths(directory=None):
        '''Lists all absolute paths of .fits files' in vels_directory'''
        return [
            os.path.join(directory, filename)
            for filename in VelsLoader.list_fits_filenames(directory)
        ]
        
    @staticmethod
    def list_fits_filenames(directory=None):
        '''Lists all .fits filenames in vels_directory'''
        return [
            filename
            for filename in os.listdir(directory)
            if filename.lower().endswith(".fits")
        ]


class VelsLoader_CCF(VelsLoader):
    
    def __init__(self, object_name, ccf_directory=None):
        super().__init__()
        if ccf_directory is not None:
            self.ccf_directory = ccf_directory
        else:
            self.ccf_directory = os.path.join(self.vels_directory, 'ccf/' + object_name)
            if not os.path.exists(self.ccf_directory): self.ccf_directory = self.current_directory

        self.object_name = object_name
        self.filenames = self.list_fits_filenames(self.ccf_directory)
        self.filepaths = self.list_fits_filepaths(self.ccf_directory)
    
    def velocity_from_fits(self, filename, hdu_index=0):
        hdul = fits.open(os.path.join(self.ccf_directory, filename))
        return {'v': hdul[hdu_index].header['v'], 'e_v': hdul[hdu_index].header['e_v'], 'unit': 'cm/s'}

