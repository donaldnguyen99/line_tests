import numpy as np

@np.vectorize
def air2vac(wavelength_air): 
    '''
    Wavelength in vacuum from wavelength in air in units of Angstroms (from VALD)
    '''
    s = 1e4 / wavelength_air
    n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s * s) + 0.0001599740894897 / (38.92568793293 - s * s)
    return wavelength_air * n

@np.vectorize
def vac2air(wavelength):
    '''
    Wavelength in air from wavelength in vacuum in units of Angstroms (from VALD)
    '''
    s = 1e4 / wavelength
    return wavelength / (1 + 0.0000834254 + 0.02406147 / (130 - s * s) + 0.00015998 / (38.9 - s * s))

@np.vectorize
def relativistic_rv_shift(rv_cm_per_s):
    '''
    Wavelength in rest frame from wavelength in observed frame and relativistic radial velocity in units of cm/s
    '''
    return np.sqrt(1 - rv_cm_per_s / 29979245800) / np.sqrt(1 + rv_cm_per_s / 29979245800)