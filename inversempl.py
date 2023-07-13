import numpy as np
import scipy
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
from matplotlib.colors import ListedColormap
from pympl import PyMPL

molecular_lidar_ratio = 8*np.pi/3 #Extinction to backscatter ratio for molecular scatterers

_atmosphere_data = np.loadtxt('us_standard.csv', skiprows=1, delimiter=',') # Read US Standard Atmosphere data
_height   = _atmosphere_data[:,0]
_pressure = _atmosphere_data[:,1] * 100 #pa (file uses hpa)
_temp_K   = _atmosphere_data[:,2]

US_standard_pressure    = scipy.interpolate.interp1d(_height, _pressure) # Interpolate US standard Atmosphere data
US_standard_temperature = scipy.interpolate.interp1d(_height, _temp_K)

def rayleigh_backscatter_coeff_at_height(LAM, height):
    '''
    LAM in um, height in km
    '''
    return rayleigh_backscatter_coeff(LAM, US_standard_pressure(height), US_standard_temperature(height))

def rayleigh_backscatter_coeff(LAM, pressure, temperature):
    '''
    Using equations from Gray G. Gimmestad and David W. Roberts Page 40
    input: wavelength (um), pressure (pa), temperature (K) 
    return: volume bsckscatter coefficient (km-1 sr-1)
    '''
    stp_backscatter_coeff = 1.39 * (0.55/LAM)**4 * 1E-6
    return stp_backscatter_coeff * (pressure * 288.15) / (101325 * temperature) * 1000

def Fernald_inversion_inwards(nrb, lidar_range, S1, S2 = molecular_lidar_ratio, calibration_beta_1 = 0, calibration_range = 20, LAM = 0.532):
    '''
    Inwards 2-components Fernald method for backscatter coefficient profile
    This function is adapted from Equation (6) of Frederick G. Fernald, Analysis of atmospheric lidar observations: some comments. 1984. Applied Optics
    input:
        nrb: 1d or 2d attenuated backscatter data. dimension is [timestamp (optional), range]; unit is (counts km2 us-1 uJ-1)
        lidar_range: height in km
        S1: Lidar ratio (extinction to backscatter ratio) for particulate scatterers
        S2: Lidar ratio for molecular scatterers
        calibration_beta_1: Aerosol backscatter coefficient calibration value at calibration range
        calibration_range: Calibration range
        LAM: lidar wavelength in um. Default to 532 nm for micropulse lidar

    return:
        beta1: retrieved backscatter coefficient (km-1 sr-1). dimension is [timestamp (optional), range]
        inv_range: range used for the inversion
    '''
    nrb = np.transpose(np.array(nrb)) # transpose so the dimension is now [range, timestamp] since timestamp is optional
    average_nrb  = nrb[lidar_range<=calibration_range]  #counts * km2 / micro s micro j 
    inv_range = lidar_range[lidar_range<=calibration_range]

    delta_r = np.mean(inv_range[1:]-inv_range[:-1]) #km
    beta2  = rayleigh_backscatter_coeff_at_height(LAM, inv_range) # standard atmosphere molecular backscatter coefficient km-1 sr-1
    result = list(range(inv_range.size))

    calibration_backscatter = calibration_beta_1+beta2[-1]
    if nrb.ndim == 1:
        result[-1] = calibration_backscatter # calibration data
    elif nrb.ndim == 2:
        result[-1] = np.array([calibration_backscatter]*nrb.shape[1])
    for i in reversed(range(inv_range.size)):
        if i>0:
            A_term = (S1-S2) * (beta2[i-1]+beta2[i]) * delta_r
            numera = average_nrb[i-1] * np.exp(A_term)
            C_term = average_nrb[i]/result[i]
            denomi = C_term + S1*(average_nrb[i]+average_nrb[i-1]*np.exp(A_term)) * delta_r
            result[i-1] = numera / denomi  # Each loop is calculating for i-1
        else:
            pass
    
    return np.transpose(np.array(result)), inv_range