import numpy as np
import scipy
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
from matplotlib.colors import ListedColormap
from pympl import PyMPL
import pywt

#Extinction to backscatter ratio for molecular scatterers
molecular_lidar_ratio = 8*np.pi/3 

# Read U.S. Standard Atmosphere
_atmosphere_data = np.loadtxt('us_standard.csv', skiprows=1, delimiter=',') # Read US Standard Atmosphere data
_height   = _atmosphere_data[:,0]
_pressure = _atmosphere_data[:,1] * 100 #pa (file uses hpa)
_temp_K   = _atmosphere_data[:,2]

US_standard_pressure    = scipy.interpolate.interp1d(_height, _pressure) # Interpolate US standard Atmosphere data
US_standard_temperature = scipy.interpolate.interp1d(_height, _temp_K)

def calculate_backscatter_coefficient(LAM, pressure, temperature):
    '''
    Using equations from Gray G. Gimmestad and David W. Roberts Page 40
    input: wavelength (um), pressure (pa), temperature (K) 
    return: volume bsckscatter coefficient (km-1 sr-1)
    '''
    stp_backscatter_coeff = 1.39 * (0.55/LAM)**4 * 1E-6
    return stp_backscatter_coeff * (pressure * 288.15) / (101325 * temperature) * 1000

def rayleigh_bcksca_coeff(LAM, height, pressure_interpolater = US_standard_pressure, temperature_interpolater = US_standard_temperature):
    '''
    molecular backscatter of atmosphere at height. 
    Default is US standard atmosphere defined in John H. Seinfeld, Atmospheric Chemistry and Physics
    LAM in um, height in km
    '''
    return calculate_backscatter_coefficient(LAM, pressure_interpolater(height), temperature_interpolater(height))

def Fernald_inversion_inwards(nrb, mpl_range, beta2, S1, S2 = molecular_lidar_ratio, calibration_beta1 = 0, calibration_range = 20, LAM = 0.532):
    '''
    Inwards 2-Components Fixed Lidar-ratio Fernald Method Inversion for Backscatter Coefficient Profile
    This function is adapted from Equation (6) of Frederick G. Fernald, Analysis of atmospheric lidar observations: some comments. 1984. Applied Optics
    input:
        nrb: 1d or 2d attenuated backscatter data. dimension is [timestamp (optional), range]; unit is (counts km2 us-1 uJ-1)
        lidar_range: height in km
        beta2: rayleigh backscatter coefficient profile. Dimension can be 1d or same to nrb dimension
        S1: Lidar ratio (extinction to backscatter ratio) for particulate scatterers
        S2: Lidar ratio for molecular scatterers
        calibration_beta1: Aerosol backscatter coefficient calibration value at calibration range. A number or a 1D array
        calibration_range: Calibration range
        LAM: lidar wavelength in um. Default to 532 nm for micropulse lidar

    return:
        beta1: retrieved backscatter coefficient (km-1 sr-1). dimension is [timestamp (optional), range]
        inv_range: range used for the inversion
    '''
    nrb = np.transpose(np.array(nrb)) # transpose so the dimension is now [range, timestamp] since timestamp is optional
    beta2 = np.transpose(np.array(beta2))
    selected_range = mpl_range[mpl_range<=calibration_range]
    selected_nrb   = nrb[mpl_range<=calibration_range]  #counts * km2 / micro s micro j 
    selected_beta2 = beta2[mpl_range<=calibration_range]

    delta_r = np.mean(selected_range[1:]-selected_range[:-1]) #km
    total_backscatter = list(range(selected_range.size))
    calibration_beta2 = selected_beta2[-1]
    calibration_backscatter = calibration_beta1+calibration_beta2

    if (nrb.ndim - calibration_backscatter.ndim == 1):
        beta2_output = selected_beta2
        total_backscatter[-1] = calibration_backscatter # calibration data
    elif nrb.ndim - calibration_backscatter.ndim == 2:
        beta2_output = np.array([selected_beta2]*nrb.shape[1])
        total_backscatter[-1] = np.array([calibration_backscatter]*nrb.shape[1])
    else:
        raise ValueError('nrb can have 1 or 2 more dimensions than the boundary condition.\
                         beta2 can be a 2D array or an 1D array. calibration_beta1 can be a 1D array or a number')
    
    for i in reversed(range(selected_range.size)):
        if i>0:
            A_term = (S1-S2) * (selected_beta2[i-1]+selected_beta2[i]) * delta_r
            numera = selected_nrb[i-1] * np.exp(A_term)
            C_term = selected_nrb[i]/total_backscatter[i]
            denomi = C_term + S1*(selected_nrb[i]+selected_nrb[i-1]*np.exp(A_term)) * delta_r
            total_backscatter[i-1] = numera / denomi  # Each loop is calculating for i-1
        else:
            pass

    return np.transpose(np.array(total_backscatter)), beta2_output, selected_range

def Klett_inversion_inwards(nrb, mpl_range, S1, calibration_beta = 0, calibration_range = 20):
    '''
    Inwards 1-Component Fixed Lidar-ratio Klett Method Inversion for Bakcscatter Coefficient Profile
    This function is adapted from Equation (9) of Frederick G. Fernald, Analysis of atmospheric lidar observations: some comments. 1984. Applied Optics
    input:
        calibration_beta can be a number or an 1D array with time dimension
    '''

    nrb = np.transpose(np.array(nrb))
    selected_range = mpl_range[mpl_range<=calibration_range]
    selected_nrb   = nrb[mpl_range<=calibration_range]

    delta_r = np.mean(selected_range[1:]-selected_range[:-1]) #km
    total_backscatter = list(range(selected_range.size))
    
    if nrb.ndim == 1:
        if isinstance(calibration_beta, (int, float)): # nrb is 1D array and calibration_beta is a number
            total_backscatter[-1] = calibration_beta
        else:
            raise ValueError('calibration_beta has to be a number when nrb is an 1D array')
    elif nrb.ndim == 2:
        if isinstance(calibration_beta, (int, float)):
            total_backscatter[-1] = np.array([calibration_beta]*nrb.shape[1]) # nrb is 1D array and calibration_beta is a 1D array
        elif isinstance(np.array(calibration_beta), np.ndarray) and np.array(calibration_beta).ndim == 1:
            total_backscatter[-1] = np.array(calibration_beta)
        else:
            raise ValueError('calibration_beta has to be a number or an 1D array when nrb is a 2D array')
    else:
        raise ValueError('nrb has to be an 1D array or a 2D array')

    for i in reversed(range(selected_range.size)):
        if i>0:
            numera = selected_nrb[i-1]
            C_term = selected_nrb[i]/total_backscatter[i]
            denomi = C_term + (selected_nrb[i] + selected_nrb[i-1])*delta_r
            total_backscatter[i-1] = numera / denomi  # Each loop is calculating for i-1
        else:
            pass

    return np.transpose(np.array(total_backscatter)), selected_range

# Adapted from https://notebook.community/CSchoel/learn-wavelets/wavelet-denoising
# Plot the whole range of coefficients for a full decomposition with the DWT
def plot_dwt(details, approx, xlim=(-300,300), **line_kwargs):
    for i in range(len(details)):
        plt.subplot(len(details)+1,1,i+1)
        d = details[len(details)-1-i]
        half = len(d)//2
        xvals = np.arange(-half,-half+len(d))* 2**i
        plt.plot(xvals, d, **line_kwargs)
        plt.xlim(xlim)
        plt.title("detail[{}]".format(i))
    plt.subplot(len(details)+1,1,len(details)+1)
    plt.title("approx")
    plt.plot(xvals, approx, **line_kwargs)
    plt.xlim(xlim)

# Adapted from https://notebook.community/CSchoel/learn-wavelets/wavelet-denoising
# Based on Incorporating Information on Neighbouring Coefficients into Wavelet Estimation
def neigh_block(details, n, sigma):
    res = []
    L0 = int(np.log2(n) // 2)
    L1 = max(1, L0 // 2)
    L = L0 + 2 * L1
    def nb_beta(sigma, L, detail):
        S2 = np.sum(detail ** 2)
        lmbd = 4.50524 # solution of lmbd - log(lmbd) = 3
        beta = (1 - lmbd * L * sigma**2 / S2)
        return max(0, beta)
    for d in details:
        d2 = d.copy()
        for start_b in range(0, len(d2), L0):
            end_b = min(len(d2), start_b + L0)
            start_B = start_b - L1
            end_B = start_B + L
            if start_B < 0:
                end_B -= start_B
                start_B = 0
            elif end_B > len(d2):
                start_B -= end_B - len(d2)
                end_B = len(d2)
            assert end_B - start_B == L
            d2[start_b:end_b] *= nb_beta(sigma, L, d2[start_B:end_B])
        res.append(d2)
    return res


def _signal_smoothing(mpl_signal, mpl_range, nb_factor=0.1, wavelet = "bior3.5", merge_range = 5):
    # apply discrete wavelet smoothing to 
    mpl_signal = np.nan_to_num(mpl_signal, nan=0.0)
    coeffs_n = pywt.wavedec(mpl_signal, wavelet, mode = 'constant')
    approx_n = coeffs_n[0]
    details_n = coeffs_n[1:]

    new_detail = neigh_block(details_n, len(mpl_signal), nb_factor)

    smoothed_signal= pywt.waverec([approx_n] + new_detail, wavelet)[1:]

    smoothed_multiplier   =     scipy.stats.norm.cdf(mpl_range, merge_range, 1)
    unsmoothed_multiplier = 1 - smoothed_multiplier
    combined_signal = smoothed_signal * smoothed_multiplier + mpl_signal * unsmoothed_multiplier

    return combined_signal


def signal_smoothing(mpl_signal_2d, mpl_range, nb_factor=0.1, wavelet = "bior3.5", merge_range = 5):
    mpl_signal_2d = np.array(mpl_signal_2d)
    if mpl_signal_2d.ndim == 1:
        return _signal_smoothing(mpl_signal_2d, mpl_range, nb_factor=nb_factor, wavelet = wavelet, merge_range = merge_range)
    elif mpl_signal_2d.ndim == 2:
        combined_signal_list = []
        for i in range(mpl_signal_2d.shape[0]):
            combined_signal_list.append(_signal_smoothing(mpl_signal_2d[i], mpl_range, nb_factor=nb_factor, wavelet = wavelet, merge_range = merge_range))
        return np.array(combined_signal_list)
    else:
        raise ValueError('mpl_signal_2d has to be an 1D array or a 2D array')


def _movingaverage(values, window):
    if window%2 == 0:
        raise ValueError('Use odd moving average window')
    weights = np.repeat(1.0, window)/window
    sma = np.convolve(values, weights, mode='valid')
    pre_sma = np.empty(window//2)
    pre_sma[:] = np.nan
    sma = np.hstack([pre_sma, sma, pre_sma])
    return sma

def movingaverage(values, window):
    # values could either be an 1d array or a 2d array.
    values = np.array(values)
    if values.ndim == 1:
        return _movingaverage(values, window)
    elif values.ndim == 2:
        result = []
        for i in range(values.shape[0]):
            result.append(_movingaverage(values[i], window))
        result = np.array(result)
        return result
    else:
        raise ValueError('values has to be an 1D array or a 2D array')

def slope_inversion(nrb, mpl_range, slope_start, slope_end, lidar_ratio = 30):
    lbNRB = np.log(nrb)
    lbNRB = np.nan_to_num(lbNRB, nan=0.0)
    selected_lnNRB = lbNRB[np.logical_and(mpl_range >= slope_start, mpl_range < slope_end)]
    selected_range = mpl_range[np.logical_and(mpl_range >= slope_start, mpl_range < slope_end)]

    slope, intercept, r, p_value, std_err = scipy.stats.linregress(selected_range, selected_lnNRB)
    slope_ext_coeff = -0.5 * slope
    slope_bcksca_coeff = slope_ext_coeff/lidar_ratio

    return slope_ext_coeff, slope_bcksca_coeff, slope, intercept, r, p_value, std_err

