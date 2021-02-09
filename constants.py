__author__ = 'Awab Asher'

import scipy.constants
import numpy as np

from simulation_parameters import(
    TOTAL_BANDWIDTH,
    NOISE_FIGURE
)

# GLOBAL CONSTANTS!
SPEED_OF_LIGHT = scipy.constants.c
BOLTZMANN_CONSTANT = scipy.constants.k
temp = 293   # Temperature in Kelvin!


# Spectrum
M_2 = 2     # Number of operators!
M_3 = 3     # Number of operators!
M_4 = 4     # Number of operators!
M_5 = 5     # Number of operators!
exclusive_bandwidth = (TOTAL_BANDWIDTH / M_2) * 10**9
exclusive_bandwidth_iii = (TOTAL_BANDWIDTH / M_3) * 10**9
exclusive_bandwidth_iv = (TOTAL_BANDWIDTH / M_4) * 10**9
exclusive_bandwidth_v = (TOTAL_BANDWIDTH / M_5) * 10**9
pooled_bandwidth = TOTAL_BANDWIDTH * 10**9

# Noise
noise_spectral_density = BOLTZMANN_CONSTANT * temp
noise_power_exclusive = (10 * np.log10(noise_spectral_density * exclusive_bandwidth)) + NOISE_FIGURE
noise_power_pooled = (10 * np.log10(noise_spectral_density * pooled_bandwidth)) + NOISE_FIGURE




