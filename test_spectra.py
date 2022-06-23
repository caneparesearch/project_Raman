from crystal_out import crystalOut
from scipy.special import voigt_profile
import numpy as np
import matplotlib.pyplot as plt

Na3PS4 = crystalOut("crystal17_output_files\\calc-Na3PS4_114_icsd121566_tzvp_PBE0\\Na3PS4_114.out")
intensities = Na3PS4.intRaman

def raman_convolution(intensities, sigma, gamma, padding=30, resolution=1000):
    data = np.array(list(intensities.items()))
    frequencies = np.linspace(min(data[:,0])-padding, max(data[:,0])+padding, resolution)
    convoluted_intensities = np.zeros((resolution))
    for i in range(data.shape[0]):
        freq = data[i,0]
        intensity = data[i,1]
        convolution_single_peak = []
        for x in frequencies:
            convoluted_amplitude = voigt_profile(x-freq, sigma, gamma)*intensity
            convolution_single_peak.append(convoluted_amplitude)
        convoluted_intensities += np.array(convolution_single_peak)

    return frequencies, convoluted_intensities

f, i = raman_convolution(intensities, 5, 5, resolution=1000)
plt.plot(f, i)
plt.show()