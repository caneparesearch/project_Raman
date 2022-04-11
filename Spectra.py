import numpy as np
import matplotlib.pyplot as plt

class Spectra:
    """
    Raman spectra object
    """
    def __init__(self, frequencies, intensities):
        self.frequencies = frequencies
        self.intensities = intensities

    def lorentzian_transform(self, sigma):
        """
        Transforms a calculated spectra into lorentzian shape for each of the calculated intensities
        """
        frequency_range = np.linspace(0,self.frequencies[-1]+100,2000)
        lorentzian_intensities = np.zeros(2000)
        for i in range(len(self.frequencies)):
            center = self.frequencies[i]
            intensity = self.intensities[i]
            lorentz_intensity = np.array([lorentzian(frequency,center,sigma)*intensity for frequency in frequency_range])
            lorentzian_intensities += lorentz_intensity
        return Spectra(frequency_range,lorentzian_intensities)

    def plot_spectra(self):
        plt.plot(self.frequencies,self.intensities)
        return plt.figure()

    def __str__(self):
        return np.array([self.frequencies,self.intensities]).T.__str__()

def lorentzian(x, center, sigma):        
    """
    returns value of lorentzian probability at x
    """
    return 1 / np.pi * sigma / ((x - center)**2 + sigma**2)