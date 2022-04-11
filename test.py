from crystal_out import crystalOut
from Spectra import Spectra
import numpy as np
import matplotlib.pyplot as plt

Na3PS4 = crystalOut("crystal17_output_files/calc-Na3PS4_114_icsd121566_tzvp_PBE0/Na3PS4_114.out")
intensities = Na3PS4.intensities
#print(intensities)
normal_spectra = Spectra(intensities[:,0], intensities[:,1])
fig = normal_spectra.plot_spectra()
print(normal_spectra)

l_spectra = normal_spectra.lorentzian_transform(0.5)

print(l_spectra)
fig = l_spectra.plot_spectra()
plt.show()