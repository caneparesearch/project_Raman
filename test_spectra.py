from crystal_out import crystalOut
import numpy as np
from scipy.special import voigt_profile

Na3PS4 = crystalOut("crystal17_output_files\\calc-Na3PS4_114_icsd121566_tzvp_PBE0\\Na3PS4_114.out")
intensities = list(Na3PS4.intRaman.values())
print(intensities)
print(voigt_profile(intensities,1,1))