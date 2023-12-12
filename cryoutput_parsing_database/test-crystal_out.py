from crystal_out import crystalOut
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import voigt_profile

Na3PS4 = crystalOut("cryoutput_parsing_database\crystal_output_files\calc-CoS2_205_icsd86351_tzvp_PBE0-checkSPIN\CoS2_205.out")
print("structure:\n", Na3PS4.structure)
print("space group", Na3PS4.space_group)
print("atomic masses", Na3PS4.atomic_masses)
print("vibrational contributions {mode:tensor}:\n", Na3PS4.vib_contributions_dielectric)
print("vibrational contributions sum:\n", Na3PS4.vib_contributions_dielectric_sum)
print("dielectric tensor:\n", Na3PS4.dielectric_tensor)
print("second electric susceptibility:\n", Na3PS4.second_electric_susceptibility)
print("third electric susceptibility:\n", Na3PS4.third_electric_susceptibiliy)
print("born charge:\n", Na3PS4.born_charge)
print("born charge on basis of normal modes:\n",Na3PS4.born_charge_normal_mode)
print("Raman and IR intensities \n", Na3PS4.raman_intensities)
print("thermodynamic terms:\n", Na3PS4.thermodynamic_terms)
print("temp, wavelength:\n", Na3PS4.raman_temp, Na3PS4.raman_wavelength)
print("IR tensor:\n", Na3PS4.get_IR_tensor().shape)
print("Raman tensor:\n", Na3PS4.get_raman_tensor().shape)
print("dynamical matrix:\n", Na3PS4.get_dynamical_matrix())

