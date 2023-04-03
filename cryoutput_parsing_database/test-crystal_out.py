from crystal_out import crystalOut
import numpy as np
import matplotlib.pyplot as plt

Na3PS4 = crystalOut("cryoutput_parsing_database\crystal17_output_files\calc-AsS2_4_icsd424590_tzvp_PBE0\AsS2_4.out")
print("structure:\n", Na3PS4.structure)
print("space group", Na3PS4.space_group)
print("atomic masses", Na3PS4.atomicMasses)
print("vibrational contributions {mode:tensor}:\n", Na3PS4.vibContributionsDielectric)
print("vibrational contributions sum:\n", Na3PS4.vibContributionsDielectricSum)
print("dielectric tensor:\n", Na3PS4.dielectricTensor)
print("second electric susceptibility:\n", Na3PS4.secondElectricSusceptibility)
print("third electric susceptibility:\n", Na3PS4.thirdElectricSusceptibility)
print("born charge:\n", Na3PS4.bornCharge)
print("born charge on basis of normal modes:\n",Na3PS4.bornChargeNormalModeBasis)
print("polycrystalline intensities", Na3PS4.intRaman)
print("thermodynamic terms:\n", Na3PS4.thermodynamicTerms)
print("temp, wavelength:\n", Na3PS4.raman_temp, Na3PS4.raman_wavelength)
#print("IR tensor:\n", Na3PS4.get_IR_tensor().shape)
#print("Raman tensor:\n", Na3PS4.get_raman_tensor().shape)
#print("dynamical matrix:\n", Na3PS4.get_dynamical_matrix())
#f, i = Na3PS4.get_convoluted_spectra(5, 5)
#plt.plot(f, i)
#plt.show()
