from crystal_out import crystalOut

Na3PS4 = crystalOut("crystal17_output_files\\calc-Te3As2_12_icsd18208_tzvp_PBE0\\Te3As2_12.out")
print("structure:\n", Na3PS4.structure)
print("space group", Na3PS4.parsed_data["space_group"])
print("atomic masses", Na3PS4.atomicMasses)
print("vibrational contributions {mode:tensor}:\n", Na3PS4.vibContributionsDielectric)
print("dielectric tensor:\n", Na3PS4.dielectricTensor)
print("second electric susceptibility:\n", Na3PS4.secondElectricSusceptibility)
print("third electric susceptibility:\n", Na3PS4.thirdElectricSusceptibility)
print("born charge:\n", Na3PS4.bornCharge)
print("born charge on basis of normal modes:\n",Na3PS4.bornChargeNormalModeBasis)
print("polycrystalline intensities", Na3PS4.intRaman)
print("thermodynamic terms:\n", Na3PS4.thermodynamicTerms)