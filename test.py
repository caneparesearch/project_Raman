from crystal_out import crystalOut

Na3PS4 = crystalOut(".\crystal17_output_files\calc-Te3As2_12_icsd18208_tzvp_PBE0\Te3As2_12.out")
print("structure:\n", Na3PS4.structure)
print("space group", Na3PS4.parsed_data["space_group"])
print("atomic masses", Na3PS4.atomic_masses)
print("vibrational contributions {mode:tensor}:\n", Na3PS4.vibrational_contributions)
print("dielectric tensor:\n", Na3PS4.dielectric_tensor)
print("second electric susceptibility:\n", Na3PS4.second_electric_susceptibility)
print("third electric susceptibility:\n", Na3PS4.third_electric_susceptibility)
print("born charge:\n", Na3PS4.born_charge)
print("thermodynamic terms:\n", Na3PS4.thermodynamic_terms)