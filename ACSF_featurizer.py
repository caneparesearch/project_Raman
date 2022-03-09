from crystal_out import crystalOut
from dscribe.descriptors import ACSF
from pymatgen.io.ase import AseAtomsAdaptor

Na3PO4_114 = crystalOut("crystal17_output_files/calc-Na3PO4_114_icsd97205_tzvp_PBE0/Na3PO4_114.out")
Na3PS4_114 = crystalOut("crystal17_output_files/calc-Na3PS4_114_icsd121566_tzvp_PBE0/Na3PS4_114.out")

O_ASE = AseAtomsAdaptor.get_atoms(Na3PO4_114.structure)
S_ASE = AseAtomsAdaptor.get_atoms(Na3PS4_114.structure)

acsf = ACSF(species = ["Na", "P", "S", "O", "Li"],
    rcut=6.0,
    periodic = True,
    g2_params=[[1, 1], [1, 2], [1, 3]],
    g4_params=[[1, 1, 1], [1, 2, 1], [1, 1, -1], [1, 2, -1]],
)
ACSF_features = acsf.create([S_ASE, O_ASE])
print(ACSF_features, ACSF_features[0].shape, ACSF_features[1].shape)