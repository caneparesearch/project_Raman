from crystal_out import crystalOut
from dscribe.descriptors import SOAP
from pymatgen.io.ase import AseAtomsAdaptor

Na3PO4_114 = crystalOut("crystal17_output_files/calc-Na3PO4_114_icsd97205_tzvp_PBE0/Na3PO4_114.out")
Na3PS4_114 = crystalOut("crystal17_output_files/calc-Na3PS4_114_icsd121566_tzvp_PBE0/Na3PS4_114.out")

O_ASE = AseAtomsAdaptor.get_atoms(Na3PO4_114.structure)
S_ASE = AseAtomsAdaptor.get_atoms(Na3PS4_114.structure)

soap = SOAP(periodic=True, 
    rcut=6,
    nmax=8,
    lmax=6,
    sigma=1,
    species = ["Na", "P", "S", "O", "Li"],
    average="inner")
SOAP_features = soap.create([S_ASE, O_ASE]) 
print(SOAP_features, SOAP_features.shape)
