from crystal_out import crystalOut
import pandas as pd
import numpy as np
#from matminer.featurizers.site.external import SOAP
from dscribe.descriptors import SOAP, ACSF
from pymatgen.io.ase import AseAtomsAdaptor

Na3PO4_114 = crystalOut("crystal17_output_files/calc-Na3PO4_114_icsd97205_tzvp_PBE0/Na3PO4_114.out")
#print(Na3PO4_114.structure)
# soap = SOAP(periodic=True, rcut=6, nmax=8, lmax=6, sigma=1)
# soap.fit([Na3PO4_114.structure])
# features = soap.featurize(Na3PO4_114.structure,1)
# df = pd.DataFrame(np.array([soap.feature_labels(), features]).T)
# print(df)

Na3_ASE = AseAtomsAdaptor.get_atoms(Na3PO4_114.structure)
soap = SOAP(periodic=True, 
    rcut=6,
    nmax=8,
    lmax=6,
    sigma=1,
    species = ["Na", "P", "O"])
SOAP_features = soap.create(Na3_ASE) # returns each position as a row with features
print(SOAP_features.shape)

acsf = ACSF(species = ["Na", "P", "O"],
    rcut=6.0,
    periodic = True,
    g2_params=[[1, 1], [1, 2], [1, 3]],
    g4_params=[[1, 1, 1], [1, 2, 1], [1, 1, -1], [1, 2, -1]],
)
ACSF_features = acsf.create(Na3_ASE)
print(ACSF_features.shape)