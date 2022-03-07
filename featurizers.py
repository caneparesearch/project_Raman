# Structure featurizers (SOAP, ACSF, CrystalNN)

# TODO: Labels 

from crystal_out import crystalOut
from dscribe.descriptors import SOAP, ACSF
from pymatgen.io.ase import AseAtomsAdaptor
from matminer.featurizers.site import CrystalNNFingerprint

Na3PO4_114 = crystalOut("crystal17_output_files/calc-Na3PO4_114_icsd97205_tzvp_PBE0/Na3PO4_114.out")
atom_list = ["Na", "P", "O"]

# SOAP featurizer 
Na3_ASE = AseAtomsAdaptor.get_atoms(Na3PO4_114.structure)
soap = SOAP(periodic=True, 
    rcut=6,
    nmax=8,
    lmax=6,
    sigma=1,
    species = atom_list)
SOAP_features = soap.create(Na3_ASE) # returns each position as a row with features
print("SOAP features shape: ", SOAP_features.shape)

# ACSF featurizer 
acsf = ACSF(species = atom_list,
    rcut=6.0,
    periodic = True,
    g2_params=[[1, 1], [1, 2], [1, 3]],
    g4_params=[[1, 1, 1], [1, 2, 1], [1, 1, -1], [1, 2, -1]],
)
ACSF_features = acsf.create(Na3_ASE)
print("ACSF features shape: ", ACSF_features.shape)

# CrystalNN featurizer 
# initialize the featurizer to use structural order parameters (i.e., octahedral, tetrahedral) rather than coordination number.
cnnf = CrystalNNFingerprint.from_preset("ops")
# get the features for the first site in the structure.
site_1_features = cnnf.featurize(Na3PO4_114.structure, 0) 
print(cnnf.feature_labels())
print(site_1_features)

############################################
# SOAP featurizer from matminer doesnt work;

#from matminer.featurizers.site.external import SOAP
# soap = SOAP(periodic=True, rcut=6, nmax=8, lmax=6, sigma=1)
# soap.fit([Na3PO4_114.structure])
# features = soap.featurize(Na3PO4_114.structure,1)
# df = pd.DataFrame(np.array([soap.feature_labels(), features]).T)
# print(df)