import crystal_out
import pandas as pd
import numpy as np
from matminer.featurizers.site.external import SOAP

Na3PO4_114 = crystal_out.crystalOut("crystal17_output_files/calc-Na3PO4_114_icsd97205_tzvp_PBE0/Na3PO4_114.out")
#print(Na3PO4_114.structure)
soap = SOAP(periodic=False, rcut=6, nmax=8, lmax=6, sigma=1)
soap.fit([Na3PO4_114.structure])
features = soap.featurize(Na3PO4_114.structure,1)
df = pd.DataFrame(np.array([soap.feature_labels(), features]).T)
print(df)