import glob
import pandas as pd
import numpy as np
import json
from crystal_out import crystalOut

all_crystal_outputs = glob.glob("crystal17_output_files/*/*.out")
print(len(all_crystal_outputs), "structures in total")

print("Structures that failed to run:")
print("popped:", all_crystal_outputs.pop(52))
print("popped:", all_crystal_outputs.pop(52))

dict_db = {"structure":[], "spaceGroup":[],"thermodynamicTerms":[],
                   "dielectricTensor_flattened":[], "dielectricTensor_shape":[], 
                   "vibContributionsDielectric_flattened":[], "vibContributionsDielectric_shape":[],
                   "secondElectricSusceptibility_flattened":[], "secondElectricSusceptibility_shape":[],
                   "thirdElectricSusceptibility_flattened":[], "thirdElectricSusceptibility_shape":[],
                   "bornChargeArray_flattened":[], "bornChargeArray_shape":[],  
                   "bornChargeNormalModeBasis_flattened":[], "bornChargeNormalModeBasis_shape":[],
                   "intRaman":[]}

for s in all_crystal_outputs:
    print("Processing: ", s)
    selected_structure = crystalOut(s)
    dict_db["structure"].append(json.dumps(selected_structure.structure.as_dict()))
    dict_db["spaceGroup"].append(selected_structure.parsed_data["space_group"])
    dict_db["thermodynamicTerms"].append(json.dumps(selected_structure.thermodynamicTerms))
    dict_db["intRaman"].append(json.dumps(selected_structure.intRaman))
    
    dict_db["dielectricTensor_shape"].append(str(selected_structure.dielectricTensor.shape)) # shape is saved so we can recover flattened array later
    dict_db["dielectricTensor_flattened"].append(str(selected_structure.dielectricTensor.flatten().tobytes()))
    vibContributionsDielectric = np.array(list(selected_structure.vibContributionsDielectric.values()))
    dict_db["vibContributionsDielectric_shape"].append(str(vibContributionsDielectric.shape))
    dict_db["vibContributionsDielectric_flattened"].append(str(vibContributionsDielectric.flatten().tobytes()))
    
    dict_db["secondElectricSusceptibility_shape"].append(str(selected_structure.secondElectricSusceptibility.shape))
    dict_db["secondElectricSusceptibility_flattened"].append(str(selected_structure.secondElectricSusceptibility.flatten().tobytes()))
    dict_db["thirdElectricSusceptibility_shape"].append(str(selected_structure.thirdElectricSusceptibility.shape))
    dict_db["thirdElectricSusceptibility_flattened"].append(str(selected_structure.thirdElectricSusceptibility.flatten().tobytes()))

    bornChargeArrayList = []
    for val in selected_structure.bornCharge.values():
        bornChargeArrayList.append(val["Born Charge"])
    bornChargeArray = np.concatenate(bornChargeArrayList, dtype = "float64")
    dict_db["bornChargeArray_shape"].append(str(bornChargeArray.shape))
    dict_db["bornChargeArray_flattened"].append(str(bornChargeArray.flatten().tobytes()))
    
    dict_db["bornChargeNormalModeBasis_shape"].append(str(selected_structure.bornChargeNormalModeBasis.shape))
    dict_db["bornChargeNormalModeBasis_flattened"].append(str(selected_structure.bornChargeNormalModeBasis.flatten().tobytes()))

df = pd.DataFrame.from_dict(dict_db)
df.to_csv(f"All_Structures_Data_{len(all_crystal_outputs)}.csv")

## use return_structure_data from crystal_out.py to get back individual structures data from csv.