import glob
import pandas as pd
import numpy as np
import json
import pickle
from bson.binary import Binary
from pymongo import MongoClient
from pymatgen.core import Lattice, Structure
from crystal_out import crystalOut

def crystal_mongo_drone():
    uri = "mongodb+srv://RamanML:CaReRamanProject00@ramanml.utaye2e.mongodb.net/?retryWrites=true&w=majority"
    mc = MongoClient(uri)
    db = mc["raman_ml"]
    all_crystal_outputs = glob.glob("cryoutput_parsing_database/crystal17_output_files/calc-*/*.out")
    print("Found", len(all_crystal_outputs), "structures to be loaded into mongodb...")

    structures = db.structures
    for s in all_crystal_outputs:
        name = s.split("calc-")[-1].split("_tzvp")[0]
        print("Working on", name)
        if name == "CoS2_205_icsd86351":
            continue
        try:
            selected_structure = crystalOut(s)
        except ValueError as e:
            print(name, e)
        bornChargeArrayList = []
        for val in selected_structure.bornCharge.values():
            bornChargeArrayList.append(val["Born Charge"])
        bornChargeArray = np.concatenate(bornChargeArrayList, dtype = "float64")

        structure = {"structure_name": name,
                    "structure":json.dumps(selected_structure.structure.as_dict()), 
                    "spaceGroup":selected_structure.space_group,
                    "thermodynamicTerms":json.dumps(selected_structure.thermodynamicTerms),
                    "dielectricTensor":Binary(pickle.dumps(selected_structure.dielectricTensor, protocol=2), subtype=128),
                    "vibContributionsDielectric":Binary(pickle.dumps(selected_structure.vibContributionsDielectric, protocol=2), subtype=128),
                    "vibContributionsDielectricSum":Binary(pickle.dumps(selected_structure.vibContributionsDielectricSum, protocol=2), subtype=128),
                    "secondElectricSusceptibility":Binary(pickle.dumps(selected_structure.secondElectricSusceptibility, protocol=2), subtype=128),
                    "thirdElectricSusceptibility":Binary(pickle.dumps(selected_structure.thirdElectricSusceptibility, protocol=2), subtype=128),
                    "bornChargeArray":Binary(pickle.dumps(selected_structure.bornCharge, protocol=2), subtype=128),
                    "bornChargeNormalModeBasis":Binary(pickle.dumps(selected_structure.bornChargeNormalModeBasis, protocol=2), subtype=128),
                    "intRaman":json.dumps(selected_structure.intRaman),
                    "ramanTemp": selected_structure.raman_temp,
                    "ramanWavelength": selected_structure.raman_wavelength}

        structure_id = structures.insert_one(structure).inserted_id
        print("Done", structure_id)


def load_structure_from_mongo(structure_name):  
    uri = "mongodb+srv://RamanML:CaReRamanProject00@ramanml.utaye2e.mongodb.net/?retryWrites=true&w=majority"
    mc = MongoClient(uri)
    db = mc["raman_ml"]
    structures = db.structures
    myquery = { "structure_name": { "$regex": f"^{structure_name}" } }
    #structure_name_out = structure_name + ".out"
    doc = structures.find(myquery)

    data = []
    for s in doc:
        structure = Structure.from_dict(json.loads(s["structure"]))
        spaceGroup = s["spaceGroup"]
        thermodynamicTerms = json.loads(s["thermodynamicTerms"])
        intRaman = json.loads(s["intRaman"])

        dielectricTensor = pickle.loads(s["dielectricTensor"])
        vibContributionsDielectric = pickle.loads(s["vibContributionsDielectric"])
        secondElectricSusceptibility = pickle.loads(s["secondElectricSusceptibility"])
        thirdElectricSusceptibility = pickle.loads(s["thirdElectricSusceptibility"])
        bornChargeArray = pickle.loads(s["bornChargeArray"])
        bornChargeNormalModeBasis = pickle.loads(s["bornChargeNormalModeBasis"])

        data.append([structure, spaceGroup, thermodynamicTerms, intRaman, dielectricTensor, vibContributionsDielectric, secondElectricSusceptibility, thirdElectricSusceptibility, bornChargeArray, bornChargeNormalModeBasis])
    
    return data


if __name__ == '__main__':
    crystal_mongo_drone()
