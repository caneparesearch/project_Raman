import glob
import json
import pickle
from bson.binary import Binary
from pymongo import MongoClient
from pymatgen.core import Lattice, Structure
from crystal_out import crystalOut
import os
SECRET = os.environ['SECRET_URI']

def crystal_mongo_drone():
    uri = SECRET
    mc = MongoClient(uri)
    db = mc["raman_ml"]
    all_crystal_outputs = glob.glob("cryoutput_parsing_database/crystal_output_files/calc-*/*.out")
    print("Found", len(all_crystal_outputs), "structures to be loaded into mongodb...")

    structures = db.structures
    # structures.drop()
    # for doc in structures.find():
    #     print(doc["structure_name"])
    # quit()
    errors = []
    for s in all_crystal_outputs:
        name = s.split("calc-")[-1].split("_tzvp")[0]
        print("Working on", name)
        if name == "CoS2_205_icsd86351" or name == "CaF2_225_icsd60559" or name == "K4P3_63_icsd64625": # dielctric tensor got problem
            continue

        myquery = { "structure_name": { "$regex": f'^{name}' } }
        if structures.count_documents(myquery, limit = 1) == 0:
            print(f"{name} not in database. Adding...")
            try:
                selected_structure = crystalOut(s)
            except ValueError as e:
                print(name, e)
                errors.append(f"{name}: {e}")
                continue
            try:
                tens_raman = selected_structure.get_raman_tensor()
            except FileNotFoundError:
                print(name, "No TENS_RAMAN file")
                tens_raman = None
            except Exception as e:
                print(e)
                tens_raman = None
            for item in selected_structure.born_charge:
                item["Born Charge"] = item["Born Charge"].tolist()

            structure = {"structure_name": name,
                        "structure":json.dumps(selected_structure.structure.as_dict()), 
                        "space_group":selected_structure.space_group,
                        "thermodynamic_terms":json.dumps(selected_structure.thermodynamic_terms),
                        "dielectric_tensor":Binary(pickle.dumps(selected_structure.dielectric_tensor, protocol=2), subtype=128),
                        "vib_contributions_dielectric":Binary(pickle.dumps(selected_structure.vib_contributions_dielectric, protocol=2), subtype=128),
                        "vib_contributions_dielectric_sum":Binary(pickle.dumps(selected_structure.vib_contributions_dielectric_sum, protocol=2), subtype=128),
                        "second_electric_susceptibility":Binary(pickle.dumps(selected_structure.second_electric_susceptibility, protocol=2), subtype=128),
                        "third_electric_susceptibiliy":Binary(pickle.dumps(selected_structure.third_electric_susceptibiliy, protocol=2), subtype=128),
                        "born_charge":json.dumps(selected_structure.born_charge),
                        "born_charge_normal_mode":Binary(pickle.dumps(selected_structure.born_charge_normal_mode, protocol=2), subtype=128),
                        "raman_IR_intensities":json.dumps(selected_structure.raman_intensities.to_dict()),
                        "raman_temp": selected_structure.raman_temp,
                        "raman_wavelength": selected_structure.raman_wavelength,
                        "raman_tensor": Binary(pickle.dumps(tens_raman, protocol=2), subtype=128)}

            structure_id = structures.insert_one(structure).inserted_id
            print("Done", structure_id)
    print("Errors:\n", "\n".join(errors))


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
        space_group = s["space_group"]
        thermodynamic_terms = json.loads(s["thermodynamic_terms"])
        raman_intensities = json.loads(s["raman_intensities"])

        dielectric_tensor = pickle.loads(s["dielectric_tensor"])
        vib_contributions_dielectric = pickle.loads(s["vib_contributions_dielectric"])
        second_electric_susceptibility = pickle.loads(s["second_electric_susceptibility"])
        third_electric_susceptibiliy = pickle.loads(s["third_electric_susceptibiliy"])
        born_charge = json.loads(s["born_charge"])
        born_charge_normal_mode = pickle.loads(s["born_charge_normal_mode"])

        data.append([structure, space_group, thermodynamic_terms, raman_intensities, dielectric_tensor, vib_contributions_dielectric, second_electric_susceptibility, third_electric_susceptibiliy, born_charge, born_charge_normal_mode])

    return data

if __name__ == '__main__':
    crystal_mongo_drone()