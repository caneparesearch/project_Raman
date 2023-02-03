from crystal_out import crystalOut
import glob

all_crystal_outputs = glob.glob("../project_Raman/crystal17_output_files/*/*.out")
all_structures = {}
for fi in all_crystal_outputs:
    structure_name = fi.split("\\")[-1].split(".")[0]
    with open(fi, "r") as f:
        data = crystalOut(fi)
        all_structures[structure_name] = data.structure

print("Number of structures processed:", len(all_structures.keys()))
print(all_structures.keys())
