import glob
from crystal_out import crystalOut

all_crystal_outputs = glob.glob("./crystal17_output_files/calc-*/*.out")
for s in all_crystal_outputs:
    name = s.split("calc-")[-1].split("_tzvp")[0]
    if name == "CoS2_205_icsd86351":
        continue
    print("Working on", name)
    try:
        selected_structure = crystalOut(s)
    except ValueError as e:
        print(name, e)