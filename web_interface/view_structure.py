import streamlit as st
"""
Draft - need to debug 
"""
st.markdown("# Compound List")
st.sidebar.markdown("# Compound List")
import pickle
import py3Dmol
import pymatgen.io.cif as cif
from pymatgen.core import Structure


# lattice = [[3.84, 0, 0], [0, 3.84, 0], [0, 0, 3.84]]
# coords = [[0, 0, 0], [0.5, 0.5, 0.5]]
# species = ["Si", "Si"]
# structure = Structure(lattice, species, coords)

# # Save the pymatgen structure to a CIF file
# cif_file = cif.CifWriter(structure)
# cif_file.write_file("structure.cif")

# Load the CIF file into py3Dmol and display the structure
with open("test.cif", "r") as f:
    cif_data = f.read()

view = py3Dmol.view()

view.addModel(cif_data, "cif",
    {"doAssembly" : True,
    "normalizeAssembly":True,
    'duplicateAssemblyAtoms':True})
view.setStyle({'sphere':{"scale":0.25},
                'stick':{"radius":0.15}})

view.addUnitCell()
view.zoomTo()
view.render()
# Get the JavaScript code for the view and display it in Streamlit
t = view.js()
html = t.startjs + "\n" + t.endjs + "\n"
st.components.v1.html(html, width=600, height=400)



# with open("data.pkl", 'rb') as f:
#     pickle_data = pickle.load(f)

# # cif_file = cif.CifWriter(pickle_data[0][0])
# cif_file = cif.CifWriter(structure)
# cif_file.write_file("structure.cif")

# # Load the CIF file into py3Dmol and display the structure
# with open("structure.cif", "r") as f:
#     cif_data = f.read()
# view = py3Dmol.view()
# view.addModel(cif_data, "cif")
# view.setStyle({'sphere': {"scale": 0.3}})
# view.addUnitCell()
# view.zoomTo()
# view.show()
# view.addModel(cif, "cif",
#     {"doAssembly" : True,
#     "normalizeAssembly":True,
#     'duplicateAssemblyAtoms':True})
# view.setStyle({'sphere':{"scale":0.15},
#                 'stick':{"radius":0.25}})
# view.addUnitCell()
# view.zoomTo()
# view.render()
# t = view.js()
# html = t.startjs + "\n" + t.endjs + "\n"
# st.components.v1.html(html, width=600, height=400)
