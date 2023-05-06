import streamlit as st
from pymongo import MongoClient
import json
import numpy as np
from scipy.special import voigt_profile
import pandas as pd
from pymatgen.core.structure import Structure
from pymatgen.io.cif import CifWriter
import plotly.graph_objects as go
import py3Dmol

#### functions ####
@st.cache_resource
def init_connection():
    return MongoClient(**st.secrets["mongo"])

client = init_connection()

@st.cache_data(ttl=600)
def fetch_data_from_mongo(structure_name):
    db = client["raman_ml"]
    structures = db.structures
    myquery = { "structure_name": { "$regex": f"^{structure_name}_" } }
    doc = structures.find(myquery)
    data = []
    for s in doc:
        name = s["structure_name"]
        intRaman = json.loads(s["intRaman"])
        structure = Structure.from_dict(json.loads(s["structure"]))
        data.append([name,intRaman,structure])
    return data

def get_convoluted_spectra(intRaman, sigma, gamma, wavenumber_range=(0,1000), resolution=1000):
    data = np.array(list(intRaman.items()))
    frequencies = np.linspace(wavenumber_range[0], wavenumber_range[1], resolution)
    convoluted_intensities = np.zeros((resolution))
    for i in range(data.shape[0]):
        freq = float(data[i, 0])
        intensity = float(data[i, 1])
        convolution_single_peak = []
        for x in frequencies:
            convoluted_amplitude = voigt_profile(x - freq, sigma, gamma) * intensity
            convolution_single_peak.append(convoluted_amplitude)
        convoluted_intensities += np.array(convolution_single_peak)

    return frequencies, convoluted_intensities

def show_data_on_page(structure):
    cif_file = CifWriter(structure)

    view = py3Dmol.view()

    view.addModel(str(cif_file), "cif",
        {"doAssembly" : True,
        "normalizeAssembly":True,
        'duplicateAssemblyAtoms':False,
        'noComputeSecondaryStructure':False})
    view.setStyle({'sphere':{"scale":0.25},
                    'stick':{"radius":0.15}})

    view.addUnitCell()
    view.zoomTo()
    view.render()
    # Get the JavaScript code for the view and display it in Streamlit
    t = view.js()
    html = t.startjs + "\n" + t.endjs + "\n"
    st.components.v1.html(html, width=600, height=400)
    return 

def plot_raman_spectra(x, y):
    fig = go.Figure(data=go.Scatter(x=x, y=y))
    fig.update_layout(title="Raman Spectra", xaxis_title="Wavenumber", yaxis_title="Intensity (a.u.)")
    return fig

#### streamlit main page ####

st.markdown("# Main page")
st.sidebar.markdown("# Main page")

st.title("Hybrid-Functional Computational Raman Database for Inorganic Compounds")

# search bar 
structure_name = st.text_input("Search a compound: e.g. As2Se3")

if structure_name:
    structures = fetch_data_from_mongo(structure_name)
    if len(structures) > 0:
        names = [s[0] for s in structures]
        tabs = st.tabs(names)
        for i in range(len(names)):
            with tabs[i]:
                show_data_on_page(structures[i][2])
                col1, col2 = st.columns(2)
                with col1:
                    sigma = st.slider("Sigma",min_value=0.0,max_value=10.0,value=1.0,key=f"{names[i]}_sigma")
                with col2:
                    gamma = st.slider("Gamma",min_value=0.0,max_value=10.0,value=1.0,key=f"{names[i]}_gamma")
                wavenumber_range = st.slider('Wavenumber range', 0.0, 1000.0, (25.0, 600.0),step=10.0,key=f"{names[i]}_range")
                frequencies, convoluted_intensities = get_convoluted_spectra(structures[i][1],sigma,gamma,wavenumber_range)
                spectra_fig = plot_raman_spectra(frequencies, convoluted_intensities)
                st.plotly_chart(spectra_fig)
                st.subheader("Calculated values")
                df = pd.DataFrame(structures[i][1].items(), columns=["Wavenumber", "Intensity (a.u.)"])
                st.table(df)
    else:
        st.write("We don't have Raman spectra for this compound yet.")
