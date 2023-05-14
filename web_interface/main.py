import streamlit as st
from pymongo import MongoClient
import json
import numpy as np
from scipy.special import voigt_profile
import pandas as pd
from pymatgen.core.structure import Structure
from pymatgen.io.cif import CifWriter
from plotly.graph_objs.layout import YAxis,XAxis,Margin
import plotly.graph_objects as go

st.set_page_config(page_title="Raman Database")
                   
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
        raman_intensities = json.loads(s["raman_intensities"])
        structure = Structure.from_dict(json.loads(s["structure"]))
        data.append([name,raman_intensities,structure])
    return data

def get_convoluted_spectra(intRaman, sigma=1, gamma=1, wavenumber_range=(0,1000), resolution=1000):
    data_freq = list(intRaman.keys())
    data_int = [x[0] for x in intRaman.values()]

    frequencies = np.linspace(wavenumber_range[0], wavenumber_range[1], resolution)
    convoluted_intensities = np.zeros((resolution))
    for i in range(len(data_freq)):
        freq = float(data_freq[i])
        intensity = float(data_int[i])
        convolution_single_peak = []
        for x in frequencies:
            convoluted_amplitude = voigt_profile(x - freq, sigma, gamma) * intensity
            convolution_single_peak.append(convoluted_amplitude)
        convoluted_intensities += np.array(convolution_single_peak)

    return frequencies, convoluted_intensities

def show_structure(structure):
    return

def plot_raman_spectra(x, y):
    layout = go.Layout(title="Raman Spectrum from Voigt Convolution",
                       xaxis=XAxis(title="Wavenumber (cm<sup>-1</sup>)"),
                       xaxis2 = XAxis(title="THz", overlaying= 'x', side= 'top'),)

    fig = go.Figure(layout=layout)
    fig.add_trace(go.Scatter(x=x, y=y))
    fig.add_trace(go.Scatter(xaxis='x2'))
    fig.update_layout(yaxis_title="Intensity (arbitrary units)")
    start, end = x[0], x[-1]
    fig.update_layout(xaxis1=dict(range=[start, end]), xaxis2=dict(range=[start/0.3335641E+02, end/0.3335641E+02]))
    return fig

#### streamlit main page ####

col1, col2 = st.columns([3,1])
with col1:
    st.markdown("""
    # Main page""")
with col2:
    st.image("web_interface/logo-inter-color.png")
st.sidebar.markdown("# Main page")

st.title("Hybrid-Functional Computational Raman Database for Inorganic Compounds")

# search bar 
structure_name = st.text_input("Search a compound: e.g. As2Se3")

if structure_name:
    structures = fetch_data_from_mongo(structure_name)
    if len(structures) > 0:
        names = [f"ICSD: {s[0].split('icsd')[-1]}" for s in structures]
        tabs = st.tabs(names)
        for i in range(len(names)):
            with tabs[i]:
                show_structure(structures[i][2])
                container = st.container()
                col1, col2 = st.columns(2)
                with col1:
                    sigma = st.slider("Sigma",help="Standard deviation of Gaussian convolution",min_value=0.0,max_value=10.0,value=1.0,key=f"{names[i]}_sigma")
                with col2:
                    gamma = st.slider("Gamma",help="Half-width at half-maximum of Cauchy distribution",min_value=0.0,max_value=10.0,value=1.0,key=f"{names[i]}_gamma")
                wavenumber_range = st.slider('Wavenumber range', 0.0, 1000.0, (25.0, 600.0),step=10.0,key=f"{names[i]}_range")
                frequencies, convoluted_intensities = get_convoluted_spectra(structures[i][1],sigma,gamma,wavenumber_range)
                spectra_fig = plot_raman_spectra(frequencies, convoluted_intensities)
                container.plotly_chart(spectra_fig)
                st.subheader("Calculated values")
                df = pd.DataFrame(structures[i][1].keys(), columns=["Wavenumber (cm-1)"],dtype=float)
                df.insert(1,"THz", df["Wavenumber (cm-1)"]/0.3335641E+02)
                df.insert(2,"Intensity (Arbitary Units)", [x[0] for x in structures[i][1].values()])
                df.insert(3,"IRREP", [x[1] for x in structures[i][1].values()])
                st.table(df)
    else:
        st.write("We don't have Raman spectra for this compound yet.")
