import streamlit as st
from pymongo import MongoClient
import json, pickle
import numpy as np
from scipy.special import voigt_profile
import pandas as pd
from pymatgen.core.structure import Structure
from plotly.graph_objs.layout import YAxis,XAxis,Margin
import plotly.graph_objects as go
from collections import OrderedDict

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
    data = OrderedDict()
    for s in doc:
        name = s["structure_name"]
        raman_IR_intensities = pd.DataFrame.from_dict(json.loads(s["raman_IR_intensities"]))
        structure = Structure.from_dict(json.loads(s["structure"]))
        born_charge = json.loads(s["born_charge"])
        raman_tensor = pickle.loads(s["raman_tensor"])
        data[name]={"raman_IR_intensities":raman_IR_intensities,"structure":structure,"born_charge":born_charge,"raman_tensor":raman_tensor}
    return data

def get_convoluted_spectra(intensities_df, raman=True, sigma=1, gamma=1, wavenumber_range=(0,1000), resolution=1000):
    data_freq = intensities_df["FREQ(CM**-1)"].to_list()
    if raman:
        data_int = intensities_df["RAMAN INTENSITIES"].to_list()
    else:
        data_int = intensities_df["INTENS"].to_list()

    frequencies = np.linspace(wavenumber_range[0], wavenumber_range[1], resolution)
    convoluted_intensities = np.zeros(resolution)
    for i in range(len(data_freq)):
        freq = float(data_freq[i])
        intensity = float(data_int[i])
        convoluted_amplitude = voigt_profile(frequencies-freq, sigma, gamma) * intensity
        #convolution_single_peak.append(convoluted_amplitude)
        convoluted_intensities += convoluted_amplitude

    return frequencies, convoluted_intensities

def show_structure(structure):
    pass

def plot_convoluted_spectra(x, y, raman=True):
    layout = go.Layout(title="Raman Spectrum from Voigt Convolution" if raman else "IR Spectrum from Voigt Convolution",
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
        names = [f"ICSD: {s.split('icsd')[-1]}" for s in structures.keys()]
        tabs = st.tabs(names)
        for i, data in enumerate(structures.values()):
            with tabs[i]:
                #show_structure(structures[i][2])
                raman_IR_intensities = data["raman_IR_intensities"]
                raman_IR_intensities.index.name='Mode'
                raman_intensities = raman_IR_intensities[raman_IR_intensities["RAMAN"] == "A"]
                IR_intensities = raman_IR_intensities[raman_IR_intensities["IR"] == "A"]
                container = st.container()
                col1, col2 = st.columns(2)
                with col1:
                    sigma = st.slider("Sigma",help="Standard deviation of Gaussian convolution",min_value=0.0,max_value=10.0,value=1.0,key=f"{names[i]}_sigma_raman")
                with col2:
                    gamma = st.slider("Gamma",help="Half-width at half-maximum of Cauchy distribution",min_value=0.0,max_value=10.0,value=1.0,key=f"{names[i]}_gamma_raman")
                wavenumber_range = st.slider('Wavenumber range', 0.0, 1000.0, (20.0, 600.0),step=10.0,key=f"{names[i]}_range_raman")
                frequencies, convoluted_intensities = get_convoluted_spectra(raman_intensities,raman=True,sigma=sigma,gamma=gamma,wavenumber_range=wavenumber_range)
                spectra_fig = plot_convoluted_spectra(frequencies, convoluted_intensities)
                container.plotly_chart(spectra_fig)
                st.subheader("Calculated Raman-active modes")
                raman_intensities.drop(["EIGV(Ha**2)", "IR","RAMAN","INTENS"], axis="columns", inplace=True)
                raman_intensities.columns.name = raman_intensities.index.name
                raman_intensities.index.name = None
                st.write(raman_intensities.to_html(), unsafe_allow_html=True)
                st.write(" ")

                col1, col2 = st.columns(2)
                with col1:
                    raman_tens_button = st.button("See Raman Tensor", key=f"raman_tens_button{i}", help="3Nx6 Raman tensor (6 independent tensor components for each mode)")
                with col2:
                    born_charge_button = st.button("See Born Charge", key=f"born_charge_button{i}", help="Born charge trace for each atom")
                data_container = st.empty()
                if raman_tens_button:
                    data_container.write(pd.DataFrame(data["raman_tensor"],columns=["XX","XY","XZ","YY","YZ","ZZ"]))
                if born_charge_button:
                    atoms = [x["Atomic symbol"] for x in data["born_charge"]]
                    born_charge_trace = [x["Dynamic Charge"]/3 for x in data["born_charge"]]
                    data_container.write(pd.DataFrame(data=born_charge_trace,index=atoms,columns=["Born Charge Trace"]))

                container = st.container()
                col1, col2 = st.columns(2)
                with col1:
                    sigma_IR = st.slider("Sigma",help="Standard deviation of Gaussian convolution",min_value=0.0,max_value=10.0,value=1.0,key=f"{names[i]}_sigma_IR")
                with col2:
                    gamma_IR = st.slider("Gamma",help="Half-width at half-maximum of Cauchy distribution",min_value=0.0,max_value=10.0,value=1.0,key=f"{names[i]}_gamma_IR")
                wavenumber_range_IR = st.slider('Wavenumber range', 0.0, 1000.0, (20.0, 600.0),step=10.0,key=f"{names[i]}_range_IR")
                frequencies, convoluted_intensities = get_convoluted_spectra(IR_intensities,raman=False,sigma=sigma_IR,gamma=gamma_IR,wavenumber_range=wavenumber_range_IR)
                spectra_fig = plot_convoluted_spectra(frequencies, convoluted_intensities,raman=False)
                container.plotly_chart(spectra_fig)

    else:
        st.write("We don't have Raman spectra for this compound yet.")
