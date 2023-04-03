import streamlit as st
#import numpy as np
#import pandas as pd
#import plotly.express as px

import plotly.graph_objects as go

#### functions ####
def fetch_data_from_mongo(structure_name):
    """
    todo
    # use load_structure_from_mongo
    """
    return 

def show_data_on_page(structure_data):
    """
    todo
    """
    st.write(structure_data)
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
    structure_data, spectra = fetch_data_from_mongo(structure_name)
    if structure_data:
        show_data_on_page(structure_data)
        
        spectra_fig = plot_raman_spectra(frequencies, convoluted_intensities)
        st.plotly_chart(spectra_fig)
    else:
        st.write("We don't have Raman spectra for this compound yet.")
