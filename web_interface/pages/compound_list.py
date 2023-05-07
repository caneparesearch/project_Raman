import streamlit as st
from pymongo import MongoClient

st.markdown("# Compound List")
st.sidebar.markdown("# Compound List")

@st.cache_resource
def init_connection():
    return MongoClient(**st.secrets["mongo"])

client = init_connection()

@st.cache_data(ttl=600)
def get_compound_list():
    db = client["raman_ml"]
    list = []
    for structure in db.structures.find():
        list.append(structure["structure_name"])
    return list

compound_list = get_compound_list()
st.table(compound_list)
