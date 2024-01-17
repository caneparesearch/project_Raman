import streamlit as st
from pymongo import MongoClient
import pandas as pd

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
        full_name = structure["structure_name"].split("_")
        list.append([full_name[0], full_name[-1].replace("icsd","")])
    return list

compound_list = get_compound_list()
url_list = [["https://raman-db.streamlit.app/?query="+x[0], x[1]] for x in compound_list]
df = pd.DataFrame(url_list, columns=["Formula", "ICSD Number"])
st.dataframe(df, column_config={
    "Formula": st.column_config.LinkColumn(
        "Formula",
        display_text=r"https://raman-db.streamlit.app/\?query=(.*)")
}, width=500)
