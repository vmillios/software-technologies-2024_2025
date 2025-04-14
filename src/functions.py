import streamlit as st
import scanpy as sc

def print_data(adata: sc.AnnData) -> None:
    st.write(f"obs: {adata.obs_names}")