import streamlit as st
import scanpy as sc

def print_data(adata: sc.AnnData) -> None:
    st.write(f"n_obs: {adata.n_obs}")
    st.write(f"n_vars: {adata.n_vars}")
    st.write(f"obs:  {", ".join(adata.obs_keys())}")