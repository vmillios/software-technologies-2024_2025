import streamlit as st
import scanpy as sc
import pandas as pd


def print_data(adata: sc.AnnData) -> None:
    df = pd.DataFrame(
        {
            "atribute": ["n_obs", "n_vars", "obs"],
            "data": [adata.n_obs, adata.n_vars, ", ".join(adata.obs_keys())]
        }
    )
    st.dataframe(df)
    # st.write(f"n_obs: {adata.n_obs}")
    # st.write(f"n_vars: {adata.n_vars}")
    # st.write(f"obs:  {", ".join(adata.obs_keys())}")