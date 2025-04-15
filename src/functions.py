import streamlit as st
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt



def print_data_preprocessing(adata: sc.AnnData) -> None:
    df = pd.DataFrame(
        {
            "atribute": ["n_obs", "n_vars", "obs"],
            "data": [adata.n_obs, adata.n_vars, ", ".join(adata.obs_keys())]
        }
    )
    st.dataframe(df)

def print_data_postprocessing(adata: sc.AnnData) -> None:
    df = pd.DataFrame(
        {
            "atribute": ["obs", "var", "uns", "obsm", "varm", "obsp"],
            "data": [", ".join(adata.obs_keys()), ", ".join(adata.var_keys()), ", ".join(adata.uns_keys()), ", ".join(adata.obsm_keys()), ", ".join(adata.varm_keys()), ", ".join(adata.obsp)]
        }
    )
    st.dataframe(df)

def preprocess(adata: sc.AnnData, min_genes_value: int, min_cells_value: int, target_sum_value: int, min_mean_value: float, max_mean_value: int, min_disp_value: float) -> sc.AnnData:
    sc.pp.filter_cells(adata, min_genes=min_genes_value)
    sc.pp.filter_genes(adata, min_cells=min_cells_value)
    adata = adata[:, [gene for gene in adata.var_names if not str(gene).startswith(tuple(['ERCC', 'MT-', 'mt-']))]]
    sc.pp.normalize_total(adata, target_sum=target_sum_value)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=min_mean_value, max_mean=max_mean_value, min_disp=min_disp_value)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    
    return adata

def choose_plot(adata : sc.AnnData):
    option = st.selectbox(
        "Select plot color",
        ("batch", "celltype"),
        index = None,
        placeholder="plot color"
    )
    if option == "batch":
        plot_show_batch(adata)
    elif option == "celltype":
        plot_show_celltype(adata)
    else:
        pass

def plot_show_batch(adata: sc.AnnData) -> None:
    fig, ax = plt.subplots(figsize=(6, 6))
    sc.pl.umap(adata, color=['batch'], legend_fontsize=10, ax=ax, show=False)
    st.pyplot(fig)

def plot_show_celltype(adata: sc.AnnData) -> None:
    fig, ax = plt.subplots(figsize=(6, 6))
    sc.pl.umap(adata, color=['celltype'], legend_fontsize=10, ax=ax, show=False)
    st.pyplot(fig)