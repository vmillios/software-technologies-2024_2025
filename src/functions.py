import streamlit as st
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np



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
        placeholder="plot color",
        key="cpfirst"
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

def choose_plot_harmony(adata):
    option = st.selectbox(
        "Select plot color",
        ("batch", "celltype"),
        index = None,
        placeholder="plot color",
        key="cpharmony"
    )
    if option == "batch":
        plot_show_batch(adata)
    elif option == "celltype":
        plot_show_celltype(adata)
    else:
        pass

def plot_show_celltype_harmony(adata):
    fig, ax = plt.subplots(figsize=(6, 6))
    sc.pl.umap(adata, color=['batch'], legend_fontsize=10, ax=ax, show=False)
    st.pyplot(fig)

def plot_show_batch_harmony(adata):
    fig, ax = plt.subplots(figsize=(6, 6))
    sc.pl.umap(adata, color=['celltype'], legend_fontsize=10, ax=ax, show=False)
    st.pyplot(fig)

def volcano_plot(degs_df):
    degs_df["neg_log10_pval"] = -np.log10(degs_df["pvals"])

    # Add a column for differential expression classification
    degs_df["diffexpressed"] = "NS"
    degs_df.loc[(degs_df["logfoldchanges"] > 1) & (degs_df["pvals"] < 0.05), "diffexpressed"] = "UP"
    degs_df.loc[(degs_df["logfoldchanges"] < -1) & (degs_df["pvals"] < 0.05), "diffexpressed"] = "DOWN"

    # Select top downregulated genes (prioritize by highest significance, then most negative log2FC)
    top_downregulated = degs_df[degs_df["diffexpressed"] == "DOWN"]
    top_downregulated = top_downregulated.sort_values(by=["neg_log10_pval", "logfoldchanges"], ascending=[False, True]).head(20)

    # Select top upregulated genes (prioritize by highest significance, then most positive log2FC)
    top_upregulated = degs_df[degs_df["diffexpressed"] == "UP"]
    top_upregulated = top_upregulated.sort_values(by=["neg_log10_pval", "logfoldchanges"], ascending=[False, False]).head(81)

    # Combine top genes
    top_genes_combined = pd.concat([top_downregulated["genes"], top_upregulated["genes"]])
    df_annotated = degs_df[degs_df["genes"].isin(top_genes_combined)]

    fig, ax = plt.subplots(figsize=(10, 6))

# Create Volcano plot
    sns.scatterplot(
        data=degs_df,
        x="logfoldchanges",
        y="neg_log10_pval",
        hue="diffexpressed",
        palette={"UP": "#bb0c00", "DOWN": "#00AFBB", "NS": "grey"},
        alpha=0.7,
        edgecolor=None,
        ax=ax
)

    # Add threshold lines
    ax.axhline(y=-np.log10(0.05), color='gray', linestyle='dashed')
    ax.axvline(x=-1, color='gray', linestyle='dashed')
    ax.axvline(x=1, color='gray', linestyle='dashed')

    # Labels and formatting
    ax.set_xlim(-11, 11)
    ax.set_ylim(25, 175)
    ax.set_xlabel("log2 Fold Change", fontsize=14)
    ax.set_ylabel("-log10 p-value", fontsize=14)
    ax.set_title("Volcano of DEGs (Disease vs Control)", fontsize=16)
    ax.legend(title="Expression", loc="upper right")

    # Show the plot in Streamlit
    st.pyplot(fig)