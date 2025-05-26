import defines
import src.functions

import streamlit as st
import scanpy as sc
import tempfile

with open(defines.MAIN_TEMPLATE, encoding="utf-8") as f:
    template = f.read()

def main() -> None:

    st.set_page_config(page_title="Single Cell App", layout="wide")

    st.write(template)
    uploaded_file = st.file_uploader("Upload a h5ad file", type=["h5ad"])
    if uploaded_file:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp:
            tmp.write(uploaded_file.read())
            tmp_path = tmp.name
        try:
            adata = sc.read(tmp_path)
            st.success("File uploaded successfully!")
            src.functions.print_data_preprocessing(adata)

            st.header("preprocessing parameters")
            min_genes_value =st.number_input ("min_genes_value",value=600)
            min_cells_value =st.number_input ("min_cells_value",value=3)
            target_sum_value =st.number_input ("target_sum_value",value=1e4)
            min_mean_value =st.number_input ("min_mean_value",value=0.0125)
            max_mean_value =st.number_input ("max_mean_value",value=3)
            min_disp_value =st.number_input ("min_disp_value",value=0.5)
            if st.button("start preprocessing"):
                adata = src.functions.preprocess(adata, min_genes_value, min_cells_value, target_sum_value, min_mean_value, max_mean_value, min_disp_value)
                src.functions.print_data_postprocessing(adata)
                st.session_state["adata"] = adata
            if "adata" in st.session_state:
                src.functions.choose_plot(st.session_state["adata"])
            



        except Exception as e:
            st.error(f"Failed to load file: {e}")


if __name__ == "__main__":
    main()
    
