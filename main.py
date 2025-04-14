import defines
import src.functions

import streamlit as st
import scanpy as sc
import tempfile

with open(defines.MAIN_TEMPLATE, encoding="utf-8") as f:
    template = f.read()

def main() -> None:
    st.write(template)
    uploaded_file = st.file_uploader("Upload a h5ad file", type=["h5ad"])
    if uploaded_file:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp:
            tmp.write(uploaded_file.read())
            tmp_path = tmp.name
        try:
            adata = sc.read(tmp_path)
            st.success("File uploaded successfully!")
            src.functions.print_data(adata)

        except Exception as e:
            st.error(f"Failed to load file: {e}")


if __name__ == "__main__":
    main()