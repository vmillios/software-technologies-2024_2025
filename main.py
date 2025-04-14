import defines

import streamlit as st

with open("templates/main.md", encoding="utf-8") as f:
    template = f.read()

def main():
    st.write(template)
    uploaded_file = st.file_uploader("Upload a h5ad file", type=["h5ad"])
    if uploaded_file:
        st.write("file has been uploaded")

if __name__ == "__main__":
    main()