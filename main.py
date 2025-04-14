import defines

import streamlit as st

with open(defines.MAIN_TEMPLATE) as f:
    template = f.read()

if __name__ == "__main__":
    st.write(template)