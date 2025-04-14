import streamlit as st
import defines

with open(defines.ABOUT_TEMPLATE) as f:
    template = f.read()

st.write(template)