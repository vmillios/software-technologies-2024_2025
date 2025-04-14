import streamlit as st
import defines

with open(defines.ABOUT_TEMPLATE, encoding="utf-8") as f:
    template = f.read()

st.write(template)