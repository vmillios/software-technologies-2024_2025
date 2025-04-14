import defines

import streamlit as st

with open("templates/main.md", encoding="utf-8") as f:
    template = f.read()

if __name__ == "__main__":
    st.write(template)