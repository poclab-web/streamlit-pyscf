"""
入力形式の確認を行うもの
"""


import streamlit as st
import stmol
from logic.molecule_handler import MoleculeHandler
import os

import py3Dmol

# Streamlit app
st.title("Molecule Viewer and Converter")

# Input type selection
input_type = st.selectbox("Select input type", ["SMILES", "XYZ"])

handler = None
start_processing = False

if input_type == "SMILES":
    smiles = st.text_area("Enter SMILES string")
    if smiles:
        st.subheader("SMILES Output")
        st.text(smiles)
        if st.button("Process Molecule"):
            handler = MoleculeHandler(smiles, input_type="smiles")
            start_processing = True

elif input_type == "XYZ":
    xyz_input = st.text_area("Enter XYZ format(symbols)")
    if xyz_input:
        st.subheader("XYZ Content")
        st.text(xyz_input)
        if st.button("Process Molecule"):
            handler = MoleculeHandler(xyz_input, input_type="xyz")
            start_processing = True

# Display molecule 2D and 3D
if start_processing and handler and handler.mol:

    st.subheader("2D Structure")
    handler.generate_2d_image("molecule_2d.png")
    st.image("molecule_2d.png", caption="2D Structure")

    st.subheader("3D Structure")
    try:
        mol_block = handler.generate_3d_molblock()
        viewer = py3Dmol.view(width=800, height=400)
        viewer.addModel(mol_block, "mol")
        viewer.setStyle({"stick": {}})
        viewer.zoomTo()
        stmol.showmol(viewer, height=400)
    except Exception as e:
        st.warning(f"Unable to generate 3D structure: {e}")

    # SDF display and download
    sdf_file_path = "output.sdf"
    handler.save_to_sdf(sdf_file_path)
    with open(sdf_file_path, "r") as f:
        sdf_content = f.read()
    st.subheader("SDF Content")
    st.text_area("SDF Content", value=sdf_content, height=300)
    with open(sdf_file_path, "rb") as f:
        st.download_button(
            label="Download SDF File",
            data=f,
            file_name="molecule.sdf",
            mime="chemical/x-mdl-sdfile"
        )

    # XYZ display and download
    xyz_file_path = "output.xyz"
    handler.save_to_xyz(xyz_file_path)
    with open(xyz_file_path, "r") as f:
        xyz_content = f.read()
    st.subheader("XYZ Content")
    st.text_area("XYZ Content", value=xyz_content, height=300)
    with open(xyz_file_path, "rb") as f:
        st.download_button(
            label="Download XYZ File",
            data=f,
            file_name="molecule.xyz",
            mime="chemical/x-xyz"
        )

    # Z-Matrix display and download
    zmatrix_file_path = "output.zmatrix"
    handler.save_to_zmatrix(zmatrix_file_path)
    with open(zmatrix_file_path, "r") as f:
        zmatrix_content = f.read()
    st.subheader("Z-Matrix Content")
    st.text_area("Z-Matrix Content", value=zmatrix_content, height=300)
    with open(zmatrix_file_path, "rb") as f:
        st.download_button(
            label="Download Z-Matrix File",
            data=f,
            file_name="molecule.zmatrix",
            mime="text/plain"
        )
else:
    if not start_processing:
        st.info("Please provide input and click 'Process Molecule' to start.")
