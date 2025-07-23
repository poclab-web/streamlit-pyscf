"""
入力形式の確認を行うもの
SMILESやXYZから2次元や3次元の構造を確認する。
SDF, XYZ, zmatrixへの変換も行う。
"""
import os
import streamlit as st
from utils.module import load_css

from rdkit import Chem
import stmol
import py3Dmol

from logic.molecule_handler import MoleculeHandler

# カスタムCSSを適用
load_css("config/styles.css")

# Streamlit app
st.title("Molecule Viewer and Converter")

# Input type selection
input_type = st.selectbox("Select input type", ["XYZ", "SMILES"])
atom_input = st.text_area(
    "Enter Molecular Structure",
    "C     -2.385271    1.648249    0.000003\nH     -3.181824    1.061240    0.407216\nH     -2.536116    1.773363   -1.051896\nH     -2.370100    2.607079    0.474676\nO     -1.453042    1.151314    0.170017\nH     -1.466511    0.300016   -0.251422"
    if input_type == "XYZ"
    else "CO",
)

handler = None
start_processing = False

if st.button("Run Check Structure"):
    try:
        handler = MoleculeHandler(atom_input, input_type=input_type.lower())
        if not handler.mol:
            raise ValueError("Invalid molecular input. Please check your format.")

        # 化合物名を取得
        compound_name = Chem.MolToInchiKey(handler.mol)
        smiles = Chem.MolToSmiles(handler.mol)

        # ディレクトリの作成
        directory = os.path.join("data", compound_name)
        os.makedirs(directory, exist_ok=True)

        col1, col2 = st.columns(2)

        # Display 2D structure in the first column
        with col1:
            st.subheader("Input 2D Structure")
            handler.generate_2d_image(f"{directory}/molecule_2d.png")
            st.image(f"{directory}/molecule_2d.png", caption=Chem.MolToSmiles(handler.mol))

        # Display 3D structure in the second column
        with col2:
            st.subheader("Input 3D Structure")
            try:
                mol_block = handler.generate_3d_molblock()
                viewer = py3Dmol.view(width=400, height=400)  # Adjust width to fit in the column
                viewer.addModel(mol_block, "mol")
                viewer.setStyle({"stick": {}})
                viewer.zoomTo()
                stmol.showmol(viewer, height=400)
            except Exception as e:
                st.warning(f"Unable to generate 3D structure: {e}")

    except Exception as e:
        st.error(f"Error processing molecule: {e}")


    # Display molecule 2D and 3D
    if handler and handler.mol:
        st.subheader("InChIKey")
        st.code(Chem.MolToInchiKey(handler.mol))

        st.code(mol_block)

        # pyscf input
        st.subheader("PySCF Input")
        pyscf_input = handler.to_pyscf_input()
        st.code("PySCF Input", value=pyscf_input, height=300)




