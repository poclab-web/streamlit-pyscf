"""
理論手法や基底関数系を選択して、シングルポイントエネルギーを計算する。

機能:
- 計算に使用する理論手法と基底関数系を選択可能。
- 量子化学計算を実行し、結果をインタラクティブに表示。
"""

import streamlit as st
from logic.molecule_handler import MoleculeHandler
from logic.calculation import theory_options, basis_set_options
from logic.calculation import run_quantum_calculation
import stmol

from rdkit import Chem
import py3Dmol

st.title("PySCF + Streamlit: Quantum Chemistry on the Cloud")

# ユーザー入力
st.header("Molecular Input")
input_type = st.selectbox("Select Input Type", ["SMILES", "XYZ"])
atom_input = st.text_area(
    "Enter Molecular Structure",
    "H 0 0 0\nH 0 0 0.74" if input_type == "XYZ" else "CCO"
)

theory = st.selectbox("Theory", theory_options)
basis_set = st.selectbox("Basis Set", basis_set_options)

# 分子構造を処理
handler = None
if st.button("Process Molecule"):
    try:
        handler = MoleculeHandler(atom_input, input_type=input_type.lower())
        if not handler.mol:
            raise ValueError("Invalid molecular input. Please check your format.")

        st.subheader("2D Structure")
        handler.generate_2d_image("molecule_2d.png")
        st.image("molecule_2d.png", caption=Chem.MolToSmiles(handler.mol))

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

    except Exception as e:
        st.error(f"Error processing molecule: {e}")

    try:
        if not handler or not handler.mol:
            raise ValueError("Please process the molecule before running the calculation.")

        # Generate PySCF input format
        pyscf_input = handler.to_pyscf_input()

        st.write("Running quantum chemistry calculation...")
        energy = run_quantum_calculation(pyscf_input, basis_set, theory)
        st.success(f"Calculated SinglePoint Energy: {energy} Hartree")
    except Exception as e:
        st.error(f"Error: {e}")

