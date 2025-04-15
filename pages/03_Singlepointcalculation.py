"""
理論手法や基底関数系を選択して、シングルポイントエネルギーを計算する。

機能:
- 計算に使用する理論手法と基底関数系を選択可能。
- 量子化学計算を実行し、結果をインタラクティブに表示。
"""

import streamlit as st

import stmol
import pandas as pd
from rdkit import Chem
import py3Dmol
import os
import numpy as np

from utils.module import load_css
from logic.molecule_handler import MoleculeHandler
from logic.calculation import theory_options, basis_set_options
from logic.calculation import run_quantum_calculation
from logic.output_handler import extract_orbital_energies

# カスタムCSSを適用
load_css("config/styles.css")

# 本文
st.title("Single Point Calculation")

# ユーザー入力
st.header("Molecular Input")
input_type = st.selectbox("Select Input Type", ["XYZ", "SMILES"])
atom_input = st.text_area(
    "Enter Molecular Structure",
    "C     -2.385271    1.648249    0.000003\nH     -3.181824    1.061240    0.407216\nH     -2.536116    1.773363   -1.051896\nH     -2.370100    2.607079    0.474676\nO     -1.453042    1.151314    0.170017\nH     -1.466511    0.300016   -0.251422"
    if input_type == "XYZ"
    else "CO",
)
theory = st.selectbox("Theory", theory_options)
basis_set = st.selectbox("Basis Set", basis_set_options)

# Charge and Spin入力セクション
with st.expander("Other Settings"):
  charge = st.number_input("Molecular Charge", min_value=-10, max_value=10, value=0, step=1)
  multiplicity = st.number_input("Spin Multiplicity (2S + 1)", min_value=1, max_value=10, value=1, step=1)
  # スピン計算
  spin = multiplicity - 1
  # 溶媒効果の設定
  solvent_model = st.selectbox("Select Solvent Model", ["None", "PCM"])
  eps = None  # epsのデフォルト値を設定

  # symmetryの選択肢を追加
  symmetry = st.selectbox("Consider Molecular Symmetry?", ["Yes", "No"])
  symmetry = True if symmetry == "Yes" else False

  # Load solvent data
  solvents_file = "config/solvents_epsilon.csv"
  solvents_data = pd.read_csv(solvents_file)

  if solvent_model == "PCM":
      # Show solvent selection dropdown
      solvent_selection = st.selectbox(
          "Select a solvent",
          [f"{row['Solvent']} (ε={row['Epsilon']})" for _, row in solvents_data.iterrows()]
      )  
      # Extract epsilon from selection
      if solvent_selection:
        eps = float(solvent_selection.split("=", 1)[-1][:-1])

      # Additional epsilon input override
      eps_input = st.selectbox("Override epsilon value?", ["No", "Yes"])

      if eps_input == "Yes":
          eps = st.number_input("Dielectric Constant (ε)", min_value=1.0, max_value=100.0, value=eps)

# 分子構造を処理
handler = None
if st.button("Run Single Point Calculation"):
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

# 計算の実行
    try:
        if not handler or not handler.mol:
            raise ValueError("Please process the molecule before running the calculation.")

        # Generate PySCF input format
        pyscf_input = handler.to_pyscf_input()

        st.write("Running quantum chemistry calculation...")
        energy = run_quantum_calculation(
            compound_name, smiles, pyscf_input, basis_set, theory, charge=charge, spin=spin, solvent_model=solvent_model, eps=eps, symmetry=symmetry
        )
        
        # 計算結果を表示
        st.success(f"Calculated SinglePoint Energy: {energy} Hartree")
        st.info(f"Results saved in: {directory}")

        # 軌道エネルギーとHOMO/LUMOを表示
        molden_file_path = os.path.join(directory, f"{compound_name}_{theory}.molden")
        if os.path.exists(molden_file_path):
            orbital_data = extract_orbital_energies(molden_file_path)

            # データの検証
            if not orbital_data or "orbital_energies" not in orbital_data:
                st.error("Failed to extract orbital energies. Please check the MOLDEN file.")
                raise ValueError("Invalid orbital data.")

            homo_index = orbital_data.get("homo_index")
            lumo_index = orbital_data.get("lumo_index")

            if homo_index is None or lumo_index is None:
                st.error("HOMO or LUMO index could not be determined. Please check the orbital data.")
                raise ValueError("HOMO or LUMO index is missing.")

            # 軌道エネルギーをデータフレームに変換し、エネルギーでソート
            orbital_df = pd.DataFrame({
                "Orbital": [f"Orbital {i + 1}" for i in range(len(orbital_data["orbital_energies"]))],
                "Energy (Hartree)": orbital_data["orbital_energies"]
            }).sort_values(by="Energy (Hartree)", ascending=False).reset_index()

            # 軌道エネルギーを表示（Energy (Hartree) のみ）
            st.subheader("Orbital Energies")
            st.dataframe(orbital_df[["Energy (Hartree)"]])  # 必要な列だけを表示

            # ソート後のHOMOとLUMOのインデックスを再計算
            sorted_homo_index = orbital_df[orbital_df["index"] == homo_index].index[0]
            sorted_lumo_index = orbital_df[orbital_df["index"] == lumo_index].index[0] if lumo_index is not None else None

            # HOMOとLUMOの情報を表示
            st.subheader("HOMO and LUMO")
            if sorted_lumo_index is not None:
                st.write(f"LUMO: Index {sorted_lumo_index + 1}, Energy {orbital_df.loc[sorted_lumo_index, 'Energy (Hartree)']:.6f} Hartree")
            else:
                st.write("LUMO: Not found")
            if sorted_homo_index is not None:
                st.write(f"HOMO: Index {sorted_homo_index + 1}, Energy {orbital_df.loc[sorted_homo_index, 'Energy (Hartree)']:.6f} Hartree")
            else:
                st.write("HOMO: Not found")

        else:
            st.warning("MOLDEN file not found. Orbital energies cannot be displayed.")
    except Exception as e:
        st.error(f"Error: {e}")
