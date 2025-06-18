"""
振動計算とIR計算
構造最適化後のものの確認を行う。
他のもののcheckpointやXYZを呼び出して、行う。
TODO: 改修中
"""

import os
import streamlit as st
import py3Dmol  # 3D可視化用ライブラリ
import stmol

import numpy as np

import matplotlib.pyplot as plt
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from logic.molecule_handler import MoleculeHandler
from logic.calculation import calculate_vibrational_frequencies
from logic.calculation import theory_options, basis_set_options
from logic.visualization import generate_cjson, write_gaussian_log

st.title("Vibrational Frequency Calculator")

# ユーザー入力の設定
st.header("Molecule Input")
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


if st.button("Calculate"):
    # 計算の実行
    with st.spinner("Calculating vibrational frequencies..."):
        try:
            # 分子を処理
            handler = MoleculeHandler(atom_input, input_type=input_type.lower())

            # 化合物名を取得
            compound_name = Chem.MolToInchiKey(handler.mol)
            smiles = Chem.MolToSmiles(handler.mol)

            # ディレクトリの作成
            directory = os.path.join("data", compound_name)
            os.makedirs(directory, exist_ok=True)

            # 2D構造の生成
            col1, col2 = st.columns(2)
            with col1:
                st.subheader("Input 2D Structure")
                handler.generate_2d_image(f"{directory}/molecule_2d.png")
                st.image(f"{directory}/molecule_2d.png", caption=Chem.MolToSmiles(handler.mol))

            # 3D構造の生成
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

            pyscf_input = handler.to_pyscf_input()

            results = calculate_vibrational_frequencies(pyscf_input, theory, basis_set, charge, spin, compound_name=compound_name)
            frequencies = results['frequencies']
            thermo_info = results['thermo_info']


        except Exception as e:
            import traceback
            st.error(f"An error occurred: {e}")
            st.text(traceback.format_exc())
            frequencies = None


    if frequencies is not None:
        st.success("✅ Calculation completed successfully!")

        # --- Vibrational Frequencies ---
        st.subheader("🔬 Vibrational Frequencies (cm⁻¹)")
        freq_array = results['frequencies']['freq_wavenumber']
        # もしリストのリストならフラット化
        if any(isinstance(f, (list, tuple, np.ndarray)) for f in freq_array):
            import itertools
            freq_array = list(itertools.chain.from_iterable(freq_array))
        if hasattr(freq_array, "tolist"):
            freq_array = freq_array.tolist()
        st.markdown(f"**Number of modes:** {len(freq_array)}")
        # DataFrameを作成し、降順ソート
        freq_df = pd.DataFrame({
            "Frequency (cm⁻¹)": [float(f"{freq:.2f}") for freq in freq_array]
        }).sort_values("Frequency (cm⁻¹)", ascending=False).reset_index(drop=True)
        st.table(freq_df)

        # --- Thermodynamic Summary ---
        st.subheader("🌡 Thermodynamic Properties (298.15 K)")
        thermo = thermo_info

        # ndarray型をリストに変換
        for k, v in thermo.items():
            if hasattr(v, "tolist"):
                thermo[k] = v.tolist()

        # 変換定数
        EH_TO_KCAL = 627.5095

        st.markdown(f"""
        **Zero-Point Energy (ZPE):** {thermo["ZPE"][0]:.4f} Eh ({thermo["ZPE"][0]*EH_TO_KCAL:.2f} kcal/mol)  
        **Enthalpy (H_tot):** {thermo["H_tot"][0]:.4f} Eh ({thermo["H_tot"][0]*EH_TO_KCAL:.2f} kcal/mol)  
        **Gibbs Free Energy (G_tot):** {thermo["G_tot"][0]:.4f} Eh ({thermo["G_tot"][0]*EH_TO_KCAL:.2f} kcal/mol)  
        **Entropy (S_tot):** {thermo["S_tot"][0]:.6f} Eh/K ({thermo["S_tot"][0]*EH_TO_KCAL:.6f} kcal/mol·K)  
        """)

        st.info(f"Results saved in: {compound_name}")

        # --- Detailed View ---
        with st.expander("🔍 Full Thermodynamic Details"):
            df = {
                "Quantity": [
                    "ZPE", "E_tot", "H_tot", "G_tot", 
                    "S_tot", "Cv_tot", "Cp_tot"
                ],
                "Value (Eh)": [
                    thermo["ZPE"][0],
                    thermo["E_tot"][0],
                    thermo["H_tot"][0],
                    thermo["G_tot"][0],
                    thermo["S_tot"][0],
                    thermo["Cv_tot"][0],
                    thermo["Cp_tot"][0],
                ]
            }
            import pandas as pd
            st.dataframe(pd.DataFrame(df), use_container_width=True)