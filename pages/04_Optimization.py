"""
PySCFを使用して分子の幾何最適化を実行する
ユーザーは分子の構造、理論手法、基底関数系、および収束条件を入力し、幾何最適化計算を実行して
結果を可視化できます。

機能:
- XYZ形式で分子構造を入力可能。
- 幾何最適化の理論手法と基底関数系を選択可能。
- 収束条件（勾配、ステップ、エネルギーの許容値、および最大反復回数）を設定可能。
- 幾何最適化計算の進行中のエネルギー収束をプロット。
- 最終的な最適化構造を表示し、必要に応じてすべての構造をXYZファイルとして保存可能。
"""

import os

import streamlit as st
import pandas as pd
from rdkit import Chem
import py3Dmol
import stmol

from utils.module import load_css
from logic.molecule_handler import MoleculeHandler
from logic.calculation import theory_options, basis_set_options, run_geometry_optimization

# カスタムCSSを適用
load_css("config/styles.css")

# 本文
st.title("Geometry Optimization")

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

    solvents_file = "config/solvents_epsilon.csv"
    solvents_data = pd.read_csv(solvents_file)

    if solvent_model == "PCM":
        solvent_selection = st.selectbox(
            "Select a solvent",
            [f"{row['Solvent']} (ε={row['Epsilon']})" for _, row in solvents_data.iterrows()]
        )
        if solvent_selection:
            eps = float(solvent_selection.split("=", 1)[-1][:-1])

        eps_input = st.selectbox("Override epsilon value?", ["No", "Yes"])
        if eps_input == "Yes":
            eps = st.number_input("Dielectric Constant (ε)", min_value=1.0, max_value=100.0, value=eps)

# 収束条件の入力
with st.expander("Convergence Parameters"):
    convergence_energy = st.number_input(
        "Energy Tolerance (Hartree)", 
        min_value=1e-7, value=1e-6, step=1e-7, format="%.7f"
    )
    convergence_grms = st.number_input(
        "Gradient RMS Tolerance (Eh/Bohr)", 
        min_value=1e-5, value=3e-4, step=1e-5, format="%.5f"
    )
    convergence_gmax = st.number_input(
        "Gradient Max Tolerance (Eh/Bohr)", 
        min_value=1e-5, value=4.5e-4, step=1e-5, format="%.5f"
    )
    convergence_drms = st.number_input(
        "Displacement RMS Tolerance (Angstrom)", 
        min_value=1e-4, value=1.2e-3, step=1e-4, format="%.4f"
    )
    convergence_dmax = st.number_input(
        "Displacement Max Tolerance (Angstrom)", 
        min_value=1e-4, value=1.8e-3, step=1e-4, format="%.4f"
    )
    max_iterations = st.number_input(
        "Max Iterations", 
        min_value=1, value=100, step=1
    )

conv_params = {
    "convergence_energy": convergence_energy,
    "convergence_grms": convergence_grms,
    "convergence_gmax": convergence_gmax,
    "convergence_drms": convergence_drms,
    "convergence_dmax": convergence_dmax,
    "max_steps": max_iterations,
}

# 最適化の実行
if st.button("Run Geometry Optimization"):
    try:
        # 分子を処理
        handler = MoleculeHandler(atom_input, input_type=input_type.lower())
        if not handler.mol:
            raise ValueError("Invalid molecular input. Please check your format.")

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

        # PySCF用入力データを作成
        pyscf_input = handler.to_pyscf_input()

        # 最適化の実行
        st.write("Running geometry optimization...")
        final_geometry = run_geometry_optimization(
            compound_name, smiles, pyscf_input, basis_set, theory, 
            charge=charge, spin=spin, solvent_model=solvent_model, eps=eps, conv_params=conv_params
        )

        # 最適化後の構造を表示
        st.write("### Final Optimized Geometry")
        for i, coord in enumerate(final_geometry):
            st.text(f"Atom {i + 1}: {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}")

        save_path = os.path.join(directory, "optimized_geometry.xyz")
        with open(save_path, 'w') as f:
            f.write(f"{len(final_geometry)}\nOptimized structure\n")
            for i, coord in enumerate(final_geometry):
                f.write(f"Atom {i + 1} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n")
        st.success(f"Optimized geometry has been saved to: {save_path}")
    except Exception as e:
        st.error(f"Error during optimization: {e}")
