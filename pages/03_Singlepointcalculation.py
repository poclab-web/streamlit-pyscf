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
from utils.pyscf_ui import (
    display_pyscf_status, 
    require_pyscf, 
    display_basis_selector, 
    display_functional_selector,
    display_calculation_options
)
from logic.molecule_handler import MoleculeHandler
from logic.calculation import theory_options, basis_set_options, run_quantum_calculation
from logic.output_handler import extract_orbital_energies

# カスタムCSSを適用
load_css("config/styles.css")

# PySCFのステータス確認
st.title("Single Point Calculation")

# PySCFのインストール状況を確認
with st.expander("🔬 PySCF Status & Configuration", expanded=False):
    pyscf_status = display_pyscf_status(show_config_section=True)

# PySCFが利用できない場合は処理を停止
if not pyscf_status.get("pyscf_available", False):
    st.error("⚠️ PySCFが利用できません。上記のインストール手順に従ってPySCFをインストールしてください。")
    st.stop()

# 本文

# ユーザー入力
st.header("📋 Molecular Input")
input_type = st.selectbox("Select Input Type", ["XYZ", "SMILES"])
atom_input = st.text_area(
    "Enter Molecular Structure",
    "C     -2.385271    1.648249    0.000003\nH     -3.181824    1.061240    0.407216\nH     -2.536116    1.773363   -1.051896\nH     -2.370100    2.607079    0.474676\nO     -1.453042    1.151314    0.170017\nH     -1.466511    0.300016   -0.251422"
    if input_type == "XYZ"
    else "CO",
)

# 計算設定セクション
st.header("⚙️ Calculation Settings")

col1, col2 = st.columns(2)

with col1:
    st.subheader("🧮 Theory Level")
    # 従来の選択肢と新しいUIを併用
    use_advanced_ui = st.checkbox("使用高级PySCF设置", value=False, help="PySCF UIの高度な設定を使用")
    
    if use_advanced_ui:
        # PySCF UIからの設定
        theory = st.selectbox(
            "Theory Method", 
            ["HF", "DFT"], 
            help="計算に使用する理論手法を選択"
        )
        if theory == "DFT":
            functional = display_functional_selector()
            theory = functional  # DFTの場合は汎関数名を使用
        basis_set = display_basis_selector()
    else:
        # 従来の設定
        theory = st.selectbox("Theory", theory_options)
        basis_set = st.selectbox("Basis Set", basis_set_options)

with col2:
    st.subheader("🔧 Calculation Options")
    if use_advanced_ui:
        # PySCF UIからの詳細設定
        calc_options = display_calculation_options()
        charge = calc_options["charge"]
        multiplicity = calc_options["multiplicity"]
        spin = multiplicity - 1
        symmetry = calc_options["symmetry"]
    else:
        # 従来の設定（expanderから移動）
        charge = st.number_input("Molecular Charge", min_value=-10, max_value=10, value=0, step=1)
        multiplicity = st.number_input("Spin Multiplicity (2S + 1)", min_value=1, max_value=10, value=1, step=1)
        spin = multiplicity - 1
        symmetry_choice = st.selectbox("Consider Molecular Symmetry?", ["Yes", "No"])
        symmetry = True if symmetry_choice == "Yes" else False

# 溶媒効果の設定
st.header("🌊 Solvent Effects")
with st.expander("溶媒効果設定", expanded=False):
    solvent_model = st.selectbox("Select Solvent Model", ["None", "PCM", "DDCOSMO"])
    eps = None  # epsのデフォルト値を設定

    # Load solvent data
    solvents_file = "config/solvents_epsilon.csv"
    if os.path.exists(solvents_file):
        solvents_data = pd.read_csv(solvents_file)

        if solvent_model in ["PCM", "DDCOSMO"]:
            # Show solvent selection dropdown
            solvent_selection = st.selectbox(
                "Select a solvent",
                [f"{row['Solvent']} (ε={row['Epsilon']})" for _, row in solvents_data.iterrows()]
            )  
            # Extract epsilon from selection
            if solvent_selection:
                eps = float(solvent_selection.split("=", 1)[-1][:-1])
    else:
        st.warning("溶媒データファイルが見つかりません。")
        if solvent_model in ["PCM", "DDCOSMO"]:
            eps = st.number_input("誘電率 (ε)", value=78.39, min_value=1.0, step=0.1)

# 分子構造を処理
st.header("🧪 Calculation Execution")

# 計算設定の概要表示
with st.expander("📋 計算設定確認", expanded=False):
    st.write("**理論手法:**", theory)
    st.write("**基底関数系:**", basis_set)
    st.write("**電荷:**", charge)
    st.write("**多重度:**", multiplicity)
    st.write("**対称性:**", "使用" if symmetry else "不使用")
    st.write("**溶媒モデル:**", solvent_model)
    if eps:
        st.write("**誘電率:**", eps)

handler = None
if st.button("🚀 Run Single Point Calculation", type="primary"):
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
            st.subheader("📋 Input 2D Structure")
            handler.generate_2d_image(f"{directory}/molecule_2d.png")
            st.image(f"{directory}/molecule_2d.png", caption=Chem.MolToSmiles(handler.mol))

        # Display 3D structure in the second column
        with col2:
            st.subheader("🧬 Input 3D Structure")
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

    import matplotlib.pyplot as plt
    import numpy as np
    import streamlit as st

    # SCFデータ（例） TODO: 表示できるようにする
    # cycles = list(range(1, 10))
    # delta_E = [1.29, -0.0334, -0.00713, -0.000359, -1.74e-05, -1.67e-06, -1.53e-07, -1.48e-08, -7.74e-10]
    # g_norm = [0.471, 0.221, 0.0383, 0.00768, 0.00187, 0.000566, 0.000147, 3.14e-05, 5.83e-06]
    # ddm = [1.57, 0.329, 0.115, 0.034, 0.00773, 0.00296, 0.000862, 0.000329, 7.51e-05]

    # # 絶対値（logスケール対応のため）
    # delta_E_abs = np.abs(delta_E)

    # # 収束ライン
    # g_conv_threshold = 3.16228e-05

    # # プロット
    # fig, ax = plt.subplots()
    # ax.plot(cycles, delta_E_abs, marker='o', label="ΔE (delta_E)")
    # ax.plot(cycles, g_norm, marker='o', label="|g| (Fock gradient)")
    # ax.plot(cycles, ddm, marker='o', label="|ddm| (Change in density matrix)")

    # # 収束ラインの追加（|g|用）
    # ax.axhline(g_conv_threshold, color='red', linestyle='--', label='|g| threshold = 3.16e-5')

    # # 軸・ラベル
    # ax.set_yscale("log")
    # ax.set_xlabel("SCF Cycle")
    # ax.set_ylabel("Value (log scale)")
    # ax.set_title("SCF convergence index trend")
    # ax.grid(True)
    # ax.legend()

    # # Streamlitで表示
    # st.pyplot(fig)


# 計算の実行
    try:
        if not handler or not handler.mol:
            raise ValueError("Please process the molecule before running the calculation.")

        # Generate PySCF input format
        pyscf_input = handler.to_pyscf_input()

        st.write("Running quantum chemistry calculation...")
        result = run_quantum_calculation(
            compound_name, smiles, pyscf_input, basis_set, theory, charge=charge, spin=spin, solvent_model=solvent_model, eps=eps, symmetry=symmetry
        )
        energy = result["energy"]
        molden_file = result["molden_file"]

        # --- エネルギー分解の表示 ---
        st.subheader("📊 Energy Decomposition")
        energy_col1, energy_col2 = st.columns(2)
        
        with energy_col1:
            st.metric("核間反発エネルギー (E_nuc)", f"{result['E_nuc']:.6f} Hartree")
            st.metric("電子-核引力項 (E_core)", f"{result['E_core']:.6f} Hartree")
            st.metric("クーロン項 (E_J)", f"{result['E_J']:.6f} Hartree")
        
        with energy_col2:
            st.metric("交換項 (E_K)", f"{result['E_K']:.6f} Hartree")
            st.metric("電子エネルギー (E_elec)", f"{result['E_elec']:.6f} Hartree")
            st.metric("**総エネルギー**", f"**{energy:.6f} Hartree**")

        # 計算結果を表示
        st.success(f"✅ Calculated SinglePoint Energy: {energy} Hartree")
        st.info(f"📁 Results saved in: {molden_file}")

        # 軌道エネルギーとHOMO/LUMOを表示

        if os.path.exists(molden_file):
            orbital_data = extract_orbital_energies(molden_file)

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

            # ソート後のHOMOとLUMOのインデックスを再計算
            sorted_homo_index = orbital_df[orbital_df["index"] == homo_index].index[0]
            sorted_lumo_index = orbital_df[orbital_df["index"] == lumo_index].index[0] if lumo_index is not None else None

            # 軌道エネルギーを表示（Energy (Hartree) のみ）
            st.subheader("⚛️ Orbital Energies")
            
            # 軌道エネルギーテーブルの表示
            with st.expander("全軌道エネルギー表", expanded=False):
                st.dataframe(orbital_df[["Energy (Hartree)"]], use_container_width=True)

            # HOMOとLUMOの情報を表示
            from logic.output_handler import convert_energy_units  # 必要な関数をインポート

            st.subheader("🎯 HOMO and LUMO")
            homo_lumo_col1, homo_lumo_col2 = st.columns(2)
            
            with homo_lumo_col1:
                if sorted_homo_index is not None:
                    homo_energy_hartree = orbital_df.loc[sorted_homo_index, 'Energy (Hartree)']
                    homo_energy_ev = convert_energy_units(homo_energy_hartree, unit="eV")
                    st.metric(
                        "HOMO Energy", 
                        f"{homo_energy_hartree:.6f} Hartree",
                        delta=f"{homo_energy_ev:.2f} eV"
                    )
                    st.write(f"**HOMO Index:** {sorted_homo_index + 1}")
                else:
                    st.error("HOMO: Not found")
            
            with homo_lumo_col2:
                if sorted_lumo_index is not None:
                    lumo_energy_hartree = orbital_df.loc[sorted_lumo_index, 'Energy (Hartree)']
                    lumo_energy_ev = convert_energy_units(lumo_energy_hartree, unit="eV")
                    st.metric(
                        "LUMO Energy", 
                        f"{lumo_energy_hartree:.6f} Hartree",
                        delta=f"{lumo_energy_ev:.2f} eV"
                    )
                    st.write(f"**LUMO Index:** {sorted_lumo_index + 1}")
                    
                    # HOMO-LUMOギャップを計算
                    if sorted_homo_index is not None:
                        gap_hartree = lumo_energy_hartree - homo_energy_hartree
                        gap_ev = convert_energy_units(gap_hartree, unit="eV")
                        st.metric(
                            "HOMO-LUMO Gap",
                            f"{gap_hartree:.6f} Hartree",
                            delta=f"{gap_ev:.2f} eV"
                        )
                else:
                    st.warning("LUMO: Not found")

        else:
            st.warning("MOLDEN file not found. Orbital energies cannot be displayed.")
    except Exception as e:
        st.error(f"Error: {e}")
