"""
個々のフラグメント（分子 A, 分子 B）のエネルギーを計算。
フラグメントが結合した複合体（A + B）のエネルギーを計算。
相互作用エネルギーを計算：
"""

import streamlit as st
from utils.module import load_css

import os
import numpy as np

import pandas as pd
from rdkit import Chem

import py3Dmol  # 3D可視化用ライブラリ
import stmol

from config.config import conv_preset_values, solvent_models, solvents_data, solvents_file

from controllers.energydecompositionanalysis import run_ced
from logic.calculation import theory_options, basis_set_options, run_quantum_calculation

from logic.molecule_handler import MoleculeHandler
import multiprocessing

# カスタムCSSを適用
load_css("config/styles.css")

# 本文
st.title("Conformational Energy Decomposition")

input1, input2 = st.columns([1, 1], gap="large")

# ユーザー入力
with input1:
    st.header("Molecular Input1")
    input_type1 = st.selectbox("Select Input Type", ["XYZ", "SMILES"], key="input_type1")
    atom_input1 = st.text_area(
        "Enter Molecular Structure",
        "O     0.000000     0.000000     0.000000\nH     0.756950     0.585880     0.000000\nH    -0.756950     0.585880     0.000000"
        if input_type1 == "XYZ"
        else "O",
        key="atom_input1"
    )
    apply_force_field1 = st.checkbox("分子力場による構造最適化を行う", value=False, key="apply_force_field1")
    
    st.header("🔧 配座生成と計算設定")
    force_field1 = st.selectbox("Select Force Field", ["MMFF", "UFF"], key="force_field1")
    num_conformers1 = st.number_input("Number of Conformers", value=1000, key="num_conformers1")

with input2:
    st.header("Molecular Input2")
    input_type2 = st.selectbox("Select Input Type", ["XYZ", "SMILES"], key="input_type2")
    atom_input2 = st.text_area(
        "Enter Molecular Structure",
        "O     0.000000     0.000000     0.000000\nH     0.958000     0.000000     0.000000\nH    -0.958000     0.000000     0.000000"
        if input_type2 == "XYZ"
        else "O",
        key="atom_input2"
    )
    apply_force_field2 = st.checkbox("分子力場による構造最適化を行う", value=False, key="apply_force_field2")
    
    st.header("🔧 配座生成と計算設定")
    force_field2 = st.selectbox("Select Force Field", ["MMFF", "UFF"], key="force_field2")
    num_conformers2 = st.number_input("Number of Conformers", value=1000, key="num_conformers2")


# 分子ハンドラーの初期化
if 'handler1' not in st.session_state:
    st.session_state.handler1 = None
if 'handler2' not in st.session_state:
    st.session_state.handler2 = None
if 'compound_name1' not in st.session_state:
    st.session_state.compound_name1 = None
if 'compound_name2' not in st.session_state:
    st.session_state.compound_name2 = None
if 'smiles1' not in st.session_state:
    st.session_state.smiles1 = None
if 'smiles2' not in st.session_state:
    st.session_state.smiles2 = None
if 'conformer_success1' not in st.session_state:
    st.session_state.conformer_success1 = False
if 'conformer_success2' not in st.session_state:
    st.session_state.conformer_success2 = False

if st.button("Check Structure"):
    try:
        # input1
        st.session_state.handler1 = MoleculeHandler(atom_input1, input_type=input_type1.lower())
        if not st.session_state.handler1.mol:
            raise ValueError("Invalid molecular input for molecule 1. Please check your format.")

        # 化合物名を取得
        st.session_state.compound_name1 = Chem.MolToInchiKey(st.session_state.handler1.mol)
        st.session_state.smiles1 = Chem.MolToSmiles(st.session_state.handler1.mol)

        # ディレクトリの作成
        directory = os.path.join("data", st.session_state.compound_name1)
        os.makedirs(directory, exist_ok=True)

        # 分子1の配座探索実行
        st.session_state.conformer_success1 = True
        if apply_force_field1 and st.session_state.handler1.mol.GetNumAtoms() > 2:
            try:
                with st.spinner(f"分子1の配座探索を実行中... ({force_field1} force field)"):
                    st.session_state.handler1.generate_conformers(num_conformers=num_conformers1, forcefield=force_field1)
                    st.session_state.handler1.keep_lowest_energy_conformer()
                    st.success(f"分子1: {num_conformers1}個の配座から最安定配座を選択しました")
            except Exception as e:
                st.warning(f"分子1の配座探索でエラーが発生しました: {e}")
                st.info("分子1: 初期構造を使用します")
                st.session_state.conformer_success1 = False
        else:
            st.session_state.conformer_success1 = not apply_force_field1  # チェックボックスがOFFなら成功とみなす

        # input2
        st.session_state.handler2 = MoleculeHandler(atom_input2, input_type=input_type2.lower())
        if not st.session_state.handler2.mol:
            raise ValueError("Invalid molecular input for molecule 2. Please check your format.")
        # 化合物名を取得
        st.session_state.compound_name2 = Chem.MolToInchiKey(st.session_state.handler2.mol)
        st.session_state.smiles2 = Chem.MolToSmiles(st.session_state.handler2.mol)
        # ディレクトリの作成
        directory2 = os.path.join("data", st.session_state.compound_name2)
        os.makedirs(directory2, exist_ok=True)

        # 分子2の配座探索実行
        st.session_state.conformer_success2 = True
        if apply_force_field2 and st.session_state.handler2.mol.GetNumAtoms() > 2:
            try:
                with st.spinner(f"分子2の配座探索を実行中... ({force_field2} force field)"):
                    st.session_state.handler2.generate_conformers(num_conformers=num_conformers2, forcefield=force_field2)
                    st.session_state.handler2.keep_lowest_energy_conformer()
                    st.success(f"分子2: {num_conformers2}個の配座から最安定配座を選択しました")
            except Exception as e:
                st.warning(f"分子2の配座探索でエラーが発生しました: {e}")
                st.info("分子2: 初期構造を使用します")
                st.session_state.conformer_success2 = False
        else:
            st.session_state.conformer_success2 = not apply_force_field2  # チェックボックスがOFFなら成功とみなす

        col1, col2 = st.columns([1, 1], gap="large")

        # Display 3D structure in the first column
        with col1:
            if apply_force_field1 and st.session_state.handler1.mol.GetNumAtoms() > 2 and st.session_state.conformer_success1:
                structure_title = "配座最適化後の3D構造"
            else:
                structure_title = "初期3D構造"
            st.subheader(structure_title)
            try:
                mol_block1 = st.session_state.handler1.generate_3d_molblock()
                viewer = py3Dmol.view(width=400, height=400)  # Adjust width to fit in the column
                viewer.addModel(mol_block1, "mol")
                viewer.setStyle({"stick": {}})
                viewer.zoomTo()
                stmol.showmol(viewer, height=400)
                st.write("inchikey:", st.session_state.compound_name1)
                if apply_force_field1 and st.session_state.handler1.mol.GetNumAtoms() > 2 and st.session_state.conformer_success1:
                    st.info(f"✅ {force_field1}力場による配座最適化済み")
                elif apply_force_field1 and st.session_state.handler1.mol.GetNumAtoms() > 2 and not st.session_state.conformer_success1:
                    st.warning(f"⚠️ {force_field1}力場による配座最適化が失敗したため、初期構造を使用")

            except Exception as e:
                st.warning(f"Unable to generate 3D structure: {e}")

        # Display 3D structure in the second column
        with col2:
            if apply_force_field2 and st.session_state.handler2.mol.GetNumAtoms() > 2 and st.session_state.conformer_success2:
                structure_title = "配座最適化後の3D構造"
            else:
                structure_title = "初期3D構造"
            st.subheader(structure_title)
            try:
                mol_block2 = st.session_state.handler2.generate_3d_molblock()
                viewer = py3Dmol.view(width=400, height=400)  # Adjust width to fit in the column
                viewer.addModel(mol_block2, "mol")
                viewer.setStyle({"stick": {}})
                viewer.zoomTo()
                stmol.showmol(viewer, height=400)
                st.write("inchikey:", st.session_state.compound_name2)
                if apply_force_field2 and st.session_state.handler2.mol.GetNumAtoms() > 2 and st.session_state.conformer_success2:
                    st.info(f"✅ {force_field2}力場による配座最適化済み")
                elif apply_force_field2 and st.session_state.handler2.mol.GetNumAtoms() > 2 and not st.session_state.conformer_success2:
                    st.warning(f"⚠️ {force_field2}力場による配座最適化が失敗したため、初期構造を使用")
            except Exception as e:
                st.warning(f"Unable to generate 3D structure: {e}")

        # compound_name1とcompound_name2が一致していることを確認
        if st.session_state.compound_name1 == st.session_state.compound_name2:
            # 同じ構造が入力されています。
            st.info("The same structure has been input for both molecules.")
        else:
            # 異なる構造が入力されています。
            st.error("Different structures have been input for the two molecules.")

    except Exception as e:
        st.error(f"Error processing molecule: {e}")

# 並列計算の有無を選択
st.header("Parallel Calculation Options")
logical_cores = multiprocessing.cpu_count()
try:
    physical_cores = multiprocessing.cpu_count(logical=False)
except TypeError:
    physical_cores = logical_cores // 2  # Fallback

st.write(f"使用しているパソコンの論理コア数: {logical_cores}")
st.write(f"使用しているパソコンの物理コア数: {physical_cores}")

# 物理コアが2以上なら並列計算を推奨
recommend_parallel = physical_cores >= 2

parallel_option = st.checkbox(
    "並列計算を有効にする（推奨）" if recommend_parallel else "並列計算を有効にする",
    value=recommend_parallel,
    key="parallel_option_checkbox"
)

if recommend_parallel:
    st.info("物理コア数が2以上のため、並列計算が推奨されます。")
else:
    st.warning("物理コア数が少ないため、並列計算は非推奨です。")

# 理論レベルと基底関数の選択
theory = st.selectbox("Theory", theory_options, key="theory_selectbox")
basis_set = st.selectbox("Basis Set", basis_set_options, key="basis_set_selectbox")

# Charge and Spin入力セクション
with st.expander("Other Settings"):
  charge = st.number_input("Molecular Charge", min_value=-10, max_value=10, value=0, step=1, key="charge_input")
  multiplicity = st.number_input("Spin Multiplicity (2S + 1)", min_value=1, max_value=10, value=1, step=1, key="multiplicity_input")
  spin = multiplicity - 1

  symmetry = st.selectbox("Consider Molecular Symmetry?", ["Yes", "No"], key="symmetry_selectbox")
  symmetry = True if symmetry == "Yes" else False

  # 溶媒効果の設定
  solvent_model = st.selectbox("Select Solvent Model", ["None", "PCM", "DDCOSMO"], key="solvent_model_selectbox")
  eps = None  # epsのデフォルト値を設定

  # Load solvent data
  solvents_file = "config/solvents_epsilon.csv"
  solvents_data = pd.read_csv(solvents_file)

  if solvent_model in ["PCM", "DDCOSMO"]:
      # Show solvent selection dropdown
      solvent_selection = st.selectbox(
          "Select a solvent",
          [f"{row['Solvent']} (ε={row['Epsilon']})" for _, row in solvents_data.iterrows()],
          key="solvent_selection_selectbox"
      )  
      # Extract epsilon from selection
      if solvent_selection:
        eps = float(solvent_selection.split("=", 1)[-1][:-1])

# 計算の設定入力
with st.expander("Setting calculation conditions for optimization"):
    st.header("🔍 量子化学計算による構造最適化")
    optimize_with_qc = st.checkbox("量子化学計算による構造最適化を行う", value=True, key="optimize_with_qc_checkbox")
    st.subheader("Convergence Parameters")
    preset = st.radio("Choose preset", ["Loose", "Normal", "Tight"], index=1, horizontal=True, key="preset_radio")
    vals = conv_preset_values[preset]
    convergence_energy = st.number_input("Energy Tolerance (Hartree)", min_value=1e-7, value=vals["energy"], step=1e-5, format="%.7f", key="convergence_energy_input")
    convergence_grms = st.number_input("Gradient RMS Tolerance (Eh/Bohr)", min_value=1e-5, value=vals["grms"], step=1e-5, format="%.5f", key="convergence_grms_input")
    convergence_gmax = st.number_input("Gradient Max Tolerance (Eh/Bohr)", min_value=1e-5, value=vals["gmax"], step=1e-5, format="%.5f", key="convergence_gmax_input")
    convergence_drms = st.number_input("Displacement RMS Tolerance (Angstrom)", min_value=1e-4, value=vals["drms"], step=1e-4, format="%.4f", key="convergence_drms_input")
    convergence_dmax = st.number_input("Displacement Max Tolerance (Angstrom)", min_value=1e-4, value=vals["dmax"], step=1e-4, format="%.4f", key="convergence_dmax_input")
    maxsteps = st.number_input("Max Iterations", min_value=1, value=100, step=1, key="maxsteps_input")
    conv_params = {
        "convergence_energy": convergence_energy,
        "convergence_grms": convergence_grms,
        "convergence_gmax": convergence_gmax,
        "convergence_drms": convergence_drms,
        "convergence_dmax": convergence_dmax,
        "maxsteps": maxsteps,
    }

# 計算の方法を出力
with st.expander("計算方法と参考文献を表示", expanded=False):
    st.markdown("### 🧪 Method for Conformational Energy Decomposition")
    st.markdown(
        "**Computational Details**  \n"
        "Molecular structures were processed using RDKit [1] for initial 3D coordinate generation.  \n"
        f"Single-point energy calculations were performed at the **{theory}/{basis_set}** level using PySCF [2].  \n"
        f"{'No solvent model was applied.' if solvent_model == 'None' else f'The solvent effect was considered using the {solvent_model} model with ε = {eps}.'}  \n"
        "For energy decomposition analysis, the total electronic energy is decomposed into:  \n"
        "- Nuclear repulsion energy (E_nuc)  \n"
        "- Core Hamiltonian energy (E_core): electron-nuclear attraction and kinetic energy  \n"
        "- Coulomb interaction energy (E_J): electron-electron repulsion  \n"
        "- Exchange energy (E_K): quantum mechanical exchange interaction  \n"
        "Energy differences between the two molecular structures are calculated and converted to kcal/mol for comparative analysis.  \n"
        "This decomposition provides insights into the physical origins of energy differences between different molecular structures or conformations."
    )
    st.markdown("---")
    st.markdown(
        "**References**  \n"
        "[1] Landrum, G. RDKit: Open-source cheminformatics. [https://www.rdkit.org](https://www.rdkit.org)  \n"
        "[2] Sun, Q. *et al.* PySCF: The Python-based Simulations of Chemistry Framework. **WIREs Comput Mol Sci** *2018*, **8**, e1340. DOI: [10.1002/wcms.1340](https://doi.org/10.1002/wcms.1340)  \n"
        "[3] PocLab streamlit-pyscf: Quantum chemistry web interface. [https://github.com/poclab-web/streamlit-pyscf](https://github.com/poclab-web/streamlit-pyscf)  \n"
        "[4] Szabo, A.; Ostlund, N. S. *Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory*; Dover Publications: New York, 1996.  \n"
        "[5] Helgaker, T.; Jørgensen, P.; Olsen, J. *Molecular Electronic-Structure Theory*; Wiley: Chichester, 2000."
    )


if st.button("Run Conformational Energy Decomposition"):
    try:
        # セッション状態の確認
        if st.session_state.handler1 is None or st.session_state.handler2 is None:
            st.error("まず'Check Structure'ボタンを押して分子を初期化してください。")
            st.stop()
        
        # PySCF入力の準備（配座探索は既にCheck Structureで実行済み）
        pyscf_input1 = st.session_state.handler1.to_pyscf_input()
        pyscf_input2 = st.session_state.handler2.to_pyscf_input()

        # 結果の表示
        st.success("Conformational Energy Decomposition の準備が完了しました。")
        st.write("入力データ:")
        st.write(f"分子1: {st.session_state.compound_name1}")
        st.write(f"分子2: {st.session_state.compound_name2}")
        st.write(f"理論レベル: {theory}")
        st.write(f"基底関数: {basis_set}")
        
        # 配座探索の状態を表示
        if apply_force_field1 and st.session_state.handler1.mol.GetNumAtoms() > 2 and st.session_state.conformer_success1:
            st.info(f"✅ 分子1: {force_field1}力場による配座最適化済み")
        elif apply_force_field1 and st.session_state.handler1.mol.GetNumAtoms() > 2 and not st.session_state.conformer_success1:
            st.warning(f"⚠️ 分子1: {force_field1}力場による配座最適化が失敗、初期構造を使用")
        else:
            st.info("ℹ️ 分子1: 初期構造を使用")
            
        if apply_force_field2 and st.session_state.handler2.mol.GetNumAtoms() > 2 and st.session_state.conformer_success2:
            st.info(f"✅ 分子2: {force_field2}力場による配座最適化済み")
        elif apply_force_field2 and st.session_state.handler2.mol.GetNumAtoms() > 2 and not st.session_state.conformer_success2:
            st.warning(f"⚠️ 分子2: {force_field2}力場による配座最適化が失敗、初期構造を使用")
        else:
            st.info("ℹ️ 分子2: 初期構造を使用")
        
        # 計算開始の表示
        with st.spinner("量子化学計算を実行中..."):
            # run_ced関数を実行してConformational Energy Decomposition計算を実行
            result = run_ced(parallel_option, st.session_state.compound_name1, st.session_state.compound_name2, 
                             pyscf_input1, pyscf_input2, st.session_state.smiles1, st.session_state.smiles2, 
                             basis_set, theory, charge, spin, opt_theory=None, opt_basis_set=None, 
                             solvent_model=solvent_model, eps=eps, symmetry=symmetry)
        
        st.success("Conformational Energy Decomposition 計算が完了しました。")
        st.subheader("🔬 エネルギー分解結果")
        
        # 安定性の比較を表示
        if 'energy_A' in result and 'energy_B' in result:
            energy_diff = result['energy_A'] - result['energy_B']
            energy_diff_kcal = energy_diff * 627.509
            
            st.subheader("⚖️ 分子の安定性比較")
            col_stability1, col_stability2 = st.columns(2)
            
            with col_stability1:
                st.metric(
                    label="分子1のエネルギー",
                    value=f"{result['energy_A']:.6f} Ha",
                    delta=f"{energy_diff_kcal:.2f} kcal/mol" if energy_diff != 0 else None
                )
            
            with col_stability2:
                st.metric(
                    label="分子2のエネルギー",
                    value=f"{result['energy_B']:.6f} Ha",
                    delta=f"{-energy_diff_kcal:.2f} kcal/mol" if energy_diff != 0 else None
                )
            
            # どちらが安定かを判定
            if abs(energy_diff_kcal) > 0.1:  # 0.1 kcal/mol以上の差がある場合
                if energy_diff < 0:
                    st.success("🏆 **分子1がより安定** (エネルギーが低い)")
                    st.write(f"分子1は分子2より {abs(energy_diff_kcal):.2f} kcal/mol 安定です")
                else:
                    st.success("🏆 **分子2がより安定** (エネルギーが低い)")
                    st.write(f"分子2は分子1より {abs(energy_diff_kcal):.2f} kcal/mol 安定です")
            else:
                st.info("⚖️ **両分子のエネルギーはほぼ同等**")
                st.write(f"エネルギー差: {abs(energy_diff_kcal):.2f} kcal/mol (微小)")
        
        # エネルギー差分をHartree単位で表示
        energy_keys = ['energy', 'E_nuc', 'E_core', 'E_J', 'E_K', 'E_elec']
        
        st.subheader("📊 エネルギー差分分析 (分子1 - 分子2)")
        col1, col2 = st.columns(2)
        
        with col1:
            st.write("**Hartree単位での差分:**")
            for key in energy_keys:
                if key in result:
                    st.write(f"Δ{key}: {result[key]:.6f} Ha")
        
        with col2:
            st.write("**kcal/mol単位での差分:**")
            for key in energy_keys:
                kcal_key = f"{key}_kcal_mol"
                if kcal_key in result:
                    st.write(f"Δ{key}: {result[kcal_key]:.2f} kcal/mol")
        
        # エネルギー成分の寄与分析
        if all(f"{key}_kcal_mol" in result for key in ['E_nuc', 'E_core', 'E_J', 'E_K']):
            st.subheader("🧮 エネルギー成分の寄与分析")
            
            # 各成分の寄与を可視化
            components = ['E_nuc', 'E_core', 'E_J', 'E_K']
            values = [result[f"{key}_kcal_mol"] for key in components]
            labels = ['核間反発', 'コア(電子-核)', 'クーロン', '交換']
            
            # データフレームとして表示
            import pandas as pd
            df_components = pd.DataFrame({
                'エネルギー成分': labels,
                'エネルギー差 (kcal/mol)': values
            })
            st.dataframe(df_components, use_container_width=True)
        
        # 収束情報の表示
        if 'converged_A' in result and 'converged_B' in result:
            st.subheader("📊 収束情報")
            st.write(f"分子1の収束: {'✅ 収束' if result['converged_A'] else '❌ 非収束'}")
            st.write(f"分子2の収束: {'✅ 収束' if result['converged_B'] else '❌ 非収束'}")
        
        # 全結果の詳細表示（デバッグ用）
        with st.expander("🔍 詳細な計算結果"):
            st.json(result)


    except Exception as e:
        st.error(f"エラーが発生しました: {e}")
        st.stop()
