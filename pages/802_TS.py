"""
PySCFを使用して遷移状態の探索を実行する
ユーザーは反応物と生成物の構造、理論手法、基底関数系を入力し、遷移状態探索計算を実行して
結果を可視化できます。

機能:
- XYZ形式で分子構造を入力可能
- 遷移状態探索の理論手法と基底関数系を選択可能
- QSD (Quadratic Steepest Descent) 最適化手法を使用
- 初期構造から遷移状態への最適化過程を可視化
- 最終的な遷移状態構造を表示し、振動解析結果も表示
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import streamlit as st
import pandas as pd
from rdkit import Chem
import py3Dmol
import stmol

from utils.module import load_css
from logic.molecule_handler import MoleculeHandler
from logic.calculation import theory_options, basis_set_options, run_quantum_calculation, normalize_basis_set, run_ts_search, perform_frequency_analysis

# カスタムCSSを適用
load_css("config/styles.css")

# 本文
st.title("Transition State Search")

st.warning("精度は改善余地あります。検討中です")

st.markdown("""
### 遷移状態探索について
遷移状態は化学反応において反応物から生成物へと変化する過程での最高エネルギー点です。
この計算では QSD (Quadratic Steepest Descent) 手法を使用して遷移状態を探索します。

**重要な注意事項:**
- 初期構造は反応座標に沿った適切な構造を提供する必要があります
- 遷移状態探索には良い初期推定構造が必要です
- 計算後の振動解析で虚振動が1つだけあることを確認してください
""")

# ユーザー入力
st.header("Molecular Input for TS Search")

# 計算レベルの推奨設定を表示
st.info("🎯 **遷移状態探索推奨設定**: DFT-B3LYP/6-31G(d) 以上を強く推奨します")

input_type = st.selectbox("Select Input Type", ["XYZ", "SMILES"])
atom_input = st.text_area(
    "Enter Molecular Structure (Initial guess for TS)",
    "C     0.000000    0.000000    0.000000\nH     0.000000    0.000000    1.100000\nH     0.000000    1.100000   -0.300000\nH     1.100000   -0.300000   -0.300000\nH    -1.100000   -0.300000   -0.300000"
    if input_type == "XYZ"
    else "C",
    help="遷移状態の初期推定構造を入力してください"
)

theory = st.selectbox("Theory", theory_options, index=1 if "DFT-B3LYP" in theory_options else 0)  # DFT-B3LYPをデフォルトに
basis_set = st.selectbox("Basis Set", basis_set_options, index=2 if "6-31g(d)" in [bs.lower() for bs in basis_set_options] else 0)  # 6-31G(d)をデフォルトに

# 選択された設定の評価
theory_eval = theory.lower().replace("-", "").replace("_", "")
basis_eval = basis_set.lower().replace("-", "").replace("*", "(d)").replace("**", "(d,p)")

if theory_eval == "hf" or basis_eval == "sto3g":
    st.warning("⚠️ 選択された設定は遷移状態探索には適していません。上記の推奨設定をご検討ください。")
elif theory_eval in ["b3lyp", "pbe", "pbe0", "m06"] and basis_eval in ["631g(d)", "6311g(d,p)", "631+g(d)", "6311+g(d,p)"]:
    st.success("✅ 良い選択です！この設定は遷移状態探索に適しています。")
elif theory_eval in ["b3lyp", "pbe", "pbe0", "m06"]:
    st.info("💡 良い理論手法です。より大きな基底関数系（6-31G(d)以上）を使用するとさらに良いでしょう。")

# Charge and Spin入力セクション
with st.expander("Other Settings"):
    charge = st.number_input("Molecular Charge", min_value=-10, max_value=10, value=0, step=1)
    multiplicity = st.number_input("Spin Multiplicity (2S + 1)", min_value=1, max_value=10, value=1, step=1)
    spin = multiplicity - 1

    symmetry = st.selectbox("Consider Molecular Symmetry?", ["Yes", "No"], index=1)
    symmetry = True if symmetry == "Yes" else False

    # 溶媒効果の設定
    solvent_model = st.selectbox("Select Solvent Model", ["None", "PCM", "DDCOSMO"])
    eps = None

    # Load solvent data
    solvents_file = "config/solvents_epsilon.csv"
    if os.path.exists(solvents_file):
        solvents_data = pd.read_csv(solvents_file)
        
        if solvent_model in ["PCM", "DDCOSMO"]:
            solvent_selection = st.selectbox(
                "Select a solvent",
                [f"{row['Solvent']} (ε={row['Epsilon']})" for _, row in solvents_data.iterrows()]
            )  
            if solvent_selection:
                eps = float(solvent_selection.split("=", 1)[-1][:-1])

# 収束条件の入力
with st.expander("TS Search Convergence Parameters"):
    convergence_energy = st.number_input(
        "Energy Tolerance (Hartree)", 
        min_value=1e-7, value=1.0e-6, step=1e-7, format="%.7f"
    )
    convergence_grms = st.number_input(
        "Gradient RMS Tolerance (Eh/Bohr)", 
        min_value=1e-6, value=1.0e-4, step=1e-6, format="%.6f"
    )
    convergence_gmax = st.number_input(
        "Gradient Max Tolerance (Eh/Bohr)", 
        min_value=1e-6, value=3.0e-4, step=1e-6, format="%.6f"
    )
    convergence_drms = st.number_input(
        "Displacement RMS Tolerance (Angstrom)", 
        min_value=1e-5, value=1.0e-4, step=1e-5, format="%.5f"
    )
    convergence_dmax = st.number_input(
        "Displacement Max Tolerance (Angstrom)", 
        min_value=1e-5, value=3.0e-4, step=1e-5, format="%.5f"
    )
    maxsteps = st.number_input(
        "Max Iterations", 
        min_value=1, value=50, step=1
    )

conv_params = {
    "convergence_energy": convergence_energy,
    "convergence_grms": convergence_grms,
    "convergence_gmax": convergence_gmax,
    "convergence_drms": convergence_drms,
    "convergence_dmax": convergence_dmax,
}

# オプション設定
with st.expander("Additional Options"):
    perform_freq = st.checkbox("Perform frequency analysis after TS search", value=True)
    save_trajectory = st.checkbox("Save optimization trajectory", value=True)

# 遷移状態探索の実行
if st.button("Run Transition State Search"):
    # 計算設定の妥当性チェック
    theory_check = theory.lower().replace("-", "").replace("_", "")
    basis_check = basis_set.lower().replace("-", "").replace("*", "(d)").replace("**", "(d,p)")
    
    if theory_check == "hf" and basis_check == "sto3g":
        st.warning("⚠️ **計算設定の警告**: HF/STO-3Gは遷移状態探索には適していません。")
        st.info("""
        **推奨設定:**
        - 理論手法: B3LYP または PBE
        - 基底関数: 6-31G(d) 以上
        - より精密な基底関数系が遷移状態の正確な特性評価に必要です
        """)
    
    if basis_check == "sto3g":
        st.warning("⚠️ STO-3G基底関数は最小基底で、遷移状態の精密な計算には不十分です。")
    
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
            st.image(f"{directory}/molecule_2d.png", caption=smiles)

        # 3D構造の生成
        with col2:
            st.subheader("Initial 3D Structure")
            try:
                mol_block = handler.generate_3d_molblock()
                viewer = py3Dmol.view(width=400, height=400)
                viewer.addModel(mol_block, "mol")
                viewer.setStyle({"stick": {}})
                viewer.zoomTo()
                stmol.showmol(viewer, height=400)
            except Exception as e:
                st.warning(f"Unable to generate 3D structure: {e}")

        # PySCF用入力データを作成
        pyscf_input = handler.to_pyscf_input()

        # 遷移状態探索の実行
        st.write("Running transition state search...")
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        status_text.text("Initializing TS search...")
        progress_bar.progress(10)
        
        final_geometry, optimizer, mf, converged, ts_xyz_file, initial_energy = run_ts_search(
            compound_name, smiles, pyscf_input, basis_set, theory, 
            charge=charge, spin=spin, solvent_model=solvent_model, eps=eps, 
            symmetry=symmetry, conv_params=conv_params, maxsteps=maxsteps
        )
        
        progress_bar.progress(70)
        status_text.text("TS search completed, analyzing results...")

        if converged:
            st.success(f"Transition state successfully found and saved to: {ts_xyz_file}")
        else:
            st.warning("TS search reached maximum iterations. Check convergence criteria.")
    
        # 遷移状態結果の表示
        st.subheader("Transition State Geometry")
        
        # エネルギー情報の詳細表示
        ts_energy = mf.e_tot
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Initial SCF Energy", f"{initial_energy:.8f} Hartree")
        with col2:
            st.metric("Final TS Energy", f"{ts_energy:.8f} Hartree", 
                     delta=f"{ts_energy - initial_energy:.8f}")
        
        # 最適化の統計情報
        if optimizer is not None:
            st.info("QSD Optimization completed")
            if hasattr(optimizer, 'max_cycle'):
                st.info(f"Maximum cycles allowed: {optimizer.max_cycle}")
        else:
            st.info("Alternative optimization method used")
        
        # 最適化された構造を表示
        st.text("Optimized TS Structure (XYZ format):")
        st.code(final_geometry, language="text")
        
        progress_bar.progress(80)
        
        # 振動解析の実行
        if perform_freq:
            status_text.text("Performing frequency analysis...")
            frequencies, freq_file = perform_frequency_analysis(mf, compound_name)
            
            if frequencies is not None:
                st.subheader("Vibrational Frequency Analysis")
                
                # 虚振動の確認
                imaginary_freqs = [f for f in frequencies if f < 0]
                real_freqs = [f for f in frequencies if f >= 0]
                
                if len(imaginary_freqs) == 1:
                    st.success(f"✅ **Perfect!** Found exactly 1 imaginary frequency: {imaginary_freqs[0]:.1f} cm⁻¹")
                    st.success("This indicates a proper transition state (first-order saddle point)")
                elif len(imaginary_freqs) == 0:
                    st.warning("⚠️ **No imaginary frequencies found**")
                    st.warning("This suggests a minimum structure, not a transition state")
                    st.info("💡 Try adjusting the initial geometry closer to the transition state")
                elif len(imaginary_freqs) > 1:
                    st.error(f"❌ **Found {len(imaginary_freqs)} imaginary frequencies**")
                    st.error("This indicates a higher-order saddle point, not a proper transition state")
                    st.info("💡 Consider improving the calculation level or initial geometry")
                
                # 周波数の詳細表示
                st.text(f"Frequency analysis results saved to: {freq_file}")
                
                # 虚振動と実振動の分別表示
                if len(imaginary_freqs) > 0:
                    st.write("**Imaginary Frequencies:**")
                    for i, freq in enumerate(imaginary_freqs):
                        st.write(f"  {i+1}. {freq:.2f} cm⁻¹")
                
                if len(real_freqs) > 0:
                    st.write("**Real Frequencies (first 10):**")
                    for i, freq in enumerate(real_freqs[:10]):
                        st.write(f"  {i+1}. {freq:.2f} cm⁻¹")
                    if len(real_freqs) > 10:
                        st.write(f"  ... and {len(real_freqs) - 10} more real frequencies")
        
        progress_bar.progress(100)
        status_text.text("Calculation completed!")

    except Exception as e:
        st.error(f"計算中にエラーが発生しました: {e}")
        st.error("設定を確認し、より小さな分子で試すか、計算レベルを下げてみてください。")

# 情報セクション
with st.expander("ℹ️ Transition State Search Tips"):
    st.markdown("""
    **良い遷移状態探索のためのヒント:**
    
    1. **初期構造の重要性**: 遷移状態に近い構造から開始することが重要です
    2. **反応座標**: 切断/形成される結合を考慮した構造を用意してください
    3. **虚振動の確認**: 遷移状態は1つの虚振動を持つ必要があります
    4. **収束条件**: より厳しい収束条件を使用することを推奨します
    5. **基底関数**: 遷移状態計算には適切なサイズの基底関数を使用してください
    
    **推奨計算レベル:**
    - **理論手法**: DFT-B3LYP または DFT-PBE（HFよりも精度が高い）
    - **基底関数**: 6-31G(d) 以上（STO-3Gは不適）
    - **最小推奨**: DFT-B3LYP/6-31G(d)
    - **高精度**: DFT-B3LYP/6-311G(d,p)
    
    **なぜHF/STO-3Gが問題か:**
    - STO-3Gは最小基底関数系で、分子の柔軟性を正しく記述できない
    - HF法は電子相関を無視するため、遷移状態の微細な構造を捉えられない
    - 結果として複数の虚振動や不正確な構造が得られる可能性がある
    
    **結果の解釈:**
    - 1つの虚振動: 正しい遷移状態（1次の鞍点）
    - 0個の虚振動: 極小構造（遷移状態ではない）
    - 2個以上の虚振動: 高次の鞍点（通常は望ましくない）→ 計算レベルの改善が必要
    """)

with st.expander("⚙️ 計算レベルの選択指針"):
    st.markdown("""
    **遷移状態探索の計算レベル選択:**
    
    | 計算レベル | 適用性 | 推奨度 | 備考 |
    |------------|--------|--------|------|
    | HF/STO-3G | × 不適 | ⭐ | 最も基本的、遷移状態には不十分 |
    | HF/6-31G(d) | △ 限定的 | ⭐⭐ | 電子相関の欠如により不正確 |
    | DFT-B3LYP/STO-3G | △ 限定的 | ⭐⭐ | 基底関数が不十分 |
    | DFT-B3LYP/6-31G(d) | ○ 推奨 | ⭐⭐⭐⭐ | 最小推奨レベル |
    | DFT-B3LYP/6-311G(d,p) | ◎ 高精度 | ⭐⭐⭐⭐⭐ | 高精度計算 |
    
    **計算時間との兼ね合い:**
    - 小分子（～10原子）: DFT-B3LYP/6-311G(d,p)
    - 中分子（10-30原子）: DFT-B3LYP/6-31G(d)
    - 大分子（30原子以上）: DFT-B3LYP/6-31G または より小さな基底
    """)

st.markdown("---")
st.markdown("**注意**: 遷移状態探索は計算コストが高い処理です。小さな分子から始めることを推奨します。")
