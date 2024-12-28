"""
Streamlitを使用した量子化学計算アプリケーション

このスクリプトは、PySCFを使用して量子化学計算を実行するための
ユーザーフレンドリーなインターフェースを提供します。ユーザーは分子構造を入力し、
理論手法や基底関数系を選択して、シングルポイントエネルギーを計算することができます。

モジュール:
- `streamlit`: インタラクティブなウェブアプリケーションを作成するために使用。
- `logic.calculation`: 理論手法や基底関数系のオプションを提供し、計算を実行する
  `run_quantum_calculation` 関数を含むヘルパーモジュール。

機能:
- XYZ形式で分子構造を入力可能。
- 計算に使用する理論手法と基底関数系を選択可能。
- 量子化学計算を実行し、結果をインタラクティブに表示。
"""

import streamlit as st

from logic.calculation import theory_options, basis_set_options
from logic.calculation import run_quantum_calculation

st.title("PySCF + Streamlit: Quantum Chemistry on the Cloud")

# ユーザー入力
st.header("Molecular Input")
atom_input = st.text_area("Atoms (XYZ format)", "H 0 0 0\nH 0 0 0.74")
theory = st.selectbox("Theory", theory_options)
basis_set = st.selectbox("Basis Set", basis_set_options)


# 計算を実行
if st.button("Run Calculation"):
    try:
        st.write("Running quantum chemistry calculation...")
        energy = run_quantum_calculation(atom_input, basis_set, theory)
        st.success(f"Calculated SinglePoint Energy: {energy} Hartree")
    except Exception as e:
        st.error(f"Error: {e}")