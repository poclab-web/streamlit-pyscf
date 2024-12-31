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

import streamlit as st
import matplotlib.pyplot as plt

from logic.calculation import theory_options, basis_set_options
from logic.calculation import run_geometry_optimization

st.title("PySCF + Streamlit: Geometry Optimization")

# サイドバーの入力
st.header("Molecular Input")
atom_input = st.text_area("Atoms (XYZ format)", "H 0 0 0\nH 0 0 0.74")
basis_set = st.selectbox("Basis Set", basis_set_options)
theory = st.selectbox("Theory", theory_options)


# 収束条件の入力
st.header("Convergence Parameters")
gradient_tol = st.number_input("Gradient Tolerance (a.u.)", min_value=1e-6, value=1e-4, step=1e-6, format="%.6f")
step_tol = st.number_input("Step Tolerance (a.u.)", min_value=1e-6, value=1e-4, step=1e-6, format="%.6f")
energy_tol = st.number_input("Energy Tolerance (Hartree)", min_value=1e-6, value=1e-6, step=1e-6, format="%.6f")
max_iterations = st.number_input("Max Iterations", min_value=1, value=100, step=1)

# 収束条件の辞書を作成
conv_params = {
    "gradient_tolerance": gradient_tol,
    "step_tolerance": step_tol,
    "energy_tolerance": energy_tol,
    "max_iterations": max_iterations,
}

# 計算実行ボタン
if st.button("Run Geometry Optimization"):
    try:
        st.write("Running geometry optimization...")
        energies, geometries = run_geometry_optimization(atom_input, basis_set, theory, conv_params)

        # エネルギー収束のプロット
        st.write("### Energy Convergence")
        plt.figure(figsize=(10, 5))
        plt.plot(energies, marker='o')
        plt.title('Energy Convergence During Geometry Optimization')
        plt.xlabel('Optimization Step')
        plt.ylabel('Total Energy (Hartree)')
        plt.grid(True)
        st.pyplot(plt)

        # 最終構造の表示
        st.write("### Final Optimized Geometry")
        st.text(f"Final Energy: {energies[-1]:.6f} Hartree")
        for coord, symbol in zip(geometries[-1], atom_input.split("\n")):
            st.text(f"{symbol.split()[0]}: {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}")

        # 全ての構造を保存するオプション
        if st.checkbox("Save All Geometries"):
            for i, coords in enumerate(geometries):
                filename = f"geometry_step_{i+1}.xyz"
                with open(filename, 'w') as f:
                    f.write(f"{len(atom_input.splitlines())}\n")
                    f.write(f"Step {i+1} Energy: {energies[i]}\n")
                    for coord, symbol in zip(coords, atom_input.split("\n")):
                        f.write(f"{symbol.split()[0]} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n")
            st.success("All geometries have been saved as XYZ files.")
    except Exception as e:
        st.error(f"Error: {e}")

