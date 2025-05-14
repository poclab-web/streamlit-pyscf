"""
振動計算とIR計算
構造最適化後のものの確認を行う。
他のもののcheckpointやXYZを呼び出して、行う。
TODO: 改修中
"""

import streamlit as st
import matplotlib.pyplot as plt
from logic.calculation import calculate_vibrational_frequencies
from logic.calculation import theory_options, basis_set_options

st.title("Vibrational Frequency Calculator")

# ユーザー入力の設定
st.header("Molecule Input")
molecule = st.text_area(
    "Enter the molecule coordinates (XYZ format):",
    """
    O  0.000000  0.000000  0.000000
    H  0.000000  0.758602  0.504284
    H  0.000000 -0.758602  0.504284
    """,
)

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
            frequencies = calculate_vibrational_frequencies(molecule, basis)
            st.success("Calculation completed successfully!")

            # 周波数の表示
            st.subheader("Vibrational Frequencies (cm⁻¹)")
            st.write(frequencies)

            # プロットの表示
            st.subheader("Frequency Spectrum")
            fig, ax = plt.subplots()
            ax.bar(range(len(frequencies)), frequencies, color="blue")
            ax.set_xlabel("Mode index")
            ax.set_ylabel("Frequency (cm⁻¹)")
            ax.set_title("Vibrational Frequencies")
            st.pyplot(fig)
        except Exception as e:
            st.error(f"An error occurred: {e}")


