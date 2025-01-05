"""
振動計算とIR計算
構造最適化後のものの確認を行う。
他のもののcheckpointやXYZを呼び出して、行う。
TODO: 改修中
"""

import streamlit as st
import matplotlib.pyplot as plt
from logic.calculation import calculate_vibrational_frequencies

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
basis = st.selectbox("Select Basis Set", ["cc-pVDZ", "cc-pVTZ", "6-31G"])

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


