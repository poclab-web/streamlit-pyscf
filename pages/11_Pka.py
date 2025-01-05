"""
pkaの計算

"""


import streamlit as st
import numpy as np

# 定数
R = 8.314  # 気体定数 (J/mol·K)
T = 298.15  # 温度 (K)

def calculate_pka(delta_g):
    """
    自由エネルギー変化からpKaを計算
    """
    # ΔG [Hartree] を kcal/mol に変換
    delta_g_kcal = delta_g * 627.509  # 1 Hartree = 627.509 kcal/mol
    pka = delta_g_kcal / (2.303 * R * T / 4184)  # 4184 J/kcal
    return pka

def calculate_ph(pka, acid_concentration, base_concentration):
    """
    pKaと酸・塩基濃度からpHを計算
    """
    ratio = base_concentration / acid_concentration
    ph = pka + np.log10(ratio)
    return ph

# StreamlitアプリのUI
st.title("pKaとpHの計算アプリ")
st.write("PySCFで計算した自由エネルギーを使用してpKaとpHを計算します。")

# 入力フィールド
st.header("1. 自由エネルギーの入力")
g_ha = st.number_input("酸 (HA) の自由エネルギー [Hartree]", value=-76.0, format="%.6f")
g_a = st.number_input("共役塩基 (A⁻) の自由エネルギー [Hartree]", value=-75.5, format="%.6f")
g_h = st.number_input("プロトン (H⁺) の自由エネルギー [Hartree]", value=-0.001, format="%.6f")

# 自由エネルギー変化の計算
delta_g = g_a + g_h - g_ha
st.write(f"自由エネルギー変化 (ΔG): {delta_g:.6f} Hartree")

# pKa計算
pka = calculate_pka(delta_g)
st.write(f"計算されたpKa: {pka:.2f}")

# pH計算
st.header("2. 酸・塩基濃度の入力")
acid_concentration = st.number_input("酸 (HA) の濃度 [mol/L]", value=0.1, format="%.4f")
base_concentration = st.number_input("共役塩基 (A⁻) の濃度 [mol/L]", value=0.01, format="%.4f")

if acid_concentration > 0 and base_concentration > 0:
    ph = calculate_ph(pka, acid_concentration, base_concentration)
    st.write(f"推定されるpH: {ph:.2f}")
else:
    st.warning("酸と塩基の濃度を正の値に設定してください。")

# 注意事項
st.info("このアプリは理想条件に基づいています。溶液のイオン強度や実験条件により結果が異なる場合があります。")
