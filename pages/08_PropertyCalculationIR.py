"""
TODO: 改修中
IRの予測
"""

import streamlit as st
import matplotlib.pyplot as plt
import sys
import os
from rdkit import Chem

import numpy as np
import pandas as pd  

from logic.data_loader import list_chk_files, load_mo_info
from logic.molecule_handler import MoleculeHandler

from pyscf import gto, dft, lib
from pyscf.hessian import rks
from pyscf.prop import infrared
from pyscf.hessian import thermo

st.title("PySCF IR計算")

st.warning("精度は改善余地あります。検討中です")

# dataディレクトリ内のchkファイル一覧を取得
chk_files = list_chk_files("data")
chk_file_path = st.selectbox("計算に使うチェックポイントファイルを選択", chk_files)

# --- 構造確認ボタン ---
if st.button("構造を確認"):
    try:
        mol, _, _, _, _, _ = load_mo_info(chk_file_path)
        xyz_str = mol.atom
        handler = MoleculeHandler(xyz_str, input_type="xyz")
        compound_name = Chem.MolToInchiKey(handler.mol)
        directory = os.path.join("data", compound_name)
        smiles = Chem.MolToSmiles(handler.mol)
        handler.generate_2d_image(f"{directory}/molecule_2d.png")
        st.image(f"{directory}/molecule_2d.png", caption=smiles)
        st.write(f"化合物名（InChIKey）: {compound_name}")
        st.write(f"SMILES: {smiles}")
    except Exception as e:
        st.error(f"構造の表示に失敗しました: {e}")


# 基底関数セットの選択
basis_set = st.selectbox(
    "基底関数セット",
    ("sto-3g", "6-31g(d)"),
    index=1
)

# IRスペクトル計算・表示ボタン
if st.button("IRスペクトルを計算・表示"):
    try:
        path = chk_file_path.replace('.chk', '')
        mol = lib.chkfile.load_mol(path + '.chk')
        mf = dft.RKS(mol)
        mf.chkfile = path + '.chk'
        mf.xc = 'b3lyp'
        mf.kernel()
        mf_ir = infrared.rks.Infrared(mf).run()
        thermo.dump_thermo(mol, thermo.thermo(mf, mf_ir.vib_dict["freq_au"], 298.15, 101325))

        # mf_ir.vib_dict からIRスペクトルデータを取得
        if (
            hasattr(mf_ir, "vib_dict")
            and isinstance(mf_ir.vib_dict, dict)
            and "freq_wavenumber" in mf_ir.vib_dict
            and "ir_inten" in mf_ir.__dict__
        ):
            freq_raw = mf_ir.vib_dict["freq_wavenumber"] * 0.960  # スケーリング
            intensity_raw = mf_ir.__dict__["ir_inten"]

            # 実数部分のみ抽出し、正の値のみを使う
            freq = freq_raw.real
            intensity = intensity_raw.real if hasattr(intensity_raw, "real") else intensity_raw
            mask = freq > 0
            freq = freq[mask]
            intensity = intensity[mask]

            # データフレームとして表示
            df = pd.DataFrame({
                "Wavenumber (cm-1)": freq,
                "Intensity (a.u.)": intensity
            })
            st.dataframe(df)

            # --- 吸光度 vs 波数の図を追加 ---
            def lorentz(x, x0, I, w):
                return I * (0.5 * w)**2 / ((x - x0)**2 + (0.5 * w)**2)
            # x軸範囲を4000～500 cm-1に設定
            x = np.linspace(4000, 500, 2000)
            y = np.zeros_like(x)
            for f, inten in zip(freq, intensity):
                y += lorentz(x, f, inten, 100)
            y = y / y.max()  # 正規化

            # 吸光度プロット
            st.write("IR吸光度スペクトル（正規化強度）")
            fig_abs, ax_abs = plt.subplots(figsize=(8, 4))
            ax_abs.plot(x, y, color="red", label="Absorbance (normalized)")
            # 計算値を縦線で表示（スティックスペクトル）
            ax_abs.vlines(freq, 0, intensity / intensity.max(), color="blue", alpha=0.6, label="Calculated Peaks")
            ax_abs.set_xlabel("Wavenumber (cm$^{-1}$)")
            ax_abs.set_ylabel("Absorbance (a.u.)")
            ax_abs.set_title("IR Absorbance Spectrum")
            ax_abs.invert_xaxis()
            ax_abs.legend()
            fig_abs.tight_layout()
            st.pyplot(fig_abs)

            # --- 透過率プロット ---
            transmittance = 1 - y
            st.write("IRスペクトルの予測を表示します（透過率換算）。スケーリングファクターは0.960です。")

            fig, ax = plt.subplots(figsize=(8, 4))
            ax.plot(x, transmittance, color="black", label="calculated IR spectrum")
            ax.set_xlabel("Wavenumber (cm$^{-1}$)")
            ax.set_ylabel("Transmittance (a.u.)")
            ax.set_title("IRspectrum")
            ax.invert_xaxis()
            ax.legend()
            fig.tight_layout()
            st.pyplot(fig)
        else:
            st.error("IRスペクトルデータが見つかりません。PySCFの出力形式をご確認ください。")
            st.write("mf_ir.vib_dict:", getattr(mf_ir, "vib_dict", None))
            st.stop()
    except Exception as e:
        st.error(f"IRスペクトルの計算・表示に失敗しました: {e}")

