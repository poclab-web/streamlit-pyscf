"""
TODO: 改修中
NMRの予測
"""

import streamlit as st
from utils.module import load_css

from pyscf import gto, scf, tools
import py3Dmol
from rdkit import Chem

from logic.data_loader import list_chk_files, load_mo_info
from logic.calculation import calc_nmr_and_shift
from logic.visualization import simulate_nmr_spectrum
from logic.molecule_handler import MoleculeHandler
from logic.visualization import simulate_nmr_spectrum_from_csv

import streamlit as st
from pyscf import gto, scf
from pyscf.prop.nmr import rhf
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob

# カスタムCSSを適用
load_css("config/styles.css")


st.title("PySCF NMRシールドテンソル計算")

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

st.write("読み込まれたファイルの構造からHFでNMRの計算を行います。")

# 元素選択を水素のみ
# TODO: 他の元素も選択可能にする
atom_type = "H（水素）"
target_element = "H"

if chk_file_path and st.button("計算実行"):
    try:
        st.write("チェックポイントから分子情報を復元中...")
        mol, _, _, _, _, mf = load_mo_info(chk_file_path)
        xyz_str = mol.atom
        handler = MoleculeHandler(xyz_str, input_type="xyz")
        compound_name = Chem.MolToInchiKey(handler.mol)
        smiles = Chem.MolToSmiles(handler.mol)
        directory = os.path.join("data", compound_name)
        img_path = f"{directory}/molecule_2d_with_index.png"
        handler.generate_2d_image_with_atom_index(img_path)
        st.image(img_path, caption=f"{smiles}（原子インデックス付き）")

        if mol is None or mf is None:
            st.error("分子情報の読み込みに失敗しました。")
        else:
            st.write("NMR計算中です...")
            results, nmr_calc_result_name = calc_nmr_and_shift(mf, mol, target_element=target_element, basis_set=basis_set)

            st.write("### 計算結果")
            df = pd.DataFrame(results)
            st.dataframe(df)

            # --- TMS reference の自動選択 ---
            tms_dir = "data/CZDYPVPMEAXLPK-UHFFFAOYSA-N"
            tms_csv_files = glob.glob(os.path.join(tms_dir, "*.csv"))
            tms_reference_csv_path = None
            for tms_csv in tms_csv_files:
                tms_df = pd.read_csv(tms_csv)
                # "Basis Set"列が存在し、かつ一致するものを探す
                if "Basis Set" in tms_df.columns:
                    if tms_df["Basis Set"].astype(str).str.lower().str.replace(" ", "")\
                        .eq(basis_set.lower().replace(" ", "")).any():
                        tms_reference_csv_path = tms_csv
                        break
            if tms_reference_csv_path is not None:
                st.success(f"TMS reference: {tms_reference_csv_path} を使用します")
            else:
                st.warning("一致する基底関数のTMS referenceが見つかりませんでした")

            if nmr_calc_result_name:
                fig, table_df = simulate_nmr_spectrum_from_csv(
                    nmr_calc_result_name,
                    target_element=target_element,
                    TMS_reference_csv_path=tms_reference_csv_path
                )
                st.pyplot(fig)
                st.write("### Chemical Shift Table")
                st.dataframe(table_df)
            else:
                st.warning("計算結果のCSVファイルが見つかりませんでした。")

    except Exception as e:
        st.error(f"エラーが発生しました: {e}")
else:
    st.info("dataディレクトリにchkファイルを置いてください。")
