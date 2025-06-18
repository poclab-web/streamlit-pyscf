"""
TODO: 改修中
誘電率の予測
"""

import streamlit as st

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from pyscf import gto, scf

import numpy as np
import os

from logic.data_loader import list_chk_files, load_mo_info
from logic.calculation import compute_electric_properties
from logic.molecule_handler import MoleculeHandler

st.title("PySCF compute_electric_properties")

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

st.write("読み込まれたファイルの構造からHFでcompute_electric_propertiesを計算します。")

# 基底関数セットの選択
basis_set = st.selectbox(
    "基底関数セット",
    ("sto-3g", "6-31g(d)"),
    index=1
)

density_g_cm3 = st.number_input(
    "密度 (g/cm³)",
    min_value=0.01,
    max_value=10.0,
    value=1.0,
    step=0.01,
    help="分子の密度を入力してください。計算に使用されます。")

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
            # 電気的性質の計算
            results = compute_electric_properties(mol, basis_set=basis_set)
            st.success("計算が完了しました！")
            st.write(f"🔷 SMILES: {results['smiles']}")
            st.write(f"✅ Dipole Moment: {results['dipole_moment']:.3f} Debye")
            st.write(f"✅ Polarizability: {results['polarizability']:.3f} a.u.")
            st.write(f"✅ Dielectric Constant (ε)(calc): {results['dielectric_constant_calc']:.2f}")
            st.write(f"✅ Dielectric Constant (ε)(corrected ≈ exp): {results['dielectric_constant_pred']:.2f}")
        
    except Exception as e:
        st.error(f"計算中にエラーが発生しました: {e}")

