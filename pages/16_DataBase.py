"""
複数計算についてのまとめについて行う。
"""
import streamlit as st
from utils.module import load_css

import streamlit as st
import pandas as pd
import os
import json
import glob
import re
import tempfile

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, inchi

from utils.module import load_css

from logic.output_handler import (
    parse_folder_and_filename,
    parse_filename,
    extract_homo_lumo_scf_from_out,
    convert_energy_units
)

from logic.database import get_molecule_from_sqlite, get_summary_statistics, get_latest_molecules, import_molecules_from_csv, get_molecules_from_sqlite
from logic.calculation import theory_options, basis_set_options
from logic.molecule_handler import MoleculeHandler

# カスタムCSSを適用
load_css("config/styles.css")

# データベースパス
db_path = "data/energy_db.sqlite"

# dataフォルダ内のフォルダ名を取得して表示
data_path = "data"  # dataフォルダのパス

# データベースの中を検索
st.title("🧪 計算したデータのsummary")

# 統計の表示
st.subheader("📊 データベースに保存済み分子の概要")
try:
    total, methods, bases, solvents = get_summary_statistics(db_path)

    st.markdown(f"- **登録分子数**: {total}")
    st.markdown(f"- **使用理論**: {', '.join(methods) if methods else 'なし'}")
    st.markdown(f"- **基底関数**: {', '.join(bases) if bases else 'なし'}")
    if solvents:
        st.markdown(f"- **使用溶媒モデル**: {', '.join(solvents)}")
    else:
        st.markdown("- **使用溶媒モデル**: なし")

except Exception as e:
    st.error(f"データベース読み込みエラー: {e}")

st.subheader("🆕 最新登録分子（件数を選択可）")
num_latest = st.number_input("表示する件数", min_value=1, max_value=50, value=5, step=1)
try:
    latest_mols = get_latest_molecules(int(num_latest), db_path=db_path)
    if latest_mols:
        df_latest = pd.DataFrame(latest_mols)
        st.dataframe(df_latest)
    else:
        st.write("データがありません。")
except Exception as e:
    st.error(f"最新分子取得エラー: {e}")


# データベースへの情報の書き込み
st.subheader("📤 データベースの情報書込")
with st.expander("CSVファイルからのインポート", expanded=True):
    # CSVファイルのアップロード
    uploaded_file = st.file_uploader("CSVファイルをアップロード", type=["csv"])
    if uploaded_file is not None:
        try:
            with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as tmp_file:
                tmp_file.write(uploaded_file.getbuffer())
                tmp_file_path = tmp_file.name
            import_molecules_from_csv(tmp_file_path, db_path=db_path)
            st.success("✅ データベースにインポートしました。")
        except Exception as e:
            st.error(f"CSVファイルのインポートに失敗しました: {e}")

# 保存済みの分子の検索機能
st.subheader("📊 データベースに保存済み分子の概要")
with st.expander("データベース検索", expanded=True):
    # 入力フォーム
    smiles_input = st.text_input("SMILES", "C=O")
    method = st.selectbox("計算理論", theory_options)
    basis = st.selectbox("基底関数",  basis_set_options)
    spin = st.number_input("スピン多重度", value=0, step=1)
    charge = st.number_input("電荷", value=0, step=1)
    # temperature = st.number_input("温度 (K)", value=298.15)
    # pressure = st.number_input("圧力 (atm)", value=1.0)

    # 溶媒モデル（オプション）
    use_solvent = st.checkbox("溶媒効果を使用する")
    solvent = st.text_input("溶媒名", "PCM") if use_solvent else None
    dielectric = st.number_input("誘電率", value=78.4) if use_solvent else None


    # 検索実行
    if st.button("検索"):
        handler = MoleculeHandler(smiles_input, input_type="smiles")
        mol = handler.mol
        inchikey_str = inchi.MolToInchiKey(mol)

        results = get_molecules_from_sqlite(inchikey_str, method, basis, spin, charge,
                                        solvent, dielectric, db_path=db_path)
        st.info(f"🗃️ 同じ条件のデータが{len(results)}件データベースに存在します。")
        if results:
            # DataFrameで全情報を表示
            df_results = pd.DataFrame(results)
            st.dataframe(df_results)

        else:
            st.warning("❌ 一致するデータが見つかりませんでした。")

