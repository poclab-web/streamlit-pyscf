"""
イオン化ポテンシャルの計算方法
中性分子とラジカルカチオンの両方を計算して、値を出力する。
"""

import os
import streamlit as st
from utils.module import load_css

from rdkit import Chem
import pandas as pd

from utils.module import load_css

from logic.data_loader import list_chk_files, load_mo_info
from logic.molecule_handler import MoleculeHandler

from logic.calculation import theory_options, basis_set_options, calculate_ionization_potential

# カスタムCSSを適用
load_css("config/styles.css")

# 本文
st.title("IP Calculation")

# dataディレクトリ内のchkファイル一覧を取得（新しいもの順にソート）
chk_paths = list_chk_files("data")
# 存在するファイルのみ、mtimeで新しい順にソート
chk_paths = [f for f in chk_paths if os.path.isfile(f)]
chk_paths = sorted(
    chk_paths,
    key=lambda x: os.path.getmtime(x),
    reverse=True
)
# ファイル名だけのリストを作成
chk_files = [os.path.basename(f) for f in chk_paths]

# 表示用のユニークなリストを作成（親ディレクトリ名/ファイル名）
chk_display_names = [
    os.path.join(os.path.basename(os.path.dirname(f)), os.path.basename(f)) for f in chk_paths
]
chk_file_display = st.selectbox("計算に使うチェックポイントファイルを選択", chk_display_names)

# 選択されたchkファイルのフルパスを取得
if chk_display_names:
    idx = chk_display_names.index(chk_file_display)
    chk_file_path = chk_paths[idx]
else:
    chk_file_path = None

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

st.write("読み込まれたファイルの構造からIPを計算します。")

# 計算手法の選択
theory = st.selectbox("Theory", theory_options)
basis_set = st.selectbox("Basis Set", basis_set_options)

# 中性分子の構造最適化の有無
optimize_neutral = st.checkbox(
    "中性分子の構造最適化を行う",
    value=False,
    help="中性分子の構造最適化を行うかどうかを選択してください。"
)

# ラジカルカチオンの構造最適化の有無
optimize_radical_cation = st.checkbox(
    "ラジカルカチオンも構造最適化を行う",
    value=False,
    help="ラジカルカチオンの構造最適化を行うかどうかを選択してください。"
)

# 溶媒効果の設定
solvent_model = st.selectbox("Select Solvent Model", ["None", "PCM", "ddCOSMO"])
eps = None  # epsのデフォルト値を設定

# Load solvent data
solvents_file = "config/solvents_epsilon.csv"
solvents_data = pd.read_csv(solvents_file)

if solvent_model in ["PCM", "ddCOSMO"]:
    # Show solvent selection dropdown
    solvent_selection = st.selectbox(
        "Select a solvent",
        [f"{row['Solvent']} (ε={row['Epsilon']})" for _, row in solvents_data.iterrows()]
    )  
    # Extract epsilon from selection
    if solvent_selection:
        eps = float(solvent_selection.split("=", 1)[-1][:-1])


# 計算の実行
if chk_file_path and st.button("計算実行"):
    try:
        st.write("チェックポイントから分子情報を復元中...")
        idx = chk_display_names.index(chk_file_display)
        selected_chk_fullpath = chk_paths[idx]
        parent_dir = os.path.basename(os.path.dirname(selected_chk_fullpath))
        chk_filename = os.path.basename(selected_chk_fullpath)
        # 親ディレクトリ名とファイル名を表示
        st.write(f"親ディレクトリ名: {parent_dir}")
        st.write(f"chkファイル名: {chk_filename}")
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
            # IP計算の実行
            st.write("IP計算中です...")
            # ここにIP計算のロジックを追加
            result = calculate_ionization_potential(
                mol, theory, basis_set,
                optimize_neutral=optimize_neutral,
                optimize_radical_cation=optimize_radical_cation,
                solvent_model=solvent_model, eps=eps
            )

            st.success("IP計算が完了しました。結果を表示します。")
            # 結果の表示
            if result:
                st.subheader("イオン化ポテンシャル（IP）計算結果")

                st.markdown(
                    f"**IP<sub>electronic</sub>（電子エネルギー差）:** {result['VIP_electronic']:.4f} eV "
                    f"({result['VIP_electronic'] * 23.0605:.2f} kcal/mol)",
                    unsafe_allow_html=True
                )
                if result["VIP_ZPE"] is not None:
                    st.markdown(
                        f"**IP<sub>ZPE</sub>（ZPE補正込み）:** {result['VIP_ZPE']:.4f} eV "
                        f"({result['VIP_ZPE'] * 23.0605:.2f} kcal/mol)",
                        unsafe_allow_html=True
                    )
                if result["VIP_Gibbs"] is not None:
                    st.markdown(
                        f"**IP<sub>Gibbs</sub>（Gibbs補正込み）:** {result['VIP_Gibbs']:.4f} eV "
                        f"({result['VIP_Gibbs'] * 23.0605:.2f} kcal/mol)",
                        unsafe_allow_html=True
                    )

                with st.expander("中性分子の詳細"):
                    st.markdown(f"- エネルギー: {result['neutral']['energy']:.8f} Hartree")
                    st.markdown(f"- ZPE: {result['neutral']['zpe']}")
                    st.markdown(f"- Gibbs Free Energy: {result['neutral']['gibbs']}")
                    st.markdown(f"- 振動数: {result['neutral']['freqs']}")

                with st.expander("カチオンの詳細"):
                    st.markdown(f"- エネルギー: {result['cation']['energy']:.8f} Hartree")
                    st.markdown(f"- ZPE: {result['cation']['zpe']}")
                    st.markdown(f"- Gibbs Free Energy: {result['cation']['gibbs']}")
                    st.markdown(f"- 振動数: {result['cation']['freqs']}")
    except Exception as e:
        st.error(f"計算中にエラーが発生しました: {e}")
