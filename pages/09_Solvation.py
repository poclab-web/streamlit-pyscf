"""
イオン化ポテンシャルの計算方法
中性分子とラジカルカチオンの両方を計算して、値を出力する。
"""

import os
import streamlit as st
from rdkit import Chem
import pandas as pd

from utils.module import load_css

from logic.data_loader import list_chk_files, load_mo_info
from logic.molecule_handler import MoleculeHandler

from logic.calculation import theory_options, basis_set_options, calculate_solvation_energy

# カスタムCSSを適用
load_css("config/styles.css")

# 本文
st.title("Solvation")

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

# --- 構造確認ボタン ---
if st.button("構造を確認"):
    try:
        idx = chk_display_names.index(chk_file_display)
        selected_chk_fullpath = chk_paths[idx]
        mol, _, _, _, _, _ = load_mo_info(selected_chk_fullpath)
        xyz_str = mol.atom
        handler = MoleculeHandler(xyz_str, input_type="xyz")
        compound_name = Chem.MolToInchiKey(handler.mol)
        directory = os.path.dirname(selected_chk_fullpath)  # 修正
        smiles = Chem.MolToSmiles(handler.mol)
        handler.generate_2d_image(f"{directory}/molecule_2d.png")
        st.image(f"{directory}/molecule_2d.png", caption=smiles)
        st.write(f"化合物名（InChIKey）: {compound_name}")
        st.write(f"SMILES: {smiles}")
    except Exception as e:
        st.error(f"構造の表示に失敗しました: {e}")

st.write("読み込まれたファイルの構造から溶媒和エネルギーを計算します。")

# 計算手法の選択
theory = st.selectbox("Theory", theory_options)
basis_set = st.selectbox("Basis Set", basis_set_options)

# 構造最適化の有無（気相・溶媒中で個別に選択）
col1, col2 = st.columns(2)
with col1:
    optimize_gas = st.checkbox(
        "気相中で構造最適化を行う",
        value=False,
        help="気相中で分子の構造最適化を行うかどうかを選択してください。"
    )
with col2:
    optimize_solvent = st.checkbox(
        "溶媒中で構造最適化を行う",
        value=False,
        help="溶媒中で分子の構造最適化を行うかどうかを選択してください。"
    )

# 溶媒効果の設定
solvent_model = st.selectbox("溶媒モデルを選択", ["PCM", "DDCOSMO"])
eps = None  # epsのデフォルト値を設定

# 溶媒データの読み込み
solvents_file = "config/solvents_epsilon.csv"
solvents_data = pd.read_csv(solvents_file)

if solvent_model in ["PCM", "DDCOSMO"]:
    # 溶媒選択ドロップダウン
    solvent_selection = st.selectbox(
        "溶媒を選択",
        [f"{row['Solvent']} (ε={row['Epsilon']})" for _, row in solvents_data.iterrows()]
    )  
    # εの抽出
    if solvent_selection:
        eps = float(solvent_selection.split("=", 1)[-1][:-1])

# 計算の実行
if chk_display_names and chk_file_display and st.button("計算実行"):
    try:
        idx = chk_display_names.index(chk_file_display)
        selected_chk_fullpath = chk_paths[idx]
        parent_dir = os.path.basename(os.path.dirname(selected_chk_fullpath))
        chk_filename = os.path.basename(selected_chk_fullpath)
        # 親ディレクトリ名とファイル名を表示
        st.write(f"親ディレクトリ名: {parent_dir}")
        st.write(f"chkファイル名: {chk_filename}")

        mol, _, _, _, _, mf = load_mo_info(selected_chk_fullpath)
        xyz_str = mol.atom
        handler = MoleculeHandler(xyz_str, input_type="xyz")
        compound_name = Chem.MolToInchiKey(handler.mol)
        smiles = Chem.MolToSmiles(handler.mol)
        directory = os.path.dirname(selected_chk_fullpath)
        img_path = f"{directory}/molecule_2d_with_index.png"
        handler.generate_2d_image_with_atom_index(img_path)
        st.image(img_path, caption=f"{smiles}（原子インデックス付き）")

        if mol is None or mf is None:
            st.error("分子情報の読み込みに失敗しました。")
        else:
            # 溶媒名を抽出
            if solvent_model in ["PCM", "DDCOSMO"] and solvent_selection:
                solvent_name = solvent_selection.split(" (")[0].strip()
            else:
                solvent_name = None

            # デバッグ: 計算に渡す値を表示
            st.write(f"mol: {mol}")
            st.write(f"theory: {theory}")
            st.write(f"basis_set: {basis_set}")
            st.write(f"optimize_gas: {optimize_gas}")
            st.write(f"optimize_solvent: {optimize_solvent}")
            st.write(f"solvent_model: {solvent_model}")
            st.write(f"eps: {eps}")
            st.write(f"solvent_name: {solvent_name}")

            # 溶媒和エネルギー計算の実行
            st.write("溶媒和エネルギー計算中です...")
            result = calculate_solvation_energy(
                mol, theory, basis_set,
                optimize_gas=optimize_gas,
                optimize_solvent=optimize_solvent,
                solvent_model=solvent_model, eps=eps,
                solvent_name=solvent_name
            )

            st.success("溶媒和エネルギー計算が完了しました。結果を表示します。")
            # 結果の表示
            if result:
                st.subheader("溶媒和エネルギー計算結果")
                st.markdown(f"**溶媒和エネルギー（ΔG_solv, kcal/mol）:** {result['solvation_energy']:.4f}")

                with st.expander("詳細情報"):
                    st.markdown(f"- ガス相エネルギー: {result['gas_phase_energy']:.8f} Hartree")
                    st.markdown(f"- 溶媒中エネルギー: {result['solvent_phase_energy']:.8f} Hartree")
                    if result.get("zpe_gas") is not None:
                        st.markdown(f"- ガス相ZPE: {result['zpe_gas']}")
                    if result.get("zpe_solvent") is not None:
                        st.markdown(f"- 溶媒中ZPE: {result['zpe_solvent']}")
                    if result.get("gibbs_gas") is not None:
                        st.markdown(f"- ガス相Gibbs Free Energy: {result['gibbs_gas']}")
                    if result.get("gibbs_solvent") is not None:
                        st.markdown(f"- 溶媒中Gibbs Free Energy: {result['gibbs_solvent']}")
    except Exception as e:
        st.error(f"計算中にエラーが発生しました: {e}")
else:
    if not chk_display_names:
        st.warning("有効なチェックポイントファイルが見つかりません。")
