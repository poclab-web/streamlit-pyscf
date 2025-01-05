"""
電子密度解析, 結合次数, 双極子モーメントなどの計算
TODO: 改修中
計算後のcheck pointを呼び出して行う。
入れられたら、NMRの予測もここに入れる。
"""

import streamlit as st
from pyscf import gto, scf, tools
import py3Dmol

# Streamlitのタイトル
st.title("PySCF CUBEファイル生成と分子構造表示アプリ")

# 分子構造入力
st.header("分子の設定")
default_atom = '''
C      0.294808    1.728651    0.000003
H     -0.501745    1.141642    0.407216
H      0.143963    1.853765   -1.051896
H      0.309979    2.687482    0.474676
H      1.227038    1.231716    0.170017
'''
atom_input = st.text_area("原子の座標 (XYZ形式)", default_atom, height=150)

# 3D表示関数
def validate_and_format_xyz(input_data):
    try:
        # 空行を削除し、各行を分割して確認
        lines = [line.strip() for line in input_data.strip().split("\n") if line.strip()]
        formatted_data = ""
        for line in lines:
            parts = line.split()
            if len(parts) != 4:  # 原子記号 + 3座標が必要
                raise ValueError("各行は '原子記号 X Y Z' の形式で入力してください")
            formatted_data += f"{parts[0]} {float(parts[1]):.6f} {float(parts[2]):.6f} {float(parts[3]):.6f}\n"
        return formatted_data.strip()
    except Exception as e:
        raise ValueError(f"入力フォーマットエラー: {e}")

def show_3d_structure(atom_data):
    # py3Dmolで3D分子構造を表示
    view = py3Dmol.view(width=800, height=400)
    view.addModel(atom_data, "xyz")  # XYZ形式のデータを読み込む
    view.setStyle({"stick": {}})  # 棒状モデルで表示
    view.zoomTo()  # ズーム設定
    return view

# 分子構造の3D表示
st.header("分子構造の3D表示")
if st.button("分子構造を表示"):
    try:
        formatted_atom_input = validate_and_format_xyz(atom_input)
        view = show_3d_structure(formatted_atom_input)
        view_html = view.js()  # py3DmolをHTMLに変換
        st.components.v1.html(view_html, height=500, width=800)
    except ValueError as e:
        st.error(f"入力エラー: {e}")
    except Exception as e:
        st.error(f"エラーが発生しました: {e}")

# 計算ボタン
if st.button("計算を実行"):
    try:
        formatted_atom_input = validate_and_format_xyz(atom_input)

        # 分子の設定
        st.write("**計算を開始します...**")
        mol = gto.Mole()
        mol.build(atom=formatted_atom_input, basis="sto-3g", charge=0, spin=0)

        # Hartree-Fock計算
        mf = scf.RHF(mol)
        energy = mf.kernel()
        st.write(f"**計算完了: エネルギー = {energy:.6f} ハートリー**")

        # CUBEファイルの生成
        cube_file = "density.cube"
        tools.cubegen.density(mol, cube_file, mf.make_rdm1())
        st.success(f"CUBEファイル {cube_file} を生成しました！")

        # ダウンロードリンク
        with open(cube_file, "rb") as file:
            btn = st.download_button(
                label="CUBEファイルをダウンロード",
                data=file,
                file_name=cube_file,
                mime="application/octet-stream"
            )
    except ValueError as e:
        st.error(f"入力エラー: {e}")
    except Exception as e:
        st.error(f"エラーが発生しました: {e}")


"""
NMRの予測

"""


import streamlit as st
from pyscf import gto, scf
from pyscf.prop.nmr import rhf

st.title("PySCF NMRシールドテンソル計算")

# 分子構造の入力
st.sidebar.header("分子構造の設定")
default_geometry = """H 0 0 0
F 0 0 1.1"""
geometry = st.sidebar.text_area("分子構造 (XYZ形式)", default_geometry, height=200)

# 基底関数セットの選択
basis_set = st.sidebar.selectbox(
    "基底関数セット",
    ("sto-3g", "cc-pvdz", "6-31g", "6-31g(d)"),
    index=1
)

# 計算ボタン
if st.sidebar.button("計算実行"):
    try:
        st.write("計算中です...")
        mol = gto.M(atom=geometry, basis=basis_set)
        mf = scf.RHF(mol)
        mf.kernel()

        nmr_calculator = rhf.NMR(mf)
        nmr_tensors = nmr_calculator.kernel()

        st.write("### 計算結果")
        results = []
        for i, tensor in enumerate(nmr_tensors):
            results.append(
                {
                    "原子番号": i + 1,
                    "元素": mol.atom_symbol(i),
                    "シールドテンソル": tensor,
                }
            )
        import pandas as pd
        df = pd.DataFrame(results)
        st.dataframe(df)

    except Exception as e:
        st.error(f"エラーが発生しました: {e}")
