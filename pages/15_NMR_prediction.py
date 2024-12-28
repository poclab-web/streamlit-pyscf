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
