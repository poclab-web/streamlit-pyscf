"""
遷移状態の計算
IRC計算を行う。
"""

# import streamlit as st
# from utils.module import load_css

# # カスタムCSSを適用
# load_css("config/styles.css")



import streamlit as st
import numpy as np
from pyscf import gto, scf, grad, hessian, geomopt
from pyscf.hessian import thermo

st.title("🚀 PySCF: TS構造からIRC風の反応経路予測")

st.markdown("""
このアプリは、PySCFで得た遷移状態構造（TS）から、
虚振動方向に沿って構造を変位させ、反応物・生成物候補を自動で最適化します。
""")

# 入力フォーム
with st.form("ts_input"):
    atom_input = st.text_area("📌 分子構造（XYZ形式）", value='''
H 0.000 0.000 -0.400
H 0.000 0.000  0.400
H 0.000 0.000  1.800
''')
    basis = st.text_input("基底関数 (basis set)", value="6-31g")
    spin = st.number_input("スピン (mol.spin = 2S)", value=1, step=1)
    charge = st.number_input("電荷 (charge)", value=0, step=1)
    displacement = st.number_input("変位量（Bohr）", value=0.1, step=0.01, format="%.2f")
    submitted = st.form_submit_button("実行")

if submitted:
    st.info("計算中...少々お待ちください")
    
    # 分子オブジェクトの作成
    mol = gto.M(atom=atom_input, basis=basis, unit="Angstrom", spin=spin, charge=charge)
    mf = scf.RHF(mol).run()

    # ヘッセ行列計算 → 虚振動の抽出
    hess = hessian.RHF(mf).kernel()
    freqs, modes = thermo.harmonic_analysis(mol, hess)[2:4]

    imag_freqs = [f for f in freqs if f < 0]
    st.write("虚振動:", imag_freqs)

    if len(imag_freqs) == 0:
        st.error("❌ 虚振動が見つかりません。TS構造ではない可能性があります。")
    else:
        # 虚振動モードの方向に沿って構造を押し出し
        ts_coords = mol.atom_coords()
        mode = modes[0].reshape(-1, 3)
        coords_fwd = ts_coords + displacement * mode
        coords_rev = ts_coords - displacement * mode

        def optimize_from_coords(coords):
            mol_new = gto.M(atom=[(a, c) for a, c in zip(mol.atom_symbol(), coords)],
                            basis=basis, unit='Bohr', charge=charge, spin=spin, verbose=0)
            mf_new = scf.RHF(mol_new).run()
            mol_opt = geomopt.optimize(mf_new)
            energy = scf.RHF(mol_opt).run().e_tot
            return mol_opt, energy

        st.write("🔁 正方向に沿って最適化中...")
        mol_fwd, e_fwd = optimize_from_coords(coords_fwd)
        st.success(f"Forward最適化完了: エネルギー = {e_fwd:.6f} Hartree")

        st.write("🔁 逆方向に沿って最適化中...")
        mol_rev, e_rev = optimize_from_coords(coords_rev)
        st.success(f"Reverse最適化完了: エネルギー = {e_rev:.6f} Hartree")

        st.subheader("✅ 最適化後の構造（Angstrom）")
        st.text("Forward構造:")
        for a, c in zip(mol.atom_symbol(), mol_fwd.atom_coords()):
            st.text(f"{a:2s}  {c[0]: .4f}  {c[1]: .4f}  {c[2]: .4f}")
        st.text("Reverse構造:")
        for a, c in zip(mol.atom_symbol(), mol_rev.atom_coords()):
            st.text(f"{a:2s}  {c[0]: .4f}  {c[1]: .4f}  {c[2]: .4f}")
