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
from pyscf import gto, scf, dft, grad, hessian, geomopt
from pyscf.hessian import thermo

st.title("🚀 PySCF: TS構造からIRC風の反応経路予測")

st.warning("精度は改善余地あります。検討中です")

st.markdown("""
このアプリは、PySCFで得た遷移状態構造（TS）から、
虚振動方向に沿って構造を変位させ、反応物・生成物候補を自動で最適化します。

**使用例：** F + H-F → F-H + F 反応の遷移状態（F-H-F⁻線形構造）を入力してください。
- 電荷: -1 (アニオン)
- スピン: 0 (一重項)        - 計算手法: HF, B3LYP, PBE, M06, wB97X から選択可能
- 基底関数: 6-31g

**注意：** 遷移状態構造（1つの虚振動を持つ構造）を入力する必要があります。
""")

# 入力フォーム
with st.form("ts_input"):
    col1, col2 = st.columns([2, 1])
    
    with col1:
        atom_input = st.text_area("📌 分子構造（XYZ形式）", value='''
F   0.000   0.000   -1.2
H   0.000   0.000    0.0
F   0.000   0.000    1.2
''', help="遷移状態構造を入力してください。例：F + H-F → F-H + F の線形遷移状態")
    
    with col2:
        # 計算手法の選択
        method = st.selectbox("計算手法", 
                             ["HF", "B3LYP", "PBE", "M06", "wB97X"],
                             help="計算に使用する理論レベルを選択してください")
        
        basis = st.text_input("基底関数 (basis set)", value="6-31g")
    
    col3, col4, col5 = st.columns(3)
    with col3:
        spin = st.number_input("スピン (mol.spin = 2S)", value=0, step=1, help="F-H-F系では通常0（一重項）")
    with col4:
        charge = st.number_input("電荷 (charge)", value=-1, step=1, help="F-H-F⁻アニオンの場合は-1")
    with col5:
        displacement = st.number_input("変位量（Bohr）", value=0.1, step=0.01, format="%.2f")
    
    submitted = st.form_submit_button("実行")

if submitted:
    st.info("計算中...少々お待ちください")
    
    # 計算手法に応じたSCF計算を行う関数
    def get_scf_method(mol, method):
        """選択された手法に応じてSCF計算オブジェクトを返す"""
        if method == "HF":
            if mol.spin == 0:
                return scf.RHF(mol)
            else:
                return scf.UHF(mol)
        else:  # DFT methods
            if mol.spin == 0:
                mf = dft.RKS(mol)
            else:
                mf = dft.UKS(mol)
            mf.xc = method.lower()  # B3LYP, PBE, M06-2X, wB97X-D
            return mf
    
    # 分子オブジェクトの作成
    mol = gto.M(atom=atom_input, basis=basis, unit="Angstrom", spin=spin, charge=charge)
    mf = get_scf_method(mol, method).run()
    
    st.success(f"✅ {method}計算完了: エネルギー = {mf.e_tot:.6f} Hartree")

    # ヘッセ行列計算 → 虚振動の抽出
    st.info("ヘッセ行列計算中...")
    if method == "HF":
        if mol.spin == 0:
            hess = hessian.RHF(mf).kernel()
        else:
            hess = hessian.UHF(mf).kernel()
    else:  # DFT methods
        if mol.spin == 0:
            hess = hessian.RKS(mf).kernel()
        else:
            hess = hessian.UKS(mf).kernel()
    
    thermo_result = thermo.harmonic_analysis(mol, hess)
    
    # 振動数と振動モードを辞書から取得
    freqs = thermo_result['freq_au']
    modes = thermo_result['norm_mode']

    imag_freqs = [f for f in freqs if f < 0]
    st.write("虚振動:", imag_freqs)

    if len(imag_freqs) == 0:
        st.error("❌ 虚振動が見つかりません。TS構造ではない可能性があります。")
    else:
        # 虚振動モードの方向に沿って構造を押し出し
        # 最初の（最も小さい）虚振動に対応するモードを取得
        imag_idx = np.where(freqs < 0)[0][0]  # 最初の虚振動のインデックス
        ts_coords = mol.atom_coords()
        mode = modes[:, imag_idx].reshape(-1, 3)  # 該当する振動モード
        coords_fwd = ts_coords + displacement * mode
        coords_rev = ts_coords - displacement * mode

        def optimize_from_coords(coords):
            mol_new = gto.M(atom=[(a, c) for a, c in zip(mol.atom_symbol(), coords)],
                            basis=basis, unit='Bohr', charge=charge, spin=spin, verbose=0)
            mf_new = get_scf_method(mol_new, method).run()
            mol_opt = geomopt.optimize(mf_new)
            energy = get_scf_method(mol_opt, method).run().e_tot
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
