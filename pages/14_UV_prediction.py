"""
UVの予測

"""

import streamlit as st
from pyscf import gto, dft, tdscf
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import matplotlib.pyplot as plt

# Streamlitアプリ
st.title("UVスペクトル予測アプリ with PySCF")
st.write("分子のSMILES表記を入力し、UVスペクトルを予測します。")

# SMILES入力
smiles_input = st.text_input("SMILES表記を入力してください", "CCO")
if st.button("予測する"):
    try:
        # SMILESから分子構造を生成
        st.text("分子構造を生成中...")
        molecule = Chem.MolFromSmiles(smiles_input)
        if not molecule:
            st.error("無効なSMILES表記です。")
            st.stop()

        # 分子構造の可視化
        st.image(Draw.MolToImage(molecule), caption="分子構造", use_container_width=True)

        # 3次元構造の生成
        molecule = Chem.AddHs(molecule)  # 水素原子を追加
        success = AllChem.EmbedMolecule(molecule, AllChem.ETKDG())  # 3次元構造の生成
        if success != 0:
            st.error("分子の3次元構造を生成できませんでした。SMILESを確認してください。")
            st.stop()
        st.success("3次元構造の生成に成功しました。")

        # XYZ形式の生成
        xyz_block = Chem.MolToXYZBlock(molecule)
        st.text("生成されたXYZ形式:")
        st.text(xyz_block)

        # 原子の数とヘッダー情報を追加
        num_atoms = molecule.GetNumAtoms()
        xyz_full = f"{num_atoms}\\nCarotene molecule\\n" + xyz_block

        # 必要な行のみを抽出（5行目以降）
        xyz_lines = xyz_full.splitlines()
        xyz_trimmed = ";".join(xyz_lines[4:])


        # PySCFの分子オブジェクトを生成
        st.text("PySCF用の入力ファイルを生成中...")
        mol = gto.M(
            atom=xyz_trimmed,           # 原子データ
            unit="Angstrom",      # 座標単位
            basis="6-31G",        # 基底関数
        )
        st.success("PySCF用の入力ファイルが生成されました。")

        # DFT計算
        st.text("DFT計算を実行中...")
        mf = dft.RKS(mol)
        mf.xc = "B3LYP"          # 汎関数指定
        energy = mf.kernel()
        st.success(f"DFT計算が完了しました。エネルギー: {energy:.6f} Ha")

        # UVスペクトル（遷移エネルギー）の計算（TDDFT）
        st.text("UVスペクトルの計算を実行中...")
        td = tdscf.TDA(mf)
        td.nstates = 10           # 最初の10個の励起状態を計算
        excitation_energies = td.kernel()
        st.success("UVスペクトルの計算が完了しました。")

        # 配列の存在チェック
        if excitation_energies is None or len(excitation_energies[0]) == 0:
            st.error("遷移エネルギーの計算結果がありません。分子や計算条件を確認してください。")
            st.stop()

        # 振動子強度の取得
        oscillator_strengths = td.oscillator_strength()
        if oscillator_strengths is None or len(oscillator_strengths) == 0:
            st.error("振動子強度の計算結果がありません。")
            st.stop()

        # 波長と振動子強度の計算
        wavelengths = []
        for i, energy in enumerate(excitation_energies[0]):
            wavelength = 1240 / (energy * 27.2114)  # エネルギー (au) を nm に変換
            wavelengths.append((wavelength, oscillator_strengths[i]))
            st.write(f"励起状態 {i + 1}: 波長 = {wavelength:.2f} nm, 振動子強度 = {oscillator_strengths[i]:.3f}")

        # UVスペクトルの可視化
        st.subheader("UVスペクトル")
        wavelengths = sorted(wavelengths, key=lambda x: x[0])  # 波長でソート
        x = [w[0] for w in wavelengths]
        y = [w[1] for w in wavelengths]

        plt.figure(figsize=(8, 4))
        plt.bar(x, y, width=10, color='blue', alpha=0.7, edgecolor='black')
        plt.title("UV")
        plt.xlabel("Wavelength (nm)")
        plt.ylabel("oscillator strength (-)")
        plt.grid(True)
        st.pyplot(plt)

    except Exception as e:
        st.error(f"エラーが発生しました: {e}")
