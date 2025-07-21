"""
一般的な、計算について、行うもの
"""

import streamlit as st
from utils.module import load_css

import os
import pandas as pd
from rdkit import Chem

from logic.molecule_handler import MoleculeHandler
from logic.calculation import (
    theory_options, basis_set_options,
    run_geometry_optimization, run_quantum_calculation,
    calculate_vibrational_frequencies
)

# カスタムCSSを適用
load_css("config/styles.css")

st.title("General Calculation Workflow")

st.warning("調整中です。まだ使用しないでください")
st.divider()

# 1. 構造入力
st.header("1. 構造の入力")
input_type = st.selectbox("入力形式", ["XYZ", "SMILES"])
atom_input = st.text_area(
    "分子構造を入力してください",
    "C     -2.385271    1.648249    0.000003\nH     -3.181824    1.061240    0.407216\nH     -2.536116    1.773363   -1.051896\nH     -2.370100    2.607079    0.474676\nO     -1.453042    1.151314    0.170017\nH     -1.466511    0.300016   -0.251422"
    if input_type == "XYZ"
    else "CO",
)
charge = st.number_input("分子電荷", min_value=-10, max_value=10, value=0, step=1)
multiplicity = st.number_input("スピン多重度 (2S+1)", min_value=1, max_value=10, value=1, step=1)
spin = multiplicity - 1

# 2. 配座探索
st.header("2. 配座探索")
num_conformers = st.number_input("生成する配座数", value=20, min_value=1)
prune_rms_thresh = st.slider("RMS類似構造の除去閾値", min_value=0.1, max_value=2.0, value=0.5, step=0.1)
force_field = st.selectbox("力場", ["UFF", "MMFF"])
# kcal/mol差の閾値を指定
kcal_threshold = st.number_input("最小エネルギーから何kcal/mol以内の配座を計算するか", min_value=0.1, max_value=50.0, value=5.0, step=0.1)

# 3. 構造最適化・1点計算・振動計算のレベル選択
st.header("3. 構造最適化・1点計算・振動計算のレベル選択")
col1, col2 = st.columns(2)
with col1:
    opt_theory = st.selectbox("構造最適化: 理論手法", theory_options, key="opt_theory")
    opt_basis = st.selectbox("構造最適化: 基底関数", basis_set_options, key="opt_basis")
with col2:
    sp_theory = st.selectbox("1点計算・振動計算: 理論手法", theory_options, key="sp_theory")
    sp_basis = st.selectbox("1点計算・振動計算: 基底関数", basis_set_options, key="sp_basis")

if st.button("全自動計算スタート"):
    try:
        # 構造入力
        handler = MoleculeHandler(atom_input, input_type=input_type.lower())
        if not handler.mol:
            st.error("分子構造の入力に失敗しました。")
            st.stop()
        compound_name = Chem.MolToInchiKey(handler.mol)
        smiles = Chem.MolToSmiles(handler.mol)
        directory = os.path.join("data", compound_name)
        os.makedirs(directory, exist_ok=True)

        st.success("分子構造の入力に成功しました。")

        # 配座探索
        st.info("配座探索を実行中...")
        conformers = handler.generate_conformers(
            num_conformers=num_conformers,
            prune_rms_threshold=prune_rms_thresh,
            forcefield=force_field
        )
        if not conformers:
            st.error("配座探索に失敗しました。")
            st.stop()
        st.success(f"{len(conformers)}個の配座を生成しました。")

        # 配座ごとのエネルギーとインデックスを取得し、エネルギー順にソート
        conf_infos = []
        for conf in conformers:
            conf_id = conf.GetId()
            energy = float(handler.mol.GetProp(f"Energy_{conf_id}"))
            conf_infos.append((conf_id, energy))
        conf_infos.sort(key=lambda x: x[1])  # エネルギー順

        hartree_to_kcal = 627.509
        min_energy = conf_infos[0][1]
        selected_conf_infos = [(cid, e) for cid, e in conf_infos if (e - min_energy) * hartree_to_kcal <= kcal_threshold]

        st.write(f"最小エネルギー配座との差が {kcal_threshold} kcal/mol 以下の配座数: {len(selected_conf_infos)}")

        for conf_id, conf_energy in selected_conf_infos:
            xyz_coordinates = handler.get_xyz_coordinates(conf_id=conf_id)
            xyz_str = "\n".join([f"{atom:<2} {x:>10.6f} {y:>10.6f} {z:>10.6f}" for atom, x, y, z in xyz_coordinates])

            st.info(f"構造最適化を実行中... (配座ID: {conf_id})")
            pyscf_input = xyz_str
            final_geometry, mf = run_geometry_optimization(
                compound_name, smiles, pyscf_input, opt_basis, opt_theory,
                charge=charge, spin=spin, solvent_model="None", eps=None, symmetry=False, conv_params={}, maxsteps=100
            )
            st.success(f"構造最適化が完了しました。（配座ID: {conf_id}）")

            # 1点計算
            st.info(f"1点エネルギー計算を実行中...（配座ID: {conf_id}）")
            def xyz_string_to_pyscf_atom(xyz_str):
                lines = xyz_str.strip().split('\n')
                atom_lines = [line for line in lines if len(line.split()) == 4]
                return "\n".join(atom_lines)
            pyscf_input = xyz_string_to_pyscf_atom(final_geometry)
            result = run_quantum_calculation(
                compound_name, smiles, pyscf_input, sp_basis, sp_theory,
                charge=charge, spin=spin, solvent_model="None", eps=None, symmetry=False
            )
            st.success(f"1点エネルギー: {result['energy']} Hartree（配座ID: {conf_id}）")

            # 振動計算
            st.info(f"振動計算を実行中...（配座ID: {conf_id}）")
            vib_results = calculate_vibrational_frequencies(pyscf_input, sp_theory, sp_basis)
            freq_array = vib_results['frequencies']['freq_wavenumber']
            st.success(f"振動計算が完了しました。（配座ID: {conf_id}）")
            st.subheader(f"振動数 (cm⁻¹)（配座ID: {conf_id}）")
            st.write(freq_array)

            # サマリー
            st.header(f"計算サマリー（配座ID: {conf_id}）")
            st.write(f"最適化構造:\n{final_geometry}")
            st.write(f"1点エネルギー: {result['energy']} Hartree")
            st.write(f"振動数: {freq_array}")

    except Exception as e:
        import traceback
        st.error(f"エラーが発生しました: {e}")
        st.text(traceback.format_exc())


