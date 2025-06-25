"""
振動計算とIR計算
構造最適化後のものの確認を行う。
他のもののcheckpointやXYZを呼び出して、行う。
TODO: 改修中
"""

import os
import streamlit as st
from utils.module import load_css

import py3Dmol  # 3D可視化用ライブラリ
import stmol

import numpy as np

import matplotlib.pyplot as plt
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import inchi  # 追加

from config.config import conv_preset_values, solvent_models, solvents_data, solvents_file

from logic.molecule_handler import MoleculeHandler
from logic.calculation import compute_molecule_properties
from logic.calculation import theory_options, basis_set_options
from logic.visualization import generate_cjson, write_gaussian_log

from logic.database import insert_molecule_with_frequencies, get_molecule_from_sqlite, get_stable_molecule_from_sqlite

# カスタムCSSを適用
load_css("config/styles.css")

st.title("Vibrational Frequency Calculator")

# ユーザー入力の設定
st.header("Molecule Input")
input_type = st.selectbox("Select Input Type", ["XYZ", "SMILES"])
atom_input = st.text_area(
    "Enter Molecular Structure",
    "C     -2.385271    1.648249    0.000003\nH     -3.181824    1.061240    0.407216\nH     -2.536116    1.773363   -1.051896\nH     -2.370100    2.607079    0.474676\nO     -1.453042    1.151314    0.170017\nH     -1.466511    0.300016   -0.251422"
    if input_type == "XYZ"
    else "CO",
)
charge = st.number_input("Molecular Charge", min_value=-10, max_value=10, value=0, step=1)
multiplicity = st.number_input("Spin Multiplicity (2S + 1)", min_value=1, max_value=10, value=1, step=1)
spin = multiplicity - 1

# 計算の設定入力
with st.expander("Setting calculation conditions"):
    st.header("🔧 配座生成と計算設定")
    apply_force_field = st.checkbox("分子力場による構造最適化を行う", value=True)
    force_field = st.selectbox("Select Force Field", ["MMFF", "UFF"])
    num_conformers = st.number_input("Number of Conformers", value=1000)

    st.header("📐 量子化学計算設定")
    theory = st.selectbox("理論レベル", theory_options, index=theory_options.index("HF"), key="theory_selector")
    basis = st.selectbox("基底関数", basis_set_options, index=basis_set_options.index("sto-3g"), key="basis_selector")
    solvent_model = st.selectbox("Select Solvent Model", solvent_models)
    eps = None
    if solvent_model in ["PCM", "ddCOSMO"]:
        solvent_selection = st.selectbox(
            "Select a solvent",
            [f"{row['Solvent']} (ε={row['Epsilon']})" for _, row in solvents_data.iterrows()]
        )
        if solvent_selection:
            eps = float(solvent_selection.split("=", 1)[-1][:-1])

    st.header("🔍 量子化学計算による構造最適化")
    optimize_with_qc = st.checkbox("量子化学計算による構造最適化を行う", value=True)
    st.subheader("Convergence Parameters")
    preset = st.radio("Choose preset", ["Loose", "Normal", "Tight"], index=1, horizontal=True)
    vals = conv_preset_values[preset]
    convergence_energy = st.number_input("Energy Tolerance (Hartree)", min_value=1e-7, value=vals["energy"], step=1e-5, format="%.7f")
    convergence_grms = st.number_input("Gradient RMS Tolerance (Eh/Bohr)", min_value=1e-5, value=vals["grms"], step=1e-5, format="%.5f")
    convergence_gmax = st.number_input("Gradient Max Tolerance (Eh/Bohr)", min_value=1e-5, value=vals["gmax"], step=1e-5, format="%.5f")
    convergence_drms = st.number_input("Displacement RMS Tolerance (Angstrom)", min_value=1e-4, value=vals["drms"], step=1e-4, format="%.4f")
    convergence_dmax = st.number_input("Displacement Max Tolerance (Angstrom)", min_value=1e-4, value=vals["dmax"], step=1e-4, format="%.4f")
    maxsteps = st.number_input("Max Iterations", min_value=1, value=100, step=1)
    conv_params = {
        "convergence_energy": convergence_energy,
        "convergence_grms": convergence_grms,
        "convergence_gmax": convergence_gmax,
        "convergence_drms": convergence_drms,
        "convergence_dmax": convergence_dmax,
        "maxsteps": maxsteps,
    }

# DBの検索有無
# 分子情報の標準化（入力が変わるたびに都度生成）
handler = MoleculeHandler(atom_input, input_type=input_type.lower())
mol = handler.mol
smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
inchi_str = inchi.MolToInchi(mol)
inchikey_str = inchi.InchiToInchiKey(inchi_str)
compound_name = inchikey_str
directory = os.path.join("data", compound_name)
os.makedirs(directory, exist_ok=True)

# ここで既存データを検索
existing = get_molecule_from_sqlite(
    inchikey=inchikey_str,
    method=theory,
    basis=basis,
    spin=spin,
    charge=charge,
    solvent=solvent_model if solvent_model not in [None, "None"] else None,
    dielectric=eps,
    temperature=298.15,
    pressure=1.0,
    db_path="data/energy_db.sqlite"
)

stable_existing = get_stable_molecule_from_sqlite(
    inchikey=inchikey_str,
    method=theory,
    basis=basis,
    spin=spin,
    charge=charge,
    solvent=solvent_model if solvent_model not in [None, "None"] else None,
    dielectric=eps,
    temperature=298.15,
    pressure=1.0,
    db_path="data/energy_db.sqlite"
)

use_db_data = False
if existing is not None:
    st.info("🗃️ 同じ条件のデータが既にデータベースに存在します。")
    use_db_data = st.radio(
        "既存のデータベースのデータを使いますか？",
        ("使う（計算をスキップ）", "使わず新たに計算する"),
        horizontal=True,
        index=0,
        key="use_db_data_radio"
    ) == "使う（計算をスキップ）"


if st.button("Calculate"):
    # 以降は use_db_data, existing を使って分岐
    if existing is not None and use_db_data:
        # chkファイルがあれば保存
        if existing["chk_file"]:
            chk_file_path = os.path.join(directory, f"{compound_name}.chk")
            with open(chk_file_path, "wb") as f:
                f.write(existing["chk_file"])
            st.success(f"chkファイルを {chk_file_path} に復元しました。")

        # 結果の表示
        freq_array = existing["frequencies"]
        num_imaginary = existing["num_imaginary"]
        g_tot = existing["g_tot"]
        zpe = existing["zpe"]

        # --- Vibrational Frequencies ---
        st.subheader("🔬 Vibrational Frequencies (cm⁻¹)")
        if any(isinstance(f, (list, tuple, np.ndarray)) for f in freq_array):
            import itertools
            freq_array = list(itertools.chain.from_iterable(freq_array))
        if hasattr(freq_array, "tolist"):
            freq_array = freq_array.tolist()

        # freq_arrayの前処理
        def to_float_real(val):
            # 文字列ならcomplexに変換
            if isinstance(val, str):
                try:
                    val = complex(val.replace(" ", ""))
                except Exception:
                    return float("nan")
            if isinstance(val, complex):
                return float(val.real)
            return float(val)

        # --- 振動数のフラット化とfloat変換 ---
        if any(isinstance(f, (list, tuple, np.ndarray)) for f in freq_array):
            import itertools
            freq_array = list(itertools.chain.from_iterable(freq_array))
        if hasattr(freq_array, "tolist"):
            freq_array = freq_array.tolist()
        freq_array = [to_float_real(f) for f in freq_array]  # ← ここで全要素をfloat化

        st.markdown(f"**Number of modes:** {len(freq_array)}")
        freq_df = pd.DataFrame({
            "Frequency (cm⁻¹)": [float(f"{freq:.2f}") for freq in freq_array]
        }).sort_values("Frequency (cm⁻¹)", ascending=False).reset_index(drop=True)
        st.table(freq_df)

        # --- Thermodynamic Summary ---
        st.subheader("🌡 Thermodynamic Properties (298.15 K)")
        thermo = existing
        EH_TO_KCAL = 627.5095
        st.markdown(f"""
        **Zero-Point Energy (ZPE):** {thermo["zpe"]:.4f} Eh ({thermo["zpe"]*EH_TO_KCAL:.2f} kcal/mol)  
        **Gibbs Free Energy (G_tot):** {thermo["g_tot"]:.4f} Eh ({thermo["g_tot"]*EH_TO_KCAL:.2f} kcal/mol)  
        **Number of Imaginary Frequencies:** {num_imaginary}
        """)
        st.info(f"Results saved in: {compound_name}")

    else:
        # ここから下は「計算を実行」する既存の処理
        with st.spinner("Calculating vibrational frequencies..."):
            try:
                handler = MoleculeHandler(atom_input, input_type=input_type.lower())
                mol = handler.mol
                inchi_str = inchi.MolToInchi(mol)
                inchikey_str = inchi.InchiToInchiKey(inchi_str)
                smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
                compound_name = inchikey_str
                directory = os.path.join("data", compound_name)
                os.makedirs(directory, exist_ok=True)

                # 2D構造の生成
                col1, col2 = st.columns(2)
                with col1:
                    st.subheader("Input 2D Structure")
                    handler.generate_2d_image(f"{directory}/molecule_2d.png")
                    st.image(f"{directory}/molecule_2d.png", caption=Chem.MolToSmiles(handler.mol))

                # 3D構造の生成
                with col2:
                    st.subheader("Input 3D Structure")
                    try:
                        mol_block = handler.generate_3d_molblock()
                        viewer = py3Dmol.view(width=400, height=400)  # Adjust width to fit in the column
                        viewer.addModel(mol_block, "mol")
                        viewer.setStyle({"stick": {}})
                        viewer.zoomTo()
                        stmol.showmol(viewer, height=400)
                    except Exception as e:
                        st.warning(f"Unable to generate 3D structure: {e}")

                pyscf_input_G = handler.to_pyscf_input()

                results = compute_molecule_properties(
                    compound_name, smiles, pyscf_input_G, charge, spin, conv_params, theory, basis,
                    solvent_model=solvent_model, eps=eps, maxsteps=maxsteps, optimize_with_qc=optimize_with_qc
                )
                frequencies = results['frequencies']
                thermo_info = results['thermo_info']
            except Exception as e:
                import traceback
                st.error(f"An error occurred: {e}")
                st.text(traceback.format_exc())
                frequencies = None


        if frequencies is not None:
            st.success("✅ Calculation completed successfully!")

            # 虚振動の数を計算（0未満の振動数の個数）
            freq_array = results['frequencies']['freq_wavenumber']
            # フラット化
            if any(isinstance(f, (list, tuple, np.ndarray)) for f in freq_array):
                import itertools
                freq_array = list(itertools.chain.from_iterable(freq_array))
            if hasattr(freq_array, "tolist"):
                freq_array = freq_array.tolist()
            # 修正: complex型対応
            num_imaginary = sum(1 for f in freq_array if (f.real if isinstance(f, complex) else f) < 0)

            # .chkファイルのパス（例: Gaussian計算で保存した場合）
            chk_file_path = os.path.join(directory, f"{compound_name}.chk")
            if not os.path.isfile(chk_file_path):
                chk_file_path = None  # ファイルがなければNone

            # 修正: complex型をfloatに変換して保存
            freq_array_to_save = [float(f.real if isinstance(f, complex) else f) for f in freq_array]

            # データベースへ保存
            try:
                insert_molecule_with_frequencies(
                    inchi=inchi_str,
                    inchikey=inchikey_str,
                    g_tot=thermo_info["G_tot"][0],  
                    zpe=thermo_info["ZPE"][0],      
                    method=theory,
                    basis=basis,
                    charge=charge,
                    spin=spin,
                    solvent=solvent_model if solvent_model not in [None, "None"] else None,
                    dielectric=eps,
                    temperature=298.15,
                    pressure=1.0,
                    frequencies=freq_array_to_save,  # ←ここを修正
                    num_imaginary=num_imaginary,
                    chk_file_path=chk_file_path,
                    db_path="data/energy_db.sqlite"
                )
                st.info("🗃️ 結果をデータベースに保存しました。")
            except Exception as e:
                st.warning(f"データベース保存時にエラーが発生しました: {e}")

            # --- Vibrational Frequencies ---
            st.subheader("🔬 Vibrational Frequencies (cm⁻¹)")
            freq_array = results['frequencies']['freq_wavenumber']
            # もしリストのリストならフラット化
            if any(isinstance(f, (list, tuple, np.ndarray)) for f in freq_array):
                import itertools
                freq_array = list(itertools.chain.from_iterable(freq_array))
            if hasattr(freq_array, "tolist"):
                freq_array = freq_array.tolist()
            st.markdown(f"**Number of modes:** {len(freq_array)}")
            # DataFrameを作成し、降順ソート
            freq_df = pd.DataFrame({
                "Frequency (cm⁻¹)": [float(f"{freq:.2f}") for freq in freq_array]
            }).sort_values("Frequency (cm⁻¹)", ascending=False).reset_index(drop=True)
            st.table(freq_df)

            # --- Thermodynamic Summary ---
            st.subheader("🌡 Thermodynamic Properties (298.15 K)")
            thermo = thermo_info

            # ndarray型をリストに変換
            for k, v in thermo.items():
                if hasattr(v, "tolist"):
                    thermo[k] = v.tolist()

            # 変換定数
            EH_TO_KCAL = 627.5095

            st.markdown(f"""
            **Zero-Point Energy (ZPE):** {thermo["ZPE"][0]:.4f} Eh ({thermo["ZPE"][0]*EH_TO_KCAL:.2f} kcal/mol)  
            **Enthalpy (H_tot):** {thermo["H_tot"][0]:.4f} Eh ({thermo["H_tot"][0]*EH_TO_KCAL:.2f} kcal/mol)  
            **Gibbs Free Energy (G_tot):** {thermo["G_tot"][0]:.4f} Eh ({thermo["G_tot"][0]*EH_TO_KCAL:.2f} kcal/mol)  
            **Entropy (S_tot):** {thermo["S_tot"][0]:.6f} Eh/K ({thermo["S_tot"][0]*EH_TO_KCAL:.6f} kcal/mol·K)  
            """)

            st.info(f"Results saved in: {compound_name}")

            # --- Detailed View ---
            with st.expander("🔍 Full Thermodynamic Details"):
                df = {
                    "Quantity": [
                        "ZPE", "E_tot", "H_tot", "G_tot", 
                        "S_tot", "Cv_tot", "Cp_tot"
                    ],
                    "Value (Eh)": [
                        thermo["ZPE"][0],
                        thermo["E_tot"][0],
                        thermo["H_tot"][0],
                        thermo["G_tot"][0],
                        thermo["S_tot"][0],
                        thermo["Cv_tot"][0],
                        thermo["Cp_tot"][0],
                    ]
                }
                import pandas as pd
                st.dataframe(pd.DataFrame(df), use_container_width=True)
        else:
            # 既存データの確認
            existing = get_molecule_from_sqlite(
                inchikey=inchikey_str,
                method=theory,
                basis=basis,
                spin=spin,
                charge=charge,
                solvent=solvent_model if solvent_model not in [None, "None"] else None,
                dielectric=eps,
                temperature=298.15,
                pressure=1.0,
                db_path="data/energy_db.sqlite"
            )

            if existing is not None:
                st.info("🗃️ 同じ条件のデータが既にデータベースに存在します。計算をスキップし、保存済みデータを表示します。")

                # chkファイルがあれば保存
                if existing["chk_file"]:
                    chk_file_path = os.path.join(directory, f"{compound_name}.chk")
                    with open(chk_file_path, "wb") as f:
                        f.write(existing["chk_file"])
                    st.success(f"chkファイルを {chk_file_path} に復元しました。")

                # 結果の表示
                freq_array = existing["frequencies"]
                num_imaginary = existing["num_imaginary"]
                g_tot = existing["g_tot"]
                zpe = existing["zpe"]

                # --- Vibrational Frequencies ---
                st.subheader("🔬 Vibrational Frequencies (cm⁻¹)")
                # もしリストのリストならフラット化
                if any(isinstance(f, (list, tuple, np.ndarray)) for f in freq_array):
                    import itertools
                    freq_array = list(itertools.chain.from_iterable(freq_array))
                if hasattr(freq_array, "tolist"):
                    freq_array = freq_array.tolist()
                st.markdown(f"**Number of modes:** {len(freq_array)}")
                # DataFrameを作成し、降順ソート
                freq_df = pd.DataFrame({
                    "Frequency (cm⁻¹)": [float(f"{freq:.2f}") for freq in freq_array]
                }).sort_values("Frequency (cm⁻¹)", ascending=False).reset_index(drop=True)
                st.table(freq_df)

                # --- Thermodynamic Summary ---
                st.subheader("🌡 Thermodynamic Properties (298.15 K)")
                thermo = existing

                # ndarray型をリストに変換
                for k, v in thermo.items():
                    if hasattr(v, "tolist"):
                        thermo[k] = v.tolist()

                # 変換定数
                EH_TO_KCAL = 627.5095

                st.markdown(f"""
                **Zero-Point Energy (ZPE):** {thermo["ZPE"][0]:.4f} Eh ({thermo["ZPE"][0]*EH_TO_KCAL:.2f} kcal/mol)  
                **Enthalpy (H_tot):** {thermo["H_tot"][0]:.4f} Eh ({thermo["H_tot"][0]*EH_TO_KCAL:.2f} kcal/mol)  
                **Gibbs Free Energy (G_tot):** {thermo["G_tot"][0]:.4f} Eh ({thermo["G_tot"][0]*EH_TO_KCAL:.2f} kcal/mol)  
                **Entropy (S_tot):** {thermo["S_tot"][0]:.6f} Eh/K ({thermo["S_tot"][0]*EH_TO_KCAL:.6f} kcal/mol·K)  
                """)

                st.info(f"Results saved in: {compound_name}")
                # --- Detailed View ---
                with st.expander("🔍 Full Thermodynamic Details"):
                    df = {
                        "Quantity": [
                            "ZPE", "E_tot", "H_tot", "G_tot", 
                            "S_tot", "Cv_tot", "Cp_tot"
                        ],
                        "Value (Eh)": [
                            thermo["ZPE"][0],
                            thermo["E_tot"][0],
                            thermo["H_tot"][0],
                            thermo["G_tot"][0],
                            thermo["S_tot"][0],
                            thermo["Cv_tot"][0],
                            thermo["Cp_tot"][0],
                        ]
                    }
                    import pandas as pd
                    st.dataframe(pd.DataFrame(df), use_container_width=True)
            else:
                st.warning("虚振動がなく、エネルギー的に最も安定な構造はデータベースにありません。")
