"""
pKaの計算
SMILESから脱プロトン化候補を選択し、量子化学計算条件を指定してpKaを計算
"""

# ライブラリーのimport
import streamlit as st
from utils.module import load_css
import os

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import inchi
from io import BytesIO

from logic.molecule_handler import MoleculeHandler
from logic.calculation import compute_molecule_properties, theory_options, basis_set_options
from config.config import conv_preset_values, solvent_models, solvents_data, solvents_file
from logic.visualization import show_imaginary_frequency_warning
from logic.output_handler import safe_dict, extract_energy

from logic.database import insert_molecule_with_frequencies, get_molecule_from_sqlite, get_stable_molecule_from_sqlite

# カスタムCSSを適用
load_css("config/styles.css")

# 定数
R = 8.314  # J/mol·K
T = 298.15  # K
HARTREE_TO_KCAL = 627.509

def deprotonation_iterator_with_bond(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    bond_infos = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            neighbors = atom.GetNeighbors()
            if not neighbors:
                continue
            neighbor = neighbors[0]
            h_idx = atom.GetIdx()
            heavy_idx = neighbor.GetIdx()
            bond = mol.GetBondBetweenAtoms(h_idx, heavy_idx)
            if bond is None:
                continue
            bond_idx = bond.GetIdx()
            bond_atoms = f"{mol.GetAtomWithIdx(h_idx).GetSymbol()}–{mol.GetAtomWithIdx(heavy_idx).GetSymbol()}"
            rw_mol = Chem.RWMol(mol)
            rw_mol.RemoveAtom(h_idx)
            atom_heavy = rw_mol.GetAtomWithIdx(heavy_idx)
            atom_heavy.SetFormalCharge(-1)
            try:
                Chem.SanitizeMol(rw_mol)
                deprot_smiles = Chem.MolToSmiles(rw_mol)
                bond_infos.append({
                    "bond_index": bond_idx,
                    "bond_atoms": bond_atoms,
                    "deprotonated_smiles": deprot_smiles
                })
            except Exception:
                continue
    return pd.DataFrame(bond_infos)

# def calculate_pka(delta_g_hartree):
#     delta_g_kcal = delta_g_hartree * HARTREE_TO_KCAL
#     pka = delta_g_kcal / (2.303 * R * T / 4184)
#     return pka

def calculate_relative_pka(delta_g_sample, delta_g_ref, pka_ref):
    """Thermodynamic Cycle A（相対法）によるpKa計算"""
    delta_g_kcal = (delta_g_sample - delta_g_ref) * HARTREE_TO_KCAL
    pka = pka_ref + delta_g_kcal / (2.303 * R * T / 4184)
    return pka

def mol_to_image(smiles):
    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol, size=(300, 300))
    buf = BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()

st.title("pKa Calculator in water")
st.markdown("**SMILESから脱プロトン化候補を選択し、量子化学計算条件を指定してpKa(水中)を計算**")
st.warning("精度は改善余地あります。検討中です")

st.divider()

# Thermodynamic Cycle A（相対法）の基準物質
# pka listのCSV読み込み
# pka listのCSVファイルのパス
pka_csv_path = "config/pKa_reference_list.csv"

if os.path.exists(pka_csv_path):
    pka_df = pd.read_csv(pka_csv_path)   
else:
    st.warning(f"{pka_csv_path} が見つかりません。")    


# pka_dfから計算に必要な物質をユーザーが選択
if not pka_df.empty:
    st.subheader("🔍 pKa Reference List")
    st.markdown("以下のリストから、pkaの基準となる物質を選択できます。")
    # デフォルトでEthanolを選択
    default_index = pka_df["name"].tolist().index("Ethanol") if "Ethanol" in pka_df["name"].tolist() else 0
    selected_pka = st.selectbox("Select a reference compound", pka_df["name"].tolist(), index=default_index)
    with st.expander("Show details of selected compound"):
        if selected_pka:
            selected_row = pka_df[pka_df["name"] == selected_pka].iloc[0]
            st.write(f"Selected Compound: {selected_row['name']}")
            st.write(f"Experimental pKa: {selected_row['pKa_exp']}")

            # アニオンのSMILESを表示
            col1, col2 = st.columns(2)
            # Display 2D structure in the first column
            with col1:
                st.subheader("Structure")
                # SMILESを表示
                st.write(f"SMILES: {selected_row['reference_smiles']}")
                st.image(mol_to_image(selected_row['reference_smiles']), caption=f"SMILES: {selected_row['reference_smiles']}")

            # Display 3D structure in the second column
            with col2:
                st.subheader("anion Structure")
                st.write(f"Anion SMILES: {selected_row['reference_anion_smiles']}")
                st.image(mol_to_image(selected_row['reference_anion_smiles']), caption=f"Anion SMILES: {selected_row['reference_anion_smiles']}")   

# ▼ SMILES入力
smiles_input = st.text_input("計算したい分子のSMILESを入力してください", "CO")

# 計算の設定入力
with st.expander("Setting calculation conditions"):
    st.header("🔧 配座生成と計算設定")
    apply_force_field = st.checkbox("分子力場による構造最適化を行う", value=True)
    force_field = st.selectbox("Select Force Field", ["MMFF", "UFF"])
    num_conformers = st.number_input("Number of Conformers", value=1000)

    st.header("📐 量子化学計算設定")
    theory = st.selectbox("理論レベル", theory_options, index=theory_options.index("HF"), key="theory_selector")
    basis = st.selectbox("基底関数", basis_set_options, index=basis_set_options.index("sto-3g"), key="basis_selector")
    # ▼ solvent_modelのデフォルトを"ddCOSMO"に
    default_solvent_model = "ddCOSMO" if "ddCOSMO" in solvent_models else solvent_models[0]
    solvent_model = st.selectbox("Select Solvent Model", solvent_models, index=solvent_models.index(default_solvent_model))
    eps = None
    if solvent_model in ["PCM", "ddCOSMO"]:
        # デフォルト溶媒（Water, ε=78.3553）を選択
        solvent_options = [f"{row['Solvent']} (ε={row['Epsilon']})" for _, row in solvents_data.iterrows()]
        default_solvent = next((i for i, row in solvents_data.iterrows() if row['Solvent'] == "Water" and float(row['Epsilon']) == 78.3553), 0)
        solvent_selection = st.selectbox(
            "Select a solvent",
            solvent_options,
            index=default_solvent
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

# 画像表示
handler = MoleculeHandler(smiles_input, input_type="smiles")
show_image = st.checkbox("🖼️ Bond Index画像を表示する", value=True)
if show_image:
    try:
        image_buffer = handler.generate_2d_image_with_bond_index()
        if image_buffer:
            st.image(image_buffer, caption="Bond Index")
        else:
            st.warning("画像バッファが生成されませんでした。")
    except Exception as e:
        st.error(f"画像生成中にエラーが発生しました: {e}")
else:
    st.info("画像は表示されません。")

# 脱プロトン化候補生成
if st.button("脱プロトン構造を生成"):
    with st.spinner("脱プロトン化候補を生成しています..."):
        df_bonds = deprotonation_iterator_with_bond(smiles_input)
        df_bonds_unique = df_bonds.drop_duplicates(subset=["deprotonated_smiles"])
        st.session_state["deprot_candidates"] = df_bonds_unique
        if df_bonds_unique.empty:
            st.warning("脱プロトン候補が見つかりませんでした。")
        else:
            st.success("脱プロトン候補を生成しました。下から選択してください。")

# 候補選択とpKa計算
if "deprot_candidates" in st.session_state:
    st.subheader("🔍 脱プロトン化候補を選択してください（1つだけ選択可）")
    df_bonds = st.session_state["deprot_candidates"]
    selected_rows = []
    options = [
        f"Bond {row['bond_index']} ({row['bond_atoms']}): [{row['deprotonated_smiles']}]"
        for _, row in df_bonds.iterrows()
    ]
    selected_option = st.radio("候補を選択", options, key="single_deprot_candidate")
    if selected_option:
        idx = options.index(selected_option)
        selected_rows.append((idx, df_bonds.iloc[idx]))


# サンプル分子（アニオン, 脱プロトン化候補ごと）
anion_existing_dict = {}
if "deprot_candidates" in st.session_state:
    df_bonds = st.session_state["deprot_candidates"]

        # ▼ DB検索・利用選択まとめ（基準物質・サンプル分子・アニオン全て）
    # 基準物質（中性）
    ref_mol = Chem.MolFromSmiles(selected_row['reference_smiles'])
    ref_inchikey = inchi.InchiToInchiKey(inchi.MolToInchi(ref_mol))
    ref_existing = get_molecule_from_sqlite(
        inchikey=ref_inchikey,
        method=theory,
        basis=basis,
        spin=0,
        charge=0,
        solvent=solvent_model if solvent_model not in [None, "None"] else None,
        dielectric=eps,
        temperature=298.15,
        pressure=1.0,
        db_path="data/energy_db.sqlite"
    )
    # 基準物質（アニオン）
    ref_anion_mol = Chem.MolFromSmiles(selected_row['reference_anion_smiles'])
    if ref_anion_mol is None:
        st.error(f"Invalid SMILES for reference anion: {selected_row['reference_anion_smiles']}")
        ref_anion_existing = None
    else:
        ref_anion_inchikey = inchi.InchiToInchiKey(inchi.MolToInchi(ref_anion_mol))
    ref_anion_existing = get_molecule_from_sqlite(
        inchikey=ref_anion_inchikey,
        method=theory,
        basis=basis,
        spin=0,
        charge=-1,
        solvent=solvent_model if solvent_model not in [None, "None"] else None,
        dielectric=eps,
        temperature=298.15,
        pressure=1.0,
        db_path="data/energy_db.sqlite"
    )
    # サンプル分子（中性）
    handler_sample = MoleculeHandler(smiles_input, input_type="smiles")
    sample_smiles = Chem.MolToSmiles(handler_sample.mol, isomericSmiles=True, canonical=True)
    sample_inchikey = inchi.InchiToInchiKey(inchi.MolToInchi(handler_sample.mol))
    sample_existing = get_molecule_from_sqlite(
        inchikey=sample_inchikey,
        method=theory,
        basis=basis,
        spin=0,
        charge=0,
        solvent=solvent_model if solvent_model not in [None, "None"] else None,
        dielectric=eps,
        temperature=298.15,
        pressure=1.0,
        db_path="data/energy_db.sqlite"
    )

    # サンプル分子（アニオン, 脱プロトン化候補ごと）
    for idx, row in df_bonds.iterrows():
        anion_smiles = row["deprotonated_smiles"]
        anion_mol = Chem.MolFromSmiles(anion_smiles)
        if anion_mol is None:
            st.error(f"Invalid SMILES for anion: {anion_smiles}")
            continue    
         
        anion_inchikey = inchi.InchiToInchiKey(inchi.MolToInchi(anion_mol))
        anion_existing = get_molecule_from_sqlite(
            inchikey=anion_inchikey,
            method=theory,
            basis=basis,
            spin=0,
            charge=-1,
            solvent=solvent_model if solvent_model not in [None, "None"] else None,
            dielectric=eps,
            temperature=298.15,
            pressure=1.0,
            db_path="data/energy_db.sqlite"
        )
        anion_existing_dict[idx] = anion_existing
    

    # ▼ DB利用選択UIまとめて表示
    use_db_data_ref = False
    if ref_existing is not None:
        st.info("🗃️ 基準物質（中性）の同じ条件のデータが既にデータベースに存在します。")
        use_db_data_ref = st.radio(
            "基準物質（中性）の既存データを使いますか？",
            ("使う（計算をスキップ）", "使わず新たに計算する"),
            horizontal=True,
            index=0,
            key="use_db_data_radio_ref"
        ) == "使う（計算をスキップ）"


    use_db_data_ref_anion = False
    if ref_anion_existing is not None:
        st.info("🗃️ 基準物質（アニオン）の同じ条件のデータが既にデータベースに存在します。")
        use_db_data_ref_anion = st.radio(
            "基準物質（アニオン）の既存データを使いますか？",
            ("使う（計算をスキップ）", "使わず新たに計算する"),
            horizontal=True,
            index=0,
            key="use_db_data_radio_ref_anion"
        ) == "使う（計算をスキップ）"


    use_db_data_sample = False
    if sample_existing is not None:
        st.info("🗃️ サンプル分子（中性）の同じ条件のデータが既にデータベースに存在します。")
        use_db_data_sample = st.radio(
            "サンプル分子（中性）の既存データを使いますか？",
            ("使う（計算をスキップ）", "使わず新たに計算する"),
            horizontal=True,
            index=0,
            key="use_db_data_radio_sample"
        ) == "使う（計算をスキップ）"


    use_db_data_anion_dict = {}
    if "deprot_candidates" in st.session_state:
        for idx, row in df_bonds.iterrows():
            use_db_data_anion = False
            if anion_existing_dict[idx] is not None:
                st.info(f"🗃️ サンプル分子（アニオン, bond {row['bond_index']}）の同じ条件のデータが既にデータベースに存在します。")
                use_db_data_anion = st.radio(
                    f"サンプル分子（アニオン, bond {row['bond_index']}）の既存データを使いますか？",
                    ("使う（計算をスキップ）", "使わず新たに計算する"),
                    horizontal=True,
                    index=0,
                    key=f"use_db_data_radio_anion_{idx}"
                ) == "使う（計算をスキップ）"
                use_db_data_anion_dict[idx] = use_db_data_anion


    if selected_rows:
        st.markdown("### 選択された脱プロトン化候補")
        selected_df = pd.DataFrame([r for _, r in selected_rows])
        st.dataframe(selected_df)

        # 計算の方法を出力
        with st.expander("Show selected method and reference"):
            st.markdown("### 🧪 Method for pKa Calculation")
            st.markdown(
                "**Computational Details**  \n"
                "Molecular structures were generated from SMILES using RDKit [1], and 3D conformers were generated using the ETKDG method.  \n"
                f"{num_conformers} conformers were generated and optimized using the {force_field} force field.  \n"
                f"The lowest-energy conformer was selected for quantum chemical geometry optimization at the **{theory}/{basis}** level using PySCF [2].  \n"
                f"{'No solvent model was applied.' if solvent_model == 'None' else f'The solvent effect was considered using the {solvent_model} model with ε = {eps}.'}  \n"
                "Proton dissociation free energies (ΔG_deprot) were computed as the difference in Gibbs free energy between the conjugate base (anion) and the neutral acid, corrected by the free energy of a proton in solution (−270.3 kcal/mol at 298.15 K).  \n\n"
                "The pKa values were estimated using the **relative method based on an isodesmic reaction (Thermodynamic Cycle A)** [4]:  \n"
                f"- The experimental pKa of a reference compound is used: **{selected_row['name']}** (pKa = {selected_row['pKa_exp']})\n"
                "- The difference in ΔG_deprot between the sample and the reference is calculated by quantum chemical calculations.\n"
                "- The sample pKa is estimated by:\n"
                "$$\n"
                "\\mathrm{p}K_a^\\mathrm{sample} = \\mathrm{p}K_a^\\mathrm{ref} + \\frac{\\Delta G_\\mathrm{deprot}^\\mathrm{sample} - \\Delta G_\\mathrm{deprot}^\\mathrm{ref}}{2.303RT}\n"
                "$$\n"
                "where $\\Delta G_\\mathrm{deprot} = G_\\mathrm{anion} - G_\\mathrm{neutral}$.\n" \
                "All calculations were performed using a custom Streamlit-based interface [3] integrating RDKit and PySCF."
            )
            st.markdown("---")
            st.markdown(
                "**References**  \n"
                "[1] Landrum, G. RDKit: Open-source cheminformatics. [https://www.rdkit.org](https://www.rdkit.org)  \n"
                "[2] Sun, Q. *et al.* PySCF: The Python-based Simulations of Chemistry Framework. **WIREs Comput Mol Sci** *2018*, **8**, e1340. DOI: [10.1002/wcms.1340](https://doi.org/10.1002/wcms.1340)  \n"
                "[3] PocLab streamlit-pyscf: Quantum chemistry web interface. [https://github.com/poclab-web/streamlit-pyscf](https://github.com/poclab-web/streamlit-pyscf)  \n"
                "[4] Sastre, S.; Casasnovas, R.; Muñoz, F.; Frau, J. Isodesmic reaction for pKa calculations of common organic molecules. **Theor. Chem. Acc.** *2013*, **132**, 1310. DOI: [10.1007/s00214-012-1310-z](https://doi.org/10.1007/s00214-012-1310-z)"
            )


        if st.button("選択されたフラグメントでpKaを一括計算"):
            pka_results = []

            # ▼ 基準物質のGibbs自由エネルギー自動計算
            st.info("基準物質のGibbs自由エネルギーを計算中...")
            ref_smiles = selected_row['reference_smiles']
            ref_anion_smiles = selected_row['reference_anion_smiles']
            ref_name = selected_row['name']
            pka_ref = selected_row['pKa_exp']


            # DB検索（基準物質・中性）
            if use_db_data_ref:
                ref_g_neutral = ref_existing["g_tot"]
                st.success("DBから基準物質（中性）のG_totを取得しました。")
                st.write(f"基準物質（中性）G_tot: {ref_g_neutral} Eh")
            else:
                handler_ref = MoleculeHandler(ref_smiles, input_type="smiles")
                inchi_str_ref = inchi.MolToInchi(handler_ref.mol)
                inchikey_str_ref = inchi.InchiToInchiKey(inchi_str_ref)
                charge=0
                spin=0
                
                directory = os.path.join("data", inchi_str_ref)
                os.makedirs(directory, exist_ok=True)

                if apply_force_field and handler_ref.mol.GetNumAtoms() > 2:
                    handler_ref.generate_conformers(num_conformers=num_conformers, forcefield=force_field)
                    handler_ref.keep_lowest_energy_conformer()

                pyscf_input_ref = handler_ref.to_pyscf_input()
                result_ref = compute_molecule_properties(
                    name=inchikey_str_ref,
                    smiles=ref_smiles,
                    xyz=pyscf_input_ref,
                    theory=theory,
                    basis_set=basis,
                    opt_theory=theory,
                    opt_basis_set=basis,
                    charge=0,
                    spin=0,
                    solvent_model=None if solvent_model == "None" else solvent_model,
                    eps=eps,
                    conv_params=conv_params,
                    maxsteps=maxsteps
                )
                ref_g_neutral = extract_energy(result_ref, "G_tot")
                show_imaginary_frequency_warning(result_ref)
                st.write(f"基準物質（中性）G_tot: {ref_g_neutral} Eh")

                frequencies = result_ref['frequencies']
                thermo_info = result_ref['thermo_info']

                if frequencies is not None:
                    st.success("✅ Calculation completed successfully!")

                    # 虚振動の数を計算（0未満の振動数の個数）
                    freq_array = result_ref['frequencies']['freq_wavenumber']
                    # フラット化
                    if any(isinstance(f, (list, tuple, np.ndarray)) for f in freq_array):
                        import itertools
                        freq_array = list(itertools.chain.from_iterable(freq_array))
                    if hasattr(freq_array, "tolist"):
                        freq_array = freq_array.tolist()
                    # 複素数を実数に変換
                    freq_array = [f.real if isinstance(f, complex) else f for f in freq_array]
                    num_imaginary = sum(1 for f in freq_array if isinstance(f, (int, float)) and f < 0)

                    # データベースへ保存
                    try:
                        insert_molecule_with_frequencies(
                            inchi=inchi_str_ref,
                            inchikey=inchikey_str_ref,
                            g_tot=thermo_info["G_tot"][0],  
                            zpe=thermo_info["ZPE"][0],      
                            method=theory,
                            basis=basis,
                            charge=0,
                            spin=0,
                            solvent=solvent_model if solvent_model not in [None, "None"] else None,
                            dielectric=eps,
                            temperature=298.15,
                            pressure=1.0,
                            frequencies=freq_array,
                            num_imaginary=num_imaginary,
                            db_path="data/energy_db.sqlite"
                        )
                        st.info("🗃️ 結果をデータベースに保存しました。")
                    except Exception as e:
                        st.warning(f"データベース保存時にエラーが発生しました: {e}")

            st.divider()
            st.info("基準物質のアニオンのGibbs自由エネルギーを計算中...")
            # DB検索（基準物質・アニオン）
            if use_db_data_ref_anion:
                ref_g_anion = ref_anion_existing["g_tot"]
                st.success("DBから基準物質（アニオン）のG_totを取得しました。")
                st.write(f"基準物質（アニオン）G_tot: {ref_g_anion} Eh")
            else:
                handler_ref_anion = MoleculeHandler(ref_anion_smiles, input_type="smiles")
                inchi_str_ref_anion = inchi.MolToInchi(handler_ref_anion.mol)
                inchikey_str_ref_anion = inchi.InchiToInchiKey(inchi_str_ref_anion)

                if apply_force_field and handler_ref_anion.mol.GetNumAtoms() > 2:
                    handler_ref_anion.generate_conformers(num_conformers=num_conformers, forcefield=force_field)
                    handler_ref_anion.keep_lowest_energy_conformer()
                pyscf_input_ref_anion = handler_ref_anion.to_pyscf_input()
                result_ref_anion = compute_molecule_properties(
                    name=inchikey_str_ref_anion,
                    smiles=ref_anion_smiles,
                    xyz=pyscf_input_ref_anion,
                    theory=theory,
                    basis_set=basis,
                    opt_theory=theory,
                    opt_basis_set=basis,
                    charge=-1,
                    spin=0,
                    solvent_model=None if solvent_model == "None" else solvent_model,
                    eps=eps,
                    conv_params=conv_params,
                    maxsteps=maxsteps
                )
                ref_g_anion = extract_energy(result_ref_anion, "G_tot")
                show_imaginary_frequency_warning(result_ref_anion)
                st.write(f"基準物質（アニオン）G_tot: {ref_g_anion} Eh")

                frequencies = result_ref_anion['frequencies']
                thermo_info = result_ref_anion['thermo_info']

                if frequencies is not None:
                    st.success("✅ Calculation completed successfully!")

                    # 虚振動の数を計算（0未満の振動数の個数）
                    freq_array = result_ref_anion['frequencies']['freq_wavenumber']
                    # フラット化
                    if any(isinstance(f, (list, tuple, np.ndarray)) for f in freq_array):
                        import itertools
                        freq_array = list(itertools.chain.from_iterable(freq_array))
                    if hasattr(freq_array, "tolist"):
                        freq_array = freq_array.tolist()
                    # 複素数を実数に変換
                    freq_array = [f.real if isinstance(f, complex) else f for f in freq_array]
                    num_imaginary = sum(1 for f in freq_array if isinstance(f, (int, float)) and f < 0)

                    # データベースへ保存
                    try:
                        insert_molecule_with_frequencies(
                            inchi=inchi_str_ref_anion,
                            inchikey=inchikey_str_ref_anion,
                            g_tot=thermo_info["G_tot"][0],  
                            zpe=thermo_info["ZPE"][0],      
                            method=theory,
                            basis=basis,
                            charge=-1,
                            spin=0,
                            solvent=solvent_model if solvent_model not in [None, "None"] else None,
                            dielectric=eps,
                            temperature=298.15,
                            pressure=1.0,
                            frequencies=freq_array,
                            num_imaginary=num_imaginary,
                            db_path="data/energy_db.sqlite"
                        )
                        st.info("🗃️ 結果をデータベースに保存しました。")
                    except Exception as e:
                        st.warning(f"データベース保存時にエラーが発生しました: {e}")

            st.divider()

            if ref_g_neutral is None or ref_g_anion is None:
                st.error("基準物質のG_totが取得できませんでした。pKa計算を中止します。")
            else:
                delta_g_ref = ref_g_anion - ref_g_neutral

                for idx, selected_row in selected_rows:
                    try:
                        

                        # DB検索（サンプル中性分子）
                        st.write("サンプル分子のGibbs自由エネルギーを計算中...")

                        if use_db_data_sample:
                            G_mol_raw = {"G_tot": sample_existing["g_tot"]}
                            g_ha = sample_existing["g_tot"]
                            st.success(f"DBからサンプル分子（中性）のG_totを取得しました。")
                            st.write(f"ギブス自由エネルギー G_tot: {g_ha} Eh")
                        else:
                            handler_neutral = MoleculeHandler(smiles_input, input_type="smiles")
                            inchi_str_sample = inchi.MolToInchi(handler_neutral.mol)
                            inchikey_str_sample = inchi.InchiToInchiKey(inchi_str_sample)

                            st.text(f"中性分子: {smiles_input}")
                            if apply_force_field and handler_neutral.mol.GetNumAtoms() > 2:
                                handler_neutral.generate_conformers(num_conformers=num_conformers, forcefield=force_field)
                                handler_neutral.keep_lowest_energy_conformer()
                            pyscf_input_neutral = handler_neutral.to_pyscf_input()
                            G_mol_raw = compute_molecule_properties(
                                name=inchikey_str_sample,
                                smiles=smiles_input,
                                xyz=pyscf_input_neutral,
                                theory=theory, 
                                basis_set=basis,
                                opt_theory=theory,
                                opt_basis_set=basis,
                                charge=0,
                                spin=0,
                                solvent_model=None if solvent_model == "None" else solvent_model,
                                eps=eps,
                                conv_params=conv_params, 
                                maxsteps=maxsteps
                            )

                            show_imaginary_frequency_warning(G_mol_raw)
                            g_ha = extract_energy(G_mol_raw, "G_tot")
                            st.write(f"ギブス自由エネルギー G_tot: {g_ha} Eh")

                            # データベースへ保存
                            frequencies = G_mol_raw['frequencies']
                            thermo_info = G_mol_raw['thermo_info']      

                            if frequencies is not None:
                                st.success("✅ Calculation completed successfully!")

                                # 虚振動の数を計算（0未満の振動数の個数）
                                freq_array = G_mol_raw['frequencies']['freq_wavenumber']
                                # フラット化
                                if any(isinstance(f, (list, tuple, np.ndarray)) for f in freq_array):
                                    import itertools
                                    freq_array = list(itertools.chain.from_iterable(freq_array))
                                if hasattr(freq_array, "tolist"):
                                    freq_array = freq_array.tolist()
                                # 複素数を実数に変換
                                freq_array = [f.real if isinstance(f, complex) else f for f in freq_array]
                                num_imaginary = sum(1 for f in freq_array if isinstance(f, (int, float)) and f < 0)

                                # データベースへ保存
                                try:
                                    insert_molecule_with_frequencies(
                                        inchi=inchi_str_sample,
                                        inchikey=inchikey_str_sample,
                                        g_tot=thermo_info["G_tot"][0],  
                                        zpe=thermo_info["ZPE"][0],      
                                        method=theory,
                                        basis=basis,
                                        charge=0,
                                        spin=0,
                                        solvent=solvent_model if solvent_model not in [None, "None"] else None,
                                        dielectric=eps,
                                        temperature=298.15,
                                        pressure=1.0,
                                        frequencies=freq_array,
                                        num_imaginary=num_imaginary,
                                        db_path="data/energy_db.sqlite"
                                    )
                                    st.info("🗃️ 結果をデータベースに保存しました。")
                                except Exception as e:
                                    st.error(f"データベースへの保存中にエラーが発生しました: {e}")
                        st.divider()
                        st.write(f"🔄 [{idx}] {selected_row['deprotonated_smiles']} を計算中...")
                        # DB検索（サンプルアニオン）
                        anion_smiles = selected_row['deprotonated_smiles']
                        anion_name = f"{Chem.MolToInchiKey(Chem.MolFromSmiles(smiles_input))}_{selected_row['bond_index']}_anion"

                        use_db_data_sample_anion = use_db_data_anion_dict.get(idx, False)
                        sample_anion_existing = anion_existing_dict.get(idx, None)

                        if use_db_data_sample_anion and sample_anion_existing is not None:
                            result_anion_raw = {"G_tot": sample_anion_existing["g_tot"]}
                            g_a = sample_anion_existing["g_tot"]
                            st.success(f"DBからサンプル分子（アニオン）のG_totを取得しました。")
                            st.write(f"ギブス自由エネルギー G_tot: {g_a} Eh")
                        else:
                            st.text(f"アニオン: {anion_smiles}")
                            handler_anion = MoleculeHandler(anion_smiles, input_type="smiles")
                            inchi_str_sample_anion = inchi.MolToInchi(handler_anion.mol)
                            inchikey_str_sample_anion = inchi.InchiToInchiKey(inchi_str_sample_anion)

                            if apply_force_field and handler_anion.mol.GetNumAtoms() > 2:
                                handler_anion.generate_conformers(num_conformers=num_conformers, forcefield=force_field)
                                handler_anion.keep_lowest_energy_conformer()
                            pyscf_input_anion = handler_anion.to_pyscf_input()
                            result_anion_raw = compute_molecule_properties(
                                name=inchikey_str_sample_anion,
                                smiles=anion_smiles,
                                xyz=pyscf_input_anion,
                                theory=theory, 
                                basis_set=basis,
                                opt_theory=theory,
                                opt_basis_set=basis,
                                charge=-1,
                                spin=0,
                                solvent_model=None if solvent_model == "None" else solvent_model,
                                eps=eps,
                                conv_params=conv_params, 
                                maxsteps=maxsteps
                            )
                            show_imaginary_frequency_warning(result_anion_raw)
                            st.write(f"データを保存しました: data/{anion_name}")
                            g_a = extract_energy(result_anion_raw, "G_tot")
                            st.write(f"ギブス自由エネルギー G_tot: {g_a} Eh")

                            # データベースへ保存
                            frequencies = result_anion_raw['frequencies']
                            thermo_info = result_anion_raw['thermo_info']

                            if frequencies is not None:
                                st.success("✅ Calculation completed successfully!")

                                # 虚振動の数を計算（0未満の振動数の個数）
                                freq_array = result_anion_raw['frequencies']['freq_wavenumber']
                                # フラット化
                                if any(isinstance(f, (list, tuple, np.ndarray)) for f in freq_array):
                                    import itertools
                                    freq_array = list(itertools.chain.from_iterable(freq_array))
                                if hasattr(freq_array, "tolist"):
                                    freq_array = freq_array.tolist()
                                # 複素数を実数に変換
                                freq_array = [f.real if isinstance(f, complex) else f for f in freq_array]
                                num_imaginary = sum(1 for f in freq_array if isinstance(f, (int, float)) and f < 0)

                                # データベースへ保存
                                try:
                                    insert_molecule_with_frequencies(
                                        inchi=inchi_str_sample_anion,
                                        inchikey=inchikey_str_sample_anion,
                                        g_tot=thermo_info["G_tot"][0],  
                                        zpe=thermo_info["ZPE"][0],      
                                        method=theory,
                                        basis=basis,
                                        charge=-1,
                                        spin=0,
                                        solvent=solvent_model if solvent_model not in [None, "None"] else None,
                                        dielectric=eps,
                                        temperature=298.15,
                                        pressure=1.0,
                                        frequencies=freq_array,
                                        num_imaginary=num_imaginary,
                                        db_path="data/energy_db.sqlite"
                                    )
                                    st.info("🗃️ 結果をデータベースに保存しました。")
                                except Exception as e:
                                    st.error(f"データベースへの保存中にエラーが発生しました: {e}")

                            st.divider()


                        # pKa計算
                        if g_ha is not None and g_a is not None:
                            delta_g_sample = g_a - g_ha
                            pka = calculate_relative_pka(delta_g_sample, delta_g_ref, pka_ref)
                            pka_results.append({
                                "bond_index": selected_row["bond_index"],
                                "bond_atoms": selected_row["bond_atoms"],
                                "deprotonated_smiles": anion_smiles,
                                "ΔG_sample (Hartree)": delta_g_sample,
                                "ΔG_ref (Hartree)": delta_g_ref,
                                "pKa": pka
                            })
                            st.success(f"[{idx}] 計算完了: pKa = {pka:.2f}")
                            st.write(f"g_ha (中性分子 G_tot): {g_ha}")
                            st.write(f"g_a (アニオン G_tot): {g_a}")
                            st.write(f"delta_g_sample: {delta_g_sample}")
                            st.write(f"delta_g_ref: {delta_g_ref}")
                            st.write(f"pKa_ref: {pka_ref}")
                            st.write(f"pKa (計算値): {pka}")
                        else:
                            st.warning(f"[{idx}] G_totが取得できませんでした。")
                    except Exception as e:
                        st.error(f"[{idx}] の計算中にエラーが発生しました: {e}")
                if pka_results:
                    st.markdown("### 📊 pKa計算結果")
                    st.dataframe(pd.DataFrame(pka_results))
