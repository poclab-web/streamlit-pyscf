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
from io import BytesIO

from logic.molecule_handler import MoleculeHandler
from logic.calculation import compute_molecule_properties, theory_options, basis_set_options
from config.config import conv_preset_values, solvent_models, solvents_data, solvents_file
from logic.visualization import show_imaginary_frequency_warning
from logic.output_handler import safe_dict, extract_energy

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

def calculate_pka(delta_g_hartree):
    delta_g_kcal = delta_g_hartree * HARTREE_TO_KCAL
    pka = delta_g_kcal / (2.303 * R * T / 4184)
    return pka

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
st.markdown("**SMILESから脱プロトン化候補を選択し、量子化学計算条件を指定してpKaを計算**")

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
    selected_pka = st.selectbox("Select a reference compound", pka_df["name"].tolist())
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

with st.expander("Show selected method and reference"):
    st.markdown("### 🧪 Method for pKa Calculation")
    st.markdown(f"""
    **Computational Details**  
    Molecular structures were generated from SMILES using RDKit [1], and 3D conformers were generated using the ETKDG method.  
    {num_conformers} conformers were generated and optimized using the {force_field} force field.  
    The lowest-energy conformer was selected for quantum chemical geometry optimization at the **{theory}/{basis}** level using PySCF [2].  
    {"No solvent model was applied." if solvent_model == "None" else f"The solvent effect was considered using the {solvent_model} model with ε = {eps}."}  
    Proton dissociation free energies (ΔG_deprot) were computed as the difference in Gibbs free energy between the conjugate base (anion) and the neutral acid, corrected by the free energy of a proton in solution (−270.3 kcal/mol at 298.15 K).  
    The pKa values were estimated using the thermodynamic relation:  
    $$ \\mathrm{{p}}K_a = \\frac{{\\Delta G_{{\\text{{deprot}}}}}}{{2.303RT}} $$
    """)
    st.markdown("---")
    st.markdown("""
    **References**  
    [1] Landrum, G. RDKit: Open-source cheminformatics. [https://www.rdkit.org](https://www.rdkit.org)  
    [2] Sun, Q. *et al.* PySCF: The Python-based Simulations of Chemistry Framework. **WIREs Comput Mol Sci** *2018*, **8**, e1340. DOI: [10.1002/wcms.1340](https://doi.org/10.1002/wcms.1340)  
    [3] PocLab streamlit-pyscf: Quantum chemistry web interface. [https://github.com/poclab-web/streamlit-pyscf](https://github.com/poclab-web/streamlit-pyscf)
    """)

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
    st.subheader("🔍 脱プロトン化候補を選択してください（複数可）")
    df_bonds = st.session_state["deprot_candidates"]
    selected_rows = []
    for i, row in df_bonds.iterrows():
        label = f"Bond {row['bond_index']} ({row['bond_atoms']}): [{row['deprotonated_smiles']}]"
        if st.checkbox(label, key=f"bond_{i}"):
            selected_rows.append((i, row))
    if selected_rows:
        st.markdown("### 選択された脱プロトン化候補")
        selected_df = pd.DataFrame([r for _, r in selected_rows])
        st.dataframe(selected_df)
        if st.button("選択された全フラグメントでpKaを一括計算"):
            pka_results = []

            # ▼ 基準物質のGibbs自由エネルギー自動計算
            st.info("基準物質のGibbs自由エネルギーを計算中...")
            ref_smiles = selected_row['reference_smiles']
            ref_anion_smiles = selected_row['reference_anion_smiles']
            ref_name = selected_row['name']
            pka_ref = selected_row['pKa_exp']

            # 中性基準物質
            handler_ref = MoleculeHandler(ref_smiles, input_type="smiles")
            if apply_force_field and handler_ref.mol.GetNumAtoms() > 2:
                handler_ref.generate_conformers(num_conformers=num_conformers, forcefield=force_field)
                handler_ref.keep_lowest_energy_conformer()
            pyscf_input_ref = handler_ref.to_pyscf_input()
            result_ref = compute_molecule_properties(
                name=f"{ref_name}_neutral",
                smiles=ref_smiles,
                xyz=pyscf_input_ref,
                theory=theory,
                basis_set=basis,
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

            # アニオン基準物質
            handler_ref_anion = MoleculeHandler(ref_anion_smiles, input_type="smiles")
            if apply_force_field and handler_ref_anion.mol.GetNumAtoms() > 2:
                handler_ref_anion.generate_conformers(num_conformers=num_conformers, forcefield=force_field)
                handler_ref_anion.keep_lowest_energy_conformer()
            pyscf_input_ref_anion = handler_ref_anion.to_pyscf_input()
            result_ref_anion = compute_molecule_properties(
                name=f"{ref_name}_anion",
                smiles=ref_anion_smiles,
                xyz=pyscf_input_ref_anion,
                theory=theory,
                basis_set=basis,
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

            if ref_g_neutral is None or ref_g_anion is None:
                st.error("基準物質のG_totが取得できませんでした。pKa計算を中止します。")
            else:
                delta_g_ref = ref_g_anion - ref_g_neutral

                for idx, selected_row in selected_rows:
                    try:
                        st.write(f"🔄 [{idx}] {selected_row['deprotonated_smiles']} のpKaを計算中...")
                        # 中性分子
                        handler_neutral = MoleculeHandler(smiles_input, input_type="smiles")
                        compound_name = Chem.MolToInchiKey(handler.mol)
                        st.text(f"中性分子: {smiles_input}")
                        if apply_force_field and handler_neutral.mol.GetNumAtoms() > 2:
                            handler_neutral.generate_conformers(num_conformers=num_conformers, forcefield=force_field)
                            handler_neutral.keep_lowest_energy_conformer()
                        pyscf_input_neutral = handler_neutral.to_pyscf_input()
                        G_mol_raw = compute_molecule_properties(
                            name=compound_name,
                            smiles=smiles_input,
                            xyz=pyscf_input_neutral,
                            theory=theory, 
                            basis_set=basis,
                            charge=0,
                            spin=0,
                            solvent_model=None if solvent_model == "None" else solvent_model,
                            eps=eps,
                            conv_params=conv_params, 
                            maxsteps=maxsteps
                        )

                        show_imaginary_frequency_warning(G_mol_raw)
                        st.write(f"データを保存しました: data/{compound_name}")
                        st.write(f"ギブス自由エネルギー G_tot: {extract_energy(G_mol_raw, 'G_tot')} Eh")

                        # アニオン
                        anion_smiles = selected_row['deprotonated_smiles']
                        anion_name = f"{compound_name}_{selected_row['bond_index']}_anion"
                        st.text(f"アニオン: {anion_smiles}")
                        handler_anion = MoleculeHandler(anion_smiles, input_type="smiles")
                        if apply_force_field and handler_anion.mol.GetNumAtoms() > 2:
                            handler_anion.generate_conformers(num_conformers=num_conformers, forcefield=force_field)
                            handler_anion.keep_lowest_energy_conformer()
                        pyscf_input_anion = handler_anion.to_pyscf_input()
                        result_anion_raw = compute_molecule_properties(
                            name=anion_name,
                            smiles=anion_smiles,
                            xyz=pyscf_input_anion,
                            theory=theory, 
                            basis_set=basis,
                            charge=-1,
                            spin=0,
                            solvent_model=None if solvent_model == "None" else solvent_model,
                            eps=eps,
                            conv_params=conv_params, 
                            maxsteps=maxsteps
                        )

                        show_imaginary_frequency_warning(result_anion_raw)
                        st.write(f"データを保存しました: data/{anion_name}")
                        st.write(f"ギブス自由エネルギー G_tot: {extract_energy(result_anion_raw, 'G_tot')} Eh")

                        # extract_energyでfloat値を取得
                        g_ha = extract_energy(G_mol_raw, "G_tot")
                        g_a = extract_energy(result_anion_raw, "G_tot")

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
