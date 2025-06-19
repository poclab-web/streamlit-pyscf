"""
結合乖離エネルギーの計算
結合を選択して、その結合を切断し、二つのラジカルを計算して値を出力する。
"""

import textwrap
import ast

import streamlit as st
import pandas as pd
import numpy as np
from rdkit import Chem

from config.config import conv_preset_values

from logic.calculation import theory_options, basis_set_options, run_quantum_calculation
from logic.bde_calculation import (
    get_fragment_dataframe,
    compute_neutral_molecule_properties,
    compute_radical_fragment_properties,
    calculate_bde_from_gibbs
)
from logic.molecule_handler import MoleculeHandler
from logic.visualization import show_imaginary_frequency_warning
from logic.output_handler import safe_dict, safe_dict

def parse_thermo_info(obj):
    """thermo_info を dict に変換（文字列だった場合のみパース）"""
    if isinstance(obj, str):
        try:
            return ast.literal_eval(obj)
        except Exception as e:
            raise ValueError(f"thermo_info のパースに失敗: {e}")
    elif isinstance(obj, dict):
        return obj
    else:
        raise TypeError("thermo_info は dict もしくは str である必要があります。")

def extract_energy(g_dict, key="G_tot", fallback_key="E_scf_only"):
    """
    thermo 辞書またはその中の 'thermo_info' からエネルギー値を取得。
    なければ fallback_key を見る。両方なければエラー。
    """
    value = g_dict.get(key)

    # トップレベルに無ければ thermo_info を見る
    if value is None and "thermo_info" in g_dict:
        thermo = g_dict["thermo_info"]
        value = thermo.get(key)
        if value is None and fallback_key:
            value = thermo.get(fallback_key)

    elif value is None and fallback_key:
        value = g_dict.get(fallback_key)

    if isinstance(value, tuple):
        return value[0]
    elif isinstance(value, (float, int)):
        return float(value)
    else:
        raise ValueError(f"{key}（または {fallback_key}）の形式が不正です。値: {value}")

EH_TO_KCAL = 627.509  # 1 Eh = 627.509 kcal/mol

def compute_all_bdes(G_mol_raw, G_f1_raw, G_f2_raw):
    thermo_mol = parse_thermo_info(G_mol_raw["thermo_info"])
    thermo_f1 = parse_thermo_info(G_f1_raw["thermo_info"])
    thermo_f2 = parse_thermo_info(G_f2_raw["thermo_info"])

    EH_TO_KCAL = 627.509

    # 中性分子（必ず補正あり）
    # E_mol = extract_energy(thermo_mol, "E_tot")
    # ZPE_mol = extract_energy(thermo_mol, "ZPE")
    H_mol = extract_energy(thermo_mol, "H_tot")
    G_mol = extract_energy(thermo_mol, "G_tot")

    # ラジカル（ZPEなどがない場合 fallback_key を使う）
    # E_f1 = extract_energy(thermo_f1, "E_tot", "E_scf_only")
    # ZPE_f1 = extract_energy(thermo_f1, "ZPE", None) or 0.0
    H_f1 = extract_energy(thermo_f1, "H_tot", "E_scf_only")
    G_f1 = extract_energy(thermo_f1, "G_tot", "E_scf_only")

    # E_f2 = extract_energy(thermo_f2, "E_tot", "E_scf_only")
    # ZPE_f2 = extract_energy(thermo_f2, "ZPE", None) or 0.0
    H_f2 = extract_energy(thermo_f2, "H_tot", "E_scf_only")
    G_f2 = extract_energy(thermo_f2, "G_tot", "E_scf_only")

    # BDEを計算（ZPE補正は一部ゼロになる）
    # bde_zpe = (E_f1 + ZPE_f1 + E_f2 + ZPE_f2 - E_mol - ZPE_mol) * EH_TO_KCAL
    bde_h    = (H_f1 + H_f2 - H_mol) * EH_TO_KCAL
    bde_g    = (G_f1 + G_f2 - G_mol) * EH_TO_KCAL

    return {
        # "ZPE-corrected BDE": round(bde_zpe, 2),
        "Enthalpy-based BDE (ΔH)": round(bde_h, 2),
        "Bond-Dissociation Free Energy (BDFE)": round(bde_g, 2)
    }


# --- UI部分 ---
solvent_models = ["None", "PCM", "ddCOSMO"]
solvents_file = "config/solvents_epsilon.csv"
solvents_data = pd.read_csv(solvents_file)

st.title("BDE Caluculator")
st.markdown("**SMILESから結合を選択し、量子化学計算条件を指定してBDE（結合解離エネルギー）を計算**")

# ▼ SMILES入力
smiles_input = st.text_input("SMILESを入力してください", "CO")

# ▼ フラグメント生成
if st.button("フラグメントを生成"):
    with st.spinner("フラグメントを生成しています..."):
        df_frag = get_fragment_dataframe([smiles_input])

        # frag1とfrag2のSMILESが同じ行を削除（完全一致の重複）
        df_frag_unique = df_frag.drop_duplicates(subset=["fragment1", "fragment2"])

        st.session_state["fragments"] = df_frag_unique
        st.success("フラグメントを生成しました。下から選択してください。")


# ▼ 計算条件の選択
st.header("🔧 配座生成と計算設定")

apply_force_field = st.checkbox("分子力場による構造最適化を行う", value=True)

# Choose force field
force_field = st.selectbox("Select Force Field", ["MMFF", "UFF"])

# Number of conformers to generate
num_conformers = st.number_input("Number of Conformers", value=1000)

# ▼ 計算条件の選択
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
st.subheader("Convergence Parameters")

# プリセット選択（loose / normal / tight）
preset = st.radio("Choose preset", ["Loose", "Normal", "Tight"], index=1, horizontal=True)

# 各プリセットに応じた初期値を設定
conv_preset_values = {
    "Loose":  {"energy": 1.0e-4, "grms": 1.0e-3, "gmax": 3.0e-3, "drms": 4.0e-3, "dmax": 6.0e-3},
    "Normal": {"energy": 1.0e-5, "grms": 5.0e-4, "gmax": 1.5e-3, "drms": 2.0e-3, "dmax": 3.0e-3},
    "Tight":  {"energy": 1.0e-6, "grms": 3.0e-4, "gmax": 1.2e-3, "drms": 1.2e-3, "dmax": 1.8e-3},
}
vals = conv_preset_values[preset]

# 数値入力（プリセット値を初期値に）
with st.expander("Manual adjustment"):
    convergence_energy = st.number_input(
        "Energy Tolerance (Hartree)", 
        min_value=1e-7, value=vals["energy"], step=1e-5, format="%.7f"
    )
    convergence_grms = st.number_input(
        "Gradient RMS Tolerance (Eh/Bohr)", 
        min_value=1e-5, value=vals["grms"], step=1e-5, format="%.5f"
    )
    convergence_gmax = st.number_input(
        "Gradient Max Tolerance (Eh/Bohr)", 
        min_value=1e-5, value=vals["gmax"], step=1e-5, format="%.5f"
    )
    convergence_drms = st.number_input(
        "Displacement RMS Tolerance (Angstrom)", 
        min_value=1e-4, value=vals["drms"], step=1e-4, format="%.4f"
    )
    convergence_dmax = st.number_input(
        "Displacement Max Tolerance (Angstrom)", 
        min_value=1e-4, value=vals["dmax"], step=1e-4, format="%.4f"
    )
    maxsteps = st.number_input(
        "Max Iterations", 
        min_value=1, value=100, step=1
    )

# 辞書にまとめる
conv_params = {
    "convergence_energy": convergence_energy,
    "convergence_grms": convergence_grms,
    "convergence_gmax": convergence_gmax,
    "convergence_drms": convergence_drms,
    "convergence_dmax": convergence_dmax,
    "maxsteps": maxsteps,
}

# フラグメント選択とBDE計算
if "fragments" in st.session_state:

    st.subheader("🔍 分解候補を選択してください（複数可）")

    handler = MoleculeHandler(smiles_input, input_type="smiles")

    # 🔘 画像表示のオプション
    show_image = st.checkbox("🖼️ Bond Index画像を表示する", value=True)

    # 💠 画像生成と表示（オプション）
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

    df_frag = st.session_state["fragments"]   
    selected_rows = []

    for i, row in df_frag.iterrows():
        frag_label = f"Bond {row['bond_index']} ({row['bond_type']}): {row['fragment1']} + {row['fragment2']}"
        if st.checkbox(frag_label, key=f"frag_{i}"):
            selected_rows.append((i, row))

    if selected_rows:
        st.markdown("### 選択された分解候補")
        selected_df = pd.DataFrame([r for _, r in selected_rows])
        st.dataframe(selected_df)

        if st.button("選択された全フラグメントでBDEを一括計算"):
            if apply_force_field:
                method_text = f"""
                **Computational Details**  
                Molecular structures were generated from SMILES inputs using RDKit [1], and 3D conformers were generated using the ETKDG method.  
                A total of {num_conformers} conformers were generated for each molecule and optimized using the {force_field} force field.  
                The lowest-energy conformer according to the {force_field} force field was selected for subsequent quantum chemical geometry optimization.  
                Quantum chemical optimizations were performed at the **{theory}/{basis}** level using PySCF [2].  
                {"No solvent model was applied." if solvent_model == "None" else f"The solvent effect was considered using the {solvent_model} model with a dielectric constant of ε = {eps}."}  
                Geometry optimizations were conducted with a convergence threshold of {convergence_energy:.1e} Hartree (energy), {convergence_gmax:.1e} Eh/Bohr (max gradient), and a maximum of {maxsteps} optimization steps.  
                All calculations were performed using a custom Streamlit-based interface [3] integrating RDKit and PySCF.
                """
            else:
                method_text = f"""
                **Computational Details**  
                Molecular structures were generated from SMILES inputs using RDKit [1], and 3D coordinates were generated.  
                These structures were directly used as initial geometries for quantum chemical optimization at the **{theory}/{basis}** level using PySCF [2].  
                {"No solvent model was applied." if solvent_model == "None" else f"The solvent effect was considered using the {solvent_model} model with a dielectric constant of ε = {eps}."}  
                Geometry optimizations were conducted with a convergence threshold of {convergence_energy:.1e} Hartree (energy), {convergence_gmax:.1e} Eh/Bohr (max gradient), and a maximum of {maxsteps} optimization steps.  
                All calculations were performed using a custom Streamlit-based interface [3] integrating RDKit and PySCF.
                """

            references = """
            **References**  
            [1] Landrum, G. RDKit: Open-source cheminformatics. [https://www.rdkit.org](https://www.rdkit.org)  
            [2] Sun, Q. *et al.* PySCF: The Python-based Simulations of Chemistry Framework. **WIREs Comput Mol Sci** *2018*, **8**, e1340. DOI: [10.1002/wcms.1340](https://doi.org/10.1002/wcms.1340)  
            [3] PocLab streamlit-pyscf: Quantum chemistry web interface. [https://github.com/poclab-web/streamlit-pyscf](https://github.com/poclab-web/streamlit-pyscf)
            """

            # Streamlit表示部分
            st.markdown("### 📄 Method")
            st.markdown(method_text)
            st.markdown("---")
            st.markdown(references)

            bde_results = []

            for idx, selected_row in selected_rows:
                frag1 = selected_row["fragment1"]
                frag2 = selected_row["fragment2"]
                bond_type = selected_row["bond_type"]

                try:
                    st.write(f"🔄 [{idx}] {frag1} + {frag2} のBDEを計算中...")

                    # 中性分子の計算
                    compound_name = Chem.MolToInchiKey(handler.mol)
                    # 原子数が2つ以上のときだけ最適化（[H]などはスキップ）
                    if apply_force_field:
                        if handler.mol.GetNumAtoms() > 2:
                            handler.generate_conformers(
                                num_conformers=num_conformers,
                                forcefield=force_field
                            )
                            handler.keep_lowest_energy_conformer()
                    pyscf_input_G = handler.to_pyscf_input()
                    G_mol_raw = compute_neutral_molecule_properties(
                        compound_name, smiles_input, pyscf_input_G,
                        theory=theory, basis_set=basis,
                        solvent_model=None if solvent_model == "None" else solvent_model,
                        eps=eps, conv_params=conv_params, maxsteps=maxsteps
                    )
                    with st.expander(f"🧪 [{idx}] 中性分子 G_mol_raw"):
                        st.json(safe_dict(G_mol_raw))

                    show_imaginary_frequency_warning(G_mol_raw)
                    st.write(f"データを保存しました: data/{compound_name}")
                    st.write(f"エンタルピー H_tot: {extract_energy(G_mol_raw, 'H_tot')} Eh")
                    st.write(f"ギブス自由エネルギー G_tot: {extract_energy(G_mol_raw, 'G_tot')} Eh")

                    # ラジカル1
                    handler_f1 = MoleculeHandler(frag1, input_type="smiles")

                    if apply_force_field:
                        if handler_f1.mol.GetNumAtoms() > 1:
                            handler_f1.generate_conformers(
                                num_conformers=num_conformers,
                                forcefield=force_field
                            )
                            handler_f1.keep_lowest_energy_conformer()

                    pyscf_input_f1 = handler_f1.to_pyscf_input()
                    frag1_name = f"{compound_name}_{selected_row['bond_index']}_frag1"

                    G_f1_raw = compute_radical_fragment_properties(
                        frag1_name, frag1, pyscf_input_f1,
                        theory=theory, basis_set=basis, charge=0, spin=1,
                        solvent_model=None if solvent_model == "None" else solvent_model,
                        eps=eps, conv_params=conv_params, maxsteps=maxsteps
                    )
                    with st.expander(f"🧪 [{idx}] ラジカル1 G_f1_raw"):
                        st.json(safe_dict(G_f1_raw))

                    show_imaginary_frequency_warning(G_f1_raw)
                    st.write(f"データを保存しました: data/{frag1_name}")
                    st.write(f"エンタルピー H_tot: {extract_energy(G_f1_raw, 'H_tot')} Eh")
                    st.write(f"ギブス自由エネルギー G_tot: {extract_energy(G_f1_raw, 'G_tot')} Eh")

                    # ラジカル2
                    handler_f2 = MoleculeHandler(frag2, input_type="smiles")
                    
                    if apply_force_field:
                        if handler_f2.mol.GetNumAtoms() > 1:
                            handler_f2.generate_conformers(
                                num_conformers=num_conformers,
                                forcefield=force_field
                            )           
                            handler_f2.keep_lowest_energy_conformer()

                    pyscf_input_f2 = handler_f2.to_pyscf_input()
                    frag2_name = f"{compound_name}_{selected_row['bond_index']}_frag2"

                    G_f2_raw = compute_radical_fragment_properties(
                        frag2_name, frag2, pyscf_input_f2,
                        theory=theory, basis_set=basis, charge=0, spin=1,
                        solvent_model=None if solvent_model == "None" else solvent_model,
                        eps=eps, conv_params=conv_params, maxsteps=maxsteps
                    )
                    with st.expander(f"🧪 [{idx}] ラジカル2 G_f2_raw"):
                        st.json(safe_dict(G_f2_raw))

                    show_imaginary_frequency_warning(G_f2_raw)
                    st.write(f"データを保存しました: data/{frag2_name}")
                    st.write(f"エンタルピー H_tot: {extract_energy(G_f2_raw, 'H_tot')} Eh")
                    st.write(f"ギブス自由エネルギー G_tot: {extract_energy(G_f2_raw, 'G_tot')} Eh")
              
                    # BDEを3種類計算
                    bde_all = compute_all_bdes(G_mol_raw, G_f1_raw, G_f2_raw)

                    # 結果にすべて追加
                    bde_results.append({
                        "bond_index": selected_row["bond_index"],
                        "bond_type": bond_type,
                        "fragment1": frag1,
                        "fragment2": frag2,
                        **bde_all  # ← 3つのBDEが自動展開される
                    })
                    st.write("データを追加しました。")
                    st.dataframe(pd.DataFrame(bde_results))
   
                except Exception as e:
                    st.error(f"[{idx}] の計算中にエラーが発生しました: {e}")

            if bde_results:
                st.success("✅ BDE計算が完了しました。")

                bde_description = """
                **Bond Dissociation Energy (BDE) Calculation**  
                BDEs were calculated as the difference in Gibbs free energy between the parent molecule and its radical fragments:

                  **BDE = G(frag1) + G(frag2) − G(parent)**

                Closed-shell species (spin = 0) were treated with restricted methods (RHF or RKS), and open-shell species with unrestricted methods (UHF or UKS).  
                Gibbs free energies (G) were obtained from thermochemical analysis at 298.15 K and 1 atm, including electronic, ZPE, thermal, and entropy contributions.  
                All optimized structures were confirmed as true minima by verifying the absence of imaginary frequencies.

                For the hydrogen radical (H•), G was estimated using a single-point energy calculation without geometry optimization or frequency analysis.  

                All calculations were performed under the same level of theory and basis set.
                """
                st.markdown("### 📊 BDE Calculation Details")

                st.markdown(bde_description)
