"""
UVの予測

"""

import streamlit as st
from utils.module import load_css
import os

from pyscf import gto, dft, tdscf

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, inchi
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
import pandas as pd

from config.config import conv_preset_values, solvent_models, solvents_data, solvents_file
from logic.molecule_handler import MoleculeHandler
from logic.calculation import compute_molecule_properties, theory_options, basis_set_options
from logic.database import insert_molecule_with_frequencies, get_molecule_from_sqlite, get_molecules_from_sqlite, get_stable_molecule_from_sqlite
from logic.visualization import plot_uv_spectrum, find_allowed_excited_state_from_lists, prepare_excited_states_table

from controllers.excited_state_calculation import calculate_excited_state


# カスタムCSSを適用
load_css("config/styles.css")

# Streamlitアプリ
st.title("UVスペクトル予測アプリ with PySCF")
st.write("分子のSMILES表記を入力し、UVスペクトルを予測します。")
st.divider()

# ユーザー入力の設定
st.header("Molecule Input")
# ▼ SMILES入力
smiles_input = st.text_input("計算したい分子のSMILESを入力してください", "C=O")
charge = st.number_input("Molecular Charge", min_value=-10, max_value=10, value=0, step=1)
multiplicity = st.number_input("Spin Multiplicity (2S + 1)", min_value=1, max_value=10, value=1, step=1)
spin = multiplicity - 1

st.header("📐 量子化学計算設定")
theory = st.selectbox("理論レベル", theory_options, index=theory_options.index("HF"), key="theory_selector")
basis_set = st.selectbox("基底関数", basis_set_options, index=basis_set_options.index("sto-3g"), key="basis_selector")
solvent_model = st.selectbox("Select Solvent Model", solvent_models)
eps = None
if solvent_model in ["PCM", "ddCOSMO"]:
    solvent_selection = st.selectbox(
        "Select a solvent",
        [f"{row['Solvent']} (ε={row['Epsilon']})" for _, row in solvents_data.iterrows()]
    )
    if solvent_selection:
        eps = float(solvent_selection.split("=", 1)[-1][:-1])

# 計算の設定入力
with st.expander("Setting calculation conditions for optimization"):
    st.header("🔧 配座生成と計算設定")
    apply_force_field = st.checkbox("分子力場による構造最適化を行う", value=True)
    force_field = st.selectbox("Select Force Field", ["MMFF", "UFF"])
    num_conformers = st.number_input("Number of Conformers", value=1000)

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
handler = MoleculeHandler(smiles_input, input_type="smiles")
mol = handler.mol
smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
inchi_str = inchi.MolToInchi(mol)
inchikey_str = inchi.InchiToInchiKey(inchi_str)
compound_name = inchikey_str
directory = os.path.join("data", compound_name)
os.makedirs(directory, exist_ok=True)

# 虚振動なしの安定な分子の取得
stable_existing = get_stable_molecule_from_sqlite(
    inchikey=inchikey_str,
    method=theory,
    basis=basis_set,
    spin=spin,
    charge=charge,
    solvent=solvent_model if solvent_model not in [None, "None"] else None,
    dielectric=eps,
    temperature=298.15,
    pressure=1.0,
    db_path="data/energy_db.sqlite"
)

existing = stable_existing

use_db_data = False
if existing is not None:
    st.info("🗃️ 同じ条件の基底状態の安定なデータが既にデータベースに存在します。")
    use_db_data = st.radio(
        "既存のデータベースのデータを使いますか？",
        ("使う（計算をスキップ）", "使わず新たに計算する"),
        horizontal=True,
        index=0,
        key="use_db_data_radio"
    ) == "使う（計算をスキップ）"

if use_db_data:
    st.success("既存のデータを使用します。計算をスキップします。")
    existing = stable_existing

# nstatesの入力
st.header("励起状態の計算設定")
nstates = st.number_input("励起状態の数",
    min_value=1,
    max_value=100,
    value=10,
    step=1,
    help="計算する励起状態の数を指定します。励起状態の数が多いほど計算時間が長くなります。"
)

# singlet / triplet
excited_spin = st.selectbox(
    "励起状態のスピン",
    ["singlet", "triplet"],
    index=0,
    help="励起状態のスピン多重度の計算を一重項、三重項かを選択できます。"
)

# tda近侍を使用するかどうか
tda = st.checkbox(
    "TDA近似を使用する",
    value=True,
    help="TDA近似を使用して遷移エネルギーを計算します。"
)   

# 励起状態の計算がDBに存在するか確認
existing_excited = get_molecules_from_sqlite(
    inchikey=inchikey_str,
    method=theory,
    basis=basis_set,
    spin=spin,
    charge=charge,
    solvent=solvent_model if solvent_model not in [None, "None"] else None,
    dielectric=eps,
    temperature=298.15,
    pressure=1.0,
    nstates=nstates,
    excited_spin=excited_spin,
    tda=tda,
    db_path="data/energy_db.sqlite"
)

selected_excited = None
if existing_excited and len(existing_excited) > 0:
    st.info(f"🗃️ 同じ条件の励起状態データが{len(existing_excited)}件データベースに存在します。")
    # 複数件の場合は選択肢を表示
    if len(existing_excited) == 1:
        selected_excited = existing_excited[0]
    else:
        # 例: idと登録日時で選択肢を作る
        options = [
            f"ID: {d['id']} / 登録日時: {d.get('timestamp', '不明')}" for d in existing_excited
        ]
        idx = st.selectbox("使用するデータを選択してください", range(len(options)), format_func=lambda i: options[i])
        selected_excited = existing_excited[idx]

        with st.expander("DBから取得したデータの詳細", expanded=False):
            st.write("DBから取得したデータ:", selected_excited)

    use_db_data_excited = st.radio(
        "既存のデータベースのデータを使いますか？",
        ("使う（計算をスキップ）", "使わず新たに計算する"),
        horizontal=True,
        index=0,
        key="use_db_data_excited_radio"
    ) == "使う（計算をスキップ）"
else:
    use_db_data_excited = False

# 計算の方法を出力
with st.expander("UVスペクトル予測の計算方法と参考文献を表示", expanded=False):
    st.markdown("### 🧪 Method for UV Spectrum Prediction")
    st.markdown(
        "**Computational Details**  \n"
        "Molecular structures were generated from SMILES using RDKit [1], and 3D conformers were generated using the ETKDG method.  \n"
        f"{num_conformers} conformers were generated and optimized using the {force_field} force field.  \n"
        f"The lowest-energy conformer was selected for quantum chemical geometry optimization at the **{theory}/{basis_set}** level using PySCF [2].  \n"
        f"{'No solvent model was applied.' if solvent_model == 'None' else f'The solvent effect was considered using the {solvent_model} model with ε = {eps}.'}  \n"
        "Excited-state calculations were performed using time-dependent DFT (TDDFT) or HF (TDHF) with "
        f"{nstates} excited states ({excited_spin.replace('_', '/')}), "
        f"{'with' if tda else 'without'} TDA approximation.  \n"
        "Excitation energies and oscillator strengths were used to simulate the UV/Vis absorption spectrum.  \n" \
        "All calculations were performed using a custom Streamlit-based interface [3] integrating RDKit and PySCF."
    )
    st.markdown("---")
    st.markdown(
        "**References**  \n"
        "[1] Landrum, G. RDKit: Open-source cheminformatics. [https://www.rdkit.org](https://www.rdkit.org)  \n"
        "[2] Sun, Q. *et al.* PySCF: The Python-based Simulations of Chemistry Framework. **WIREs Comput Mol Sci** *2018*, **8**, e1340. DOI: [10.1002/wcms.1340](https://doi.org/10.1002/wcms.1340)  \n"
        "[3] PocLab streamlit-pyscf: Quantum chemistry web interface. [https://github.com/poclab-web/streamlit-pyscf](https://github.com/poclab-web/streamlit-pyscf)  \n"
        "[4] Bauernschmitt, R.; Ahlrichs, R. Treatment of electronic excitations within the adiabatic approximation of time dependent density functional theory. **Chem. Phys. Lett.** *1996*, **256**, 454–464. DOI: [10.1016/0009-2614(96)00440-X](https://doi.org/10.1016/0009-2614(96)00440-X)"
    )

# --- ここから追加: 既存データで可視化 ---
if use_db_data_excited and selected_excited is not None:

    excitation_energies = selected_excited.get("excited_energies", None)
    oscillator_strengths = selected_excited.get("oscillator_strengths", None)

    if excitation_energies is None or oscillator_strengths is None:
        st.error("データベースに励起エネルギーまたは振動子強度が保存されていません。")
        st.stop()

    if type(excitation_energies) == str:
        try:
            excitation_energies = eval(excitation_energies)  # 文字列をリストに変換
        except Exception as e:
            st.write(type(excitation_energies))
    
    if type(oscillator_strengths) == str:
        try:
            oscillator_strengths = eval(oscillator_strengths)  # 文字列をリ
        except Exception as e:
            st.write(type(oscillator_strengths))

    # 波長と振動子強度の計算
    wavelengths = []
    energies_list = excitation_energies

    if isinstance(excitation_energies, (list, np.ndarray)) and len(excitation_energies) > 0 and isinstance(excitation_energies[0], (list, np.ndarray)):
        energies_list = excitation_energies[0]  # 二重リストの場合は一次元に

    for i, energy in enumerate(energies_list):
        try:
            energy_float = float(energy)
            wavelength = 1240 / (energy_float * 27.2114)  # エネルギー (au) を nm に変換
        except Exception as e:
            wavelength = None
        osc = oscillator_strengths[i] if i < len(oscillator_strengths) else 0
        wavelengths.append((wavelength, osc))

    # データ処理
    df, selected_state = prepare_excited_states_table(excitation_energies, oscillator_strengths, threshold=0.01)

    # 表示
    st.write("### 励起状態リスト")
    # ▼ 表形式の準備（singlet の場合のみオシレーター強度付き）
    if excited_spin == "singlet":
        df = pd.DataFrame({
            "State": list(range(1, len(excitation_energies) + 1)),
            "Energy (eV)": excitation_energies,
            "Oscillator Strength": oscillator_strengths
        })

        with st.expander("励起状態の詳細を表示", expanded=True):
            st.dataframe(df.style.highlight_max(subset=["Oscillator Strength"], color="lightgreen", axis=0))

        # 最初に有意な遷移（f >= 0.01）を強調
        f_input = st.number_input("有意な振動子強度の閾値 (f ≥ )", min_value=0.0, max_value=1.0, value=0.01, step=0.01)
        for i, (e, f) in enumerate(zip(excitation_energies, oscillator_strengths)):
            f = float(f) if isinstance(f, (int, float)) else 0.0
            if f >= f_input:
                st.success(f"有意な最初の一重項遷移(S1 energy): State {i+1}, Energy = {e:.2f} eV, f = {f_input:.4f}")
                break
        else:
            st.warning("f ≥ 0.01 の有意な一重項遷移は見つかりませんでした。")


        # 波長範囲の設定
        st.subheader("波長範囲の設定")
        wavelength_min, wavelength_max = st.slider(
            "表示する波長範囲 (nm)",
            min_value=190,
            max_value=800,
            value=(190, 800),
            step=1
        )

        # UVスペクトルの可視化
        st.subheader("UVスペクトル（滑らか＋スティック）")
        fig, ax = plot_uv_spectrum(wavelengths, wavelength_min, wavelength_max)
        st.pyplot(fig)

        # DataFrame化・ダウンロード
        df_uv = pd.DataFrame(wavelengths, columns=["wavelength_nm", "oscillator_strength"])
        csv = df_uv.to_csv(index=False)
        st.download_button(
            label="波長・振動子強度CSVをダウンロード",
            data=csv,
            file_name="uv_spectrum.csv",
            mime="text/csv"
        )

    else:
        df = pd.DataFrame({
            "State": list(range(1, len(excitation_energies) + 1)),
            "Energy (eV)": excitation_energies,
        })

        with st.expander("励起状態の詳細を表示", expanded=False):
            st.dataframe(df)

        # Tripletは通常 f ≈ 0 のため、エネルギーのみ報告


        import ast  # ファイル冒頭のimportに追加

        # Tripletは通常 f ≈ 0 のため、エネルギーのみ報告
        energies = excitation_energies
        if isinstance(energies, str):
            try:
                energies = ast.literal_eval(energies)
            except Exception:
                st.error("励起エネルギーのデータ形式が不正です。")
                st.stop()
        t1_energy = float(energies[0])
        st.success(f"三重項の最低励起状態 T₁ のエネルギー: {t1_energy:.2f} eV")
        st.info("三重項状態は通常スピン禁制遷移のため、オシレーター強度は表示されません。")

    st.stop()  # 以降の計算処理はスキップ


# --- ここまで追加 ---
if st.button("計算を開始"):
    handler = MoleculeHandler(smiles_input, input_type="smiles")
    inchi_str = inchi.MolToInchi(handler.mol)
    inchikey_str = inchi.InchiToInchiKey(inchi_str)

    directory = os.path.join("data", inchi_str)
    os.makedirs(directory, exist_ok=True)

    if apply_force_field and handler.mol.GetNumAtoms() > 2:
        handler.generate_conformers(num_conformers=num_conformers, forcefield=force_field)
        handler.keep_lowest_energy_conformer()

        pyscf_input = handler.to_pyscf_input()

    result = calculate_excited_state(inchi_str, smiles, pyscf_input, conv_params=conv_params,
                                theory=theory, basis_set=basis_set, charge=charge, spin=spin,
                                solvent_model=solvent_model, eps=eps, maxsteps=maxsteps,
                                nstates=nstates, excited_spin=excited_spin, tda=tda)

    td = result['td']
    excitation_energies = result['excitation_energies']
    xyz_opt = result['xyz_opt']
    mf = result['mf']
    G_tot = result['G_tot']
    frequencies = result['frequencies']
    thermo_info = result['thermo_info']
    

    excitation_energies = td.kernel()
    energies = excitation_energies[0]  # これが1次元のnp.ndarray
    st.success("UVスペクトルの計算が完了しました。")

    # 配列の存在チェック
    if excitation_energies is None or len(excitation_energies[0]) == 0:
        st.error("遷移エネルギーの計算結果がありません。分子や計算条件を確認してください。")
        st.stop()

    # 振動子強度の取得
    oscillator_strengths = td.oscillator_strength()

    if oscillator_strengths is None or getattr(oscillator_strengths, "size", 0) == 0:
        st.error("振動子強度の計算結果がありません。")
        st.stop()

    # energiesとoscillator_strengthsも含めて結果に追加
    if frequencies is not None:
        st.success("✅ Calculation completed successfully!")

        # 虚振動の数を計算（0未満の振動数の個数）
        freq_array = result['frequencies']['freq_wavenumber']
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
                inchi=inchi_str,
                inchikey=inchikey_str,
                g_tot=G_tot.tolist() if isinstance(G_tot, np.ndarray) else G_tot[0],  
                zpe=thermo_info.get('ZPE', None)[0],  
                method=theory,
                basis=basis_set,
                charge=charge,
                spin=spin,
                solvent=solvent_model if solvent_model not in [None, "None"] else None,
                dielectric=eps,
                temperature=298.15,
                pressure=1.0,
                frequencies=freq_array,
                num_imaginary=num_imaginary,
                nstates=nstates,
                excited_spin=excited_spin,
                tda=tda,
                excited_energies=energies.tolist(),  
                oscillator_strengths=oscillator_strengths.tolist(),  
                db_path="data/energy_db.sqlite"
            )
            st.success("データベースに保存しました。")
        except Exception as e:
            st.error(f"データベースへの保存に失敗しました: {e}")

    excited_energies=energies.tolist()
    oscillator_strengths=oscillator_strengths.tolist()

    # 波長と振動子強度の計算
    df, selected_state = prepare_excited_states_table(excited_energies, oscillator_strengths, threshold=0.01)

    # 表示
    st.write("### 励起状態リスト")
    # ▼ 表形式の準備（singlet の場合のみオシレーター強度付き）
    if excited_spin == "singlet":
        df = pd.DataFrame({
            "State": list(range(1, len(excited_energies) + 1)),
            "Energy (eV)": excited_energies,
            "Oscillator Strength": oscillator_strengths
        })

        with st.expander("励起状態の詳細を表示", expanded=True):
            st.dataframe(df.style.highlight_max(subset=["Oscillator Strength"], color="lightgreen", axis=0))

        # 最初に有意な遷移（f >= 0.01）を強調
        f_input = st.number_input("有意な振動子強度の閾値 (f ≥ )", min_value=0.0, max_value=1.0, value=0.01, step=0.01)
        for i, (e, f) in enumerate(zip(excited_energies, oscillator_strengths)):
            if f >= f_input:
                st.success(f"有意な最初の一重項遷移(S1 energy): State {i+1}, Energy = {e:.2f} eV, f = {f_input:.4f}")
                break
        else:
            st.warning("f ≥ 0.01 の有意な一重項遷移は見つかりませんでした。")

        # 波長と振動子強度の計算
        wavelengths = []

        for i, energy in enumerate(excited_energies):
            try:
                energy_float = float(energy)
                wavelength = 1240 / (energy_float * 27.2114)  # エネルギー (au) を nm に変換
            except Exception as e:
                wavelength = None
            wavelengths.append((wavelength, oscillator_strengths[i]))

        # DataFrame化
        df_uv = pd.DataFrame(wavelengths, columns=["wavelength_nm", "oscillator_strength"])
        csv = df_uv.to_csv(index=False)
        st.download_button(
            label="波長・振動子強度CSVをダウンロード",
            data=csv,
            file_name="uv_spectrum.csv",
            mime="text/csv"
        )

    else:
        df = pd.DataFrame({
            "State": list(range(1, len(excited_energies) + 1)),
            "Energy (eV)": excited_energies,
        })

        with st.expander("励起状態の詳細を表示", expanded=False):
            st.dataframe(df)

        # Tripletは通常 f ≈ 0 のため、エネルギーのみ報告
        t1_energy = excited_energies[0]
        st.success(f"三重項の最低励起状態 T₁ のエネルギー: {t1_energy:.2f} eV")
        st.info("三重項状態は通常スピン禁制遷移のため、オシレーター強度は表示されません。")


    # 計算結果をセッションに保存
    st.session_state['uv_result'] = {
        'wavelengths': wavelengths,
        'oscillator_strengths': oscillator_strengths,
        'td': td,
        'energies': excited_energies,
        'mf': mf,
        # 必要なら他のデータも
    }

    # 波長範囲の設定（ユーザーが変更可能）
    st.subheader("波長範囲の設定")
    wavelength_min, wavelength_max = st.slider(
        "表示する波長範囲 (nm)",
        min_value=190,
        max_value=800,
        value=(190, 800),
        step=1
    )

    # UVスペクトルの可視化（滑らかな形状＋スティックスペクトル）
    st.subheader("UVスペクトル（滑らか＋スティック）")
    fig, ax = plot_uv_spectrum(wavelengths, wavelength_min, wavelength_max)
    st.pyplot(fig)

    if mf is not None:
        # --- 励起状態ごとの遷移詳細を表示 ---
        with st.expander("励起状態ごとの遷移詳細を表示", expanded=True):
            st.subheader("励起状態ごとの主な遷移（MOペアと係数）")
            pyscf_mol = mf.mol  # ← ここでPySCF分子オブジェクトを取得
            nocc = pyscf_mol.nelectron // 2
            nmo = mf.mo_coeff.shape[1]
            nvirt = nmo - nocc
            st.info(f"HOMO: MO {nocc-1}, LUMO: MO {nocc}")

            for i, (e, xy) in enumerate(zip(td.e, td.xy)):
                osc_strength = oscillator_strengths[i] if i < len(oscillator_strengths) else None
                wavelength_nm = 1240 / (e * 27.2114)  # eVに変換してから波長に
                st.markdown(
                    f"**励起状態 {i+1}**  \n"
                    f"励起エネルギー: {e:.4f} Hartree = {e * 27.2114:.2f} eV = {wavelength_nm:.1f} nm  \n"
                    f"振動子強度: {osc_strength:.4f}" if osc_strength is not None else ""
                )
                transitions = []
                # xy[0]が多次元配列の場合に備えてflatten
                for idx, c in enumerate(np.ravel(xy[0])):
                    c_real = c.real if isinstance(c, complex) else c
                    if abs(c_real) > 0.01:
                        occ = idx // nvirt
                        virt = nocc + (idx % nvirt)
                        transitions.append(f"MO {occ} → MO {virt}, 係数: {c_real:.4f}")
                if transitions:
                    st.markdown("\n".join([f"- {t}" for t in transitions]))
                else:
                    st.write("主な遷移はありません（係数が小さいため省略）")










