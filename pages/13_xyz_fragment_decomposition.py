"""
個々のフラグメント（分子 A, 分子 B）のエネルギーを計算。
フラグメントが結合した複合体（A + B）のエネルギーを計算。
相互作用エネルギーを計算：
"""

import streamlit as st
from utils.module import load_css
import numpy as np

import py3Dmol  # 3D可視化用ライブラリ
import stmol

from logic.molecule_handler import MoleculeHandler
from controllers.energydecompositionanalysis import get_hf_energy_decomposition
from controllers.xyz_fragment_decomposition import (
    separate_molecules_by_distance,
    separate_molecules_by_clustering,
    analyze_fragment_separation,
    get_fragment_interaction_distance
)

from logic.calculation import theory_options, basis_set_options, run_quantum_calculation



# カスタムCSSを適用
load_css("config/styles.css")

st.title("XYZ座標フラグメント分解解析 (XYZ Fragment Decomposition Analysis)")
st.markdown("XYZ座標から2分子を自動分解してエネルギー分解解析を実行します。")


# 分子入力セクション
st.header("分子入力")

st.subheader("XYZ座標入力")
xyz_input = st.text_area(
    "複合体のXYZ座標を入力してください",
    value="""6
Water dimer
O 0.0000 0.0000 0.0000
H 0.7570 0.0000 0.5860
H -0.7570 0.0000 0.5860
O 0.0000 0.0000 3.0000
H 0.7570 0.0000 3.5860
H -0.7570 0.0000 3.5860""",
    height=200,
    key="xyz_input"
)

st.subheader("分子分解設定")
manual_input = st.checkbox(
    "手動XYZ入力モード",
    value=False,
    help="チェックすると分子A・Bの座標を個別に入力できます",
    key="manual_input"
)

if not manual_input:
    separation_method = st.selectbox(
        "分子分離方法",
        ["距離ベース分離", "クラスタリング分離"],
        key="separation_method"
    )

    if separation_method == "距離ベース分離":
        cutoff_distance = st.slider(
            "分子間距離の閾値 (Å)",
            min_value=1.0, max_value=5.0, value=2.5, step=0.1,
            help="この距離以下の原子は同じ分子として扱われます"
        )
    elif separation_method == "クラスタリング分離":
        cluster_method = st.selectbox(
            "クラスタリング手法",
            ["重心距離ベース", "密度ベース"],
            key="cluster_method"
        )
        
elif manual_input:
    st.info("分子A、分子B、複合体ABのXYZ座標を個別に入力してください")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("分子A")
        xyz_molecule_a = st.text_area(
            "分子AのXYZ座標",
            value="""3
Water molecule A
O 0.0000 0.0000 0.0000
H 0.7570 0.0000 0.5860
H -0.7570 0.0000 0.5860""",
            height=150,
            key="xyz_molecule_a"
        )
    
    with col2:
        st.subheader("分子B")
        xyz_molecule_b = st.text_area(
            "分子BのXYZ座標",
            value="""3
Water molecule B
O 0.0000 0.0000 3.0000
H 0.7570 0.0000 3.5860
H -0.7570 0.0000 3.5860""",
            height=150,
            key="xyz_molecule_b"
        )

# 分子構造の生成と表示
st.header("分子構造の生成")


# 分子構造生成の状態管理
if 'molecules_generated' not in st.session_state:
    st.session_state.molecules_generated = False

# 分子構造生成ボタン
if st.button("分子構造を生成", type="secondary") or st.session_state.molecules_generated:
    # 入力データの確認
    if not manual_input:
        # XYZ入力の確認
        if not xyz_input or not xyz_input.strip():
            st.error("XYZ座標を入力してください")
            st.stop()
    elif manual_input:
        # 手動入力の確認
        if not xyz_molecule_a or not xyz_molecule_a.strip():
            st.error("分子AのXYZ座標を入力してください")
            st.stop()
        if not xyz_molecule_b or not xyz_molecule_b.strip():
            st.error("分子BのXYZ座標を入力してください")
            st.stop()
    
    try:
        with st.spinner("分子構造を解析中..."):
            
            if not manual_input:
                # 自動分離モード
                st.write("🔍 XYZ座標の解析を開始...")
                handler = MoleculeHandler(xyz_input, input_type="xyz")
                st.write(f"全体構造: {handler.mol.GetNumAtoms()} 原子")
                
                # 分子の分解処理
                if separation_method == "距離ベース分離":
                    st.write("🔍 距離ベース分離を実行中...")
                    fragments = separate_molecules_by_distance(handler.mol, cutoff_distance=cutoff_distance)
                    
                elif separation_method == "クラスタリング分離":
                    st.write(f"🔍 {cluster_method}による分離を実行中...")
                    fragments = separate_molecules_by_clustering(handler.mol, n_clusters=2, method="simple")
                
                # 分離結果の確認
                if len(fragments) >= 2:
                    mol1 = fragments[0]
                    mol2 = fragments[1]
                    
                    # 分離品質の分析
                    analysis = analyze_fragment_separation(fragments)
                    
                    st.success(f"2つの分子に分離しました:")
                    st.write(f"- 分子A: {mol1.GetNumAtoms()} 原子 ({analysis['fragment_formulas'][0]})")
                    st.write(f"- 分子B: {mol2.GetNumAtoms()} 原子 ({analysis['fragment_formulas'][1]})")
                    st.write(f"- 分離品質: {analysis['separation_quality']}")
                    
                    # 分子間距離を計算
                    interaction_distance = get_fragment_interaction_distance(mol1, mol2)
                    if interaction_distance:
                        st.write(f"- 分子間最短距離: {interaction_distance:.2f} Å")
                    
                elif len(fragments) == 1:
                    st.error("1つの分子しか検出されませんでした。分離パラメータを調整してください")
                    st.stop()
                else:
                    st.error("分子の分離に失敗しました")
                    st.stop()
                
                # 分離した分子の結合
                from rdkit import Chem
                combo = Chem.CombineMols(mol1, mol2)
                
            elif manual_input:
                # 手動入力モード
                st.write("🔍 手動入力XYZ座標の解析を開始...")
                
                # 各分子のMoleculeHandlerを作成
                handler_a = MoleculeHandler(xyz_molecule_a, input_type="xyz")
                handler_b = MoleculeHandler(xyz_molecule_b, input_type="xyz")
                
                mol1 = handler_a.mol
                mol2 = handler_b.mol
                
                st.write(f"分子A: {mol1.GetNumAtoms()} 原子")
                st.write(f"分子B: {mol2.GetNumAtoms()} 原子")
                
                # 分子の結合
                from rdkit import Chem
                combo = Chem.CombineMols(mol1, mol2)
                
                # 分子間距離を計算
                interaction_distance = get_fragment_interaction_distance(mol1, mol2)
                if interaction_distance:
                    st.write(f"- 分子間最短距離: {interaction_distance:.2f} Å")
                
                st.success("手動入力された分子構造を正常に読み込みました")
            
            # MoleculeHandlerの作成
            st.write("🔍 PySCF入力形式に変換中...")
            
            try:
                handler = MoleculeHandler(combo, input_type="rdkit")
                handler_1 = MoleculeHandler(mol1, input_type="rdkit")
                handler_2 = MoleculeHandler(mol2, input_type="rdkit")
                
                # PySCF形式の原子座標を取得（改行区切り → セミコロン区切りに変換）
                atom_coords_A = handler_1.to_pyscf_input().replace('\n', '; ')
                atom_coords_B = handler_2.to_pyscf_input().replace('\n', '; ')
                atom_coords_AB = handler.to_pyscf_input().replace('\n', '; ')
                
                st.success("PySCF入力形式への変換が完了しました")
                
            except Exception as e:
                st.error(f"PySCF入力形式への変換でエラーが発生しました: {e}")
                raise e
            
            # セッション状態に保存
            st.session_state.molecules_generated = True
            st.session_state.atom_coords_A = atom_coords_A
            st.session_state.atom_coords_B = atom_coords_B  
            st.session_state.atom_coords_AB = atom_coords_AB
            st.session_state.handler = handler
            st.session_state.mol1 = mol1
            st.session_state.mol2 = mol2
            
        st.success("分子構造の生成が完了しました")
        
        # 分離結果の表示
        st.subheader("分離結果")
        if not manual_input:
            st.info(f"💡 {separation_method}による分子分離が完了しました")
        else:
            st.info(f"💡 手動XYZ入力による分子構造の読み込みが完了しました")
        
        # 分子構造を3つのカラムで表示
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.subheader("分子A")
            st.write(f"原子数: {mol1.GetNumAtoms()}")
            with st.expander("原子座標 (PySCF形式)"):
                st.code(st.session_state.atom_coords_A.replace('; ', '\n'))
        
        with col2:
            st.subheader("分子B")
            st.write(f"原子数: {mol2.GetNumAtoms()}")
            with st.expander("原子座標 (PySCF形式)"):
                st.code(st.session_state.atom_coords_B.replace('; ', '\n'))
        
        with col3:
            st.subheader("複合体AB")
            st.write(f"原子数: {handler.mol.GetNumAtoms()}")
            with st.expander("原子座標 (PySCF形式)"):
                st.code(st.session_state.atom_coords_AB.replace('; ', '\n'))

    except Exception as e:
        st.error(f"分子の分解処理に失敗しました: {e}")
        st.session_state.molecules_generated = False
        st.stop()
else:
    st.info("XYZ座標を入力して、分子構造生成ボタンをクリックしてください。")

# 3D構造表示
st.header("分子の3D構造")

if st.session_state.get('molecules_generated', False):
    # 表示する分子の選択
    display_option = st.radio(
        "表示する分子を選択してください",
        ["複合体AB", "分子A", "分子B"],
        key="display_option",
        horizontal=True
    )
    
    try:
        # 選択された分子に応じてMolオブジェクトを決定
        if display_option == "複合体AB":
            display_mol = st.session_state.handler.mol
            display_title = f"複合体AB ({display_mol.GetNumAtoms()} 原子)"
        elif display_option == "分子A":
            # 分子AのMoleculeHandlerを作成
            mol_a_handler = MoleculeHandler(st.session_state.mol1, input_type="rdkit")
            display_mol = mol_a_handler.mol
            display_title = f"分子A ({display_mol.GetNumAtoms()} 原子)"
        else:  # 分子B
            # 分子BのMoleculeHandlerを作成
            mol_b_handler = MoleculeHandler(st.session_state.mol2, input_type="rdkit")
            display_mol = mol_b_handler.mol
            display_title = f"分子B ({display_mol.GetNumAtoms()} 原子)"
        
        st.subheader(display_title)
        
        # 3D構造の表示
        if display_option == "複合体AB":
            mol_block = st.session_state.handler.generate_3d_molblock()
        elif display_option == "分子A":
            mol_block = mol_a_handler.generate_3d_molblock()
        else:  # 分子B
            mol_block = mol_b_handler.generate_3d_molblock()
        
        # 分子情報の表示
        col1, col2 = st.columns(2)
        
        with col1:
            st.write("**分子情報:**")
            st.write(f"原子数: {display_mol.GetNumAtoms()}")
            
            # 分子式を計算
            try:
                from rdkit.Chem import rdMolDescriptors
                formula = rdMolDescriptors.CalcMolFormula(display_mol)
                st.write(f"分子式: {formula}")
            except:
                st.write("分子式: 計算できませんでした")
        
        with col2:
            st.write("**表示設定:**")
            # スタイルの選択
            style_option = st.selectbox(
                "表示スタイル",
                ["sphere", "stick", "ball_and_stick"],
                index=0,
                key=f"style_{display_option}"
            )
            
            # 追加設定
            if style_option == "stick":
                radius = st.slider(
                    "結合の太さ",
                    min_value=0.05, max_value=0.3, value=0.1, step=0.05,
                    key=f"radius_{display_option}"
                )
            # sphereとball_and_stickは自動設定のみ
        
        # スタイルに応じてビューアーを設定
        viewer = py3Dmol.view(width=600, height=400)
        viewer.addModel(mol_block, "mol")
        
        if style_option == "stick":
            viewer.setStyle({"stick": {"radius": radius}})
        elif style_option == "sphere":
            # ファンデルワールス半径を使用
            viewer.setStyle({"sphere": {"scale": 1.0}})  # scaleでファンデルワールス半径の倍率を指定
        elif style_option == "ball_and_stick":
            # 自動設定：ファンデルワールス半径を使用
            viewer.setStyle({"sphere": {"scale": 0.5}, "stick": {"radius": 0.1}})
        
        viewer.zoomTo()
        stmol.showmol(viewer, height=400)
        
    except Exception as e:
        st.warning(f"3D構造の表示ができませんでした: {e}")
        st.error(f"詳細なエラー: {str(e)}")
else:
    st.info("3D構造を表示するには、まず分子構造を生成してください。")

# エネルギー分解解析の実行
st.header("エネルギー分解解析")

# サイドバーで計算設定
st.header("計算設定")

# 理論レベルと基底関数の選択
theory = st.selectbox("Theory", theory_options, key="theory_selectbox")
basis_set = st.selectbox("Basis Set", basis_set_options, key="basis_set_selectbox")

charge = st.number_input("電荷", value=0, step=1)
spin = st.number_input("スピン多重度 (2S)", value=0, step=1)

# 並列計算の有無を選択
st.header("Parallel Calculation Options")
import multiprocessing

logical_cores = multiprocessing.cpu_count()
try:
    physical_cores = multiprocessing.cpu_count(logical=False)
except TypeError:
    physical_cores = logical_cores // 2  # Fallback

st.write(f"使用しているパソコンの論理コア数: {logical_cores}")
st.write(f"使用しているパソコンの物理コア数: {physical_cores}")

# 物理コアが3以上なら並列計算を推奨（分子A、分子B、複合体ABの3つの計算のため）
recommend_parallel = physical_cores >= 3

parallel_option = st.checkbox(
    "並列計算を有効にする（推奨）" if recommend_parallel else "並列計算を有効にする",
    value=recommend_parallel,
    key="parallel_option_checkbox"
)

if recommend_parallel:
    st.info("物理コア数が3以上のため、並列計算が推奨されます。（分子A、分子B、複合体ABの3つの計算を並列実行）")
elif physical_cores >= 2:
    st.info("物理コア数が2以上です。並列計算により計算時間を短縮できる可能性があります。")
else:
    st.warning("物理コア数が少ないため、並列計算は非推奨です。")

# 計算方法と参考文献の表示
with st.expander("計算方法と参考文献を表示", expanded=False):
    st.markdown("### 🧪 Method for XYZ Fragment Decomposition Analysis")
    st.markdown(
        "**Computational Details**  \n"
        "Complex molecular structures were processed from XYZ coordinates using RDKit [1] for molecular fragment decomposition.  \n"
        f"Fragment separation was performed using **{'手動XYZ入力' if manual_input else '自動分子分離'}** method"
        f"{f' ({separation_method})' if not manual_input else ''}.  \n"
        f"Single-point energy calculations were performed at the **{theory}/{basis_set}** level using PySCF [2].  \n"
        "The Energy Decomposition Analysis (EDA) provides detailed insights into intermolecular interactions by decomposing the total electronic energy into:  \n"
        "- **Nuclear repulsion energy (E_nuc)**: Electrostatic repulsion between nuclei  \n"
        "- **Core Hamiltonian energy (E_core)**: Electron-nuclear attraction and kinetic energy  \n"
        "- **Coulomb interaction energy (E_J)**: Classical electron-electron repulsion  \n"
        "- **Exchange energy (E_K)**: Quantum mechanical exchange interaction  \n"
        "- **Total electronic energy (E_elec)**: Sum of all electronic contributions  \n"
        "The interaction energy is calculated as: **ΔE_int = E_AB - (E_A + E_B)**  \n"
        "where E_AB is the energy of the complex and E_A, E_B are the energies of isolated fragments.  \n"
        "Each energy component is decomposed to understand the physical origins of intermolecular interactions.  \n"
        "Energy values are provided in Hartree (Ha) and converted to kcal/mol for comparative analysis (1 Ha = 627.509 kcal/mol).  \n"
        "This decomposition analysis is essential for understanding non-covalent interactions, hydrogen bonding, and van der Waals forces in pre-formed complexes."
    )
    st.markdown("---")
    st.markdown(
        "**References**  \n"
        "[1] Landrum, G. RDKit: Open-source cheminformatics. [https://www.rdkit.org](https://www.rdkit.org)  \n"
        "[2] Sun, Q. *et al.* PySCF: The Python-based Simulations of Chemistry Framework. **WIREs Comput Mol Sci** *2018*, **8**, e1340. DOI: [10.1002/wcms.1340](https://doi.org/10.1002/wcms.1340)  \n"
        "[3] PocLab streamlit-pyscf: Quantum chemistry web interface. [https://github.com/poclab-web/streamlit-pyscf](https://github.com/poclab-web/streamlit-pyscf)  \n"
        "[4] Szabo, A.; Ostlund, N. S. *Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory*; Dover Publications: New York, 1996.  \n"
        "[5] Helgaker, T.; Jørgensen, P.; Olsen, J. *Molecular Electronic-Structure Theory*; Wiley: Chichester, 2000.  \n"
        "[6] Kitaura, K.; Morokuma, K. A new energy decomposition scheme for molecular interactions within the Hartree-Fock approximation. **Int. J. Quantum Chem.** *1976*, **10**, 325-340. DOI: [10.1002/qua.560100211](https://doi.org/10.1002/qua.560100211)"
    )

if st.button("計算実行", type="primary"):
    if not st.session_state.get('molecules_generated', False):
        st.warning("まず分子構造を生成してください。")
    elif xyz_input and xyz_input.strip():
        try:
            # セッション状態から座標を取得
            atom_coords_A = st.session_state.atom_coords_A
            atom_coords_B = st.session_state.atom_coords_B
            atom_coords_AB = st.session_state.atom_coords_AB
            
            # 分子からSMILES文字列を生成
            from rdkit import Chem
            mol1 = st.session_state.mol1
            mol2 = st.session_state.mol2
            
            try:
                smiles_A = Chem.MolToSmiles(mol1)
                smiles_B = Chem.MolToSmiles(mol2)
            except:
                smiles_A = f"Fragment_A_{mol1.GetNumAtoms()}atoms"
                smiles_B = f"Fragment_B_{mol2.GetNumAtoms()}atoms"
            
            # 計算パラメータの準備
            compound_names = [f"Mol_A_{smiles_A}", f"Mol_B_{smiles_B}", f"Complex_{smiles_A}_{smiles_B}"]
            smiles_list = [smiles_A, smiles_B, f"{smiles_A}.{smiles_B}"]
            atom_inputs = [atom_coords_A, atom_coords_B, atom_coords_AB]
            
            # プログレスバーと状態表示
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            status_text.text("エネルギー分解解析を開始しています...")
            progress_bar.progress(10)
            
            # エネルギー分解解析の実行
            with st.spinner("量子化学計算を実行中..."):
                eda_results = get_hf_energy_decomposition(
                    compound_names=compound_names,
                    smiles_list=smiles_list,
                    atom_inputs=atom_inputs,
                    basis=basis_set,
                    theory=theory,
                    charge=charge,
                    spin=spin,
                    parallel_option=parallel_option
                )
            
            progress_bar.progress(90)
            status_text.text("結果を処理中...")
            
            # 結果の表示
            st.success("計算が完了しました！")
            progress_bar.progress(100)
            status_text.text("完了")
            
            # 計算情報の表示
            st.subheader("計算情報")
            calc_info = eda_results.get("calculation_info", {})
            info_col1, info_col2, info_col3 = st.columns(3)
            
            with info_col1:
                st.metric("理論レベル", calc_info.get('theory', 'Unknown'))
                st.metric("基底関数", calc_info.get('basis', 'Unknown'))
                
            with info_col2:
                st.metric("電荷", calc_info.get('charge', 'Unknown'))
                st.metric("スピン", calc_info.get('spin', 'Unknown'))
                
            with info_col3:
                st.metric("並列処理", "有効" if parallel_option else "無効")
                st.metric("計算ID", calc_info.get('timestamp', 'Unknown')[:10] + "...")
            
            # エラーチェック
            energy_a = eda_results["molecule_A"]
            energy_b = eda_results["molecule_B"]
            energy_ab = eda_results["complex_AB"]
            
            if "error" in energy_a or "error" in energy_b or "error" in energy_ab:
                st.error("計算中にエラーが発生しました:")
                if "error" in energy_a:
                    st.error(f"分子A: {energy_a['error']}")
                if "error" in energy_b:
                    st.error(f"分子B: {energy_b['error']}")
                if "error" in energy_ab:
                    st.error(f"複合体AB: {energy_ab['error']}")
            else:
                # 収束チェック
                converged_a = energy_a.get("converged", False)
                converged_b = energy_b.get("converged", False)
                converged_ab = energy_ab.get("converged", False)
                
                if not (converged_a and converged_b and converged_ab):
                    st.warning("一部の計算が収束していません")
                    conv_col1, conv_col2, conv_col3 = st.columns(3)
                    with conv_col1:
                        st.metric("分子A収束", "✓" if converged_a else "✗")
                    with conv_col2:
                        st.metric("分子B収束", "✓" if converged_b else "✗")
                    with conv_col3:
                        st.metric("複合体AB収束", "✓" if converged_ab else "✗")
                
                # エネルギー分解結果の表示
                if "energy_decomposition" in eda_results:
                    decomp = eda_results["energy_decomposition"]
                    
                    st.subheader("個別エネルギー")
                    energy_keys = ["energy", "E_nuc", "E_core", "E_J", "E_K", "E_elec"]
                    
                    # テーブル形式で表示
                    import pandas as pd
                    
                    individual_data = []
                    for key in energy_keys:
                        if f"{key}_A" in decomp and f"{key}_B" in decomp and f"{key}_AB" in decomp:
                            display_key = "E_total" if key == "energy" else key
                            individual_data.append({
                                "エネルギー成分": display_key,
                                "分子A (Ha)": f"{decomp[f'{key}_A']:+.6f}",
                                "分子B (Ha)": f"{decomp[f'{key}_B']:+.6f}",
                                "複合体AB (Ha)": f"{decomp[f'{key}_AB']:+.6f}"
                            })
                    
                    if individual_data:
                        df_individual = pd.DataFrame(individual_data)
                        st.dataframe(df_individual, use_container_width=True)
                    
                    st.subheader("相互作用エネルギー")
                    interaction_data = []
                    for key in energy_keys:
                        delta_key = f"Δ{key}"
                        kcal_key = f"Δ{key}_kcal_mol"
                        if delta_key in decomp and kcal_key in decomp:
                            display_key = "ΔE_total" if key == "energy" else f"Δ{key}"
                            if isinstance(decomp[delta_key], (int, float, np.number)):
                                interaction_data.append({
                                    "エネルギー成分": display_key,
                                    "相互作用エネルギー (Ha)": f"{decomp[delta_key]:+.6f}",
                                    "相互作用エネルギー (kcal/mol)": f"{decomp[kcal_key]:+.2f}"
                                })
                            else:
                                interaction_data.append({
                                    "エネルギー成分": display_key,
                                    "相互作用エネルギー (Ha)": str(decomp[delta_key]),
                                    "相互作用エネルギー (kcal/mol)": "N/A"
                                })
                    
                    if interaction_data:
                        df_interaction = pd.DataFrame(interaction_data)
                        st.dataframe(df_interaction, use_container_width=True)
                        
                        # 重要な結果をハイライト
                        if len(interaction_data) > 0:
                            total_interaction = interaction_data[0]  # E_totalの相互作用エネルギー
                            st.info(f"**全相互作用エネルギー**: {total_interaction['相互作用エネルギー (Ha)']} Ha = {total_interaction['相互作用エネルギー (kcal/mol)']} kcal/mol")
                    
                    # 結果のダウンロード
                    st.subheader("結果のダウンロード")
                    
                    # JSONファイルとして結果をダウンロード
                    import json
                    result_json = json.dumps(eda_results, indent=2, default=str)
                    st.download_button(
                        label="結果をJSONでダウンロード",
                        data=result_json,
                        file_name=f"eda_results_{calc_info.get('timestamp', 'unknown')}.json",
                        mime="application/json"
                    )
                    
                    # CSVファイルとして相互作用エネルギーをダウンロード
                    if interaction_data:
                        csv_data = pd.DataFrame(interaction_data).to_csv(index=False)
                        st.download_button(
                            label="相互作用エネルギーをCSVでダウンロード",
                            data=csv_data,
                            file_name=f"interaction_energies_{calc_info.get('timestamp', 'unknown')}.csv",
                            mime="text/csv"
                        )
                        
        except Exception as e:
            st.error(f"計算中にエラーが発生しました: {str(e)}")
            st.exception(e)
    else:
        st.warning("XYZ座標を入力してください")