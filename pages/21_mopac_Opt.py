"""
MOPAC半経験的量子化学計算ソフトウェアを使用して構造最適化を行う。

機能:
- 計算に使用する理論手法を選択可能（PM7, PM6, AM1, MNDO）。
- 構造最適化計算を実行し、結果をインタラクティブに表示。
- 最適化前後の構造比較、エネルギー変化、および幾何学的パラメータを表示。
"""

import streamlit as st
import stmol
import pandas as pd
from rdkit import Chem
import py3Dmol
import os
import numpy as np
import subprocess
import shutil

from utils.module import load_css
from utils.mopac_ui import require_mopac
from logic.molecule_handler import MoleculeHandler
from logic.mopac_calculation import MopacCalculator, theory_options

# カスタムCSSを適用
load_css("config/styles.css")

# 本文
st.title("MOPAC Geometry Optimization")

# MOPACインストール状況の確認
require_mopac()

# ユーザー入力
st.header("Molecular Input")
input_type = st.selectbox("Select Input Type", ["XYZ", "SMILES"])
atom_input = st.text_area(
    "Enter Molecular Structure",
    "C     -2.385271    1.648249    0.000003\nH     -3.181824    1.061240    0.407216\nH     -2.536116    1.773363   -1.051896\nH     -2.370100    2.607079    0.474676\nO     -1.453042    1.151314    0.170017\nH     -1.466511    0.300016   -0.251422"
    if input_type == "XYZ"
    else "CO",
)

# MOPAC理論手法の選択
theory = st.selectbox("Theory", theory_options)

# その他の設定
with st.expander("Other Settings"):
    charge = st.number_input("Molecular Charge", min_value=-10, max_value=10, value=0, step=1)
    multiplicity = st.number_input("Spin Multiplicity (2S + 1)", min_value=1, max_value=10, value=1, step=1)
    
    # 配座探索の設定
    st.subheader("Initial Conformer Search")
    use_conformer_search = st.checkbox("Perform initial conformer search", value=False, 
                                      help="最適化前に分子力場による配座探索を実行")
    
    # 複数配座計算の設定
    multi_conformer_calculation = st.checkbox("Calculate multiple conformers", value=False,
                                             help="複数配座を同時にMOPAC計算して比較")
    
    if use_conformer_search:
        col1, col2 = st.columns(2)
        with col1:
            num_conformers = st.number_input("Number of conformers", min_value=1, max_value=100, value=10, step=1,
                                           help="生成する配座数")
            forcefield = st.selectbox("Force field", ["UFF", "MMFF"], index=1,
                                    help="使用する分子力場")
        with col2:
            max_iterations = st.number_input("Max iterations", min_value=50, max_value=1000, value=200, step=50,
                                           help="最適化の最大反復回数")
            prune_threshold = st.number_input("RMS pruning threshold", min_value=0.1, max_value=2.0, value=0.5, step=0.1,
                                            help="類似配座の除去閾値（Å）")
    
    # 複数配座計算の詳細設定
    if multi_conformer_calculation:
        st.subheader("Multi-Conformer Calculation Settings")
        col1, col2 = st.columns(2)
        with col1:
            max_conformers_to_calculate = st.number_input("Max conformers to calculate", 
                                                        min_value=1, max_value=50, value=5, step=1,
                                                        help="MOPAC計算する最大配座数")
            show_all_results = st.checkbox("Show all conformer results", value=True,
                                         help="全配座の結果を表示（OFFの場合は最安定配座のみ）")
        with col2:
            energy_threshold = st.number_input("Energy threshold (kcal/mol)", 
                                             min_value=0.0, max_value=50.0, value=10.0, step=1.0,
                                             help="表示する配座のエネルギー閾値")
            parallel_calculation = st.checkbox("Show calculation progress", value=True,
                                             help="計算進行状況を表示")
    
    # MOPACの最適化オプション
    st.subheader("Optimization Options")
    precise = st.checkbox("Use PRECISE mode", value=True, help="高精度モード（厳密な収束）")
    gnorm = st.number_input("Gradient norm convergence", min_value=0.1, max_value=10.0, value=1.0, step=0.1,
                           help="勾配収束判定値（kcal/mol/Å）")
    
    # 最適化の詳細設定
    col1, col2 = st.columns(2)
    with col1:
        optimize_cartesian = st.checkbox("Cartesian coordinates", value=False, 
                                       help="デカルト座標系で最適化（通常は内部座標系）")
        ts_search = st.checkbox("Transition state search", value=False,
                               help="遷移状態構造探索（TS）")
    with col2:
        symmetry = st.checkbox("Use molecular symmetry", value=False,
                              help="分子対称性を利用した最適化")
        isotope = st.checkbox("Include isotope effects", value=False,
                             help="同位体効果を考慮")
    
    # 溶媒効果の設定
    st.subheader("Solvent Effects")
    use_solvent = st.checkbox("Include solvent effects", value=False, 
                             help="COSMO溶媒効果モデルを使用")
    
    solvent = None
    if use_solvent:
        # 溶媒データを読み込み
        import pandas as pd
        try:
            solvents_df = pd.read_csv("config/solvents_epsilon.csv")
            solvent_options = ["Custom"] + solvents_df["Solvent"].tolist()
            
            col1, col2 = st.columns(2)
            with col1:
                solvent_choice = st.selectbox("Solvent", solvent_options, 
                                            help="溶媒を選択してください")
                
            with col2:
                if solvent_choice == "Custom":
                    epsilon = st.number_input("Dielectric constant (ε)", 
                                            min_value=1.0, max_value=100.0, 
                                            value=78.36, step=0.1,
                                            help="溶媒の誘電率")
                    solvent = {"name": "Custom", "epsilon": epsilon}
                else:
                    selected_solvent = solvents_df[solvents_df["Solvent"] == solvent_choice]
                    if not selected_solvent.empty:
                        epsilon = selected_solvent["Epsilon"].iloc[0]
                        st.info(f"ε = {epsilon}")
                        solvent = {"name": solvent_choice, "epsilon": epsilon}
                        
        except Exception as e:
            st.error(f"Error loading solvent data: {e}")
            st.write("Using manual epsilon input")
            epsilon = st.number_input("Dielectric constant (ε)", 
                                    min_value=1.0, max_value=100.0, 
                                    value=78.36, step=0.1)
            solvent = {"name": "Custom", "epsilon": epsilon}

# 分子構造を処理
handler = None
if st.button("Run MOPAC Optimization"):
    try:
        handler = MoleculeHandler(atom_input, input_type=input_type.lower())
        if not handler.mol:
            raise ValueError("Invalid molecular input. Please check your format.")

        # 化合物名を取得
        compound_name = Chem.MolToInchiKey(handler.mol)
        smiles = Chem.MolToSmiles(handler.mol)

        # 初期構造の保存
        initial_mol = Chem.Mol(handler.mol)
        
        # 配座探索の実行
        if use_conformer_search:
            st.write("Performing initial conformer search...")
            try:
                with st.spinner("Generating conformers..."):
                    conformers = handler.generate_conformers(
                        num_conformers=num_conformers,
                        max_iterations=max_iterations,
                        prune_rms_threshold=prune_threshold,
                        forcefield=forcefield
                    )
                
                if conformers:
                    st.success(f"Generated {len(conformers)} conformers")
                    
                    # エネルギー情報を表示
                    st.subheader("Initial Conformer Energies")
                    energy_data = []
                    for i, conf in enumerate(conformers):
                        conf_id = conf.GetId()
                        if handler.mol.HasProp(f"Energy_{conf_id}"):
                            energy = float(handler.mol.GetProp(f"Energy_{conf_id}"))
                            energy_data.append({"Conformer": conf_id, "Energy (kcal/mol)": energy})
                    
                    if energy_data:
                        df = pd.DataFrame(energy_data)
                        st.dataframe(df)
                        
                        # 最低エネルギー配座を選択
                        if st.checkbox("Use lowest energy conformer for optimization", value=True):
                            handler.keep_lowest_energy_conformer()
                            min_energy = min(energy_data, key=lambda x: x["Energy (kcal/mol)"])
                            st.info(f"Using conformer {min_energy['Conformer']} with energy {min_energy['Energy (kcal/mol)']:.6f} kcal/mol")
                    else:
                        st.warning("No energy information available for conformers")
                else:
                    st.warning("No conformers were generated. Using original structure.")
            except Exception as e:
                st.error(f"Error during conformer generation: {e}")
                st.info("Continuing with original structure...")

        # ディレクトリの作成
        directory = os.path.join("data", compound_name)
        os.makedirs(directory, exist_ok=True)

        # 初期構造の表示
        st.subheader("Initial Structure")
        col1, col2 = st.columns(2)

        # Display 2D structure in the first column
        with col1:
            st.write("**2D Structure**")
            handler.generate_2d_image(f"{directory}/molecule_2d.png")
            st.image(f"{directory}/molecule_2d.png", caption=Chem.MolToSmiles(handler.mol))

        # Display initial XYZ coordinates in the second column
        with col2:
            st.write("**Initial 3D Structure**")
            try:
                # 3D構造表示
                mol_block = handler.generate_3d_molblock()
                viewer = py3Dmol.view(width=400, height=400)
                viewer.addModel(mol_block, "mol")
                viewer.setStyle({"stick": {}})
                viewer.zoomTo()
                stmol.showmol(viewer, height=400)
                
                # 初期構造のXYZ座標を生成
                mol = handler.mol
                initial_xyz = f"{mol.GetNumAtoms()}\n"
                initial_xyz += f"Initial Structure - {Chem.MolToSmiles(mol)}\n"
                
                if mol.GetNumConformers() > 0:
                    conf = mol.GetConformer()
                    for i in range(mol.GetNumAtoms()):
                        atom = mol.GetAtomWithIdx(i)
                        pos = conf.GetAtomPosition(i)
                        initial_xyz += f"{atom.GetSymbol():<2} {pos.x:>12.6f} {pos.y:>12.6f} {pos.z:>12.6f}\n"
                else:
                    initial_xyz += "No coordinates available\n"
                
                st.text_area("Initial XYZ", initial_xyz, height=200, key="preview_initial_xyz")
                
                # 配座探索が実行された場合は情報を表示
                if use_conformer_search and handler.mol.GetNumConformers() > 0:
                    conf_id = 0  # 現在表示されている配座ID
                    if handler.mol.HasProp(f"Energy_{conf_id}"):
                        energy = float(handler.mol.GetProp(f"Energy_{conf_id}"))
                        st.caption(f"Initial Energy: {energy:.6f} kcal/mol ({forcefield})")
                    else:
                        st.caption(f"Optimized with {forcefield}")
                        
            except Exception as e:
                st.warning(f"Unable to generate initial 3D structure or XYZ coordinates: {e}")

        # MOPAC構造最適化の実行
        if multi_conformer_calculation:
            st.write("Running MOPAC multiple conformer optimization...")
            
            # 作業ディレクトリを指定
            work_dir = os.path.join(directory, "mopac_work")
            
            try:
                # MopacCalculatorのインスタンス作成
                calculator = MopacCalculator(handler, work_dir=work_dir)
                
                # 最適化計算のキーワード設定
                keywords = []
                if precise:
                    keywords.append("PRECISE")
                if gnorm != 1.0:
                    keywords.append(f"GNORM={gnorm:.1f}")
                if optimize_cartesian:
                    keywords.append("CARTESIAN")
                if ts_search:
                    keywords.append("TS")
                if symmetry:
                    keywords.append("SYMMETRY")
                if isotope:
                    keywords.append("ISOTOPE")
                
                # 溶媒効果の追加
                if use_solvent and solvent:
                    keywords.append(f"COSMO EPS={solvent['epsilon']:.4f}")
                
                # 複数配座の最適化実行
                with st.spinner("Generating conformers and optimizing with MOPAC..."):
                    multi_result = calculator.optimize_multiple_conformers(
                        theory=theory,
                        charge=charge,
                        multiplicity=multiplicity,
                        keywords=keywords,
                        title="Multi-conformer MOPAC Optimization",
                        max_conformers=max_conformers_to_calculate,
                        num_conformers=num_conformers,
                        max_iterations=max_iterations,
                        prune_rms_threshold=prune_threshold,
                        forcefield=forcefield
                    )
                
                # 複数配座の結果表示
                st.subheader("🧬 Multi-Conformer Optimization Results")
                
                # 設定情報の表示
                if 'generation_parameters' in multi_result:
                    params = multi_result['generation_parameters']
                    st.info(f"🔬 **Generated with {params['forcefield']} force field** | "
                           f"Requested: {params['num_conformers_requested']} conformers | "
                           f"Max iterations: {params['max_iterations']} | "
                           f"RMS threshold: {params['prune_rms_threshold']} Å")
                
                # 計算サマリー
                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    st.metric("Total Conformers", multi_result['total_conformers'])
                with col2:
                    st.metric("Successful", multi_result['successful_calculations'])
                with col3:
                    st.metric("Failed", multi_result['failed_calculations'])
                with col4:
                    if multi_result['energy_range']:
                        st.metric("Energy Range", f"{multi_result['energy_range']['span']:.3f} kcal/mol")
                
                # エネルギーランキング表
                if multi_result['energy_ranking']:
                    st.subheader("📊 Conformer Energy Ranking")
                    
                    ranking_data = []
                    for i, conf_energy in enumerate(multi_result['energy_ranking']):
                        conf_id = conf_energy['conformer_id']
                        relative_energy = conf_energy['heat_of_formation'] - multi_result['energy_ranking'][0]['heat_of_formation']
                        
                        # 初期エネルギーとの差分計算
                        energy_change = None
                        if conf_energy['initial_energy'] is not None:
                            energy_change = conf_energy['heat_of_formation'] - conf_energy['initial_energy']
                        
                        ranking_data.append({
                            "Rank": i + 1,
                            "Conformer ID": conf_id,
                            "Heat of Formation (kcal/mol)": f"{conf_energy['heat_of_formation']:.6f}",
                            "Relative Energy (kcal/mol)": f"{relative_energy:.6f}",
                            "Initial Energy (kcal/mol)": f"{conf_energy['initial_energy']:.6f}" if conf_energy['initial_energy'] is not None else "N/A",
                            "Energy Change (kcal/mol)": f"{energy_change:.6f}" if energy_change is not None else "N/A",
                            "Force Field": conf_energy.get('forcefield', 'Unknown')
                        })
                    
                    df_ranking = pd.DataFrame(ranking_data)
                    st.dataframe(df_ranking, use_container_width=True)
                
                # 最安定配座の詳細表示
                if multi_result['best_conformer']:
                    st.subheader("🏆 Best Conformer Details")
                    best_conf_id = multi_result['best_conformer']['conformer_id']
                    best_result = multi_result['conformer_results'][best_conf_id]
                    
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Best Conformer ID", best_conf_id)
                    with col2:
                        st.metric("Heat of Formation", f"{multi_result['best_conformer']['heat_of_formation']:.6f} kcal/mol")
                    with col3:
                        if multi_result['best_conformer']['initial_energy'] is not None:
                            energy_improvement = multi_result['best_conformer']['heat_of_formation'] - multi_result['best_conformer']['initial_energy']
                            st.metric("Energy Improvement", f"{energy_improvement:.6f} kcal/mol")
                    
                    # 最安定配座の構造表示
                    try:
                        best_optimized_mol = calculator.read_optimized_structure(best_result['arc_file']) if best_result.get('arc_file') else None
                        if best_optimized_mol:
                            st.subheader("🧬 Best Conformer Structure")
                            
                            col1, col2 = st.columns(2)
                            with col1:
                                st.write("**Initial Best Conformer**")
                                # 初期構造の表示
                                initial_mol_copy = Chem.Mol(initial_mol)
                                initial_mol_copy.RemoveAllConformers()
                                original_conf = initial_mol.GetConformer(best_conf_id)
                                new_conf = Chem.Conformer(initial_mol_copy.GetNumAtoms())
                                for i in range(initial_mol_copy.GetNumAtoms()):
                                    pos = original_conf.GetAtomPosition(i)
                                    new_conf.SetAtomPosition(i, pos)
                                initial_mol_copy.AddConformer(new_conf, assignId=True)
                                
                                initial_handler_best = MoleculeHandler(initial_mol_copy, input_type="rdkit")
                                initial_molblock_best = initial_handler_best.generate_3d_molblock()
                                initial_viewer_best = py3Dmol.view(width=400, height=400)
                                initial_viewer_best.addModel(initial_molblock_best, "mol")
                                initial_viewer_best.setStyle({"stick": {}})
                                initial_viewer_best.zoomTo()
                                stmol.showmol(initial_viewer_best, height=400)
                            
                            with col2:
                                st.write("**Optimized Best Conformer**")
                                # 最適化構造の表示
                                optimized_handler_best = MoleculeHandler(best_optimized_mol, input_type="rdkit")
                                optimized_molblock_best = optimized_handler_best.generate_3d_molblock()
                                optimized_viewer_best = py3Dmol.view(width=400, height=400)
                                optimized_viewer_best.addModel(optimized_molblock_best, "mol")
                                optimized_viewer_best.setStyle({"stick": {}})
                                optimized_viewer_best.zoomTo()
                                stmol.showmol(optimized_viewer_best, height=400)
                    except Exception as e:
                        st.warning(f"Could not display best conformer structure: {e}")
                
                # 全配座の結果表示（オプション）
                if show_all_results and multi_result['conformer_results']:
                    st.subheader("📋 All Conformer Results")
                    
                    # エネルギー閾値でフィルタリング
                    if multi_result['best_conformer']:
                        min_energy = multi_result['best_conformer']['heat_of_formation']
                        filtered_results = {}
                        for conf_id, result in multi_result['conformer_results'].items():
                            if (result.get('success') and 
                                'heat_of_formation' in result and
                                result['heat_of_formation'] - min_energy <= energy_threshold):
                                filtered_results[conf_id] = result
                        
                        if filtered_results:
                            for conf_id, result in filtered_results.items():
                                with st.expander(f"Conformer {conf_id} (ΔE = {result['heat_of_formation'] - min_energy:.3f} kcal/mol)"):
                                    col1, col2 = st.columns(2)
                                    with col1:
                                        if 'summary' in result and 'key_results' in result['summary']:
                                            key_results = result['summary']['key_results']
                                            for key, value in key_results.items():
                                                st.write(f"**{key.replace('_', ' ').title()}**: {value}")
                                    with col2:
                                        st.write(f"**Success**: {result.get('success', False)}")
                                        st.write(f"**Return Code**: {result.get('return_code', 'N/A')}")
                                        if result.get('arc_file'):
                                            st.write(f"**ARC File**: {os.path.basename(result['arc_file'])}")
                        else:
                            st.info("No conformers within the energy threshold")
                
                # 複数配座計算では、従来の単一配座結果表示をスキップ
                result = {"success": True, "multi_conformer": True, "best_result": multi_result}
                
            except Exception as e:
                st.error(f"Error in multi-conformer calculation: {e}")
                import traceback
                st.text_area("Error Details", traceback.format_exc(), height=200, key="multi_error_details")
                st.stop()
        
        else:
            st.write("Running MOPAC geometry optimization...")
            
            # 作業ディレクトリを指定
            work_dir = os.path.join(directory, "mopac_work")
            
            try:
                # MopacCalculatorのインスタンス作成
                calculator = MopacCalculator(handler, work_dir=work_dir)
                
                # 最適化計算の実行
                keywords = []
                if precise:
                    keywords.append("PRECISE")
                if gnorm != 1.0:
                    keywords.append(f"GNORM={gnorm:.1f}")
                if optimize_cartesian:
                    keywords.append("CARTESIAN")
                if ts_search:
                    keywords.append("TS")
                if symmetry:
                    keywords.append("SYMMETRY")
                if isotope:
                    keywords.append("ISOTOPE")
                
                # 溶媒効果の追加
                if use_solvent and solvent:
                    keywords.append(f"COSMO EPS={solvent['epsilon']:.4f}")
                
                result = calculator.optimize_geometry(
                    theory=theory,
                    charge=charge,
                    multiplicity=multiplicity,
                    keywords=keywords,
                    title="MOPAC Geometry Optimization"
                )
                
                # 計算結果の詳細確認と表示
                st.subheader("Optimization Results")
                
                # 実際の成功/失敗を判断
                actual_success = (
                    result.get('success', False) and 
                    result.get('return_code', -1) == 0 and
                    result.get('output_file') and 
                    os.path.exists(result.get('output_file', ''))
                )
                
            except Exception as e:
                st.error(f"Error running MOPAC optimization: {e}")
                import traceback
                st.text_area("Error Details", traceback.format_exc(), height=200, key="single_error_details")
                st.stop()
        
        # 結果の表示処理（複数配座 vs 単一配座の分岐処理）
        if result.get("multi_conformer"):
            # 複数配座計算の場合は既に表示済みなので何もしない
            pass
        else:
            # 単一配座計算の結果表示
            
            if actual_success:
                st.success("✅ MOPAC geometry optimization completed successfully!")
                
                # 最適化された構造の読み込み
                optimized_mol = None
                structure_source = None
                
                # ARCファイルから構造を読み込み
                if result.get('arc_file') and os.path.exists(result['arc_file']):
                    try:
                        st.write("📁 Reading optimized structure from ARC file...")
                        arc_file_size = os.path.getsize(result['arc_file'])
                        st.write(f"ARC file size: {arc_file_size} bytes")
                        
                        optimized_mol = calculator.read_optimized_structure(result['arc_file'])
                        if optimized_mol:
                            structure_source = "ARC file"
                            st.success(f"✅ Optimized structure loaded successfully from ARC file! ({optimized_mol.GetNumAtoms()} atoms)")
                            
                            # 分子の基本情報を表示
                            if optimized_mol.GetNumConformers() > 0:
                                conf = optimized_mol.GetConformer()
                                positions = []
                                for i in range(optimized_mol.GetNumAtoms()):
                                    pos = conf.GetAtomPosition(i)
                                    positions.append((pos.x, pos.y, pos.z))
                                st.write(f"Coordinate range: X: {min(p[0] for p in positions):.3f} to {max(p[0] for p in positions):.3f}")
                                st.write(f"                  Y: {min(p[1] for p in positions):.3f} to {max(p[1] for p in positions):.3f}")
                                st.write(f"                  Z: {min(p[2] for p in positions):.3f} to {max(p[2] for p in positions):.3f}")
                        else:
                            st.warning("⚠️ Could not parse structure from ARC file")
                    except Exception as e:
                        st.warning(f"⚠️ Could not load optimized structure from ARC file: {e}")
                        import traceback
                        st.text(traceback.format_exc())
                
                # ARCファイルから読み込めない場合、OUTファイルから試す
                if not optimized_mol and result.get('output_file') and os.path.exists(result['output_file']):
                    try:
                        st.write("📁 Attempting to read structure from output file...")
                        output_file_size = os.path.getsize(result['output_file'])
                        st.write(f"Output file size: {output_file_size} bytes")
                        
                        optimized_mol = calculator.read_structure_from_output(result['output_file'])
                        if optimized_mol:
                            structure_source = "Output file"
                            st.success(f"✅ Optimized structure loaded from output file! ({optimized_mol.GetNumAtoms()} atoms)")
                            
                            # 分子の基本情報を表示
                            if optimized_mol.GetNumConformers() > 0:
                                conf = optimized_mol.GetConformer()
                                positions = []
                                for i in range(optimized_mol.GetNumAtoms()):
                                    pos = conf.GetAtomPosition(i)
                                    positions.append((pos.x, pos.y, pos.z))
                                st.write(f"Coordinate range: X: {min(p[0] for p in positions):.3f} to {max(p[0] for p in positions):.3f}")
                                st.write(f"                  Y: {min(p[1] for p in positions):.3f} to {max(p[1] for p in positions):.3f}")
                                st.write(f"                  Z: {min(p[2] for p in positions):.3f} to {max(p[2] for p in positions):.3f}")
                        else:
                            st.warning("⚠️ Could not parse structure from output file")
                    except Exception as e:
                        st.warning(f"⚠️ Could not load structure from output file: {e}")
                        import traceback
                        st.text(traceback.format_exc())
                
                # どちらからも読み込めない場合
                if not optimized_mol:
                    st.error("❌ Could not load optimized structure from any source")
                    st.write("Available files:")
                    if result.get('arc_file'):
                        arc_exists = os.path.exists(result['arc_file'])
                        arc_size = os.path.getsize(result['arc_file']) if arc_exists else 0
                        st.write(f"- ARC file: {result['arc_file']} (exists: {arc_exists}, size: {arc_size} bytes)")
                    if result.get('output_file'):
                        out_exists = os.path.exists(result['output_file'])
                        out_size = os.path.getsize(result['output_file']) if out_exists else 0
                        st.write(f"- Output file: {result['output_file']} (exists: {out_exists}, size: {out_size} bytes)")
                else:
                    st.info(f"📊 Structure loaded from: {structure_source}")
                
                # 主要な計算結果を目立つように表示
                if 'summary' in result and 'key_results' in result['summary']:
                    key_results = result['summary']['key_results']
                    
                    st.subheader("🔬 Optimization Results")
                    
                    # メトリクス形式で主要結果を表示
                    col1, col2, col3, col4 = st.columns(4)
                    
                    with col1:
                        if 'heat_of_formation' in key_results:
                            st.metric("生成熱", key_results['heat_of_formation'], help="Heat of Formation")
                    
                    with col2:
                        if 'homo_energy' in key_results:
                            st.metric("HOMO", key_results['homo_energy'], help="Highest Occupied Molecular Orbital Energy")
                    
                    with col3:
                        if 'lumo_energy' in key_results:
                            st.metric("LUMO", key_results['lumo_energy'], help="Lowest Unoccupied Molecular Orbital Energy")
                    
                    with col4:
                        if 'homo_lumo_gap' in key_results:
                            st.metric("HOMO-LUMO Gap", key_results['homo_lumo_gap'], help="Energy gap between HOMO and LUMO")
                    
                    # 追加の結果（双極子モーメント、電子エネルギーなど）
                    if 'dipole_moment' in key_results or 'electronic_energy' in key_results:
                        st.subheader("📊 Additional Properties")
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            if 'dipole_moment' in key_results:
                                st.metric("双極子モーメント", key_results['dipole_moment'], help="Dipole Moment")
                        
                        with col2:
                            if 'electronic_energy' in key_results:
                                st.metric("電子エネルギー", key_results['electronic_energy'], help="Electronic Energy")
                
                # XYZ座標の表示
                if optimized_mol:
                    st.subheader("📊 Structure Coordinates")
                    
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        st.write("**Initial Structure (XYZ)**")
                        try:
                            # 初期構造のXYZ座標を生成
                            initial_xyz = f"{initial_mol.GetNumAtoms()}\n"
                            initial_xyz += "Initial Structure\n"
                            
                            if initial_mol.GetNumConformers() > 0:
                                conf = initial_mol.GetConformer()
                                for i in range(initial_mol.GetNumAtoms()):
                                    atom = initial_mol.GetAtomWithIdx(i)
                                    pos = conf.GetAtomPosition(i)
                                    initial_xyz += f"{atom.GetSymbol():<2} {pos.x:>12.6f} {pos.y:>12.6f} {pos.z:>12.6f}\n"
                            else:
                                initial_xyz += "No coordinates available\n"
                            
                            st.text_area("Initial XYZ", initial_xyz, height=300, key="initial_xyz")
                            
                            # ダウンロードボタン
                            st.download_button(
                                label="Download Initial XYZ",
                                data=initial_xyz,
                                file_name=f"{compound_name}_initial.xyz",
                                mime="text/plain"
                            )
                        except Exception as e:
                            st.error(f"Error generating initial XYZ: {e}")
                    
                    with col2:
                        st.write("**Optimized Structure (XYZ)**")
                        try:
                            # 最適化構造のXYZ座標を生成
                            optimized_xyz = f"{optimized_mol.GetNumAtoms()}\n"
                            optimized_xyz += "Optimized Structure\n"
                            
                            if optimized_mol.GetNumConformers() > 0:
                                conf = optimized_mol.GetConformer()
                                for i in range(optimized_mol.GetNumAtoms()):
                                    atom = optimized_mol.GetAtomWithIdx(i)
                                    pos = conf.GetAtomPosition(i)
                                    optimized_xyz += f"{atom.GetSymbol():<2} {pos.x:>12.6f} {pos.y:>12.6f} {pos.z:>12.6f}\n"
                            else:
                                optimized_xyz += "No coordinates available\n"
                            
                            st.text_area("Optimized XYZ", optimized_xyz, height=300, key="optimized_xyz")
                            
                            # ダウンロードボタン
                            st.download_button(
                                label="Download Optimized XYZ",
                                data=optimized_xyz,
                                file_name=f"{compound_name}_optimized.xyz",
                                mime="text/plain"
                            )
                        except Exception as e:
                            st.error(f"Error generating optimized XYZ: {e}")
                
                # 3D構造表示
                st.subheader("🧬 3D Structure Comparison")
                
                try:
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        st.write("**Initial Structure (3D)**")
                        if initial_mol and initial_mol.GetNumConformers() > 0:
                            initial_handler = MoleculeHandler(initial_mol, input_type="rdkit")
                            initial_molblock = initial_handler.generate_3d_molblock()
                            initial_viewer = py3Dmol.view(width=400, height=400)
                            initial_viewer.addModel(initial_molblock, "mol")
                            initial_viewer.setStyle({"stick": {}})
                            initial_viewer.zoomTo()
                            stmol.showmol(initial_viewer, height=400)
                        else:
                            st.warning("Initial structure has no 3D coordinates")
                    
                    with col2:
                        st.write("**Optimized Structure (3D)**")
                        if optimized_mol and optimized_mol.GetNumConformers() > 0:
                            optimized_handler = MoleculeHandler(optimized_mol, input_type="rdkit")
                            optimized_molblock = optimized_handler.generate_3d_molblock()
                            optimized_viewer = py3Dmol.view(width=400, height=400)
                            optimized_viewer.addModel(optimized_molblock, "mol")
                            optimized_viewer.setStyle({"stick": {}})
                            optimized_viewer.zoomTo()
                            stmol.showmol(optimized_viewer, height=400)
                        else:
                            st.warning("Optimized structure has no 3D coordinates")
                    
                except Exception as e:
                    st.error(f"Error displaying 3D structures: {e}")
                    import traceback
                    st.text(traceback.format_exc())
                
                # エネルギー比較の解析
                try:
                    st.subheader("⚡ Energy Comparison")
                    
                    # 初期構造の単点計算を実行
                    with st.spinner("初期構造の単点エネルギー計算を実行中..."):
                        initial_energy_result = calculator.single_point_energy(
                            theory=theory,
                            charge=charge,
                            multiplicity=multiplicity,
                            title=f"{compound_name}_initial_single_point"
                        )
                    
                    # エネルギー比較を表示
                    col1, col2, col3 = st.columns(3)
                    
                    initial_hof = None
                    optimized_hof = None
                    
                    with col1:
                        st.write("**初期構造**")
                        if (initial_energy_result.get('success') and 
                            'heat_of_formation' in initial_energy_result):
                            initial_hof = initial_energy_result['heat_of_formation']
                            st.metric("生成熱", f"{initial_hof:.3f} kcal/mol", help="Initial Structure Heat of Formation")
                        else:
                            st.error("初期構造のエネルギー計算に失敗")
                    
                    with col2:
                        st.write("**最適化構造**")
                        if 'heat_of_formation' in result:
                            optimized_hof = result['heat_of_formation']
                            st.metric("生成熱", f"{optimized_hof:.3f} kcal/mol", help="Optimized Structure Heat of Formation")
                        elif ('summary' in result and 'key_results' in result['summary'] and 
                              'heat_of_formation' in result['summary']['key_results']):
                            # summary から取得（文字列形式の可能性があるため解析）
                            hof_str = result['summary']['key_results']['heat_of_formation']
                            try:
                                optimized_hof = float(hof_str.split()[0])
                                st.metric("生成熱", f"{optimized_hof:.3f} kcal/mol", help="Optimized Structure Heat of Formation")
                            except:
                                st.metric("生成熱", hof_str, help="Optimized Structure Heat of Formation")
                        else:
                            st.error("最適化構造のエネルギー情報が見つからない")
                    
                    with col3:
                        st.write("**エネルギー差**")
                        if initial_hof is not None and optimized_hof is not None:
                            energy_diff = optimized_hof - initial_hof
                            delta_symbol = "🔻" if energy_diff < 0 else "🔺" if energy_diff > 0 else "➡️"
                            st.metric("ΔHf", f"{energy_diff:.3f} kcal/mol", 
                                    delta=f"{delta_symbol} {'安定化' if energy_diff < 0 else '不安定化' if energy_diff > 0 else '変化なし'}",
                                    help="最適化後の生成熱変化（負の値は安定化を示す）")
                            
                            # エネルギー変化の解釈
                            if abs(energy_diff) < 0.1:
                                st.info("💡 エネルギー変化は僅かです。初期構造が既に安定でした。")
                            elif energy_diff < -1.0:
                                st.success("💡 大きな安定化が見られます。構造最適化により分子が大幅に安定になりました。")
                            elif energy_diff < 0:
                                st.success("💡 適度な安定化が見られます。構造最適化により分子が安定になりました。")
                            else:
                                st.warning("💡 エネルギーが上昇しています。計算設定を確認してください。")
                        else:
                            st.error("エネルギー差を計算できませんでした")
                
                except Exception as energy_error:
                    st.error(f"エネルギー比較エラー: {energy_error}")

                # 構造変化の解析
                try:
                    st.subheader("📏 Structural Changes")
                    
                    # 分子の基本情報を比較
                    col1, col2 = st.columns(2)
                    with col1:
                        st.write("**Initial Structure Info:**")
                        st.write(f"原子数: {initial_mol.GetNumAtoms()}")
                        st.write(f"結合数: {initial_mol.GetNumBonds()}")
                        st.write(f"環数: {initial_mol.GetRingInfo().NumRings()}")
                    
                    with col2:
                        st.write("**Optimized Structure Info:**")
                        st.write(f"原子数: {optimized_mol.GetNumAtoms()}")
                        st.write(f"結合数: {optimized_mol.GetNumBonds()}")
                        st.write(f"環数: {optimized_mol.GetRingInfo().NumRings()}")
                    
                    # RMSD計算（エラーハンドリング付き）
                    rmsd_calculated = False
                    if initial_mol.GetNumAtoms() == optimized_mol.GetNumAtoms():
                        try:
                            from rdkit.Chem import rdMolAlign
                            
                            # 分子のアライメントを試行
                            # まず、同じ分子かどうかをチェック
                            initial_smiles = Chem.MolToSmiles(initial_mol)
                            optimized_smiles = Chem.MolToSmiles(optimized_mol)
                            
                            if initial_smiles == optimized_smiles:
                                # 同じ分子の場合、直接RMSD計算
                                rmsd = rdMolAlign.AlignMol(optimized_mol, initial_mol)
                                st.metric("RMSD", f"{rmsd:.4f} Å", help="Root Mean Square Deviation between initial and optimized structures")
                                rmsd_calculated = True
                            else:
                                # 異なる分子構造の場合、座標ベースでの距離計算
                                st.warning("⚠️ Molecular structures differ - calculating coordinate-based displacement")
                                
                                if (initial_mol.GetNumConformers() > 0 and 
                                    optimized_mol.GetNumConformers() > 0 and
                                    initial_mol.GetNumAtoms() == optimized_mol.GetNumAtoms()):
                                    
                                    initial_conf = initial_mol.GetConformer()
                                    optimized_conf = optimized_mol.GetConformer()
                                    
                                    total_displacement = 0.0
                                    for i in range(initial_mol.GetNumAtoms()):
                                        initial_pos = initial_conf.GetAtomPosition(i)
                                        optimized_pos = optimized_conf.GetAtomPosition(i)
                                        
                                        dx = optimized_pos.x - initial_pos.x
                                        dy = optimized_pos.y - initial_pos.y
                                        dz = optimized_pos.z - initial_pos.z
                                        
                                        displacement = (dx*dx + dy*dy + dz*dz)**0.5
                                        total_displacement += displacement
                                    
                                    avg_displacement = total_displacement / initial_mol.GetNumAtoms()
                                    st.metric("Average Atomic Displacement", f"{avg_displacement:.4f} Å", 
                                            help="Average displacement of atoms between initial and optimized structures")
                                    rmsd_calculated = True
                        
                        except Exception as rmsd_error:
                            st.warning(f"Could not calculate RMSD: {rmsd_error}")
                            
                            # フォールバック: 座標ベースの比較
                            try:
                                if (initial_mol.GetNumConformers() > 0 and 
                                    optimized_mol.GetNumConformers() > 0):
                                    
                                    initial_conf = initial_mol.GetConformer()
                                    optimized_conf = optimized_mol.GetConformer()
                                    
                                    # 重心の変化
                                    initial_center = [0.0, 0.0, 0.0]
                                    optimized_center = [0.0, 0.0, 0.0]
                                    
                                    for i in range(initial_mol.GetNumAtoms()):
                                        initial_pos = initial_conf.GetAtomPosition(i)
                                        optimized_pos = optimized_conf.GetAtomPosition(i)
                                        
                                        initial_center[0] += initial_pos.x
                                        initial_center[1] += initial_pos.y
                                        initial_center[2] += initial_pos.z
                                        
                                        optimized_center[0] += optimized_pos.x
                                        optimized_center[1] += optimized_pos.y
                                        optimized_center[2] += optimized_pos.z
                                    
                                    n_atoms = initial_mol.GetNumAtoms()
                                    initial_center = [x/n_atoms for x in initial_center]
                                    optimized_center = [x/n_atoms for x in optimized_center]
                                    
                                    center_displacement = ((optimized_center[0] - initial_center[0])**2 + 
                                                         (optimized_center[1] - initial_center[1])**2 + 
                                                         (optimized_center[2] - initial_center[2])**2)**0.5
                                    
                                    st.metric("Center of Mass Displacement", f"{center_displacement:.4f} Å",
                                            help="Displacement of molecular center of mass")
                                    rmsd_calculated = True
                                    
                            except Exception as fallback_error:
                                st.error(f"Could not perform any structural comparison: {fallback_error}")
                    else:
                        st.warning(f"⚠️ Cannot calculate RMSD: Different number of atoms (Initial: {initial_mol.GetNumAtoms()}, Optimized: {optimized_mol.GetNumAtoms()})")
                    
                    # 追加の幾何学的情報
                    if not rmsd_calculated:
                        st.write("**Individual Structure Analysis:**")
                        
                    # 結合長の統計（可能であれば）
                    try:
                        if optimized_mol.GetNumBonds() > 0:
                            bond_lengths = []
                            conf = optimized_mol.GetConformer()
                            for bond in optimized_mol.GetBonds():
                                atom1_idx = bond.GetBeginAtomIdx()
                                atom2_idx = bond.GetEndAtomIdx()
                                
                                pos1 = conf.GetAtomPosition(atom1_idx)
                                pos2 = conf.GetAtomPosition(atom2_idx)
                                
                                distance = ((pos1.x - pos2.x)**2 + 
                                           (pos1.y - pos2.y)**2 + 
                                           (pos1.z - pos2.z)**2)**0.5
                                bond_lengths.append(distance)
                            
                            if bond_lengths:
                                avg_bond_length = sum(bond_lengths) / len(bond_lengths)
                                min_bond_length = min(bond_lengths)
                                max_bond_length = max(bond_lengths)
                                
                                st.write("**Bond Length Statistics (Optimized):**")
                                col1, col2, col3 = st.columns(3)
                                with col1:
                                    st.metric("Average", f"{avg_bond_length:.3f} Å")
                                with col2:
                                    st.metric("Minimum", f"{min_bond_length:.3f} Å")
                                with col3:
                                    st.metric("Maximum", f"{max_bond_length:.3f} Å")
                    except Exception as bond_error:
                        st.info("Could not analyze bond lengths")
                    
                except Exception as e:
                    st.warning(f"Could not perform structural analysis: {e}")
                    import traceback
                    st.text(traceback.format_exc())
                
                # デバッグ情報（折りたたみ状態で表示）
                with st.expander("🔧 Debug Information", expanded=False):
                    st.write("**Calculation Result Keys:**")
                    st.write(list(result.keys()))
                    st.write("**Return Code:**", result.get('return_code', 'Unknown'))
                    st.write("**Success Flag:**", result.get('success', False))
                    
                    # ファイルの存在確認
                    input_file = result.get('input_file', '')
                    output_file = result.get('output_file', '')
                    arc_file = result.get('arc_file', '')
                    
                    st.write("**File Information:**")
                    st.write(f"- Input file: {input_file}")
                    st.write(f"- Input file exists: {os.path.exists(input_file) if input_file else False}")
                    st.write(f"- Output file: {output_file}")
                    st.write(f"- Output file exists: {os.path.exists(output_file) if output_file else False}")
                    st.write(f"- ARC file: {arc_file}")
                    st.write(f"- ARC file exists: {os.path.exists(arc_file) if arc_file else False}")
                    
                    # 作業ディレクトリの内容を表示
                    work_directory = result.get('work_directory', work_dir)
                    st.write(f"- Work directory: {work_directory}")
                    if os.path.exists(work_directory):
                        files = os.listdir(work_directory)
                        st.write(f"- Files in work directory: {files}")
                    else:
                        st.write("- Work directory does not exist")
                    
                    # 標準出力・エラー出力
                    if result.get('stdout'):
                        st.text_area("STDOUT", result['stdout'], height=150, key="debug_stdout")
                    if result.get('stderr'):
                        st.text_area("STDERR", result['stderr'], height=150, key="debug_stderr")
                    
                    # サマリー情報の表示
                    if 'summary' in result:
                        st.write("**Summary:**")
                        st.json(result['summary'])
                    
                    # 出力ファイル情報の表示
                    if 'output_files' in result:
                        st.write("**Generated Files:**")
                        output_files = result['output_files']
                        generated_files = []
                        for file_type, file_info in output_files.items():
                            if file_info.get('exists', False):
                                generated_files.append(f"✅ {file_type.upper()}: {file_info['size_human']} - {file_info['path']}")
                        
                        if generated_files:
                            for file_info in generated_files:
                                st.write(file_info)
                        else:
                            st.write("No output files were generated")
                    
                    # 実行されたコマンドの表示
                    if 'command_executed' in result:
                        st.write(f"**Command executed:** `{result['command_executed']}`")
            
            else:
                st.error("❌ MOPAC optimization failed!")
                
                # 失敗の詳細分析
                st.write("**Failure Analysis:**")
                if not result.get('success', False):
                    st.write("- Success flag is False")
                if result.get('return_code', -1) != 0:
                    st.write(f"- Non-zero return code: {result.get('return_code')}")
                if not result.get('output_file'):
                    st.write("- No output file specified")
                elif not os.path.exists(result.get('output_file', '')):
                    st.write("- Output file does not exist")
                
                # エラーメッセージの表示
                if 'error' in result:
                    st.error(f"Error: {result['error']}")
                
                # STDERRの表示
                if result.get('stderr'):
                    st.text_area("Error Output", result['stderr'], height=200, key="error_output")
                
                # 失敗の場合は以降の処理をスキップ
                st.stop()
            
            # 成功時の追加詳細情報
            with st.expander("📄 Calculation Details"):
                st.write(f"**Input File**: {os.path.basename(result['input_file'])}")
                if result.get('output_file'):
                    st.write(f"**Output File**: {os.path.basename(result['output_file'])}")
                if result.get('arc_file'):
                    st.write(f"**ARC File**: {os.path.basename(result['arc_file'])}")
                
                # 計算設定の表示
                st.write(f"**Theory**: {theory}")
                st.write(f"**Charge**: {charge}")
                st.write(f"**Multiplicity**: {multiplicity}")
                st.write(f"**Gradient norm**: {gnorm} kcal/mol/Å")
                if optimize_cartesian:
                    st.write("**Coordinate system**: Cartesian")
                else:
                    st.write("**Coordinate system**: Internal")
                if ts_search:
                    st.write("**Calculation type**: Transition State Search")
                else:
                    st.write("**Calculation type**: Energy Minimization")
                if use_solvent and solvent:
                    st.write(f"**Solvent**: {solvent['name']} (ε = {solvent['epsilon']:.4f})")
                else:
                    st.write("**Solvent**: Gas phase")
                
                # 正常終了の確認
                if result.get('normal_termination'):
                    st.success("Calculation terminated normally")
                else:
                    st.warning("Calculation may not have terminated normally")
                
                # エラーメッセージがある場合
                if 'errors' in result:
                    st.error("Errors found in calculation:")
                    for error in result['errors']:
                        st.write(f"- {error}")
            
            # 出力ファイルの内容表示（オプション）
            with st.expander("📋 Output File Content"):
                if result.get('output_file') and os.path.exists(result['output_file']):
                    try:
                        with open(result['output_file'], 'r', encoding='utf-8', errors='ignore') as f:
                            output_content = f.read()
                        st.text_area("MOPAC Output", output_content, height=300, key="mopac_output_content")
                    except Exception as e:
                        st.error(f"Error reading output file: {e}")
                else:
                    st.write("Output file not available")

    except Exception as e:
        st.error(f"Error processing molecule: {e}")

# 使用方法の説明
with st.expander("Usage Information"):
    st.markdown("""
    ### MOPAC Geometry Optimization
    
    このページではMOPAC半経験的量子化学計算ソフトウェアを使用して分子構造の最適化を行います。
    
    **理論手法:**
    - **PM7**: 最新の半経験的法。高精度で汎用的
    - **PM6**: 以前の主力モデル。PM7より少し粗い
    - **AM1**: 古典的なモデル（古いが軽量）
    - **MNDO**: 最も基本的な手法（教材向き）
    
    **最適化オプション:**
    - **PRECISE**: 高精度モード（より厳密な収束判定）
    - **GNORM**: 勾配収束判定値（デフォルト: 1.0 kcal/mol/Å）
    - **CARTESIAN**: デカルト座標系での最適化（通常は内部座標系）
    - **TS**: 遷移状態構造探索（エネルギー極大点を探索）
    - **SYMMETRY**: 分子対称性を利用した最適化
    - **ISOTOPE**: 同位体効果を考慮
    
    **配座探索:**
    - **UFF**: Universal Force Field（軽量、汎用的）
    - **MMFF**: Merck Molecular Force Field（高精度、有機分子向け）
    - 最適化前に複数の配座を生成し、最適な開始点を選択
    - 局所最小値への収束を回避
    
    **複数配座同時計算:**
    - **Calculate multiple conformers**: 複数配座を同時にMOPAC計算
    - **Max conformers to calculate**: 計算する最大配座数（計算時間との兼ね合い）
    - **Energy threshold**: 表示する配座のエネルギー閾値（最安定配座からの相対エネルギー）
    - **Show all conformer results**: 全配座の詳細結果を表示
    - 最安定配座の自動選択と詳細表示
    - 配座間のエネルギー比較とランキング表示
    - 初期構造（分子力場）と最適化後（MOPAC）のエネルギー変化追跡
    
    **溶媒効果:**
    - **COSMO**: COnductor-like Screening MOdel（連続誘電体モデル）
    - 60種類以上の溶媒から選択可能
    - カスタム誘電率（ε）の指定も可能
    - 溶媒中での構造最適化
    
    **出力される結果:**
    - 最適化前後の構造比較（3D表示）
    - 構造変化の解析（RMSD）
    - エネルギー情報（生成熱、電子エネルギー）
    - 分子軌道エネルギー（HOMO/LUMO）
    - 双極子モーメント
    - 幾何学的パラメータの変化
    - 複数配座のエネルギーランキング
    - 最安定配座の詳細情報
    - 各配座の構造最適化結果比較
    
    **注意:**
    - MOPACがシステムにインストールされている必要があります
    - 遷移状態探索（TS）は通常の最適化より困難で時間がかかります
    - 大きな分子や複雑な系では計算時間が長くなる場合があります
    - 配座探索を使用すると計算時間が延長されますが、より適切な最適化結果が得られます
    - 複数配座計算では各配座を独立して最適化するため、計算時間が配座数に比例して増加します
    - カルテシアン座標系は内部座標系で問題が生じた場合の代替手段です
    - 複数配座計算では、最安定配座が初期の分子力場計算の結果と異なる場合があります
    """)
