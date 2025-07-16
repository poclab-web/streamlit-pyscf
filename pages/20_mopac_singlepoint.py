"""
MOPAC半経験的量子化学計算ソフトウェアを使用してシングルポイントエネルギーを計算する。

機能:
- 計算に使用する理論手法を選択可能（PM7, PM6, AM1, MNDO）。
- MOPAC計算を実行し、結果をインタラクティブに表示。
- 生成熱、電子エネルギー、HOMO/LUMOエネルギーを表示。
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
from logic.molecule_handler import MoleculeHandler
from logic.mopac_calculation import MopacCalculator, theory_options, check_mopac_installation

# カスタムCSSを適用
load_css("config/styles.css")

# 本文
st.title("MOPAC Single Point Calculation")

# MOPACインストール状況の表示
st.subheader("MOPAC Installation Status")
mopac_status = check_mopac_installation()

if mopac_status["installed"]:
    st.success(f"✅ MOPAC is installed at: `{mopac_status['path']}`")
        
    # エラーがある場合は警告として表示
    if mopac_status["error"]:
        st.warning(f"Note: {mopac_status['error']}")
else:
    st.error("❌ MOPAC is not installed or not found in PATH")
    st.write(f"Error: {mopac_status['error']}")
    st.markdown("""
    **MOPACをインストールするには:**
    1. [MOPAC公式サイト](http://openmopac.net/)からダウンロード
    2. `mopac`コマンドがターミナルで実行できることを確認
    """)

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
    st.subheader("Conformer Search")
    use_conformer_search = st.checkbox("Perform conformer search", value=False, 
                                      help="分子力場による配座探索を実行してから計算")
    
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
    
    # MOPACの追加オプション
    st.subheader("MOPAC Options")
    precise = st.checkbox("Use PRECISE mode", value=True, help="高精度モード（厳密な収束）")
    
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
    
    # 計算タイプの選択
    calculation_type = st.selectbox(
        "Calculation Type", 
        ["Single Point"],
        help="シングルポイント計算"
    )

# 分子構造を処理
handler = None
if st.button("Run MOPAC Calculation"):
    try:
        handler = MoleculeHandler(atom_input, input_type=input_type.lower())
        if not handler.mol:
            raise ValueError("Invalid molecular input. Please check your format.")

        # 化合物名を取得
        compound_name = Chem.MolToInchiKey(handler.mol)
        smiles = Chem.MolToSmiles(handler.mol)

        # 配座探索の実行
        if use_conformer_search:
            st.write("Performing conformer search...")
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
                    st.subheader("Conformer Energies")
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
                        if st.checkbox("Use lowest energy conformer for calculation", value=True):
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

        col1, col2 = st.columns(2)

        # Display 2D structure in the first column
        with col1:
            st.subheader("Input 2D Structure")
            handler.generate_2d_image(f"{directory}/molecule_2d.png")
            st.image(f"{directory}/molecule_2d.png", caption=Chem.MolToSmiles(handler.mol))

        # Display 3D structure in the second column
        with col2:
            if use_conformer_search and handler.mol.GetNumConformers() > 0:
                st.subheader("Optimized 3D Structure")
            else:
                st.subheader("Input 3D Structure")
            try:
                mol_block = handler.generate_3d_molblock()
                viewer = py3Dmol.view(width=400, height=400)
                viewer.addModel(mol_block, "mol")
                viewer.setStyle({"stick": {}})
                viewer.zoomTo()
                stmol.showmol(viewer, height=400)
                
                # 配座探索が実行された場合は情報を表示
                if use_conformer_search and handler.mol.GetNumConformers() > 0:
                    conf_id = 0  # 現在表示されている配座ID
                    if handler.mol.HasProp(f"Energy_{conf_id}"):
                        energy = float(handler.mol.GetProp(f"Energy_{conf_id}"))
                        st.caption(f"Energy: {energy:.6f} kcal/mol ({forcefield})")
                    else:
                        st.caption(f"Optimized with {forcefield}")
            except Exception as e:
                st.warning(f"Unable to generate 3D structure: {e}")

        # MOPAC計算の実行
        st.write("Running MOPAC calculation...")
        
        # 作業ディレクトリを指定
        work_dir = os.path.join(directory, "mopac_work")
        
        try:
            # MopacCalculatorのインスタンス作成
            calculator = MopacCalculator(handler, work_dir=work_dir)
            
            # 計算タイプに応じて実行
            if calculation_type == "Single Point":
                # シングルポイント計算
                keywords = ["1SCF"]  # 構造最適化なし
                if precise:
                    keywords.append("PRECISE")
                
                # 溶媒効果の追加
                if use_solvent and solvent:
                    keywords.append(f"COSMO EPS={solvent['epsilon']:.4f}")
                
                result = calculator.single_point_energy(
                    theory=theory,
                    charge=charge,
                    multiplicity=multiplicity,
                    keywords=keywords,
                    title="MOPAC Single Point Calculation"
                )
            else:
                # 構造最適化
                keywords = []
                if precise:
                    keywords.append("PRECISE")
                
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
            st.subheader("Calculation Results")
            
            # 実際の成功/失敗を判断
            actual_success = (
                result.get('success', False) and 
                result.get('return_code', -1) == 0 and
                result.get('output_file') and 
                os.path.exists(result.get('output_file', ''))
            )
            
            if actual_success:
                st.success("✅ MOPAC calculation completed successfully!")
                
                # 主要な計算結果を目立つように表示
                if 'summary' in result and 'key_results' in result['summary']:
                    key_results = result['summary']['key_results']
                    
                    st.subheader("🔬 Key Results")
                    
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
                
                # デバッグ情報（折りたたみ状態で表示）
                with st.expander("🔧 Debug Information", expanded=False):
                    st.write("**Calculation Result Keys:**")
                    st.write(list(result.keys()))
                    st.write("**Return Code:**", result.get('return_code', 'Unknown'))
                    st.write("**Success Flag:**", result.get('success', False))
                    
                    # ファイルの存在確認
                    input_file = result.get('input_file', '')
                    output_file = result.get('output_file', '')
                    
                    st.write("**File Information:**")
                    st.write(f"- Input file: {input_file}")
                    st.write(f"- Input file exists: {os.path.exists(input_file) if input_file else False}")
                    st.write(f"- Output file: {output_file}")
                    st.write(f"- Output file exists: {os.path.exists(output_file) if output_file else False}")
                    
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
                st.error("❌ MOPAC calculation failed!")
                
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
            st.error(f"Error running MOPAC calculation: {e}")
            import traceback
            st.text_area("Error Details", traceback.format_exc(), height=200, key="error_details")

    except Exception as e:
        st.error(f"Error processing molecule: {e}")

# 使用方法の説明
with st.expander("Usage Information"):
    st.markdown("""
    ### MOPAC Single Point Calculation
    
    このページではMOPAC半経験的量子化学計算ソフトウェアを使用してシングルポイントエネルギー計算を行います。
    
    **理論手法:**
    - **PM7**: 最新の半経験的法。高精度で汎用的
    - **PM6**: 以前の主力モデル。PM7より少し粗い
    - **AM1**: 古典的なモデル（古いが軽量）
    - **MNDO**: 最も基本的な手法（教材向き）
    
    **配座探索:**
    - **UFF**: Universal Force Field（軽量、汎用的）
    - **MMFF**: Merck Molecular Force Field（高精度、有機分子向け）
    - 複数の配座を生成し、最低エネルギー配座を自動選択
    - RMS閾値により類似配座を除去
    
    **溶媒効果:**
    - **COSMO**: COnductor-like Screening MOdel（連続誘電体モデル）
    - 60種類以上の溶媒から選択可能（水、アセトニトリル、メタノール等）
    - カスタム誘電率（ε）の指定も可能
    - 溶媒分子との相互作用を近似的に考慮
    
    **計算タイプ:**
    - **Single Point**: 入力された構造でのエネルギー計算のみ
    
    **出力される結果:**
    - 生成熱 (Heat of Formation) [kcal/mol]
    - 電子エネルギー [eV]
    - 双極子モーメント [Debye]
    - HOMO/LUMOエネルギー [eV]
    - HOMO-LUMOギャップ [eV]
    
    **注意:**
    - MOPACがシステムにインストールされている必要があります
    - 配座探索は計算時間を延長しますが、より安定な構造を得られます
    - 溶媒効果を含む計算は気相計算より時間がかかる場合があります
    - 大きな分子の場合、計算に時間がかかる場合があります
    """)
