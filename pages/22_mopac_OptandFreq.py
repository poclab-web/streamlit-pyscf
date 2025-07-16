"""
MOPAC半経験的量子化学計算ソフトウェアを使用して構造最適化と振動数計算を順次実行する。

機能:
- 計算に使用する理論手法を選択可能（PM7, PM6, AM1, MNDO）。
- 構造最適化後に振動数計算を自動実行。
- 最適化前後の構造比較、エネルギー変化、振動解析結果を表示。
- IR スペクトルの可視化、熱化学量（エンタルピー、Gibbs自由エネルギー）の計算。
- ゼロ点振動エネルギー（ZPVE）の補正を含む詳細な熱力学データ。
"""

import streamlit as st
import pandas as pd
import numpy as np
import os
import traceback
import matplotlib.pyplot as plt
import py3Dmol
import stmol
from rdkit import Chem

from utils.module import load_css
from logic.molecule_handler import MoleculeHandler
from logic.mopac_calculation import MopacCalculator, theory_options, check_mopac_installation


def _get_vibration_assignment(frequency):
    """
    振動数から振動の帰属を推定する簡易関数
    
    Args:
        frequency: 振動数 (cm⁻¹)
        
    Returns:
        str: 振動の帰属
    """
    if frequency < 0:
        return "Imaginary"
    elif frequency < 500:
        return "Bending/Rocking"
    elif frequency < 1000:
        return "C-C stretch/bend"
    elif frequency < 1300:
        return "C-O stretch"
    elif frequency < 1700:
        return "C-H bend"
    elif frequency < 2000:
        return "C=C/C=O stretch"
    elif frequency < 2500:
        return "C≡C/C≡N stretch"
    elif frequency < 3100:
        return "C-H stretch"
    elif frequency < 3700:
        return "O-H/N-H stretch"
    else:
        return "High frequency"


def _extract_thermochemical_data(output_content, temperature):
    """
    MOPAC出力ファイルから熱化学データを抽出する
    
    Args:
        output_content: 出力ファイルの内容
        temperature: 計算温度
        
    Returns:
        dict: 熱化学データ
    """
    thermo_data = {}
    
    try:
        lines = output_content.split('\n')
        
        for i, line in enumerate(lines):
            line_upper = line.upper()
            
            # ゼロ点振動エネルギー（複数の形式に対応）
            if any(phrase in line_upper for phrase in ["ZERO POINT ENERGY", "ZPE", "ZERO-POINT"]):
                try:
                    # 数値を抽出
                    import re
                    numbers = re.findall(r'[-+]?\d*\.?\d+', line)
                    if numbers:
                        zpve = float(numbers[-1])  # 最後の数値を使用
                        thermo_data["Zero Point Energy"] = f"{zpve:.3f} kcal/mol"
                except:
                    pass
            
            # エンタルピー
            elif any(phrase in line_upper for phrase in ["TOTAL ENTHALPY", "ENTHALPY"]) and "KCAL" in line_upper:
                try:
                    import re
                    numbers = re.findall(r'[-+]?\d*\.?\d+', line)
                    if numbers:
                        enthalpy = float(numbers[-1])
                        thermo_data["Total Enthalpy"] = f"{enthalpy:.3f} kcal/mol"
                except:
                    pass
            
            # Gibbs自由エネルギー
            elif any(phrase in line_upper for phrase in ["TOTAL FREE ENERGY", "GIBBS", "FREE ENERGY"]) and "KCAL" in line_upper:
                try:
                    import re
                    numbers = re.findall(r'[-+]?\d*\.?\d+', line)
                    if numbers:
                        gibbs = float(numbers[-1])
                        thermo_data["Gibbs Free Energy"] = f"{gibbs:.3f} kcal/mol"
                except:
                    pass
            
            # エントロピー
            elif any(phrase in line_upper for phrase in ["TOTAL ENTROPY", "ENTROPY"]) and ("CAL" in line_upper or "J" in line_upper):
                try:
                    import re
                    numbers = re.findall(r'[-+]?\d*\.?\d+', line)
                    if numbers:
                        entropy = float(numbers[-1])
                        unit = "cal/mol/K" if "CAL" in line_upper else "J/mol/K"
                        thermo_data["Total Entropy"] = f"{entropy:.3f} {unit}"
                except:
                    pass
            
            # 熱容量
            elif any(phrase in line_upper for phrase in ["HEAT CAPACITY", "CP"]) and ("CAL" in line_upper or "J" in line_upper):
                try:
                    import re
                    numbers = re.findall(r'[-+]?\d*\.?\d+', line)
                    if numbers:
                        heat_capacity = float(numbers[-1])
                        unit = "cal/mol/K" if "CAL" in line_upper else "J/mol/K"
                        thermo_data["Heat Capacity"] = f"{heat_capacity:.3f} {unit}"
                except:
                    pass
            
            # 生成熱
            elif any(phrase in line_upper for phrase in ["HEAT OF FORMATION", "HOF", "FORMATION"]) and "KCAL" in line_upper:
                try:
                    import re
                    numbers = re.findall(r'[-+]?\d*\.?\d+', line)
                    if numbers:
                        hof = float(numbers[-1])
                        thermo_data["Heat of Formation"] = f"{hof:.3f} kcal/mol"
                except:
                    pass
        
        # 振動数から熱化学量を計算（出力ファイルに含まれていない場合）
        if not thermo_data:
            # 振動数リストがあれば概算計算を試行
            frequencies = []
            for line in lines:
                if "CM-1" in line.upper() and any(word in line.upper() for word in ["FREQ", "VIBRATION"]):
                    try:
                        import re
                        numbers = re.findall(r'[-+]?\d*\.?\d+', line)
                        for num_str in numbers:
                            try:
                                freq_val = float(num_str)
                                if 0 < freq_val < 5000:
                                    frequencies.append(freq_val)
                            except ValueError:
                                continue
                    except:
                        continue
            
            if frequencies:
                # 簡易的なZPVE計算 (hc*ν/2 の和)
                h = 6.626e-34  # Planck constant
                c = 2.998e10   # Speed of light in cm/s
                na = 6.022e23  # Avogadro's number
                kcal_to_j = 4184
                
                zpve_j = sum(h * c * freq / 2 for freq in frequencies) * na
                zpve_kcal = zpve_j / kcal_to_j
                thermo_data["Zero Point Energy (calc)"] = f"{zpve_kcal:.3f} kcal/mol"
        
        return thermo_data
    
    except Exception as e:
        return {}


# カスタムCSSを適用
load_css("config/styles.css")

# 本文
st.title("MOPAC Geometry Optimization + Frequency Analysis")

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

# 計算設定
with st.expander("Calculation Settings"):
    charge = st.number_input("Molecular Charge", min_value=-10, max_value=10, value=0, step=1)
    multiplicity = st.number_input("Spin Multiplicity (2S + 1)", min_value=1, max_value=10, value=1, step=1)
    
    # 最適化の設定
    st.subheader("Optimization Settings")
    precise = st.checkbox("Use PRECISE mode", value=True, help="高精度モード（厳密な収束）")
    gnorm = st.number_input("Gradient norm convergence", min_value=0.1, max_value=10.0, value=1.0, step=0.1,
                           help="勾配収束判定値（kcal/mol/Å）")
    
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

    # 振動数計算の設定
    st.subheader("Frequency Calculation Settings")
    col1, col2 = st.columns(2)
    with col1:
        calculate_thermo = st.checkbox("Calculate thermochemical data", value=True,
                                     help="熱化学量（エンタルピー、Gibbs自由エネルギーなど）を計算")
        temperature = st.number_input("Temperature (K)", min_value=0.0, max_value=2000.0, value=298.15, step=1.0,
                                    help="熱化学量計算の温度")
    with col2:
        show_ir_spectrum = st.checkbox("Show IR spectrum", value=True,
                                     help="IRスペクトルを表示")
        ir_intensity_threshold = st.number_input("IR intensity threshold", min_value=0.0, max_value=100.0, value=1.0, step=0.1,
                                                help="表示するIR振動のピーク強度閾値")
    
    # 詳細解析オプション（新機能）
    st.subheader("🔬 Advanced Analysis Options")
    
    analysis_col1, analysis_col2, analysis_col3 = st.columns(3)
    
    with analysis_col1:
        st.write("**電子構造解析**")
        include_vectors = st.checkbox("分子軌道ベクトル", value=False, help="VECTORSキーワード: 分子軌道係数を出力")
        include_localize = st.checkbox("局在化軌道計算", value=False, help="LOCALIZEキーワード: 局在化軌道を計算")
        include_mullik = st.checkbox("Mulliken電荷解析", value=False, help="MULLIKキーワード: Mulliken電荷解析を実行")
        include_bonds = st.checkbox("結合次数解析", value=False, help="BONDSキーワード: 結合次数行列を出力")
    
    with analysis_col2:
        st.write("**静電特性**")
        include_esp = st.checkbox("静電ポテンシャル", value=False, help="ESPキーワード: 静電ポテンシャルを計算")
        include_polar = st.checkbox("分極率計算", value=False, help="POLARキーワード: 分極率と超分極率を計算")
        include_charges = st.checkbox("詳細電荷解析", value=False, help="CHARGESキーワード: 全電荷と各原子電荷を出力")
        include_denout = st.checkbox("電子密度出力", value=False, help="DENOUTキーワード: 電子密度データを出力")
    
    with analysis_col3:
        st.write("**相互作用・出力**")
        include_disp = st.checkbox("分散・水素結合", value=False, help="DISPキーワード: 分散エネルギーと水素結合寄与")
        include_super = st.checkbox("反応性指標", value=False, help="SUPERキーワード: 親核・親電子非局在化能")
        include_large = st.checkbox("拡張出力", value=False, help="LARGEキーワード: 詳細情報の拡張出力")
        include_pdbout = st.checkbox("PDB形式出力", value=False, help="PDBOUTキーワード: PDB形式での構造出力")
    
    # 溶媒効果の設定
    st.subheader("Solvent Effects")
    use_solvent = st.checkbox("Include solvent effects", value=False, 
                             help="COSMO溶媒効果モデルを使用")
    
    solvent = None
    if use_solvent:
        # 溶媒データを読み込み
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
                        solvent = {"name": solvent_choice, "epsilon": epsilon}
                        st.info(f"Dielectric constant: {epsilon:.4f}")
                        
        except Exception as e:
            st.error(f"Error loading solvent data: {e}")
            st.write("Using manual epsilon input")
            epsilon = st.number_input("Dielectric constant (ε)", 
                                    min_value=1.0, max_value=100.0, 
                                    value=78.36, step=0.1)
            solvent = {"name": "Custom", "epsilon": epsilon}

# 計算の実行
if st.button("Run MOPAC Optimization + Frequency"):
    try:
        handler = MoleculeHandler(atom_input, input_type=input_type.lower())
        if not handler.mol:
            raise ValueError("Invalid molecular input. Please check your format.")

        # 化合物名を取得
        compound_name = Chem.MolToInchiKey(handler.mol)
        smiles = Chem.MolToSmiles(handler.mol)

        # 初期構造の保存
        initial_mol = Chem.Mol(handler.mol)
        
        # ディレクトリの作成
        directory = os.path.join("data", compound_name)
        os.makedirs(directory, exist_ok=True)

        # 作業ディレクトリを指定
        work_dir = os.path.join(directory, "mopac_work")
        
        try:
            # MopacCalculatorのインスタンス作成
            calculator = MopacCalculator(handler, work_dir=work_dir)
            
            # 最適化計算のキーワード設定
            opt_keywords = []
            if precise:
                opt_keywords.append("PRECISE")
            if gnorm != 1.0:
                opt_keywords.append(f"GNORM={gnorm:.1f}")
            if optimize_cartesian:
                opt_keywords.append("CARTESIAN")
            if ts_search:
                opt_keywords.append("TS")
            if symmetry:
                opt_keywords.append("SYMMETRY")
            if isotope:
                opt_keywords.append("ISOTOPE")
            
            # 詳細解析オプションを最適化にも追加
            if include_vectors:
                opt_keywords.append("VECTORS")
            if include_localize:
                opt_keywords.append("LOCALIZE")
            if include_mullik:
                opt_keywords.append("MULLIK")
            if include_bonds:
                opt_keywords.append("BONDS")
            if include_esp:
                opt_keywords.append("ESP")
            if include_polar:
                opt_keywords.append("POLAR")
            if include_charges:
                opt_keywords.append("CHARGES")
            if include_denout:
                opt_keywords.append("DENOUT")
            if include_disp:
                opt_keywords.append("DISP")
            if include_super:
                opt_keywords.append("SUPER")
            if include_large:
                opt_keywords.append("LARGE")
            if include_pdbout:
                opt_keywords.append("PDBOUT")
            
            # 溶媒効果の追加
            if use_solvent and solvent:
                opt_keywords.append(f"COSMO EPS={solvent['epsilon']:.4f}")
            
            # ステップ1: 構造最適化の実行
            st.subheader("🔧 Step 1: Geometry Optimization")
            with st.spinner("Running geometry optimization..."):
                opt_result = calculator.optimize_geometry(
                    theory=theory,
                    charge=charge,
                    multiplicity=multiplicity,
                    keywords=opt_keywords,
                    title="Geometry Optimization for Frequency Analysis"
                )
            
            # 最適化結果の確認
            opt_success = (
                opt_result.get('success', False) and 
                opt_result.get('return_code', -1) == 0 and
                opt_result.get('output_file') and 
                os.path.exists(opt_result.get('output_file', ''))
            )
            
            if not opt_success:
                st.error("❌ Geometry optimization failed!")
                st.write("**Failure details:**")
                if 'error' in opt_result:
                    st.error(f"Error: {opt_result['error']}")
                if opt_result.get('stderr'):
                    st.text_area("Error Output", opt_result['stderr'], height=200)
                st.stop()
            
            st.success("✅ Geometry optimization completed successfully!")
            
            # 最適化結果の表示
            if 'summary' in opt_result and 'key_results' in opt_result['summary']:
                key_results = opt_result['summary']['key_results']
                
                col1, col2, col3, col4 = st.columns(4)
                
                with col1:
                    if 'heat_of_formation' in key_results:
                        st.metric("Heat of Formation", key_results['heat_of_formation'])
                    else:
                        st.metric("Heat of Formation", "N/A")
                
                with col2:
                    if 'electronic_energy' in key_results:
                        st.metric("Electronic Energy", key_results['electronic_energy'])
                    else:
                        st.metric("Electronic Energy", "N/A")
                
                with col3:
                    if 'homo_energy' in key_results:
                        st.metric("HOMO Energy", key_results['homo_energy'])
                    else:
                        st.metric("HOMO Energy", "N/A")
                
                with col4:
                    if 'lumo_energy' in key_results:
                        st.metric("LUMO Energy", key_results['lumo_energy'])
                    else:
                        st.metric("LUMO Energy", "N/A")
                        
                # 詳細解析データの表示
                if any([include_mullik, include_charges, include_polar, include_disp]):
                    st.subheader("🔬 Advanced Analysis Results")
                    
                    # Mulliken電荷解析
                    if include_mullik and 'mulliken_charges' in key_results:
                        with st.expander("Mulliken Charge Analysis"):
                            st.write(key_results['mulliken_charges'])
                    
                    # 電荷解析
                    if include_charges and 'atomic_charges' in key_results:
                        with st.expander("Detailed Charge Analysis"):
                            st.write(key_results['atomic_charges'])
                    
                    # 分極率データ
                    if include_polar and 'polarizability' in key_results:
                        with st.expander("Polarizability Data"):
                            st.write(key_results['polarizability'])
                    
                    # 分散・水素結合データ
                    if include_disp and 'dispersion_data' in key_results:
                        with st.expander("Dispersion and Hydrogen Bond Analysis"):
                            st.write(key_results['dispersion_data'])
            
            # 最適化された構造の読み込み
            optimized_mol = None
            if opt_result.get('arc_file') and os.path.exists(opt_result['arc_file']):
                try:
                    optimized_mol = calculator.read_optimized_structure(opt_result['arc_file'])
                    if optimized_mol:
                        st.info("📊 Optimized structure loaded from ARC file")
                except Exception as e:
                    st.warning(f"Could not load optimized structure: {e}")
            
            # ステップ2: 振動数計算の実行
            st.subheader("🎵 Step 2: Frequency Analysis")
            
            # 振動数計算用のキーワード設定
            freq_keywords = ["FORCE"]
            if calculate_thermo:
                freq_keywords.append("THERMO")
            if temperature != 298.15:
                freq_keywords.append(f"T={temperature:.2f}")
            
            # 詳細解析オプションを振動数計算にも追加
            if include_vectors:
                freq_keywords.append("VECTORS")
            if include_localize:
                freq_keywords.append("LOCALIZE")
            if include_mullik:
                freq_keywords.append("MULLIK")
            if include_bonds:
                freq_keywords.append("BONDS")
            if include_esp:
                freq_keywords.append("ESP")
            if include_polar:
                freq_keywords.append("POLAR")
            if include_charges:
                freq_keywords.append("CHARGES")
            if include_denout:
                freq_keywords.append("DENOUT")
            if include_disp:
                freq_keywords.append("DISP")
            if include_super:
                freq_keywords.append("SUPER")
            if include_large:
                freq_keywords.append("LARGE")
            if include_pdbout:
                freq_keywords.append("PDBOUT")
            
            # 溶媒効果の追加（最適化と同じ設定）
            if use_solvent and solvent:
                freq_keywords.append(f"COSMO EPS={solvent['epsilon']:.4f}")
            
            # 最適化された構造を振動数計算の初期構造として使用
            temp_handler = handler
            if optimized_mol and optimized_mol.GetNumConformers() > 0:
                try:
                    # 最適化された構造で新しいハンドラーを作成
                    from logic.molecule_handler import MoleculeHandler
                    
                    # 最適化された構造のXYZ形式を作成
                    optimized_xyz = f"{optimized_mol.GetNumAtoms()}\n"
                    optimized_xyz += "Optimized structure for frequency calculation\n"
                    
                    conf = optimized_mol.GetConformer()
                    for i in range(optimized_mol.GetNumAtoms()):
                        atom = optimized_mol.GetAtomWithIdx(i)
                        pos = conf.GetAtomPosition(i)
                        optimized_xyz += f"{atom.GetSymbol():2s}  {pos.x:12.6f}  {pos.y:12.6f}  {pos.z:12.6f}\n"
                    
                    # 最適化された構造で新しいMoleculeHandlerを作成
                    temp_handler = MoleculeHandler(optimized_xyz, input_type="xyz")
                    
                    # 新しいCalculatorを作成
                    freq_calculator = MopacCalculator(temp_handler, work_dir=work_dir)
                    
                    st.info("🔄 Using optimized structure for frequency calculation")
                    
                except Exception as e:
                    st.warning(f"Could not use optimized structure, using original: {e}")
                    freq_calculator = calculator
            else:
                freq_calculator = calculator
            
            with st.spinner("Running frequency calculation..."):
                # 振動数計算の入力ファイルを作成
                freq_input_file = freq_calculator.create_input_file(
                    theory=theory,
                    keywords=freq_keywords,
                    charge=charge,
                    multiplicity=multiplicity,
                    title="Frequency Analysis",
                    filename="frequency"
                )
                
                # 振動数計算を実行
                freq_result = freq_calculator.run_calculation(freq_input_file)
                freq_result['summary'] = freq_calculator.get_calculation_summary(freq_result)
            
            # 振動数計算結果の確認
            freq_success = (
                freq_result.get('success', False) and 
                freq_result.get('return_code', -1) == 0 and
                freq_result.get('output_file') and 
                os.path.exists(freq_result.get('output_file', ''))
            )
            
            # 振動数計算が失敗した場合の追加チェック
            if not freq_success:
                # より詳細なエラー解析
                error_msg = freq_result.get('stderr', '')
                output_content = ""
                
                if freq_result.get('output_file') and os.path.exists(freq_result.get('output_file')):
                    try:
                        with open(freq_result.get('output_file'), 'r') as f:
                            output_content = f.read()
                    except:
                        pass
                
                # グラディエントが大きすぎる場合の対処
                if ("GRADIENT IS TOO LARGE" in output_content.upper() or 
                    "GEOMETRY IS NOT AT A STATIONARY POINT" in output_content.upper()):
                    
                    st.warning("⚠️ Initial frequency calculation failed due to non-converged geometry")
                    st.info("🔧 Attempting more precise optimization before frequency calculation...")
                    
                    # より厳密な最適化を実行
                    precise_opt_keywords = opt_keywords.copy()
                    if "PRECISE" not in precise_opt_keywords:
                        precise_opt_keywords.append("PRECISE")
                    if not any("GNORM" in keyword for keyword in precise_opt_keywords):
                        precise_opt_keywords.append("GNORM=0.1")
                    
                    # 厳密な最適化では詳細解析オプションは一時的に除外（安定性を優先）
                    analysis_keywords = ["VECTORS", "LOCALIZE", "MULLIK", "BONDS", "ESP", "POLAR", 
                                       "CHARGES", "DENOUT", "DISP", "SUPER", "LARGE", "PDBOUT"]
                    precise_opt_keywords = [kw for kw in precise_opt_keywords if not any(ak in kw for ak in analysis_keywords)]
                    
                    with st.spinner("Running more precise optimization..."):
                        precise_opt_result = calculator.optimize_geometry(
                            theory=theory,
                            charge=charge,
                            multiplicity=multiplicity,
                            keywords=precise_opt_keywords,
                            title="Precise Optimization for Frequency Analysis"
                        )
                    
                    # 精密最適化が成功した場合、再度振動数計算を試行
                    if (precise_opt_result.get('success', False) and 
                        precise_opt_result.get('return_code', -1) == 0):
                        
                        st.success("✅ Precise optimization completed")
                        
                        # 精密最適化された構造を読み込み
                        if precise_opt_result.get('arc_file') and os.path.exists(precise_opt_result['arc_file']):
                            try:
                                precise_optimized_mol = calculator.read_optimized_structure(precise_opt_result['arc_file'])
                                if precise_optimized_mol:
                                    # 精密最適化された構造で振動数計算
                                    precise_xyz = f"{precise_optimized_mol.GetNumAtoms()}\n"
                                    precise_xyz += "Precisely optimized structure for frequency calculation\n"
                                    
                                    conf = precise_optimized_mol.GetConformer()
                                    for i in range(precise_optimized_mol.GetNumAtoms()):
                                        atom = precise_optimized_mol.GetAtomWithIdx(i)
                                        pos = conf.GetAtomPosition(i)
                                        precise_xyz += f"{atom.GetSymbol():2s}  {pos.x:12.6f}  {pos.y:12.6f}  {pos.z:12.6f}\n"
                                    
                                    precise_handler = MoleculeHandler(precise_xyz, input_type="xyz")
                                    precise_calculator = MopacCalculator(precise_handler, work_dir=work_dir)
                                    
                                    with st.spinner("Running frequency calculation on precisely optimized structure..."):
                                        precise_freq_input = precise_calculator.create_input_file(
                                            theory=theory,
                                            keywords=freq_keywords,  # 詳細解析オプション付きの振動数計算
                                            charge=charge,
                                            multiplicity=multiplicity,
                                            title="Frequency Analysis (Precise)",
                                            filename="frequency_precise"
                                        )
                                        
                                        freq_result = precise_calculator.run_calculation(precise_freq_input)
                                        freq_result['summary'] = precise_calculator.get_calculation_summary(freq_result)
                                        
                                        # 最適化された構造を更新
                                        optimized_mol = precise_optimized_mol
                            except Exception as e:
                                st.error(f"Error in precise optimization: {e}")
                    
                    # 再度成功チェック
                    freq_success = (
                        freq_result.get('success', False) and 
                        freq_result.get('return_code', -1) == 0 and
                        freq_result.get('output_file') and 
                        os.path.exists(freq_result.get('output_file', ''))
                    )
            
            if not freq_success:
                st.error("❌ Frequency calculation failed!")
                st.write("**Failure details:**")
                if 'error' in freq_result:
                    st.error(f"Error: {freq_result['error']}")
                if freq_result.get('stderr'):
                    st.text_area("Frequency Error Output", freq_result['stderr'], height=200)
                st.stop()
            
            st.success("✅ Frequency calculation completed successfully!")
            
            # 振動数結果の解析と表示
            st.subheader("📊 Frequency Analysis Results")
            
            # 出力ファイルから振動数データを読み込み
            freq_output_file = freq_result.get('output_file')
            if freq_output_file and os.path.exists(freq_output_file):
                try:
                    with open(freq_output_file, 'r') as f:
                        freq_output_content = f.read()
                    
                    # 振動数の抽出
                    frequencies = []
                    intensities = []
                    
                    # MOPACの出力から振動数とIR強度を抽出
                    lines = freq_output_content.split('\n')
                    
                    # 振動数データを探す（複数の可能な形式に対応）
                    for i, line in enumerate(lines):
                        # 方法1: DESCRIPTION OF VIBRATIONSセクション
                        if "DESCRIPTION OF VIBRATIONS" in line:
                            # このセクションの後の行を解析
                            for j in range(i+1, min(i+300, len(lines))):  # より多くの行をチェック
                                next_line = lines[j].strip()
                                if not next_line:
                                    continue
                                
                                # 振動数の行を探す（よりフレキシブルな検索）
                                if "FREQUENCY" in next_line and any(char.isdigit() for char in next_line):
                                    try:
                                        # 数値を抽出
                                        import re
                                        numbers = re.findall(r'[-+]?\d*\.?\d+', next_line)
                                        for num_str in numbers:
                                            try:
                                                freq_val = float(num_str)
                                                # 妥当な振動数範囲をチェック（負の値も含む虚振動数対応）
                                                if -1000 < freq_val < 5000 and abs(freq_val) > 1:
                                                    frequencies.append(freq_val)
                                                    intensities.append(1.0)  # デフォルト値
                                                    break  # 最初の妥当な数値のみ使用
                                            except ValueError:
                                                continue
                                    except (ValueError, IndexError):
                                        continue
                                
                                # セクションの終了をチェック
                                if ("THERMODYNAMIC" in next_line or 
                                    "TOTAL ENERGY" in next_line or
                                    "SCF CALCULATION" in next_line):
                                    break
                        
                        # 方法2: 直接的な振動数行の検出
                        elif ("VIBRATION" in line and "FREQUENCY" in line and 
                              any(char.isdigit() for char in line)):
                            try:
                                # 次の行に振動数がある場合
                                if i+1 < len(lines):
                                    freq_line = lines[i+1].strip()
                                    if "FREQUENCY" in freq_line:
                                        import re
                                        numbers = re.findall(r'[-+]?\d*\.?\d+', freq_line)
                                        for num_str in numbers:
                                            try:
                                                freq_val = float(num_str)
                                                if -1000 < freq_val < 5000 and abs(freq_val) > 1:
                                                    frequencies.append(freq_val)
                                                    intensities.append(1.0)
                                                    break
                                            except ValueError:
                                                continue
                            except:
                                continue
                        
                        # 方法3: シンプルな振動数リスト（行内に FREQUENCY と数値）
                        elif ("FREQUENCY" in line and 
                              any(char.isdigit() for char in line) and
                              "CM-1" not in line.upper()):  # CM-1が含まれていない場合のチェック
                            try:
                                # 行から数値を抽出
                                import re
                                numbers = re.findall(r'[-+]?\d*\.?\d+', line)
                                for num_str in numbers:
                                    try:
                                        freq_val = float(num_str)
                                        if -1000 < freq_val < 5000 and abs(freq_val) > 1:
                                            frequencies.append(freq_val)
                                            intensities.append(1.0)
                                            break  # 最初の妥当な数値のみ使用
                                    except ValueError:
                                        continue
                            except:
                                continue
                    
                    # 重複削除と並び替え
                    if frequencies:
                        freq_intensity_pairs = list(zip(frequencies, intensities))
                        # 重複削除（近い値は同じとみなす）
                        unique_pairs = []
                        for freq, intensity in freq_intensity_pairs:
                            is_duplicate = False
                            for existing_freq, _ in unique_pairs:
                                if abs(freq - existing_freq) < 0.1:  # 0.1 cm-1の差は同じとみなす
                                    is_duplicate = True
                                    break
                            if not is_duplicate:
                                unique_pairs.append((freq, intensity))
                        
                        frequencies, intensities = zip(*unique_pairs) if unique_pairs else ([], [])
                        frequencies = list(frequencies)
                        intensities = list(intensities)
                    
                    # 振動数データの表示
                    if frequencies:
                        st.write(f"**Found {len(frequencies)} vibrational modes**")
                        
                        # 振動数統計
                        col1, col2, col3, col4 = st.columns(4)
                        with col1:
                            st.metric("Number of Modes", len(frequencies))
                        with col2:
                            st.metric("Lowest Frequency", f"{min(frequencies):.1f} cm⁻¹")
                        with col3:
                            st.metric("Highest Frequency", f"{max(frequencies):.1f} cm⁻¹")
                        with col4:
                            # 虚振動数のチェック
                            imaginary_freqs = [f for f in frequencies if f < 0]
                            if imaginary_freqs:
                                st.metric("Imaginary Frequencies", len(imaginary_freqs), delta="⚠️")
                            else:
                                st.metric("Imaginary Frequencies", "0", delta="✅")
                        
                        # IRスペクトルの表示
                        if show_ir_spectrum and len(frequencies) == len(intensities):
                            st.subheader("📈 IR Spectrum")
                            
                            # IRスペクトルのプロット
                            fig, ax = plt.subplots(figsize=(12, 6))
                            
                            # 強度閾値でフィルタリング
                            filtered_data = [(f, i) for f, i in zip(frequencies, intensities) if i >= ir_intensity_threshold]
                            
                            if filtered_data:
                                filtered_freqs, filtered_intensities = zip(*filtered_data)
                                
                                # 棒グラフでスペクトルを描画
                                ax.bar(filtered_freqs, filtered_intensities, width=20, alpha=0.7, color='blue')
                                ax.set_xlabel('Frequency (cm⁻¹)')
                                ax.set_ylabel('IR Intensity')
                                ax.set_title(f'IR Spectrum (threshold: {ir_intensity_threshold})')
                                ax.grid(True, alpha=0.3)
                                
                                # X軸を逆順にする（一般的なIRスペクトルの表示方法）
                                ax.invert_xaxis()
                                
                                st.pyplot(fig)
                                
                                # 主要なピークの表示
                                strong_peaks = [(f, i) for f, i in filtered_data if i > max(filtered_intensities) * 0.1]
                                if strong_peaks:
                                    st.write("**Major IR peaks:**")
                                    peak_data = []
                                    for freq, intensity in sorted(strong_peaks, key=lambda x: x[1], reverse=True):
                                        peak_data.append({
                                            "Frequency (cm⁻¹)": f"{freq:.1f}",
                                            "Intensity": f"{intensity:.2f}",
                                            "Assignment": _get_vibration_assignment(freq)
                                        })
                                    
                                    df_peaks = pd.DataFrame(peak_data)
                                    st.dataframe(df_peaks, use_container_width=True)
                            else:
                                st.info(f"No IR peaks above threshold ({ir_intensity_threshold})")
                        
                        # 振動数テーブルの表示
                        with st.expander("📋 Complete Frequency Table"):
                            freq_data = []
                            for i, (freq, intensity) in enumerate(zip(frequencies, intensities)):
                                freq_data.append({
                                    "Mode": i + 1,
                                    "Frequency (cm⁻¹)": f"{freq:.2f}",
                                    "IR Intensity": f"{intensity:.4f}",
                                    "Type": "Imaginary" if freq < 0 else "Real",
                                    "Assignment": _get_vibration_assignment(freq)
                                })
                            
                            df_frequencies = pd.DataFrame(freq_data)
                            st.dataframe(df_frequencies, use_container_width=True)
                    
                    else:
                        st.warning("No frequency data found in the output file")
                        
                        # デバッグ情報を表示
                        with st.expander("🔍 Debug Information"):
                            st.write("**Searching for frequency data in output file...**")
                            
                            # 出力ファイルの一部を表示
                            debug_lines = []
                            for i, line in enumerate(lines):
                                if any(keyword in line.upper() for keyword in ["FREQ", "VIBRATION", "CM-1", "FORCE", "THERMO"]):
                                    debug_lines.append(f"Line {i+1}: {line.strip()}")
                            
                            if debug_lines:
                                st.write("**Lines containing frequency-related keywords:**")
                                for debug_line in debug_lines[:20]:  # 最初の20行のみ表示
                                    st.text(debug_line)
                                if len(debug_lines) > 20:
                                    st.write(f"... and {len(debug_lines) - 20} more lines")
                            else:
                                st.write("No frequency-related keywords found in output")
                            
                            # 振動数抽出の詳細情報
                            st.write(f"**Raw frequencies found**: {len(frequencies)}")
                            if frequencies:
                                st.write("**Raw frequency values:**")
                                st.text(", ".join([f"{f:.2f}" for f in frequencies[:20]]))
                                if len(frequencies) > 20:
                                    st.write(f"... and {len(frequencies) - 20} more")
                            
                            # 出力ファイルのサイズと場所
                            st.write(f"**Output file**: {freq_output_file}")
                            try:
                                file_size = os.path.getsize(freq_output_file)
                                st.write(f"**File size**: {file_size} bytes")
                            except:
                                st.write("**File size**: Could not determine")
                            
                            # 計算が正常終了したかチェック
                            normal_termination = any("NORMAL TERMINATION" in line.upper() or 
                                                   "CALCULATION COMPLETED" in line.upper() or
                                                   "JOB FINISHED" in line.upper() 
                                                   for line in lines)
                            st.write(f"**Normal termination**: {normal_termination}")
                            
                            # FORCE キーワードがあったかチェック
                            force_keyword_found = any("FORCE" in line.upper() for line in lines)
                            st.write(f"**FORCE keyword found**: {force_keyword_found}")
                            
                            # DESCRIPTION OF VIBRATIONSセクションの確認
                            vibration_section_found = any("DESCRIPTION OF VIBRATIONS" in line.upper() for line in lines)
                            st.write(f"**DESCRIPTION OF VIBRATIONS found**: {vibration_section_found}")
                            
                            # グラディエントエラーのチェック
                            gradient_errors = [line.strip() for line in lines if 
                                             "GRADIENT IS TOO LARGE" in line.upper() or
                                             "GEOMETRY IS NOT AT A STATIONARY POINT" in line.upper()]
                            if gradient_errors:
                                st.write("**Gradient errors found:**")
                                for error_line in gradient_errors:
                                    st.text(error_line)
                            
                            # 最初の50行を表示（診断用）
                            with st.expander("📄 First 50 lines of output file"):
                                for i, line in enumerate(lines[:50]):
                                    st.text(f"{i+1:3d}: {line}")
                            
                            # 最後の50行を表示（診断用）
                            with st.expander("📄 Last 50 lines of output file"):
                                for i, line in enumerate(lines[-50:], len(lines)-49):
                                    if i > 0:
                                        st.text(f"{i:3d}: {line}")
                    
                    # 熱化学データの抽出と表示
                    if calculate_thermo:
                        st.subheader("🌡️ Thermochemical Data")
                        
                        thermo_data = _extract_thermochemical_data(freq_output_content, temperature)
                        
                        if thermo_data:
                            col1, col2 = st.columns(2)
                            
                            with col1:
                                st.write("**Energy Components:**")
                                for key, value in thermo_data.items():
                                    if "energy" in key.lower() or "enthalpy" in key.lower():
                                        st.metric(key, value)
                            
                            with col2:
                                st.write("**Thermodynamic Properties:**")
                                for key, value in thermo_data.items():
                                    if "entropy" in key.lower() or "gibbs" in key.lower() or "heat_capacity" in key.lower():
                                        st.metric(key, value)
                        else:
                            st.info("Thermochemical data not found or could not be parsed")
                    
                except Exception as e:
                    st.error(f"Error analyzing frequency results: {e}")
                    st.text_area("Error Details", traceback.format_exc(), height=200)
            

            # 計算サマリーの表示
            with st.expander("📄 Calculation Summary"):
                st.write("### Optimization Results")
                st.json(opt_result.get('summary', {}))
                
                st.write("### Frequency Results")
                st.json(freq_result.get('summary', {}))
                
                # 設定情報の表示
                st.write("### Calculation Settings")
                settings = {
                    "Theory": theory,
                    "Charge": charge,
                    "Multiplicity": multiplicity,
                    "Temperature": f"{temperature} K",
                    "Precise mode": precise,
                    "Gradient norm": f"{gnorm} kcal/mol/Å",
                    "Coordinate system": "Cartesian" if optimize_cartesian else "Internal",
                    "Calculation type": "Transition State Search" if ts_search else "Energy Minimization",
                    "Solvent": solvent['name'] if use_solvent and solvent else "Gas phase",
                    "Thermochemistry": calculate_thermo,
                    "Optimization keywords": " ".join(opt_keywords),
                    "Frequency keywords": " ".join(freq_keywords)
                }
                
                # 使用された詳細解析オプション
                analysis_options = []
                if include_vectors: analysis_options.append("VECTORS")
                if include_localize: analysis_options.append("LOCALIZE")
                if include_mullik: analysis_options.append("MULLIK")
                if include_bonds: analysis_options.append("BONDS")
                if include_esp: analysis_options.append("ESP")
                if include_polar: analysis_options.append("POLAR")
                if include_charges: analysis_options.append("CHARGES")
                if include_denout: analysis_options.append("DENOUT")
                if include_disp: analysis_options.append("DISP")
                if include_super: analysis_options.append("SUPER")
                if include_large: analysis_options.append("LARGE")
                if include_pdbout: analysis_options.append("PDBOUT")
                
                if analysis_options:
                    settings["Advanced analysis"] = " ".join(analysis_options)
                else:
                    settings["Advanced analysis"] = "None"
                
                st.json(settings)
                
        except Exception as e:
            st.error(f"Error in calculation: {e}")
            st.text_area("Error Details", traceback.format_exc(), height=200)

        except Exception as e:
            st.error(f"Error in calculation: {e}")
            st.text_area("Error Details", traceback.format_exc(), height=200)

    except Exception as e:
        st.error(f"Error processing molecule: {e}")
        st.text_area("Error Details", traceback.format_exc(), height=200)


# 使用方法の説明
with st.expander("Usage Information"):
    st.markdown("""
    ### MOPAC Geometry Optimization + Frequency Analysis
    
    このページではMOPAC半経験的量子化学計算ソフトウェアを使用して分子構造の最適化と振動数計算を順次実行します。
    
    **計算の流れ:**
    1. **構造最適化**: 分子構造をエネルギー最小構造に最適化
    2. **振動数計算**: 最適化された構造で振動解析を実行
    3. **結果表示**: 構造変化、振動モード、IRスペクトル、熱化学データを表示
    
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
    
    **詳細解析オプション（新機能）:**
    - **VECTORS**: 分子軌道ベクトルの出力
    - **LOCALIZE**: 局在化軌道の計算
    - **MULLIK**: Mulliken電荷解析
    - **BONDS**: 結合次数行列の出力
    - **ESP**: 静電ポテンシャル計算
    - **POLAR**: 分極率・超分極率計算
    - **CHARGES**: 詳細電荷解析
    - **DENOUT**: 電子密度データ出力
    - **DISP**: 分散エネルギー・水素結合解析
    - **SUPER**: 親核・親電子非局在化能
    - **LARGE**: 詳細情報の拡張出力
    - **PDBOUT**: PDB形式での構造出力
    
    **振動数計算オプション:**
    - **FORCE**: 振動解析の実行（IRスペクトル計算）
    - **THERMO**: 熱化学量の計算（エンタルピー、Gibbs自由エネルギーなど）
    - **Temperature**: 熱化学量計算の温度設定（デフォルト: 298.15 K）
    
    **出力される結果:**
    - 最適化前後の構造比較（3D表示）
    - 構造変化の解析（RMSD）
    - エネルギー情報（生成熱、電子エネルギー）
    - 分子軌道エネルギー（HOMO/LUMO）
    - 振動数一覧（実振動数/虚振動数の判定）
    - IRスペクトル（可視化とピーク帰属）
    - 熱化学データ（ゼロ点エネルギー、エンタルピー、Gibbs自由エネルギー、エントロピー）
    - 主要なIRピークの帰属表
    
    **詳細解析結果（新機能）:**
    - **分子軌道データ**: 軌道係数と形状
    - **局在化軌道**: 化学結合の詳細解析
    - **Mulliken電荷**: 各原子の電荷分布
    - **結合次数**: 原子間結合の強さ
    - **静電ポテンシャル**: 分子表面の静電特性
    - **分極率**: 外部電場に対する応答特性
    - **分散・水素結合**: 非共有結合相互作用の詳細
    - **反応性指標**: 親核・親電子攻撃部位の予測
    - **電子密度データ**: 可視化用電子密度ファイル
    - **PDB構造**: タンパク質解析ソフト対応形式
    
    **振動数計算の意味:**
    - **実振動数（正の値）**: 安定な構造を示す
    - **虚振動数（負の値）**: 遷移状態や不安定構造を示す
    - **ゼロ点振動エネルギー（ZPVE）**: 絶対零度での残留振動エネルギー
    - **IRスペクトル**: 赤外分光法で観測される振動吸収
    
    **溶媒効果:**
    - **COSMO**: COnductor-like Screening MOdel（連続誘電体モデル）
    - 最適化と振動数計算の両方に同じ溶媒効果を適用
    - 溶媒中での構造最適化と振動解析
     **注意:**
    - MOPACがシステムにインストールされている必要があります
    - 振動数計算は最適化後の構造で実行されるため、計算時間が2倍程度かかります
    - 遷移状態探索（TS）の場合、1つの虚振動数が期待されます
    - 大きな分子では振動モード数が多くなり、計算時間が延長されます
    - IRスペクトルの強度閾値を調整して、表示するピークを制御できます
    - **詳細解析オプション**を多数選択すると計算時間が大幅に増加します
    - **ESP、POLAR、DISP**などの計算は特に時間がかかります
    - **大きな分子**では詳細解析オプションの使用は慎重に検討してください
    - **PDB出力**はタンパク質解析ソフト（PyMOL、ChimeraX等）での可視化に便利です
    """)
                        