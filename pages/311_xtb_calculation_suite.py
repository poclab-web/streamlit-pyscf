"""
xTBを使用した包括的な量子化学計算スイート

機能:
- xTB（extended tight-binding）による高速な半経験的量子化学計算
- シングルポイント計算、構造最適化、振動計算
- 溶媒効果の考慮（ALPB溶媒モデル）
- 熱力学パラメータの計算
"""

import streamlit as st
import stmol
import py3Dmol
import os
import time
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Draw

from utils.module import load_css
from utils.xtb_ui import display_xtb_status, require_xtb, display_gfn_selector, display_solvent_selector
from logic.molecule_handler import MoleculeHandler
from logic.xtb_calculation import run_xtb_calculation

# カスタムCSSを適用
load_css("config/styles.css")

st.title("🧪 xTB計算")

st.markdown("extended tight-binding（xTB）による高速量子化学計算")

# xTBのインストール状況を確認
# require_xtb()

# # xTBの状況表示
# with st.expander("🔧 xTB Status", expanded=False):
#     display_xtb_status(key_suffix="_main")

st.divider()

# ユーザー入力
st.header("🧮 分子の入力")
input_type = st.selectbox("入力タイプを選択", ["SMILES", "XYZ"])
atom_input = st.text_area(
    "分子構造を入力",
    "CCO" if input_type == "SMILES" else "C     -2.385271    1.648249    0.000003\nH     -3.181824    1.061240    0.407216\nH     -2.536116    1.773363   -1.051896\nH     -2.370100    2.607079    0.474676\nO     -1.453042    1.151314    0.170017\nH     -1.466511    0.300016   -0.251422",
    key="molecular_input"
)

# SDFファイルアップロード
with st.expander("📁 SDFファイルアップロード", expanded=False):
    uploaded_file = st.file_uploader(
        "SDFファイルをアップロード", 
        type=['sdf'], 
        help="分子構造データファイル（SDF形式）をアップロードできます",
        key="sdf_uploader"
    )
    
    if uploaded_file is not None:
        atom_input = uploaded_file.getvalue().decode('utf-8')
        input_type = "SDF"
        st.success("✅ SDFファイルを読み込みました！")

# 分子構造を処理
handler = None
if st.button("🔄 分子構造を生成", type="primary", key="generate_structure_btn"):
    try:
        handler = MoleculeHandler(atom_input, input_type=input_type.lower())
        if not handler.mol:
            raise ValueError("Invalid molecular input. Please check your format.")
        
        # セッションステートに保存
        st.session_state.handler = handler
        st.success("✅ 分子構造を生成しました！")
        st.rerun()
    except Exception as e:
        st.error(f"❌ 分子解析エラー: {e}")

# セッションステートから分子ハンドラを取得
if "handler" in st.session_state:
    handler = st.session_state.handler
    
    # 分子が読み込まれている場合の表示と計算設定
    if handler and handler.mol is not None:
        st.divider()
        
        # 分子情報の表示
        with st.container():
            st.subheader("📊 分子情報")
            
            col1, col2, col3 = st.columns(3)
            
            with col1:
                # RDKitから直接分子式を取得
                try:
                    molecular_formula = rdMolDescriptors.CalcMolFormula(handler.mol)
                    st.metric("分子式", molecular_formula)
                except:
                    st.metric("分子式", "N/A")
            
            with col2:
                try:
                    molecular_weight = rdMolDescriptors.CalcExactMolWt(handler.mol)
                    st.metric("分子量", f"{molecular_weight:.2f}")
                except:
                    st.metric("分子量", "N/A")
            
            with col3:
                num_atoms = handler.mol.GetNumAtoms()
                st.metric("原子数", num_atoms)

        # 分子構造の2D/3D表示
        with st.expander("🔍 分子構造の確認", expanded=True):
            view_col1, view_col2 = st.columns(2)
            
            with view_col1:
                st.subheader("2D構造")
                try:
                    # 2D画像を生成して表示
                    img = Draw.MolToImage(handler.mol, size=(300, 300))
                    st.image(img, caption="2D分子構造", use_container_width=True)
                except Exception as e:
                    st.error(f"2D構造表示エラー: {e}")
            
            with view_col2:
                st.subheader("3D構造")
                try:
                    # 3D構造をMOLブロック形式で取得
                    mol_block = handler.generate_3d_molblock()
                    if mol_block:
                        # py3Dmolを使用した3D表示
                        viewer = py3Dmol.view(width=400, height=300)
                        viewer.addModel(mol_block, "mol")
                        viewer.setStyle({"stick": {}})
                        viewer.zoomTo()
                        stmol.showmol(viewer, height=300)
                    else:
                        st.warning("3D構造を生成できませんでした")
                except Exception as e:
                    st.error(f"3D構造表示エラー: {e}")

        st.divider()

        # 計算設定
        with st.container():
            st.subheader("⚙️ xTB計算設定")
            
            # 計算タイプの選択（ボタン式）
            st.write("**計算タイプ**")
            calc_type_cols = st.columns(4)
            
            calc_types = [
                ("Single Point", "sp", "エネルギー計算のみ"),
                ("Optimization", "opt", "構造最適化"),
                ("Frequency", "freq", "振動解析"),
                ("Opt + Freq", "opt+freq", "最適化 + 振動解析")
            ]
            
            # デフォルト選択
            if "selected_calc_type" not in st.session_state:
                st.session_state.selected_calc_type = "Single Point"
            
            for i, (display_name, value, description) in enumerate(calc_types):
                with calc_type_cols[i]:
                    if st.button(
                        display_name,
                        key=f"calc_type_{value}",
                        help=description,
                        type="primary" if st.session_state.selected_calc_type == display_name else "secondary"
                    ):
                        st.session_state.selected_calc_type = display_name
            
            calculation_type = st.session_state.selected_calc_type
            st.info(f"選択中: **{calculation_type}**")
            
            st.divider()
            
            # GFNモデルの選択（ボタン式）
            st.write("**GFNモデル**")
            gfn_cols = st.columns(3)
            
            gfn_models = [
                ("GFN0", 0, "最も高速、低精度"),
                ("GFN1", 1, "バランス良好"),
                ("GFN2", 2, "高精度、重い")
            ]
            
            # デフォルト選択
            if "selected_gfn" not in st.session_state:
                st.session_state.selected_gfn = 1
            
            for i, (display_name, value, description) in enumerate(gfn_models):
                with gfn_cols[i]:
                    if st.button(
                        display_name,
                        key=f"gfn_{value}",
                        help=description,
                        type="primary" if st.session_state.selected_gfn == value else "secondary"
                    ):
                        st.session_state.selected_gfn = value
            
            gfn = st.session_state.selected_gfn
            st.info(f"選択中: **GFN{gfn}** - {gfn_models[gfn][2]}")
            
            st.divider()
            
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                # 溶媒の選択
                solvent = display_solvent_selector(key_suffix="_xtb_calc")
            
            with col2:
                # 電荷とスピン多重度
                charge = st.number_input("電荷", value=0, step=1, help="分子の電荷を指定", key="charge_input")
            
            with col3:
                uhf = st.number_input("不対電子数", value=0, min_value=0, step=1, help="不対電子の数（スピン多重度-1）", key="uhf_input")
            
            with col4:
                # 計算精度の設定
                accuracy = st.selectbox(
                    "計算精度", 
                    options=[0.1, 0.5, 1.0, 2.0, 5.0, 10.0],
                    index=2,  # デフォルトは1.0
                    help="計算精度（低い値ほど高精度・時間がかかる）",
                    key="accuracy_input"
                )

        st.divider()

        # 計算実行
        if st.button("🚀 計算を開始", type="primary", key="start_calculation_btn"):
            # 分子データから適切なディレクトリ名を生成
            try:
                inchikey = Chem.MolToInchiKey(handler.mol)
                directory = f"data/{inchikey}"
            except:
                import time
                timestamp = str(int(time.time()))
                directory = f"data/molecule_{timestamp}"
            
            # 計算設定の表示
            with st.expander("📝 計算設定の確認", expanded=True):
                st.write(f"• 計算タイプ: {calculation_type}")
                st.write(f"• GFNモデル: GFN{gfn}")
                st.write(f"• 電荷: {charge}")
                st.write(f"• 不対電子数: {uhf}")
                if solvent:
                    st.write(f"• 溶媒: {solvent}")
                else:
                    st.write("• 溶媒: None (Gas phase)")
                st.write(f"• 入力形式: SDF → XYZ (if needed)")
            
            with st.spinner("Running xTB calculation..."):
                # 溶媒設定
                s = solvent  # display_solvent_selectorはNoneまたは溶媒名を返す
                
                # 計算タイプの変換
                if calculation_type == "Optimization":
                    calc_type = "opt"
                elif calculation_type == "Frequency":
                    calc_type = "freq"
                elif calculation_type == "Opt + Freq":
                    calc_type = "opt+freq"
                else:
                    calc_type = "sp"
                
                # セグメンテーションフォルトの原因となりそうな条件をチェック
                potential_issues = []
                if num_atoms > 200:
                    potential_issues.append("Very large molecule (>200 atoms)")
                if molecular_weight > 1000:
                    potential_issues.append("Very high molecular weight (>1000)")
                
                # ラジカルや不安定な構造をチェック
                try:
                    if handler.is_radical():
                        potential_issues.append("Radical species detected")
                except:
                    pass
                    
                if potential_issues:
                    st.warning(f"⚠️ Potential stability issues detected: {', '.join(potential_issues)}")
                    st.info("Switching to most conservative settings...")
                    # 最も保守的な設定に強制変更
                    calc_type = "sp"
                    gfn = 0
                    s = None
                
                # 直接最終作業ディレクトリを使用
                final_work_dir = os.path.join(directory, "xtb_results")
                os.makedirs(final_work_dir, exist_ok=True)
                
                st.info(f"🗂️ Using work directory: {final_work_dir}")
                
                # xTB計算の実行（SDF形式を使用、失敗時はXYZ形式でリトライ）
                result = run_xtb_calculation(
                    handler, 
                    final_work_dir, 
                    calculation_type=calc_type,
                    gfn=gfn, 
                    charge=charge, 
                    uhf=uhf, 
                    solvent=s,
                    accuracy=accuracy,
                    input_format="sdf"  # SDF形式を指定
                )
                
                # SDF形式でセグメンテーションフォルトが発生した場合、XYZ形式でリトライ
                if not result['success'] and 'SIGSEGV' in result.get('error', ''):
                    st.warning("SDF format caused segmentation fault. Retrying with XYZ format...")
                    result = run_xtb_calculation(
                        handler, 
                        final_work_dir, 
                        calculation_type=calc_type,
                        gfn=gfn, 
                        charge=charge, 
                        uhf=uhf, 
                        solvent=s,
                        accuracy=accuracy,
                        input_format="xyz"  # XYZ形式でリトライ
                    )
                    
                    # XYZ形式でもセグメンテーションフォルトが発生した場合の追加対策
                    if not result['success'] and 'SIGSEGV' in result.get('error', ''):
                        st.warning("XYZ format also failed. Trying simplified settings...")
                        
                        # より安全な設定でリトライ（溶媒なし、GFN0、精度緩和）
                        result = run_xtb_calculation(
                            handler, 
                            final_work_dir, 
                            calculation_type="sp",  # シングルポイント計算のみ
                            gfn=0,  # 最もシンプルなGFN0
                            charge=charge, 
                            uhf=uhf, 
                            solvent=None,  # 溶媒効果なし
                            accuracy=10.0,  # 精度を下げる
                            input_format="xyz"
                        )
                        
                        if result['success']:
                            st.info("✅ Calculation succeeded with simplified settings (no solvent, GFN0, single point)")
                        else:
                            st.error("⚠️ All retry attempts failed. This may be due to:")
                            st.error("- xTB version compatibility issues")
                            st.error("- Problematic molecular structure")
                            st.error("- System memory issues")
                            st.error("- xTB installation problems")
                
                # 計算結果はすでに最終ディレクトリに保存されている
                if result['success']:
                    st.success(f"✅ Results saved to: {final_work_dir}")
                else:
                    st.error("❌ Calculation failed - no results to save.")

        # 結果の表示
        # 入力ファイル名の表示（成功・失敗に関わらず表示）
        if 'result' in locals():
            st.info(f"📁 Input file used for xTB: `{result.get('input_file', 'N/A')}`")
            
            # 実行されたコマンドの表示
            if 'command_executed' in result:
                st.code(f"$ {result['command_executed']}", language="bash")
            
            if result['success']:
                st.success("✅ Calculation completed successfully!")
                
                # エネルギー結果の表示
                if 'energy' in result and result['energy'] is not None:
                    energy = result['energy']
                    col1, col2, col3 = st.columns(3)
                    
                    with col1:
                        st.metric("エネルギー (Hartree)", f"{energy:.6f}")
                    
                    with col2:
                        # kcal/molに変換
                        energy_kcal = energy * 627.509474
                        st.metric("エネルギー (kcal/mol)", f"{energy_kcal:.2f}")
                    
                    with col3:
                        # kJ/molに変換
                        energy_kj = energy * 2625.4996394
                        st.metric("エネルギー (kJ/mol)", f"{energy_kj:.2f}")
                
                # Moldenファイルの表示
                if 'molden_file' in result and result['molden_file']:
                    st.subheader("📁 出力ファイル")
                    st.success(f"✅ Moldenファイル: `{result['molden_file']}`")
                    
                    # ファイルサイズの表示
                    try:
                        import os
                        file_size = os.path.getsize(result['molden_file'])
                        if file_size > 1024 * 1024:
                            size_str = f"{file_size / (1024 * 1024):.2f} MB"
                        elif file_size > 1024:
                            size_str = f"{file_size / 1024:.2f} KB"
                        else:
                            size_str = f"{file_size} bytes"
                        st.info(f"📊 ファイルサイズ: {size_str}")
                    except:
                        pass
                    
                    # Moldenファイルの内容の一部を表示（オプション）
                    with st.expander("👀 Moldenファイルの内容を確認", expanded=False):
                        try:
                            with open(result['molden_file'], 'r') as f:
                                molden_content = f.read()
                                # 最初の50行程度を表示
                                lines = molden_content.split('\n')
                                preview_lines = lines[:50]
                                preview_content = '\n'.join(preview_lines)
                                if len(lines) > 50:
                                    preview_content += f"\n... (残り {len(lines) - 50} 行)"
                                st.text(preview_content)
                        except Exception as e:
                            st.error(f"ファイル読み込みエラー: {e}")
                elif 'molden_file' in result and result['molden_file'] is None:
                    st.warning("⚠️ Moldenファイルが生成されませんでした")
                
                # 追加の計算結果（双極子モーメント、HOMO-LUMOギャップなど）
                additional_results = []
                if 'dipole_moment' in result and result['dipole_moment'] is not None:
                    additional_results.append(f"双極子モーメント: {result['dipole_moment']:.3f} Debye")
                
                if 'homo_lumo_gap' in result and result['homo_lumo_gap'] is not None:
                    additional_results.append(f"HOMO-LUMOギャップ: {result['homo_lumo_gap']:.3f} eV")
                
                if additional_results:
                    st.subheader("📋 追加の計算結果")
                    for res in additional_results:
                        st.write(f"• {res}")
                
                # 振動計算結果の表示
                if 'frequencies' in result and result['frequencies']:
                    st.subheader("🎵 振動解析結果")
                    
                    frequencies = result['frequencies']
                    
                    # 振動数の統計表示
                    col1, col2, col3, col4 = st.columns(4)
                    with col1:
                        st.metric("振動モード数", len(frequencies))
                    with col2:
                        if frequencies:
                            st.metric("最低振動数", f"{min(frequencies):.1f} cm⁻¹")
                    with col3:
                        if frequencies:
                            st.metric("最高振動数", f"{max(frequencies):.1f} cm⁻¹")
                    with col4:
                        # 実振動数の数をカウント
                        real_frequencies = [f for f in frequencies if f > 0]
                        st.metric("実振動数", len(real_frequencies))
                    
                    # 虚振動の詳細解析
                    imaginary_frequencies = [f for f in frequencies if f < 0]
                    if imaginary_frequencies:
                        st.error(f"⚠️ 虚振動が検出されました: {len(imaginary_frequencies)}個")
                        
                        # 虚振動の詳細表示
                        with st.expander("🔍 虚振動の詳細", expanded=True):
                            st.markdown("**検出された虚振動:**")
                            for i, freq in enumerate(imaginary_frequencies):
                                st.write(f"• モード {i+1}: {freq:.1f} cm⁻¹")
                            
                            st.warning("""
                            **虚振動の意味:**
                            - **1個の虚振動**: 遷移状態（TS）の可能性
                            - **複数の虚振動**: 不安定な構造、または最適化不足
                            - **大きな虚振動** (|ν| > 200 cm⁻¹): 構造に重大な問題
                            """)
                            
                            st.info("""
                            **対処方法:**
                            - 構造最適化をより厳密に実行
                            - より良い初期構造から開始
                            - 異なるGFNレベルで計算
                            - 溶媒効果の有無を変更
                            """)
                    else:
                        st.success("✅ 虚振動は検出されませんでした（安定構造）")
                    
                    # 振動数の分布表示
                    with st.expander("📊 振動数の分布", expanded=False):
                        import pandas as pd
                        
                        # 振動数を範囲別に分類
                        freq_ranges = {
                            "虚振動 (< 0)": len([f for f in frequencies if f < 0]),
                            "低振動 (0-500)": len([f for f in frequencies if 0 <= f < 500]),
                            "中振動 (500-1500)": len([f for f in frequencies if 500 <= f < 1500]),
                            "伸縮振動 (1500-4000)": len([f for f in frequencies if 1500 <= f < 4000]),
                            "高振動 (≥4000)": len([f for f in frequencies if f >= 4000])
                        }
                        
                        # データフレームとして表示
                        df = pd.DataFrame(list(freq_ranges.items()), columns=['振動範囲', 'モード数'])
                        st.dataframe(df, use_container_width=True)
                        
                        # 全振動数のリスト表示
                        with st.expander("📋 全振動数リスト", expanded=False):
                            # 振動数を昇順にソート
                            sorted_frequencies = sorted(frequencies)
                            
                            # 5列で表示
                            cols = st.columns(5)
                            for i, freq in enumerate(sorted_frequencies):
                                with cols[i % 5]:
                                    color = "🔴" if freq < 0 else "🟢"
                                    st.write(f"{color} {freq:.1f} cm⁻¹")
                    
                    # 特徴的な振動の同定
                    with st.expander("🔬 特徴的な振動の同定", expanded=False):
                        characteristic_vibrations = []
                        
                        for freq in frequencies:
                            if freq < 0:
                                characteristic_vibrations.append(f"{freq:.1f} cm⁻¹: 虚振動")
                            elif 0 <= freq < 200:
                                characteristic_vibrations.append(f"{freq:.1f} cm⁻¹: 低振動（並進・回転様）")
                            elif 200 <= freq < 800:
                                characteristic_vibrations.append(f"{freq:.1f} cm⁻¹: 変角振動")
                            elif 800 <= freq < 1300:
                                characteristic_vibrations.append(f"{freq:.1f} cm⁻¹: C-C, C-N, C-O伸縮")
                            elif 1300 <= freq < 1800:
                                characteristic_vibrations.append(f"{freq:.1f} cm⁻¹: 変形振動")
                            elif 1800 <= freq < 2000:
                                characteristic_vibrations.append(f"{freq:.1f} cm⁻¹: C=C, C=N伸縮")
                            elif 2000 <= freq < 2300:
                                characteristic_vibrations.append(f"{freq:.1f} cm⁻¹: C≡C, C≡N伸縮")
                            elif 2500 <= freq < 3100:
                                characteristic_vibrations.append(f"{freq:.1f} cm⁻¹: C-H伸縮")
                            elif 3100 <= freq < 3700:
                                characteristic_vibrations.append(f"{freq:.1f} cm⁻¹: N-H, O-H伸縮")
                            elif freq >= 3700:
                                characteristic_vibrations.append(f"{freq:.1f} cm⁻¹: 高周波O-H, N-H伸縮")
                        
                        for vib in characteristic_vibrations[:20]:  # 最初の20個まで表示
                            st.write(f"• {vib}")
                        
                        if len(characteristic_vibrations) > 20:
                            st.write(f"... および他 {len(characteristic_vibrations) - 20} 個のモード")
                    
                    # 熱力学パラメータの表示
                    thermo_data = []
                    if 'zero_point_energy' in result and result['zero_point_energy'] is not None:
                        thermo_data.append(("ゼロ点エネルギー", result['zero_point_energy'], "Hartree"))
                    if 'thermal_correction' in result and result['thermal_correction'] is not None:
                        thermo_data.append(("熱補正", result['thermal_correction'], "Hartree"))
                    if 'enthalpy_correction' in result and result['enthalpy_correction'] is not None:
                        thermo_data.append(("エンタルピー補正", result['enthalpy_correction'], "Hartree"))
                    if 'entropy' in result and result['entropy'] is not None:
                        thermo_data.append(("エントロピー", result['entropy'], "cal/(mol·K)"))
                    if 'free_energy_correction' in result and result['free_energy_correction'] is not None:
                        thermo_data.append(("自由エネルギー補正", result['free_energy_correction'], "Hartree"))
                    
                    if thermo_data:
                        st.subheader("🌡️ 熱力学パラメータ")
                        for name, value, unit in thermo_data:
                            st.write(f"• {name}: {value:.6f} {unit}")
                        
                        # 温度とエネルギーの詳細表示
                        if 'zero_point_energy' in result and result['zero_point_energy'] is not None:
                            zpe_kcal = result['zero_point_energy'] * 627.509474
                            st.info(f"💡 ゼロ点エネルギー: {zpe_kcal:.2f} kcal/mol")
                    
                    # 構造の安定性評価
                    st.subheader("🎯 構造安定性の評価")
                    if not imaginary_frequencies:
                        st.success("✅ **安定な極小構造**: 虚振動なし")
                        st.info("この構造は局所的または全体的極小に対応します")
                    elif len(imaginary_frequencies) == 1:
                        st.warning("⚠️ **遷移状態候補**: 1つの虚振動")
                        st.info("この構造は化学反応の遷移状態である可能性があります")
                    else:
                        st.error("❌ **不安定構造**: 複数の虚振動")
                        st.info("さらなる構造最適化が必要です")
                
                # 振動計算のデバッグ情報表示（Frequency計算が実行された場合）
                elif calculation_type in ["Frequency", "Opt + Freq"] or calc_type in ["freq", "opt+freq"]:
                    st.warning("⚠️ 振動計算が実行されましたが、振動数データが見つかりませんでした")
                    
                    # デバッグ情報を表示
                    with st.expander("🔍 振動計算デバッグ情報", expanded=True):
                        st.write("**計算設定:**")
                        st.write(f"• 計算タイプ: {calculation_type}")
                        st.write(f"• 内部計算タイプ: {calc_type}")
                        
                        # 実行されたコマンドの詳細情報
                        st.write("**実行されたコマンド:**")
                        if 'command_executed' in result:
                            st.code(f"$ {result['command_executed']}", language="bash")
                            
                            # コマンドの解析
                            cmd_parts = result['command_executed'].split()
                            st.write("**コマンド解析:**")
                            
                            if "--hess" in cmd_parts:
                                st.write("✅ --hess: 振動計算が含まれています")
                            else:
                                st.write("❌ --hess: 振動計算が含まれていません")
                            
                            if "--thermo" in cmd_parts:
                                st.write("✅ --thermo: 熱化学補正が含まれています")
                            else:
                                st.write("❌ --thermo: 熱化学補正が含まれていません")
                            
                            if "--temp" in cmd_parts:
                                temp_idx = cmd_parts.index("--temp")
                                if temp_idx + 1 < len(cmd_parts):
                                    st.write(f"✅ --temp: 温度設定 {cmd_parts[temp_idx + 1]} K")
                            else:
                                st.write("❌ --temp: 温度設定がありません")
                            
                            if "--press" in cmd_parts:
                                press_idx = cmd_parts.index("--press")
                                if press_idx + 1 < len(cmd_parts):
                                    st.write(f"✅ --press: 圧力設定 {cmd_parts[press_idx + 1]} atm")
                            else:
                                st.write("❌ --press: 圧力設定がありません")
                        else:
                            st.write("コマンド情報が見つかりません")
                        
                        # 結果辞書の内容を確認
                        st.write("**結果辞書に含まれるキー:**")
                        result_keys = list(result.keys())
                        for key in result_keys:
                            st.write(f"• {key}: {type(result[key])}")
                        
                        # 振動関連のキーを特別に表示
                        freq_related_keys = [k for k in result_keys if 'freq' in k.lower() or 'vibr' in k.lower()]
                        if freq_related_keys:
                            st.write("**振動関連のキー:**")
                            for key in freq_related_keys:
                                st.write(f"• {key}: {result[key]}")
                        
                        # 作業ディレクトリのファイル一覧
                        if 'work_directory' in result:
                            work_dir = result['work_directory']
                            st.write(f"**作業ディレクトリ ({work_dir}) のファイル:**")
                            try:
                                import os
                                if os.path.exists(work_dir):
                                    files = os.listdir(work_dir)
                                    for file in sorted(files):
                                        st.write(f"• {file}")
                                else:
                                    st.write("作業ディレクトリが存在しません")
                            except Exception as e:
                                st.write(f"ファイル一覧取得エラー: {e}")
                        
                        st.info("""
                        **考えられる原因:**
                        - 振動計算が完了していない
                        - 出力ファイルの形式が予期したものと異なる
                        - xTBのバージョンによる出力形式の違い
                        - 振動計算でエラーが発生した
                        """)
                        
                        st.markdown("**推奨対処法:**")
                        st.markdown("1. より単純な分子で振動計算をテスト")
                        st.markdown("2. 最適化なしの振動計算（すでに最適化された構造から）")
                        st.markdown("3. 異なるGFNレベルで計算")
                        st.markdown("4. 詳細ログの確認")
                
                # 最適化された構造の表示（最適化計算の場合）
                if 'optimized_xyz' in result and result['optimized_xyz']:
                    st.subheader("🎯 最適化された構造")
                    
                    # 最適化されたXYZ座標を表示
                    with st.expander("📄 最適化されたXYZ座標", expanded=False):
                        st.text(result['optimized_xyz'])
                    
                    # 最適化された構造の3D表示を試行
                    try:
                        # XYZ形式から3D構造を生成して表示
                        xyz_content = result['optimized_xyz']
                        if xyz_content:
                            st.subheader("🔍 最適化された3D構造")
                            # XYZ内容をMOLブロック形式に変換して表示
                            # 注：完全な変換には追加のライブラリが必要な場合があります
                            with st.expander("3D構造ビューア", expanded=True):
                                st.code(xyz_content, language="text")
                                st.info("💡 より高度な3D表示には追加の可視化ツールが推奨されます")
                    except Exception as e:
                        st.warning(f"最適化構造の3D表示でエラーが発生しました: {e}")
                
                # 計算時間と作業ディレクトリの情報
                st.subheader("ℹ️ 計算情報")
                if 'work_directory' in result:
                    st.info(f"📁 作業ディレクトリ: `{result['work_directory']}`")
                
                # ファイル情報の表示
                info_items = []
                if 'molden_file' in result and result['molden_file']:
                    info_items.append(f"📄 Moldenファイル: `{result['molden_file']}`")
                if 'optimized_xyz_file' in result and result['optimized_xyz_file']:
                    info_items.append(f"📄 最適化構造ファイル: `{result['optimized_xyz_file']}`")
                
                if info_items:
                    st.write("**出力ファイル:**")
                    for item in info_items:
                        st.info(item)
                
                # 詳細ログの表示
                with st.expander("📝 詳細ログ", expanded=False):
                    if 'stdout' in result and result['stdout']:
                        st.subheader("標準出力")
                        st.text(result['stdout'])
                    
                    if 'stderr' in result and result['stderr']:
                        st.subheader("エラー出力")
                        st.text(result['stderr'])
            
            else:
                st.error("❌ 計算が失敗しました")
                
                # エラー情報の表示
                if 'error' in result:
                    st.error(f"エラー詳細: {result['error']}")
                
                # デバッグ情報の表示
                with st.expander("🔍 デバッグ情報", expanded=True):
                    if 'stdout' in result and result['stdout']:
                        st.subheader("標準出力")
                        st.text(result['stdout'])
                    
                    if 'stderr' in result and result['stderr']:
                        st.subheader("エラー出力") 
                        st.text(result['stderr'])
                    
                    if 'return_code' in result:
                        st.write(f"終了コード: {result['return_code']}")
                    
                    # opt+freq計算特有のデバッグ情報を表示
                    if calculation_type == "Opt + Freq" or calc_type == "opt+freq":
                        st.subheader("🔍 opt+freq計算のデバッグ情報")
                        
                        # 最適化結果の詳細
                        if 'optimization_result' in result:
                            opt_result = result['optimization_result']
                            st.write("**最適化ステップの結果:**")
                            st.write(f"• 成功: {opt_result.get('success', 'N/A')}")
                            if opt_result.get('energy'):
                                st.write(f"• エネルギー: {opt_result['energy']:.6f} Hartree")
                        
                        # 作業ディレクトリの内容
                        if 'work_directory_files' in result:
                            st.write("**作業ディレクトリ内のファイル:**")
                            for file in result['work_directory_files']:
                                st.write(f"• {file}")
                        
                        if 'debug_info' in result:
                            debug_info = result['debug_info']
                            st.write("**詳細デバッグ情報:**")
                            st.write(f"• 作業ディレクトリ: {debug_info.get('work_dir', 'N/A')}")
                            st.write(f"• ファイル総数: {debug_info.get('file_count', 'N/A')}")
                            st.write(f"• XYZファイル: {debug_info.get('xyz_files', [])}")
                            st.write(f"• coordファイル: {debug_info.get('coord_files', [])}")
                            st.write(f"• 検出された拡張子: {debug_info.get('all_extensions', [])}")
                    
                    # 実行されたコマンドの分析
                    if 'command_executed' in result:
                        st.subheader("実行コマンドの分析")
                        cmd = result['command_executed']
                        st.code(cmd, language="bash")
                        
                        # コマンドオプションの確認
                        if "opt+freq" in calculation_type or "opt+freq" in calc_type:
                            st.write("**注意**: opt+freq計算は2段階で実行されます")
                            st.write("1. 構造最適化 (--opt)")
                            st.write("2. 振動解析 (--hess)")
                            st.write("エラーがどの段階で発生したかを確認してください")
                
                # トラブルシューティングの提案
                st.subheader("🛠️ トラブルシューティング")
                st.markdown("""
                **一般的な解決方法:**
                - 分子構造を確認（不正な結合や原子配置）
                - より単純なGFNモデル（GFN0またはGFN1）を使用
                - 溶媒効果を無効にして気相計算を試行
                - シングルポイント計算から開始
                - 分子サイズを小さくする
                """)

else:
    st.info("👆 上記で分子構造を入力して「分子構造を生成」ボタンを押してください")
