"""
MOPAC計算結果の分子軌道可視化

機能:
- arcファイルまたはoutファイルから分子軌道データを読み取り
- HOMO/LUMO軌道の3D可視化
- 軌道エネルギー図の表示
- 電子密度表面の可視化
- 分子軌道係数の解析
"""

import streamlit as st
import pandas as pd
import os
import tempfile
import stmol

from utils.module import load_css
from utils.mopac_ui import require_mopac
from logic.mopac_visualization import (
    list_out_files,
    display_file_info,
    parse_mopac_out_file,
    create_molecular_orbital_plot,
    visualize_molecule_3d,
    create_orbital_dataframe
)

# カスタムCSSを適用
load_css("config/styles.css")

# ページタイトル
st.title("MOPAC分子軌道可視化")
st.markdown("MOPACのoutファイルから分子軌道情報を読み取り、可視化します。")

# MOPACインストール確認
require_mopac()

st.divider()

# ファイル選択方法
st.header("1. MOPAC計算結果ファイルの選択")

input_method = st.radio(
    "ファイルの選択方法",
    ["dataフォルダから選択", "ファイルをアップロード"],
    help="dataフォルダに保存されたファイルを選択するか、新しいファイルをアップロードできます"
)

data_file_path = None

if input_method == "dataフォルダから選択":
    # dataフォルダ内のMOPACファイルを探索
    mopac_files = list_out_files("data")
    total_files = len(mopac_files["out"])
    
    if total_files > 0:
        st.subheader(f"利用可能なMOPACファイル ({total_files}件)")
        st.info("mopac_workディレクトリ内のoutファイルのみを表示しています。")
        
        available_files = mopac_files["out"]
        
        if available_files:
            selected_file = st.selectbox(
                "ファイルを選択してください",
                available_files,
                format_func=lambda x: f"[MOPAC-OUT] {os.path.basename(x)} ({os.path.dirname(x).replace('data/', '')})"
            )
            
            if selected_file:
                data_file_path = selected_file
                
                # ファイル情報表示
                with st.expander("ファイル詳細情報"):
                    display_file_info(selected_file)
        else:
            st.warning("MOPACのoutファイルがありません。")
    else:
        st.warning("dataフォルダ内のmopac_workディレクトリにMOPACのoutファイルが見つかりませんでした。")
        st.info("MOPACの計算を実行してファイルを生成するか、ファイルアップロード機能を使用してください。")

else:  # ファイルをアップロード
    uploaded_file = st.file_uploader(
        "outファイルをアップロードしてください",
        type=['out'],
        help="MOPAC計算で生成されたoutファイルを選択してください"
    )
    
    if uploaded_file is not None:
        # 一時ファイルに保存
        with tempfile.NamedTemporaryFile(delete=False, suffix=f".{uploaded_file.name.split('.')[-1]}") as tmp_file:
            tmp_file.write(uploaded_file.getvalue())
            data_file_path = tmp_file.name

# ファイルが選択された場合の処理
if data_file_path is not None:
    try:
        # ファイル解析（outファイルのみ）
        data = parse_mopac_out_file(data_file_path)
        
        if not data['atoms']:
            st.error("ファイルから分子データを読み取れませんでした。")
        else:
            st.success(f"ファイルの読み取りが完了しました。原子数: {len(data['atoms'])}")
            
            # 基本情報表示
            st.header("2. 計算結果概要")
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("分子情報")
                st.write(f"原子数: {len(data['atoms'])}")
                st.write(f"原子構成: {', '.join(set(data['atoms']))}")
                
                if data['total_energy']:
                    st.write(f"全エネルギー: {data['total_energy']:.6f} eV")
                if data['heat_of_formation']:
                    st.write(f"生成熱: {data['heat_of_formation']:.6f} kcal/mol")
                if data['dipole_moment']:
                    st.write(f"双極子モーメント: {data['dipole_moment']:.6f} Debye")
            
            with col2:
                if data['orbital_energies']:
                    st.subheader("軌道情報")
                    st.write(f"軌道数: {len(data['orbital_energies'])}")
                    
                    if data['homo_index'] is not None:
                        homo_energy = data['orbital_energies'][data['homo_index']]
                        st.write(f"HOMO軌道: #{data['homo_index']+1} ({homo_energy:.3f} eV)")
                        
                        lumo_candidates = [e for e in data['orbital_energies'] if e > homo_energy]
                        if lumo_candidates:
                            lumo_energy = min(lumo_candidates)
                            lumo_index = data['orbital_energies'].index(lumo_energy)
                            st.write(f"LUMO軌道: #{lumo_index+1} ({lumo_energy:.3f} eV)")
                            
                            gap = lumo_energy - homo_energy
                            st.write(f"HOMO-LUMOギャップ: {gap:.3f} eV")
            
            # 3D分子構造表示
            st.header("3. 分子構造")
            
            if data['coordinates']:
                view = visualize_molecule_3d(data['atoms'], data['coordinates'])
                stmol.showmol(view, height=400, width=600)
            
            # 軌道エネルギー図
            st.header("4. Molecular Orbital Energy Diagram")
            
            if data['orbital_energies']:
                fig = create_molecular_orbital_plot(data['orbital_energies'], data['homo_index'])
                st.pyplot(fig)
                
                # 軌道エネルギー表
                st.subheader("軌道エネルギー詳細")
                orbital_df = create_orbital_dataframe(data['orbital_energies'], data['homo_index'])
                st.dataframe(orbital_df, use_container_width=True)
            
            # ダウンロードオプション
            st.header("5. 結果のダウンロード")
            
            if data['orbital_energies']:
                # CSVファイルとしてダウンロード
                orbital_df = create_orbital_dataframe(data['orbital_energies'], data['homo_index'])
                csv_data = orbital_df.to_csv(index=False, encoding='utf-8-sig')
                st.download_button(
                    label="軌道エネルギーデータをCSVでダウンロード",
                    data=csv_data,
                    file_name="orbital_energies.csv",
                    mime="text/csv"
                )
        
    except Exception as e:
        st.error(f"ファイル処理中にエラーが発生しました: {e}")
        st.exception(e)
    
    finally:
        # 一時ファイルの場合のみ削除（アップロードファイル）
        if input_method == "ファイルをアップロード" and hasattr(data_file_path, 'startswith') and '/tmp' in data_file_path:
            if os.path.exists(data_file_path):
                os.unlink(data_file_path)

else:
    if input_method == "dataフォルダから選択":
        st.info("利用可能なMOPACファイルがありません。計算を実行してからお試しください。")
    else:
        st.info("MOPACのoutファイルをアップロードしてください。")
    
    # 使用例の表示
    with st.expander("使用方法"):
        st.markdown("""
        ### 対応ファイル形式
        
        **outファイル**: MOPAC計算の詳細ログ（推奨）
        
        ### 必要なMOPACキーワード
        
        分子軌道情報を取得するために、MOPAC計算時に以下のキーワードを含めてください：
        
        **基本キーワード:**
        - `VECTORS`: 分子軌道ベクトルを出力
        - `GRAPH MO`: 分子軌道を3Dグラフィック形式で出力
        - `PRECISE`: 高精度計算（推奨）
        
        **詳細解析キーワード（新機能）:**
        - `ALLVEC`: すべての固有ベクトルを出力
        - `GRAPHF`: Jmol用の分子軌道可視化ファイル
        - `HTML`: JSmolによるWeb表示用HTMLファイル
        - `LOCALIZE`: 局在化軌道の計算
        - `MULLIK`: Mulliken電荷解析
        - `BONDS`: 結合次数行列
        - `CHARGES`: 詳細電荷解析
        - `ESP`: 静電ポテンシャル計算
        - `POLAR`: 分極率・超分極率
        - `DISP`: 分散エネルギー・水素結合解析
        - `SUPER`: 親核・親電子非局在化能
        - `HYPERFINE`: 超微細結合定数（ESR用）
        
        ### 計算例
        
        **基本的な分子軌道計算:**
        ```
        PM7 PRECISE VECTORS GRAPH MO
        Example molecule
        
        C   0.0000   0.0000   0.0000
        H   1.0900   0.0000   0.0000
        H  -0.3633   1.0281   0.0000
        H  -0.3633  -0.5140   0.8900
        H  -0.3633  -0.5140  -0.8900
        ```
        
        **詳細解析付きの計算:**
        ```
        PM7 PRECISE VECTORS ALLVEC GRAPHF HTML LOCALIZE MULLIK BONDS CHARGES ESP POLAR
        Advanced analysis example
        
        C   0.0000   0.0000   0.0000
        H   1.0900   0.0000   0.0000
        H  -0.3633   1.0281   0.0000
        H  -0.3633  -0.5140   0.8900
        H  -0.3633  -0.5140  -0.8900
        ```
        
        ### 表示される情報
        
        **基本情報:**
        - 分子の3D構造
        - 軌道エネルギー図
        - HOMO/LUMOエネルギーとギャップ
        - 軌道エネルギー一覧表
        
        **詳細解析情報（新機能）:**
        - **HTML可視化**: JSmolによるWeb表示用ファイル
        - **全分子軌道**: すべての固有ベクトル情報
        - **局在化軌道**: 化学結合の詳細解析
        - **Mulliken電荷**: 各原子の電荷分布
        - **結合次数**: 原子間結合の強さ
        - **静電ポテンシャル**: 分子表面の静電特性
        - **分極率**: 外部電場に対する応答特性
        - **分散・水素結合**: 非共有結合相互作用
        - **反応性指標**: 親核・親電子攻撃部位の予測
        - **超微細結合**: ESR/EPRスペクトル解析
        
        ### 分子軌道計算機能
        
        **可視化特化機能:**
        - 既存の分子構造からMOPAC分子軌道計算を実行
        - GRAPH MOキーワードによる3D分子軌道出力
        - 計算結果とファイル情報の詳細表示
        - 複数の理論手法とオプションをサポート
        
        **可視化オプション:**
        - **HTML可視化**: JSmolによるインタラクティブ表示
        - **全軌道出力**: すべての分子軌道の可視化データ
        - **軌道係数**: 詳細な分子軌道係数情報
        
        **詳細解析は22_mopac_OptandFreqで実行**
        - このページは可視化に特化
        - 詳細な電子構造解析は構造最適化+振動数計算ページで実行
        - Mulliken解析、ESP、分極率などの高度な解析機能
        
        ### dataフォルダからの選択
        
        - MOPACのoutファイルのみを表示します
        - mopac_workディレクトリ内のファイルを対象とします
        - ファイルの詳細情報とプレビューを確認できます
        - 計算履歴を管理しやすくなります
        """)
