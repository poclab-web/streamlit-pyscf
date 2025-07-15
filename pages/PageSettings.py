"""
ページ表示設定の管理ページ

このページでは、サイドバーに表示するページの設定を一括で管理できます。
"""

import streamlit as st
import os
import sys

# プロジェクトルートを追加
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from config.user_preferences import UserPreferences
from utils.module import clean_module_name

st.title("⚙️ ページ表示設定")
st.markdown("サイドバーに表示するページの設定を管理します。")

# ユーザー設定管理
user_prefs = UserPreferences()

# 現在の設定を読み込み
pages_directory = "../pages"  # pagesフォルダからの相対パス
if not os.path.exists(pages_directory):
    pages_directory = "pages"  # ルートからの絶対パス

if os.path.exists(pages_directory):
    page_files = [f for f in os.listdir(pages_directory) 
                 if f.endswith('.py') and f != '__init__.py']
    
    if page_files:
        current_settings = user_prefs.load_page_visibility()
        
        st.header("📋 ページ表示設定")
        st.markdown("チェックボックスでサイドバーに表示するページを選択してください。")
        
        # タブで機能を分割
        tab1, tab2, tab3 = st.tabs(["個別設定", "一括操作", "現在の設定"])
        
        with tab1:
            st.subheader("個別にページを設定")
            
            # 設定変更用のフォーム
            with st.form("page_visibility_form"):
                new_settings = {}
                
                # 2列レイアウトで表示
                col1, col2 = st.columns(2)
                
                for i, page_file in enumerate(sorted(page_files)):
                    clean_name = clean_module_name(page_file)
                    current_value = current_settings.get(page_file, True)
                    
                    # 偶数番目は左列、奇数番目は右列
                    with col1 if i % 2 == 0 else col2:
                        new_settings[page_file] = st.checkbox(
                            clean_name, 
                            value=current_value,
                            key=f"setting_{page_file}",
                            help=f"ファイル名: {page_file}"
                        )
                
                # 保存ボタン
                save_settings = st.form_submit_button("💾 設定を保存", type="primary")
                
                if save_settings:
                    user_prefs.save_page_visibility(new_settings)
                    st.success("✅ 設定を保存しました！")
                    st.balloons()
                    
                    # セッション状態も更新
                    for page_file, visibility in new_settings.items():
                        st.session_state[f"visibility_{page_file}"] = visibility
                    
                    # 少し待ってからページを再読み込み
                    st.rerun()
        
        with tab2:
            st.subheader("一括操作")
            
            col1, col2, col3 = st.columns(3)
            
            with col1:
                if st.button("✅ 全て表示", use_container_width=True):
                    all_visible = {page: True for page in page_files}
                    user_prefs.save_page_visibility(all_visible)
                    st.success("全てのページを表示に設定しました")
                    # セッション状態も更新
                    for page_file in page_files:
                        st.session_state[f"visibility_{page_file}"] = True
                    st.rerun()
            
            with col2:
                if st.button("❌ 全て非表示", use_container_width=True):
                    all_hidden = {page: False for page in page_files}
                    user_prefs.save_page_visibility(all_hidden)
                    st.warning("全てのページを非表示に設定しました")
                    # セッション状態も更新
                    for page_file in page_files:
                        st.session_state[f"visibility_{page_file}"] = False
                    st.rerun()
            
            with col3:
                if st.button("🔄 デフォルトに戻す", use_container_width=True):
                    default_settings = {page: True for page in page_files}
                    user_prefs.save_page_visibility(default_settings)
                    st.info("デフォルト設定に戻しました")
                    # セッション状態も更新
                    for page_file in page_files:
                        st.session_state[f"visibility_{page_file}"] = True
                    st.rerun()
            
            # カテゴリ別設定（推奨ページセット）
            st.markdown("---")
            st.subheader("🎯 推奨設定")
            
            col1, col2 = st.columns(2)
            
            with col1:
                if st.button("⚗️ 基本計算セット", use_container_width=True):
                    basic_calc_pages = [
                        '01_GeneralCalculation.py',
                        '04_Optimization.py', 
                        '02_StructureCheck.py',
                        '06_OPTandFreq.py',
                        'PageSettings.py'
                    ]
                    basic_settings = {page: page in basic_calc_pages for page in page_files}
                    user_prefs.save_page_visibility(basic_settings)
                    st.success("基本計算セットを適用しました")
                    for page_file in page_files:
                        st.session_state[f"visibility_{page_file}"] = basic_settings[page_file]
                    st.rerun()
            
            with col2:
                if st.button("📊 物性計算セット", use_container_width=True):
                    property_calc_pages = [
                        '08_PropertyCalculationIR.py',
                        '08_PropertyCalculationNMR.py',
                        '08_PropertyCalculationPolarizability.py',
                        '12_UV_SpectrumPrediction.py',
                        '07_Visualization.py',
                        'PageSettings.py'
                    ]
                    property_settings = {page: page in property_calc_pages for page in page_files}
                    user_prefs.save_page_visibility(property_settings)
                    st.success("物性計算セットを適用しました")
                    for page_file in page_files:
                        st.session_state[f"visibility_{page_file}"] = property_settings[page_file]
                    st.rerun()
        
        with tab3:
            st.subheader("📈 現在の設定状況")
            
            # 統計情報
            visible_count = sum(1 for v in current_settings.values() if v)
            total_count = len(page_files)
            hidden_count = total_count - visible_count
            
            # メトリクス表示
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("総ページ数", total_count)
            with col2:
                st.metric("表示中", visible_count, delta=f"{visible_count/total_count*100:.1f}%")
            with col3:
                st.metric("非表示", hidden_count, delta=f"{hidden_count/total_count*100:.1f}%")
            
            # 詳細リスト
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("🟢 表示中のページ")
                visible_pages = [page for page, visible in current_settings.items() 
                               if visible and page in page_files]
                if visible_pages:
                    for page in sorted(visible_pages):
                        st.write(f"✅ {clean_module_name(page)}")
                else:
                    st.info("表示中のページはありません")
            
            with col2:
                st.subheader("🔴 非表示のページ")
                hidden_pages = [page for page, visible in current_settings.items() 
                              if not visible and page in page_files]
                if hidden_pages:
                    for page in sorted(hidden_pages):
                        st.write(f"❌ {clean_module_name(page)}")
                else:
                    st.info("非表示のページはありません")
        
        # 設定ファイルの場所を表示
        st.markdown("---")
        st.caption(f"📁 設定ファイル: `{user_prefs.config_file}`")
        
        # エクスポート/インポート機能
        with st.expander("⚙️ 詳細設定"):
            st.subheader("設定のエクスポート/インポート")
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.subheader("📤 エクスポート")
                if st.button("設定をダウンロード"):
                    import json
                    config_json = json.dumps(current_settings, indent=2, ensure_ascii=False)
                    st.download_button(
                        label="💾 JSONファイルをダウンロード",
                        data=config_json,
                        file_name="page_visibility_config.json",
                        mime="application/json"
                    )
            
            with col2:
                st.subheader("📥 インポート")
                uploaded_file = st.file_uploader("設定ファイルをアップロード", type=['json'])
                if uploaded_file is not None:
                    try:
                        import json
                        config_data = json.load(uploaded_file)
                        if st.button("設定を適用"):
                            user_prefs.save_page_visibility(config_data)
                            st.success("設定をインポートしました！")
                            st.rerun()
                    except Exception as e:
                        st.error(f"設定ファイルの読み込みに失敗しました: {e}")
    
    else:
        st.warning("📂 Pagesフォルダ内にPythonファイルが見つかりません。")

else:
    st.error("📁 Pagesフォルダが見つかりません。")

# ヘルプセクション
with st.expander("❓ ヘルプ"):
    st.markdown("""
    ### 使い方
    
    1. **個別設定**: 各ページのチェックボックスで表示/非表示を個別に設定
    2. **一括操作**: 全て表示/非表示やデフォルトに戻すなどの一括操作
    3. **推奨設定**: 用途別にあらかじめ設定されたページセット
    
    ### 推奨設定について
    
    - **基本計算セット**: 一般的な量子化学計算に必要最小限のページ
    - **物性計算セット**: 分子の物性計算に特化したページ
    
    ### 注意事項
    
    - 設定はすぐに保存され、次回起動時にも反映されます
    - 全て非表示にするとサイドバーに何も表示されなくなります
    - 設定ファイルをバックアップしておくことをお勧めします
    """)
