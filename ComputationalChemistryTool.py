import streamlit as st
import os
from config.user_preferences import UserPreferences
from utils.module import clean_module_name, load_readme, get_module_docstrings

# ページ設定
st.set_page_config(
    page_title="Computational Chemistry Tool",
    page_icon="⚗️",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ユーザー設定管理
user_prefs = UserPreferences()

# CSS for better styling
st.markdown("""
<style>
.main-header {
    font-size: 2.5rem;
    font-weight: bold;
    text-align: center;
    margin-bottom: 2rem;
    color: #1f77b4;
}
.category-header {
    font-size: 1.2rem;
    font-weight: bold;
    margin-top: 1rem;
    margin-bottom: 0.5rem;
    color: #666;
}
</style>
""", unsafe_allow_html=True)

# セッション状態の初期化
if "page_visibility_settings" not in st.session_state:
    st.session_state.page_visibility_settings = None

# Streamlit App
st.markdown('<h1 class="main-header">⚗️ Computational Chemistry Tool</h1>', unsafe_allow_html=True)

def home_page():
    """ホームページ"""
    st.markdown("# 🧪 Computational Chemistry Tool")
    st.markdown("### PySCFを活用した量子化学計算プラットフォーム")
    
    # イントロダクション
    with st.container():
        st.markdown("""
        このアプリケーションは、**PySCF**（Python-based Simulations of Chemistry Framework）を使用して、
        様々な量子化学計算を簡単に実行できるWebインターフェースです。
        """)
    
    st.markdown("---")
    
    # 機能概要
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("### 🎯 主な機能")
        st.markdown("""
        - 🧬 **分子構造最適化**
        - ⚡ **一点エネルギー計算**
        - 📊 **IR/NMR/UV-Visスペクトル**
        - 🔍 **Energy Decomposition Analysis**
        - 💧 **溶媒効果計算**
        - 📈 **分子軌道可視化**
        - 🎛️ **カスタマイズ可能な設定**
        """)
    
    with col2:
        st.markdown("### 🚀 特徴")
        st.markdown("""
        - 👥 **ユーザーフレンドリー**な操作
        - 🔧 **柔軟な計算設定**
        - 💾 **計算結果の保存・管理**
        - 📱 **レスポンシブ**デザイン
        - ⚙️ **ページ表示のカスタマイズ**
        - 🎨 **美しいUI/UX**
        - 📊 **詳細な結果表示**
        """)
    
    st.markdown("---")
    
    # 利用可能なページカテゴリ
    st.markdown("### 📚 利用可能な機能カテゴリ")
    
    # ページ統計を取得
    pages_directory = "pages"
    page_files = []
    if os.path.exists(pages_directory):
        page_files = [f for f in os.listdir(pages_directory) if f.endswith('.py') and f != '__init__.py']
    
    current_settings = user_prefs.load_page_visibility()
    visible_pages = [f for f in page_files if current_settings.get(f, True)]
    
    categories = {
        "🧪 基本計算": {
            "keywords": ['general', 'optimization', 'structure', 'singlepoint', 'conformational', 'opt'],
            "description": "分子構造最適化、一点エネルギー計算、配座解析"
        },
        "� 可視化と解析": {
            "keywords": ['visualization', 'energy', 'decomposition', 'fragment', 'analysis'],
            "description": "分子軌道可視化、エネルギー分解解析、フラグメント解析"
        },
        "⚡ 物性計算": {
            "keywords": ['ionization', 'solvation', 'bond', 'pka', 'property'],
            "description": "イオン化ポテンシャル、溶媒効果、結合解離エネルギー、pKa計算"
        },
        "📊 スペクトル計算": {
            "keywords": ['spectrum', 'ir', 'nmr', 'uv', 'polarizability'],
            "description": "IR、NMR、UV-Visスペクトル予測、分極率計算"
        },
        "🔄 遷移状態計算": {
            "keywords": ['transition', 'neb', 'ts', 'irc', 'reaction'],
            "description": "遷移状態探索、反応経路計算、IRC解析"
        },
        "⚙️ システム・設定": {
            "keywords": ['settings', 'database', 'summarization', 'system'],
            "description": "設定管理、データベース、結果集計"
        }
    }
    
    for category, info in categories.items():
        category_files = [f for f in page_files if any(keyword in f.lower() for keyword in info['keywords'])]
        visible_in_category = [f for f in category_files if f in visible_pages]
        
        with st.expander(f"{category} ({len(visible_in_category)}/{len(category_files)} ページ表示中)", expanded=False):
            st.markdown(f"**説明**: {info['description']}")
            
            if visible_in_category:
                st.markdown("**利用可能なページ:**")
                for page_file in sorted(visible_in_category):
                    clean_name = clean_module_name(page_file)
                    st.markdown(f"• {clean_name}")
            else:
                st.info("このカテゴリで表示中のページはありません。")
            
            if len(category_files) > len(visible_in_category):
                hidden_count = len(category_files) - len(visible_in_category)
                st.markdown(f"💡 {hidden_count}個のページが非表示になっています。設定ページで表示/非表示を変更できます。")
    
    st.markdown("---")
    
    # クイックアクション
    st.markdown("### 🚀 クイックアクション")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if st.button("⚙️ ページ設定", use_container_width=True, type="primary"):
            # st.navigationを使用している場合は、直接のページ切り替えはできない
            st.info("左側のナビゲーションから「システム > ページ設定」を選択してください")
    
    with col2:
        if visible_pages:
            if st.button("🧪 計算開始", use_container_width=True, type="secondary"):
                st.info("左側のナビゲーションから計算ページを選択してください")
    
    with col3:
        st.markdown("ℹ️ **ヘルプ**")
        st.markdown("各ページにはガイドが含まれています")
    
    st.markdown("---")
    
    # システム情報
    st.markdown("### 💻 システム情報")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        try:
            import pyscf
            pyscf_version = pyscf.__version__
        except:
            pyscf_version = "未インストール"
        st.info(f"🐍 PySCF: {pyscf_version}")
    
    with col2:
        import streamlit as st_module
        st.info(f"🌊 Streamlit: {st_module.__version__}")
    
    with col3:
        st.info(f"📄 総ページ数: {len(page_files)}")
    
    # フッター
    st.markdown("---")
    st.markdown("""
    <div style='text-align: center; color: #666; font-size: 0.9em;'>
        Made with ❤️ using Streamlit and PySCF<br>
        量子化学計算をもっと身近に
    </div>
    """, unsafe_allow_html=True)

def create_page_from_file(file_path, file_name):
    """ファイルパスからページ関数を作成"""
    def page_function():
        try:
            import importlib.util
            spec = importlib.util.spec_from_file_location(file_name[:-3], file_path)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
        except Exception as e:
            st.error(f"ページの読み込み中にエラーが発生しました: {e}")
    
    return page_function

def get_dynamic_pages():
    """設定に基づいて動的にページを作成"""
    pages_directory = "pages"
    page_dict = {}
    
    if not os.path.exists(pages_directory) or not os.path.isdir(pages_directory):
        return page_dict
    
    # Pythonファイルのリストを取得
    page_files = [f for f in os.listdir(pages_directory) if f.endswith('.py') and f != '__init__.py']
    
    if not page_files:
        return page_dict
    
    # 設定ファイルを初期化
    user_prefs.initialize_page_settings(page_files)
    
    # 設定ファイルから表示するページを取得
    config_settings = user_prefs.load_page_visibility()
    
    # カテゴリ別にページを分類
    calculation_pages = []
    visualization_pages = []
    property_pages = []
    spectrum_pages = []
    transition_pages = []
    system_pages = []
    
    for file_name in sorted(page_files):
        if not config_settings.get(file_name, True):
            continue  # 非表示設定のページはスキップ
        
        # PageSettings.pyは設定ページとして別途処理するのでスキップ
        if file_name == "PageSettings.py":
            continue
        
        file_path = os.path.join(pages_directory, file_name)
        clean_name = clean_module_name(file_name)
        
        # アイコンとカテゴリを決定
        if file_name == "13_ConformationalEnergyDecomposition.py":
            # 特別な処理：ConformationalEnergyDecompositionは可視化と解析に分類
            icon = ":material/analytics:"
            page = st.Page(
                create_page_from_file(file_path, file_name),
                title=clean_name,
                icon=icon,
                url_path=f"vis_{file_name[:-3]}"
            )
            visualization_pages.append(page)
        elif any(keyword in file_name.lower() for keyword in ['general', 'optimization', 'structure', 'singlepoint', 'opt']) or (file_name.lower().startswith('05_conformational') and 'energy' not in file_name.lower()):
            icon = ":material/science:"
            page = st.Page(
                create_page_from_file(file_path, file_name),
                title=clean_name,
                icon=icon,
                url_path=f"calc_{file_name[:-3]}"
            )
            calculation_pages.append(page)
        elif any(keyword in file_name.lower() for keyword in ['ionization', 'solvation', 'bond', 'pka']) and 'property' not in file_name.lower():
            icon = ":material/bolt:"
            page = st.Page(
                create_page_from_file(file_path, file_name),
                title=clean_name,
                icon=icon,
                url_path=f"prop_{file_name[:-3]}"
            )
            property_pages.append(page)
        elif any(keyword in file_name.lower() for keyword in ['visualization', 'energydecomposition', 'fragment']):
            icon = ":material/analytics:"
            page = st.Page(
                create_page_from_file(file_path, file_name),
                title=clean_name,
                icon=icon,
                url_path=f"vis_{file_name[:-3]}"
            )
            visualization_pages.append(page)
        elif any(keyword in file_name.lower() for keyword in ['neb', 'ts', 'irc', 'transition']):
            icon = ":material/timeline:"
            page = st.Page(
                create_page_from_file(file_path, file_name),
                title=clean_name,
                icon=icon,
                url_path=f"trans_{file_name[:-3]}"
            )
            transition_pages.append(page)
        elif any(keyword in file_name.lower() for keyword in ['propertycalculationir', 'propertycalculationnmr', 'uv_spectrum', 'polarizability', 'spectrum']) or file_name.startswith('08_PropertyCalculation'):
            icon = ":material/graphic_eq:"
            page = st.Page(
                create_page_from_file(file_path, file_name),
                title=clean_name,
                icon=icon,
                url_path=f"spec_{file_name[:-3]}"
            )
            spectrum_pages.append(page)
        else:
            icon = ":material/settings:"
            page = st.Page(
                create_page_from_file(file_path, file_name),
                title=clean_name,
                icon=icon,
                url_path=f"sys_{file_name[:-3]}"
            )
            system_pages.append(page)
    
    # カテゴリ別にページ辞書を構築
    if calculation_pages:
        page_dict["基本計算"] = calculation_pages
    if visualization_pages:
        page_dict["可視化と解析"] = visualization_pages
    if property_pages:
        page_dict["物性計算"] = property_pages
    if spectrum_pages:
        page_dict["スペクトル計算"] = spectrum_pages
    if transition_pages:
        page_dict["遷移状態計算"] = transition_pages
    if system_pages:
        page_dict["システム・設定"] = system_pages
    
    return page_dict

def settings_page():
    """設定ページ"""
    st.markdown("## ⚙️ ページ表示設定")
    st.markdown("ナビゲーションに表示するページのカテゴリと個別ページを管理します。")
    
    pages_directory = "pages"
    if not os.path.exists(pages_directory):
        st.error("Pagesフォルダが見つかりません。")
        return
    
    page_files = [f for f in os.listdir(pages_directory) if f.endswith('.py') and f != '__init__.py']
    
    if not page_files:
        st.warning("Pagesフォルダ内にPythonファイルが見つかりません。")
        return
    
    current_settings = user_prefs.load_page_visibility()
    
    # 統計情報を表示
    visible_count = sum(1 for v in current_settings.values() if v)
    total_count = len(page_files)
    
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("📋 総ページ数", total_count)
    with col2:
        st.metric("✅ 表示中", visible_count, delta=f"{visible_count/total_count*100:.0f}%")
    with col3:
        st.metric("❌ 非表示", total_count - visible_count, delta=f"{(total_count-visible_count)/total_count*100:.0f}%")
    
    st.markdown("---")
    
    # タブで機能を分割
    tab1, tab2, tab3 = st.tabs(["🎯 個別設定", "🔧 一括操作", "📊 設定状況"])
    
    with tab1:
        st.markdown("### カテゴリ別ページ設定")
        
        with st.form("page_visibility_form"):
            new_settings = {}
            processed_files = set()  # 処理済みファイルを追跡
            
            # カテゴリ別に分けて表示
            categories = {
                "🧪 基本計算": ['general', 'optimization', 'structure', 'singlepoint', 'conformational', 'opt'],
                "� 可視化と解析": ['visualization', 'energy', 'decomposition', 'fragment', 'analysis'],
                "⚡ 物性計算": ['ionization', 'solvation', 'bond', 'pka'],
                "📊 スペクトル計算": ['spectrum', 'ir', 'nmr', 'uv', 'polarizability'],
                "🔄 遷移状態計算": ['transition', 'neb', 'ts', 'irc', 'reaction'],
                "⚙️ システム・設定": ['settings', 'database', 'summarization', 'system']
            }
            
            for category, keywords in categories.items():
                st.markdown(f"#### {category}")
                category_files = []
                for page_file in sorted(page_files):
                    if page_file in processed_files:
                        continue  # 既に処理済みのファイルはスキップ
                    if any(keyword in page_file.lower() for keyword in keywords):
                        category_files.append(page_file)
                        processed_files.add(page_file)  # 処理済みとしてマーク
                
                if category_files:
                    cols = st.columns(min(3, len(category_files)))
                    for i, page_file in enumerate(category_files):
                        clean_name = clean_module_name(page_file)
                        current_value = current_settings.get(page_file, True)
                        with cols[i % len(cols)]:
                            new_settings[page_file] = st.checkbox(
                                clean_name,
                                value=current_value,
                                key=f"setting_{category}_{i}_{page_file}",  # カテゴリとインデックスを含む
                                help=f"ファイル: {page_file}"
                            )
                else:
                    st.info(f"このカテゴリにはページがありません。")
                
                st.markdown("")
            
            # その他のファイル
            other_files = []
            for page_file in sorted(page_files):
                if page_file not in processed_files:
                    other_files.append(page_file)
            
            if other_files:
                st.markdown("#### 📂 その他")
                cols = st.columns(min(3, len(other_files)))
                for i, page_file in enumerate(other_files):
                    clean_name = clean_module_name(page_file)
                    current_value = current_settings.get(page_file, True)
                    with cols[i % len(cols)]:
                        new_settings[page_file] = st.checkbox(
                            clean_name,
                            value=current_value,
                            key=f"setting_other_{i}_{page_file}",  # その他カテゴリ用のキー
                            help=f"ファイル: {page_file}"
                        )
            
            # 保存ボタン
            col1, col2, col3 = st.columns([1, 1, 1])
            with col2:
                save_settings = st.form_submit_button(
                    "💾 設定を保存", 
                    type="primary", 
                    use_container_width=True
                )
            
            if save_settings:
                user_prefs.save_page_visibility(new_settings)
                st.success("✅ 設定を保存しました！")
                st.balloons()
                st.rerun()
    
    with tab2:
        st.markdown("### 一括操作")
        st.markdown("複数のページの表示設定を一度に変更できます。")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            if st.button("✅ 全て表示", use_container_width=True, type="primary"):
                all_visible = {page: True for page in page_files}
                user_prefs.save_page_visibility(all_visible)
                st.success("全てのページを表示に設定しました")
                st.rerun()
        
        with col2:
            if st.button("❌ 全て非表示", use_container_width=True):
                all_hidden = {page: False for page in page_files}
                user_prefs.save_page_visibility(all_hidden)
                st.warning("全てのページを非表示に設定しました")
                st.rerun()
        
        with col3:
            if st.button("🔄 デフォルト", use_container_width=True):
                user_prefs.reset_user_settings()
                st.info("デフォルト設定に戻しました")
                st.rerun()
        
        st.markdown("---")
        st.markdown("### 🎯 プリセット設定")
        
        col1, col2 = st.columns(2)
        
        with col1:
            if st.button("🧪 基本計算セット", use_container_width=True):
                basic_calc_pages = [f for f in page_files if any(keyword in f.lower() for keyword in ['general', 'optimization', 'structure', 'opt']) or f == 'PageSettings.py']
                basic_settings = {page: page in basic_calc_pages for page in page_files}
                user_prefs.save_page_visibility(basic_settings)
                st.success("基本計算セットを適用しました")
                st.rerun()
        
        with col2:
            if st.button("📊 スペクトル計算セット", use_container_width=True):
                spectrum_calc_pages = [f for f in page_files if any(keyword in f.lower() for keyword in ['spectrum', 'ir', 'nmr', 'uv', 'polarizability']) or f == 'PageSettings.py']
                spectrum_settings = {page: page in spectrum_calc_pages for page in page_files}
                user_prefs.save_page_visibility(spectrum_settings)
                st.success("スペクトル計算セットを適用しました")
                st.rerun()
    
    with tab3:
        st.markdown("### 📊 現在の設定詳細")
        
        # カテゴリ別の設定状況
        categories = {
            "🧪 基本計算": ['general', 'optimization', 'structure', 'singlepoint', 'conformational', 'opt'],
            "� 可視化と解析": ['visualization', 'energy', 'decomposition', 'fragment', 'analysis'],
            "⚡ 物性計算": ['ionization', 'solvation', 'bond', 'pka'],
            "📊 スペクトル計算": ['spectrum', 'ir', 'nmr', 'uv', 'polarizability'],
            "🔄 遷移状態計算": ['transition', 'neb', 'ts', 'irc', 'reaction'],
            "⚙️ システム・設定": ['settings', 'database', 'summarization', 'system']
        }
        
        for category, keywords in categories.items():
            st.markdown(f"#### {category}")
            category_files = [f for f in page_files if any(keyword in f.lower() for keyword in keywords)]
            
            if category_files:
                visible_in_category = [f for f in category_files if current_settings.get(f, True)]
                hidden_in_category = [f for f in category_files if not current_settings.get(f, True)]
                
                col1, col2 = st.columns(2)
                
                with col1:
                    st.markdown("**✅ 表示中:**")
                    if visible_in_category:
                        for f in visible_in_category:
                            st.write(f"• {clean_module_name(f)}")
                    else:
                        st.info("表示中のページはありません")
                
                with col2:
                    st.markdown("**❌ 非表示:**")
                    if hidden_in_category:
                        for f in hidden_in_category:
                            st.write(f"• {clean_module_name(f)}")
                    else:
                        st.info("非表示のページはありません")
            else:
                st.info("このカテゴリにはページがありません。")
            
            st.markdown("")
        
        # 設定ファイル情報
        st.markdown("---")
        st.markdown("### 📁 設定ファイル情報")
        
        config_info = user_prefs.get_config_info()
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("**🔧 デフォルト設定（Git管理）**")
            if config_info["default_config_exists"]:
                st.info(f"📍 パス: `{config_info['default_config_path']}`")
                if "default_config_modified" in config_info:
                    st.info(f"🕒 更新: {config_info['default_config_modified']}")
            else:
                st.warning("デフォルト設定ファイルが見つかりません")
        
        with col2:
            st.markdown("**👤 ユーザー設定（Git管理外）**")
            if config_info["user_config_exists"]:
                st.info(f"📍 パス: `{config_info['user_config_path']}`")
                if "user_config_modified" in config_info:
                    st.info(f"🕒 更新: {config_info['user_config_modified']}")
                
                # ユーザー設定のリセットボタン
                if st.button("🔄 ユーザー設定をリセット", key="reset_user_settings"):
                    user_prefs.reset_user_settings()
                    st.success("ユーザー設定をリセットしました")
                    st.rerun()
            else:
                st.info("ユーザー設定はデフォルトを使用中")
        
        # 設定の詳細情報
        st.markdown("---")
        st.markdown("### � 設定詳細")
        
        user_changes = user_prefs.get_user_changes()
        if user_changes:
            st.markdown(f"**デフォルトから変更された設定: {len(user_changes)}件**")
            
            with st.expander("変更された設定を表示", expanded=False):
                for page, visible in user_changes.items():
                    status = "✅ 表示" if visible else "❌ 非表示"
                    clean_name = clean_module_name(page)
                    st.write(f"• {clean_name}: {status}")
        else:
            st.info("全ての設定がデフォルト値を使用しています")

# メインナビゲーション
def main():
    # ホームページを作成
    home = st.Page(home_page, title="ホーム", icon=":material/home:", url_path="home")
    
    # 設定ページを作成
    settings = st.Page(settings_page, title="ページ設定", icon=":material/settings:", url_path="page_settings")
    
    # 動的ページを取得
    dynamic_pages = get_dynamic_pages()
    
    # ナビゲーション辞書を構築
    page_dict = {"メイン": [home]}
    
    # 動的ページを追加
    page_dict.update(dynamic_pages)
    
    # 設定ページを追加
    page_dict["システム"] = [settings]
    
    # ナビゲーションを作成
    pg = st.navigation(page_dict)
    
    # ページを実行
    pg.run()

if __name__ == "__main__":
    main()
