"""
Streamlit用MOPAC設定UI

MOPACのインストール状況確認とパス設定のためのStreamlitコンポーネント
"""

import streamlit as st
from logic.mopac_calculation import check_mopac_installation
from config.external_software_config import save_mopac_path, validate_mopac_path, clear_mopac_config

def display_mopac_status(show_config_section=True):
    """
    MOPACのインストール状況を表示し、必要に応じて設定UIを表示する
    
    Args:
        show_config_section (bool): 設定セクションを表示するかどうか
        
    Returns:
        bool: MOPACが利用可能な場合True、そうでなければFalse
    """
    st.subheader("MOPAC Installation Status")
    mopac_status = check_mopac_installation()

    if mopac_status["installed"]:
        # インストール済みの場合
        source_text = {
            'config': '(設定ファイルから)',
            'path': '(PATHから自動検出)',
            'not_found': ''
        }.get(mopac_status.get('source', 'not_found'), '')
        
        st.success(f"✅ MOPAC is installed at: `{mopac_status['path']}` {source_text}")
            
        # エラーがある場合は警告として表示
        if mopac_status["error"]:
            st.warning(f"Note: {mopac_status['error']}")
        
        # 設定をクリアするオプション
        if show_config_section and mopac_status.get('source') == 'config':
            col1, col2 = st.columns([3, 1])
            with col2:
                if st.button("設定をクリア", help="保存されたMOPACパスをクリアして自動検出に戻す"):
                    if clear_mopac_config():
                        st.success("設定をクリアしました")
                        st.rerun()
                    else:
                        st.error("設定のクリアに失敗しました")
        
        return True

    else:
        # インストールされていない場合
        st.error("❌ MOPAC is not installed or not found in PATH")
        st.write(f"Error: {mopac_status['error']}")
        
        if show_config_section:
            # パス入力セクション
            st.subheader("🔧 MOPAC Path Configuration")
            st.markdown("""
            MOPACが見つからない場合、手動でパスを指定できます。
            MOPACバイナリの完全なファイルパスを入力してください。
            """)
            
            # パス入力フォーム
            with st.form("mopac_path_form"):
                mopac_path_input = st.text_input(
                    "MOPAC Binary Path",
                    placeholder="/home/username/tools/mopac2023/bin/MOPAC",
                    help="MOPACバイナリの完全なファイルパス"
                )
                
                col1, col2 = st.columns(2)
                with col1:
                    submit_button = st.form_submit_button("パスを保存", type="primary")
                with col2:
                    test_button = st.form_submit_button("パスをテスト")
                
                if test_button and mopac_path_input:
                    # パスの検証のみ実行
                    validation = validate_mopac_path(mopac_path_input)
                    if validation['valid']:
                        st.success(f"✅ Valid MOPAC binary: `{validation['absolute_path']}`")
                    else:
                        st.error(f"❌ Invalid path: {validation['error']}")
                
                if submit_button and mopac_path_input:
                    # パスの検証と保存
                    validation = validate_mopac_path(mopac_path_input)
                    if validation['valid']:
                        if save_mopac_path(validation['absolute_path']):
                            st.success(f"✅ MOPAC path saved: `{validation['absolute_path']}`")
                            st.info("ページを再読み込みして設定を反映してください")
                            if st.button("ページを再読み込み"):
                                st.rerun()
                        else:
                            st.error("Failed to save MOPAC path")
                    else:
                        st.error(f"❌ Invalid path: {validation['error']}")
            
            # インストール手順
            st.markdown("""
            **MOPACをインストールするには:**
            1. [MOPAC公式サイト](http://openmopac.net/)からダウンロード
            2. インストール後、上記のフォームでバイナリのパスを指定
            3. または`mopac`コマンドがターミナルで実行できるように環境変数PATHを設定
            """)
        
        return False

def require_mopac():
    """
    MOPACが必要なページで使用する関数
    MOPACが利用できない場合はページの実行を停止する
    
    Returns:
        bool: MOPACが利用可能な場合True（実際にはFalseの場合は停止する）
    """
    if not display_mopac_status():
        st.stop()
    return True
