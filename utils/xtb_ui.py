"""
Streamlit用xTB設定UI

xTBのインストール状況確認と基本機能テストのためのStreamlitコンポーネント
"""

import streamlit as st
import subprocess
import os
import platform
import tempfile
import shutil
from pathlib import Path


def check_xtb_installation():
    """
    xTBのインストール状況をチェックする
    
    Returns:
        dict: インストール状況の情報
            - installed (bool): インストールされているかどうか
            - version (str): バージョン情報
            - path (str): 実行ファイルのパス
            - error (str): エラーメッセージ
    """
    try:
        result = subprocess.run(["xtb", "--version"], capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            version_lines = result.stdout.splitlines()
            version_info = ""
            for line in version_lines[:3]:  # 最初の3行を取得
                if line.strip():
                    version_info += line.strip() + "\n"
            
            return {
                'installed': True,
                'version': version_info.strip(),
                'path': "xtb",  # PATHから実行可能
                'error': ""
            }
        else:
            return {
                'installed': False,
                'version': "",
                'path': "",
                'error': f"xTB returned error code: {result.returncode}"
            }
    except FileNotFoundError:
        return {
            'installed': False,
            'version': "",
            'path': "",
            'error': "xTB executable not found in PATH"
        }
    except subprocess.TimeoutExpired:
        return {
            'installed': False,
            'version': "",
            'path': "",
            'error': "xTB version check timed out"
        }
    except Exception as e:
        return {
            'installed': False,
            'version': "",
            'path': "",
            'error': f"Unexpected error: {str(e)}"
        }


def check_xtb_status():
    """xTBのインストール状況と利用可能性をチェック"""
    status = {
        "xtb_available": False,
        "xtb_version": None,
        "xtb_path": None,
        "available_features": {},
        "error_messages": [],
        "diagnostic_info": {}
    }
    
    # xTBコアのチェック
    try:
        xtb_status = check_xtb_installation()
        if xtb_status['installed']:
            status["xtb_available"] = True
            status["xtb_version"] = xtb_status['version']
            status["xtb_path"] = xtb_status['path']
            
            # 追加の診断情報を収集
            try:
                # 実行パスの確認
                which_result = subprocess.run(["which", "xtb"], capture_output=True, text=True, timeout=5)
                if which_result.returncode == 0:
                    status["diagnostic_info"]["executable_path"] = which_result.stdout.strip()
                
                # 基本的なヘルプの確認
                help_result = subprocess.run(["xtb", "--help"], capture_output=True, text=True, timeout=5)
                status["diagnostic_info"]["help_accessible"] = help_result.returncode == 0
                
                # 一時ディレクトリへの書き込み権限確認
                temp_dir = tempfile.gettempdir()
                status["diagnostic_info"]["temp_writable"] = os.access(temp_dir, os.W_OK)
                status["diagnostic_info"]["temp_dir"] = temp_dir
                
                # プラットフォーム情報
                status["diagnostic_info"]["platform"] = platform.platform()
                
            except Exception as e:
                status["diagnostic_info"]["check_error"] = str(e)
                
        else:
            status["error_messages"].append(f"xTB installation error: {xtb_status['error']}")
    except Exception as e:
        status["error_messages"].append(f"xTB check error: {e}")
    
    # xTB機能のチェック
    if status["xtb_available"]:
        features_to_check = {
            "gfn0": "--gfn 0",
            "gfn1": "--gfn 1", 
            "gfn2": "--gfn 2",
            "optimization": "--opt",
            "alpb_solvent": "--alpb",
            "frequencies": "--hess",
            "properties": "--prop"
        }
        
        for feature_name, test_flag in features_to_check.items():
            try:
                # ヘルプに該当フラグが含まれているかチェック
                help_result = subprocess.run(["xtb", "--help"], capture_output=True, text=True, timeout=5)
                if help_result.returncode == 0 and test_flag in help_result.stdout:
                    status["available_features"][feature_name] = True
                else:
                    status["available_features"][feature_name] = False
            except:
                status["available_features"][feature_name] = False
    
    return status


def test_xtb_basic_functionality():
    """xTBの基本機能をテスト（改良版）"""
    try:
        # 簡単なテスト分子（水分子）のXYZ座標
        h2o_xyz = """3
Water molecule
O  0.0000  0.0000  0.0000
H  0.0000  0.0000  0.9600
H  0.9270  0.0000 -0.2400
"""
        
        # 一時ファイルを作成してテスト
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
            f.write(h2o_xyz)
            temp_xyz = f.name
        
        try:
            # より安全なディレクトリ設定
            work_dir = tempfile.mkdtemp(prefix="xtb_basic_")
            
            # xTBシングルポイント計算のテスト（最も安全なオプション）
            result = subprocess.run(
                ["xtb", temp_xyz, "--gfn", "1", "--acc", "2.0"],
                capture_output=True,
                text=True,
                timeout=30,
                cwd=work_dir
            )
            
            if result.returncode == 0:
                return True, "基本機能テスト成功: GFN1レベルでの水分子計算が正常に完了しました"
            else:
                return False, f"基本機能テスト失敗: xTB計算でエラーが発生しました\n{result.stderr}"
                
        finally:
            # ファイルとディレクトリのクリーンアップ
            if os.path.exists(temp_xyz):
                os.unlink(temp_xyz)
            
            # 作業ディレクトリが作成されている場合は削除
            try:
                if os.path.exists(work_dir):
                    shutil.rmtree(work_dir)
            except:
                pass
                
    except subprocess.TimeoutExpired:
        return False, "基本機能テスト失敗: xTB計算がタイムアウトしました（30秒）"
    except Exception as e:
        return False, f"基本機能テスト失敗: {str(e)}"


def test_xtb_optimization():
    """xTB構造最適化機能をテスト"""
    try:
        # メタン分子のXYZ座標（わざと歪んだ構造）
        ch4_xyz = """5
Methane molecule (distorted)
C  0.0000  0.0000  0.0000
H  0.0000  0.0000  1.5000
H  1.5000  0.0000 -0.5000
H -0.7500  1.2990 -0.5000
H -0.7500 -1.2990 -0.5000
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
            f.write(ch4_xyz)
            temp_xyz = f.name
        
        try:
            work_dir = tempfile.mkdtemp(prefix="xtb_opt_")
            
            result = subprocess.run(
                ["xtb", temp_xyz, "--opt", "--gfn", "1"],
                capture_output=True,
                text=True,
                timeout=60,
                cwd=work_dir
            )
            
            if result.returncode == 0:
                return True, "構造最適化テスト成功: メタン分子の最適化が正常に完了しました"
            else:
                return False, f"構造最適化テスト失敗: {result.stderr}"
                
        finally:
            if os.path.exists(temp_xyz):
                os.unlink(temp_xyz)
            try:
                if os.path.exists(work_dir):
                    shutil.rmtree(work_dir)
            except:
                pass
                
    except Exception as e:
        return False, f"構造最適化テスト失敗: {str(e)}"


def test_xtb_solvent_calculation():
    """xTB溶媒効果計算をテスト"""
    try:
        # 水分子の溶媒効果計算
        h2o_xyz = """3
Water molecule
O  0.0000  0.0000  0.0000
H  0.0000  0.0000  0.9600
H  0.9270  0.0000 -0.2400
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xyz', delete=False) as f:
            f.write(h2o_xyz)
            temp_xyz = f.name
        
        try:
            work_dir = tempfile.mkdtemp(prefix="xtb_solv_")
            
            result = subprocess.run(
                ["xtb", temp_xyz, "--alpb", "water", "--gfn", "1"],
                capture_output=True,
                text=True,
                timeout=60,
                cwd=work_dir
            )
            
            if result.returncode == 0:
                return True, "溶媒効果計算テスト成功: ALPB水溶媒モデルでの計算が正常に完了しました"
            else:
                return False, f"溶媒効果計算テスト失敗: {result.stderr}"
                
        finally:
            if os.path.exists(temp_xyz):
                os.unlink(temp_xyz)
            try:
                if os.path.exists(work_dir):
                    shutil.rmtree(work_dir)
            except:
                pass
                
    except Exception as e:
        return False, f"溶媒効果計算テスト失敗: {str(e)}"


def display_xtb_status(show_config_section=True, key_suffix=""):
    """
    xTBのインストール状況を表示し、必要に応じて設定UIを表示する
    
    Args:
        show_config_section (bool): 設定セクションを表示するかどうか
        key_suffix (str): ウィジェットキーに追加するサフィックス
        
    Returns:
        dict: xTBの状況を示す辞書
    """
    st.subheader("🔬 xTB Installation Status")
    status = check_xtb_status()

    # 基本インストール状況の表示
    col1, col2 = st.columns(2)
    
    with col1:
        # xTBコア
        if status["xtb_available"]:
            st.success(f"✅ xTB v{status['xtb_version'].split()[0] if status['xtb_version'] else 'unknown'}")
        else:
            st.error("❌ xTB 利用不可")
    
    with col2:
        # 利用可能機能の概要
        if status["xtb_available"]:
            available_features = sum(1 for feat in status["available_features"].values() if feat)
            total_features = len(status["available_features"])
            st.info(f"🔧 利用可能機能: {available_features}/{total_features}")

    # 機能詳細の表示
    if status["xtb_available"]:
        with st.expander("🔬 xTB機能詳細"):
            col1, col2 = st.columns(2)
            
            feature_items = list(status["available_features"].items())
            mid_point = len(feature_items) // 2
            
            with col1:
                for feature_name, available in feature_items[:mid_point]:
                    display_name = feature_name.replace("_", " ").title()
                    if available:
                        st.success(f"✅ {display_name}")
                    else:
                        st.warning(f"⚠️ {display_name}")
            
            with col2:
                for feature_name, available in feature_items[mid_point:]:
                    display_name = feature_name.replace("_", " ").title()
                    if available:
                        st.success(f"✅ {display_name}")
                    else:
                        st.warning(f"⚠️ {display_name}")
        
        # 診断情報の表示
        if status.get("diagnostic_info"):
            with st.expander("🔍 診断情報"):
                diag_info = status["diagnostic_info"]
                if "executable_path" in diag_info:
                    st.info(f"実行パス: {diag_info['executable_path']}")
                if "platform" in diag_info:
                    st.info(f"プラットフォーム: {diag_info['platform']}")
                if "temp_dir" in diag_info:
                    st.info(f"一時ディレクトリ: {diag_info['temp_dir']}")
                    if diag_info.get("temp_writable"):
                        st.success("✅ 一時ディレクトリへの書き込み権限: OK")
                    else:
                        st.error("❌ 一時ディレクトリへの書き込み権限: NG")

    # エラーメッセージがある場合
    if status["error_messages"]:
        st.subheader("⚠️ エラー詳細")
        for msg in status["error_messages"]:
            st.error(msg)

    # xTBが利用可能な場合のテスト
    if status["xtb_available"] and show_config_section:
        # 動作テストセクション
        st.subheader("🧪 xTB 動作テスト")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            if st.button("基本計算テスト", help="xTBの基本的なシングルポイント計算をテスト", key=f"xtb_basic_test_btn{key_suffix}"):
                with st.spinner("基本計算テスト中..."):
                    test_success, test_message = test_xtb_basic_functionality()
                    if test_success:
                        st.success(test_message)
                    else:
                        st.error(test_message)
        
        with col2:
            if st.button("構造最適化テスト", help="xTBの構造最適化機能をテスト", key=f"xtb_opt_test_btn{key_suffix}"):
                with st.spinner("構造最適化テスト中..."):
                    test_success, test_message = test_xtb_optimization()
                    if test_success:
                        st.success(test_message)
                    else:
                        st.error(test_message)
        
        with col3:
            if st.button("溶媒効果テスト", help="xTBの溶媒効果計算をテスト", key=f"xtb_solv_test_btn{key_suffix}"):
                with st.spinner("溶媒効果テスト中..."):
                    test_success, test_message = test_xtb_solvent_calculation()
                    if test_success:
                        st.success(test_message)
                    else:
                        st.error(test_message)

    elif not status["xtb_available"] and show_config_section:
        # インストール手順の表示
        st.subheader("🔧 インストール手順")
        st.markdown("""
        **xTBをインストールするには:**
        """)
        
        st.code("conda install -c conda-forge xtb", language="bash")
        st.markdown("- xTB: 拡張タイトバインディング量子化学計算プログラム")
        
        st.markdown("""
        **その他のインストール方法:**
        """)
        st.code("""# Homebrewを使用する場合（macOS）
brew install xtb

# pipを使用する場合（非推奨）
pip install xtb-python

# ソースからコンパイル（上級者向け）
git clone https://github.com/grimme-lab/xtb.git
cd xtb
meson setup build --prefix=$HOME/.local
ninja -C build install""", language="bash")

        st.markdown("""
        **推奨インストール手順:**
        1. conda環境でxTBをインストール（最も安定）
        2. ページを再読み込みして状況を確認
        3. 動作テストを実行して正常性を確認
        
        **注意:**
        - conda環境での使用を強く推奨
        - macOSの場合、Xcode Command Line Toolsが必要な場合があります
        - 計算実行時は十分なディスク容量を確保してください
        """)

    return status


def require_xtb():
    """
    xTBが必要なページで使用する関数
    xTBが利用できない場合はページの実行を停止する
    
    Returns:
        dict: xTBが利用可能な場合のステータス（実際には利用不可の場合は停止する）
    """
    status = display_xtb_status(key_suffix="_require")
    
    if not status["xtb_available"]:
        st.stop()
    
    return status


def get_available_gfn_models():
    """利用可能なGFNモデルのリストを取得"""
    gfn_models = [
        {"value": 0, "name": "GFN0-xTB", "description": "最軽量で高速。大きな分子に適用可能"},
        {"value": 1, "name": "GFN1-xTB", "description": "バランスの取れた精度と速度"},
        {"value": 2, "name": "GFN2-xTB", "description": "最高精度。有機分子に最適"}
    ]
    return gfn_models


def get_available_solvents():
    """利用可能な溶媒のリストを取得"""
    solvents = [
        {"value": None, "name": "なし（気相）", "description": "溶媒効果なし"},
        {"value": "water", "name": "水", "description": "ε = 78.39"},
        {"value": "methanol", "name": "メタノール", "description": "ε = 32.66"},
        {"value": "ethanol", "name": "エタノール", "description": "ε = 24.55"},
        {"value": "acetonitrile", "name": "アセトニトリル", "description": "ε = 37.5"},
        {"value": "dmso", "name": "DMSO", "description": "ε = 46.7"},
        {"value": "chloroform", "name": "クロロホルム", "description": "ε = 4.81"},
        {"value": "toluene", "name": "トルエン", "description": "ε = 2.38"}
    ]
    return solvents


def display_gfn_selector(default_gfn=1, key_suffix=""):
    """GFNモデル選択ウィジェットを表示"""
    available_gfn = get_available_gfn_models()
    
    gfn_options = [f"GFN{model['value']} - {model['description']}" for model in available_gfn]
    gfn_values = [model['value'] for model in available_gfn]
    
    selected_index = st.selectbox(
        "GFN Model",
        range(len(gfn_options)),
        format_func=lambda x: gfn_options[x],
        index=default_gfn,
        help="計算に使用するGFNモデルを選択してください",
        key=f"gfn_selector{key_suffix}"
    )
    
    return gfn_values[selected_index]


def display_solvent_selector(default_solvent=None, key_suffix=""):
    """溶媒選択ウィジェットを表示"""
    available_solvents = get_available_solvents()
    
    solvent_options = [f"{solv['name']} - {solv['description']}" for solv in available_solvents]
    solvent_values = [solv['value'] for solv in available_solvents]
    
    default_index = 0
    if default_solvent:
        try:
            default_index = solvent_values.index(default_solvent)
        except ValueError:
            default_index = 0
    
    selected_index = st.selectbox(
        "溶媒",
        range(len(solvent_options)),
        format_func=lambda x: solvent_options[x],
        index=default_index,
        help="ALPB溶媒モデルを使用する溶媒を選択してください",
        key=f"solvent_selector{key_suffix}"
    )
    
    return solvent_values[selected_index]


def display_calculation_options():
    """xTB計算オプションの設定UIを表示"""
    st.subheader("⚙️ 計算設定")
    
    col1, col2 = st.columns(2)
    
    with col1:
        charge = st.number_input(
            "分子電荷",
            value=0,
            step=1,
            help="分子の電荷を指定してください"
        )
        
        calculation_type = st.selectbox(
            "計算タイプ",
            ["Single Point", "Optimization"],
            help="Single Point: エネルギー計算のみ, Optimization: 構造最適化"
        )
    
    with col2:
        uhf = st.number_input(
            "不対電子数 (UHF)",
            value=0,
            min_value=0,
            step=1,
            help="不対電子の数を指定してください"
        )
        
        accuracy = st.selectbox(
            "計算精度",
            ["normal", "crude", "tight"],
            index=0,
            help="計算の精度を選択してください"
        )
    
    return {
        "charge": charge,
        "uhf": uhf,
        "calculation_type": calculation_type,
        "accuracy": accuracy
    }
