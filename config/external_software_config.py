"""
外部ソフトウェア設定管理モジュール

外部ソフトウェア（MOPAC、Gaussian、ORCAなど）のパス設定を保存・読み込みする機能を提供します。
設定はexternal_software_paths.jsonファイルに保存されます（このファイルはgitで管理されません）。
"""

import json
import os
from pathlib import Path
import shutil

# 設定ファイルのパス
CONFIG_FILE = Path("config/external_software_paths.json")

def save_software_path(software_name, software_path):
    """
    外部ソフトウェアのパスを設定ファイルに保存する
    
    Args:
        software_name (str): ソフトウェア名（例: 'mopac', 'gaussian', 'orca'）
        software_path (str): ソフトウェアバイナリのフルパス
    """
    try:
        # 設定ディレクトリを作成（存在しない場合）
        CONFIG_FILE.parent.mkdir(exist_ok=True)
        
        # 既存の設定を読み込み（存在する場合）
        config = {}
        if CONFIG_FILE.exists():
            try:
                with open(CONFIG_FILE, 'r', encoding='utf-8') as f:
                    config = json.load(f)
            except (json.JSONDecodeError, UnicodeDecodeError):
                config = {}
        
        # ソフトウェアパスを更新
        config[software_name] = str(software_path)
        
        # 設定を保存
        with open(CONFIG_FILE, 'w', encoding='utf-8') as f:
            json.dump(config, f, indent=2, ensure_ascii=False)
        
        return True
    except Exception as e:
        print(f"Error saving {software_name} path: {e}")
        return False

def load_software_path(software_name):
    """
    設定ファイルから指定したソフトウェアのパスを読み込む
    
    Args:
        software_name (str): ソフトウェア名
        
    Returns:
        str or None: ソフトウェアパス（設定されていない場合はNone）
    """
    try:
        if not CONFIG_FILE.exists():
            return None
        
        with open(CONFIG_FILE, 'r', encoding='utf-8') as f:
            config = json.load(f)
        
        return config.get(software_name)
    except Exception as e:
        print(f"Error loading {software_name} path: {e}")
        return None

def validate_software_path(software_path):
    """
    ソフトウェアパスが有効かどうかを検証する
    
    Args:
        software_path (str): 検証するソフトウェアパス
        
    Returns:
        dict: 検証結果
            - valid (bool): パスが有効かどうか
            - error (str): エラーメッセージ（エラーがない場合は空文字列）
            - absolute_path (str): 絶対パス
    """
    result = {
        'valid': False,
        'error': '',
        'absolute_path': ''
    }
    
    try:
        if not software_path:
            result['error'] = 'Path is empty'
            return result
        
        # パスをPathオブジェクトに変換
        path_obj = Path(software_path)
        
        # ファイルが存在するかチェック
        if not path_obj.exists():
            result['error'] = f'File does not exist: {software_path}'
            return result
        
        # ファイルかどうかチェック
        if not path_obj.is_file():
            result['error'] = f'Path is not a file: {software_path}'
            return result
        
        # 実行可能かチェック
        if not os.access(path_obj, os.X_OK):
            result['error'] = f'File is not executable: {software_path}'
            return result
        
        # 絶対パスに変換
        absolute_path = str(path_obj.absolute())
        result['valid'] = True
        result['absolute_path'] = absolute_path
        
        return result
        
    except Exception as e:
        result['error'] = f'Error validating path: {e}'
        return result

def get_software_executable(software_name, possible_names=None):
    """
    指定したソフトウェアの実行可能ファイルのパスを取得する
    
    優先順位:
    1. 設定ファイルに保存されたパス
    2. PATHから自動検索
    
    Args:
        software_name (str): ソフトウェア名
        possible_names (list): 検索する実行ファイル名のリスト（省略時は software_name を使用）
        
    Returns:
        dict: 実行可能ファイル情報
            - path (str): 実行可能ファイルのパス（見つからない場合はNone）
            - source (str): パスの取得元（'config', 'path', 'not_found'）
            - error (str): エラーメッセージ
    """
    result = {
        'path': None,
        'source': 'not_found',
        'error': ''
    }
    
    # 1. 設定ファイルから読み込み
    saved_path = load_software_path(software_name)
    if saved_path:
        validation = validate_software_path(saved_path)
        if validation['valid']:
            result['path'] = validation['absolute_path']
            result['source'] = 'config'
            return result
        else:
            result['error'] = f"Saved path is invalid: {validation['error']}"
    
    # 2. PATHから自動検索
    if possible_names is None:
        possible_names = [software_name]
    
    for name in possible_names:
        path = shutil.which(name)
        if path:
            result['path'] = path
            result['source'] = 'path'
            return result
    
    # 見つからない場合
    if not result['error']:
        result['error'] = f'{software_name} executable not found in PATH or config'
    
    return result

def clear_software_config(software_name=None):
    """
    ソフトウェア設定をクリアする
    
    Args:
        software_name (str): 削除するソフトウェア名（省略時は全設定を削除）
        
    Returns:
        bool: 成功した場合True
    """
    try:
        if not CONFIG_FILE.exists():
            return True
            
        if software_name is None:
            # 全設定を削除
            CONFIG_FILE.unlink()
        else:
            # 特定のソフトウェア設定のみ削除
            with open(CONFIG_FILE, 'r', encoding='utf-8') as f:
                config = json.load(f)
            
            if software_name in config:
                del config[software_name]
                
                # 設定が空になった場合はファイルを削除
                if not config:
                    CONFIG_FILE.unlink()
                else:
                    with open(CONFIG_FILE, 'w', encoding='utf-8') as f:
                        json.dump(config, f, indent=2, ensure_ascii=False)
        
        return True
    except Exception as e:
        print(f"Error clearing software config: {e}")
        return False

def list_configured_software():
    """
    設定済みのソフトウェア一覧を取得する
    
    Returns:
        dict: 設定済みソフトウェアの辞書
    """
    try:
        if not CONFIG_FILE.exists():
            return {}
        
        with open(CONFIG_FILE, 'r', encoding='utf-8') as f:
            config = json.load(f)
        
        return config
    except Exception as e:
        print(f"Error listing configured software: {e}")
        return {}

# MOPAC用の便利関数（後方互換性のため）
def save_mopac_path(mopac_path):
    """MOPACのパスを設定ファイルに保存する（後方互換性用）"""
    return save_software_path('mopac', mopac_path)

def load_mopac_path():
    """設定ファイルからMOPACのパスを読み込む（後方互換性用）"""
    return load_software_path('mopac')

def validate_mopac_path(mopac_path):
    """MOPACパスが有効かどうかを検証する（後方互換性用）"""
    return validate_software_path(mopac_path)

def get_mopac_executable():
    """MOPACの実行可能ファイルのパスを取得する（後方互換性用）"""
    mopac_names = ['mopac', 'MOPAC2016.exe', 'MOPAC2023.exe', 'MOPAC']
    return get_software_executable('mopac', mopac_names)

def clear_mopac_config():
    """MOPAC設定をクリアする（後方互換性用）"""
    return clear_software_config('mopac')
