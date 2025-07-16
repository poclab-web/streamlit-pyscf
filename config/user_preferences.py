import json
import os
from typing import Dict, List

class UserPreferences:
    def __init__(self, default_config_file="config/page_visibility.json", user_config_file="config/user_page_visibility.json"):
        self.default_config_file = default_config_file
        self.user_config_file = user_config_file
        self.config_dir = os.path.dirname(user_config_file)
        self._ensure_config_dir()
        
    def _ensure_config_dir(self):
        """設定ディレクトリが存在しない場合は作成"""
        if not os.path.exists(self.config_dir):
            os.makedirs(self.config_dir)
    
    def _load_default_config(self) -> Dict[str, bool]:
        """デフォルト設定を読み込み（Git管理対象）"""
        if os.path.exists(self.default_config_file):
            try:
                with open(self.default_config_file, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                
                # 新しい構造化フォーマットの場合（各ページにvisibleとcategoryが含まれる）
                if isinstance(next(iter(data.values()), None), dict):
                    flat_settings = {}
                    for page_name, page_info in data.items():
                        # メタデータキーは除外
                        if page_name.startswith('_'):
                            continue
                        if isinstance(page_info, dict) and "visible" in page_info:
                            flat_settings[page_name] = page_info["visible"]
                        else:
                            # 後方互換性のため、dictだがvisibleキーがない場合はTrueとする
                            flat_settings[page_name] = True
                    return flat_settings
                
                # 古いフラット形式の場合（後方互換性）
                else:
                    return data
                    
            except (json.JSONDecodeError, FileNotFoundError):
                return {}
        return {}
    
    def _load_user_config(self) -> Dict[str, bool]:
        """ユーザー設定を読み込み（Git管理対象外）"""
        if os.path.exists(self.user_config_file):
            try:
                with open(self.user_config_file, 'r', encoding='utf-8') as f:
                    return json.load(f)
            except (json.JSONDecodeError, FileNotFoundError):
                return {}
        return {}
    
    def load_page_visibility(self) -> Dict[str, bool]:
        """ページ表示設定を読み込み（デフォルト設定にユーザー設定をマージ）"""
        default_settings = self._load_default_config()
        user_settings = self._load_user_config()
        
        # デフォルト設定をベースにユーザー設定でオーバーライド
        merged_settings = default_settings.copy()
        merged_settings.update(user_settings)
        
        return merged_settings
    
    def save_page_visibility(self, visibility_settings: Dict[str, bool]):
        """ページ表示設定をユーザー設定ファイルに保存"""
        # デフォルト設定と比較して、変更があったもののみ保存
        default_settings = self._load_default_config()
        user_changes = {}
        
        for page, visible in visibility_settings.items():
            default_visible = default_settings.get(page, True)  # デフォルトは表示
            if visible != default_visible:
                user_changes[page] = visible
        
        # ユーザー設定ファイルに保存
        with open(self.user_config_file, 'w', encoding='utf-8') as f:
            json.dump(user_changes, f, indent=2, ensure_ascii=False)
    
    def update_page_visibility(self, page_name: str, visible: bool):
        """特定ページの表示設定を更新"""
        settings = self.load_page_visibility()
        settings[page_name] = visible
        self.save_page_visibility(settings)
    
    def get_visible_pages(self, all_pages: List[str]) -> List[str]:
        """表示すべきページのリストを取得"""
        visibility_settings = self.load_page_visibility()
        visible_pages = []
        
        for page in all_pages:
            # デフォルトは表示（True）
            if visibility_settings.get(page, True):
                visible_pages.append(page)
        
        return visible_pages
    
    def initialize_page_settings(self, page_files: List[str]):
        """初回起動時に設定ファイルを初期化"""
        # デフォルト設定ファイルが存在しない場合は作成
        if not os.path.exists(self.default_config_file):
            default_settings = {page_file: True for page_file in page_files}
            with open(self.default_config_file, 'w', encoding='utf-8') as f:
                json.dump(default_settings, f, indent=2, ensure_ascii=False)
        
        # 元のJSONデータを読み込み（構造を保持するため）
        original_data = {}
        if os.path.exists(self.default_config_file):
            try:
                with open(self.default_config_file, 'r', encoding='utf-8') as f:
                    original_data = json.load(f)
            except (json.JSONDecodeError, FileNotFoundError):
                original_data = {}
        
        # フラット形式の設定を取得（比較用）
        flat_settings = self._load_default_config()
        
        # 新しいページファイルや削除されたページファイルをチェック
        updated = False
        
        # 新しい構造化フォーマットの場合
        if isinstance(next(iter(original_data.values()), None), dict):
            # 新しいページファイルのデフォルト設定を追加
            for page_file in page_files:
                if page_file not in original_data:
                    # 新しいページには基本計算カテゴリを割り当て
                    original_data[page_file] = {
                        "visible": True,
                        "category": "基本計算"
                    }
                    updated = True
            
            # 存在しないページファイルの設定を削除
            files_to_remove = []
            for page_file in original_data.keys():
                if not page_file.startswith('_') and page_file not in page_files:  # メタデータキーは除外
                    files_to_remove.append(page_file)
                    updated = True
            
            for page_file in files_to_remove:
                del original_data[page_file]
        
        # 古いフラット形式の場合
        else:
            # 新しいページファイルのデフォルト設定を追加
            for page_file in page_files:
                if page_file not in flat_settings:
                    original_data[page_file] = True
                    updated = True
            
            # 存在しないページファイルの設定を削除
            files_to_remove = []
            for page_file in flat_settings.keys():
                if page_file not in page_files:
                    files_to_remove.append(page_file)
                    updated = True
            
            for page_file in files_to_remove:
                if page_file in original_data:
                    del original_data[page_file]
        
        # デフォルト設定ファイルを更新（元の構造を保持）
        if updated:
            with open(self.default_config_file, 'w', encoding='utf-8') as f:
                json.dump(original_data, f, indent=2, ensure_ascii=False)
        
        # ユーザー設定から存在しないページの設定を削除
        user_settings = self._load_user_config()
        user_updated = False
        user_files_to_remove = []
        
        for page_file in user_settings.keys():
            if page_file not in page_files:
                user_files_to_remove.append(page_file)
                user_updated = True
        
        for page_file in user_files_to_remove:
            del user_settings[page_file]
        
        if user_updated:
            with open(self.user_config_file, 'w', encoding='utf-8') as f:
                json.dump(user_settings, f, indent=2, ensure_ascii=False)
    
    def reset_user_settings(self):
        """ユーザー設定をリセット（デフォルト設定に戻す）"""
        if os.path.exists(self.user_config_file):
            os.remove(self.user_config_file)
    
    def get_user_changes(self) -> Dict[str, bool]:
        """デフォルトから変更されたユーザー設定のみを取得"""
        return self._load_user_config()
    
    def has_user_changes(self) -> bool:
        """ユーザー設定に変更があるかどうかを確認"""
        return bool(self._load_user_config())
    
    def get_config_info(self) -> Dict[str, str]:
        """設定ファイルの情報を取得"""
        info = {
            "default_config_exists": os.path.exists(self.default_config_file),
            "user_config_exists": os.path.exists(self.user_config_file),
            "default_config_path": self.default_config_file,
            "user_config_path": self.user_config_file
        }
        
        if info["default_config_exists"]:
            import datetime
            mtime = os.path.getmtime(self.default_config_file)
            info["default_config_modified"] = datetime.datetime.fromtimestamp(mtime).strftime("%Y年%m月%d日 %H:%M:%S")
        
        if info["user_config_exists"]:
            import datetime
            mtime = os.path.getmtime(self.user_config_file)
            info["user_config_modified"] = datetime.datetime.fromtimestamp(mtime).strftime("%Y年%m月%d日 %H:%M:%S")
        
        return info
    
    def get_categories_info(self) -> Dict:
        """カテゴリー情報を取得"""
        if os.path.exists(self.default_config_file):
            try:
                with open(self.default_config_file, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                    return data.get("categories", {})
            except (json.JSONDecodeError, FileNotFoundError):
                return {}
        return {}
    
    def get_metadata(self) -> Dict:
        """メタデータを取得"""
        if os.path.exists(self.default_config_file):
            try:
                with open(self.default_config_file, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                    return data.get("_metadata", {})
            except (json.JSONDecodeError, FileNotFoundError):
                return {}
        return {}
    
    def get_page_category(self, page_file: str) -> str:
        """ページファイルが属するカテゴリーを取得"""
        if os.path.exists(self.default_config_file):
            try:
                with open(self.default_config_file, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                    
                # 新しい構造化フォーマットの場合
                if isinstance(next(iter(data.values()), None), dict):
                    if page_file in data and not page_file.startswith('_') and "category" in data[page_file]:
                        return data[page_file]["category"]
                    
            except (json.JSONDecodeError, FileNotFoundError):
                pass
        return "その他"  # デフォルトカテゴリ
    
    def get_pages_by_category(self) -> Dict[str, List[str]]:
        """カテゴリ別にページをグループ化して取得"""
        category_pages = {}
        
        if os.path.exists(self.default_config_file):
            try:
                with open(self.default_config_file, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                    
                # 新しい構造化フォーマットの場合
                if isinstance(next(iter(data.values()), None), dict):
                    for page_file, page_info in data.items():
                        # メタデータキーは除外
                        if page_file.startswith('_'):
                            continue
                        if isinstance(page_info, dict) and "category" in page_info:
                            category = page_info["category"]
                            if category not in category_pages:
                                category_pages[category] = []
                            category_pages[category].append(page_file)
                    
            except (json.JSONDecodeError, FileNotFoundError):
                pass
        
        return category_pages
