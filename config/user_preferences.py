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
                
                # 新しい構造化フォーマットの場合
                if "categories" in data:
                    flat_settings = {}
                    for category_info in data["categories"].values():
                        if "pages" in category_info:
                            flat_settings.update(category_info["pages"])
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
        
        # デフォルト設定を読み込み
        default_settings = self._load_default_config()
        
        # 新しいページファイルのデフォルト設定を追加
        updated = False
        for page_file in page_files:
            if page_file not in default_settings:
                default_settings[page_file] = True
                updated = True
        
        # 存在しないページファイルの設定を削除
        files_to_remove = []
        for page_file in default_settings.keys():
            if page_file not in page_files:
                files_to_remove.append(page_file)
                updated = True
        
        for page_file in files_to_remove:
            del default_settings[page_file]
        
        # デフォルト設定ファイルを更新
        if updated:
            with open(self.default_config_file, 'w', encoding='utf-8') as f:
                json.dump(default_settings, f, indent=2, ensure_ascii=False)
        
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
                    return data.get("metadata", {})
            except (json.JSONDecodeError, FileNotFoundError):
                return {}
        return {}
    
    def get_category_for_page(self, page_file: str) -> str:
        """ページファイルが属するカテゴリーを取得"""
        categories = self.get_categories_info()
        for category_key, category_info in categories.items():
            if "pages" in category_info and page_file in category_info["pages"]:
                return category_key
        return "other"
    
    def get_pages_by_category(self, category_key: str) -> Dict[str, bool]:
        """指定したカテゴリーのページ設定を取得"""
        categories = self.get_categories_info()
        if category_key in categories and "pages" in categories[category_key]:
            return categories[category_key]["pages"]
        return {}
