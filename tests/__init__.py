"""
streamlit-pyscf テストスイート

このパッケージには、streamlit-pyscfプロジェクトの各種テストが含まれています。

テストモジュール:
- xtb: xTB計算のテスト
"""

__version__ = "1.0.0"

# テストモジュールの便利関数をエクスポート
def run_xtb_tests():
    """xTBテストを実行"""
    from .xtb import run_all_tests
    return run_all_tests()

def run_basic_xtb_test():
    """基本的なxTBテストを実行"""
    from .xtb import run_basic_test
    return run_basic_test()
