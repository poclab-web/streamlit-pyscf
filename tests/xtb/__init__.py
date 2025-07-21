"""
xTB テストパッケージ

このパッケージには、xTB (extended tight-binding) 計算の動作確認用テストが含まれています。

パッケージ構成:
- minimal_xtb_test.py - 最小限のxTB動作確認テスト
- safe_xtb_test.py - 異なるGFNレベルでの安全性テスト  
- final_xtb_test.py - 包括的なベンチマークテスト (推奨)
- test_xtb_simple.py - プロジェクト統合テスト
- production_xtb_test.py - 本格的統合テスト
"""

__version__ = "1.0.0"

# パッケージからインポート可能な関数
__all__ = [
    "run_all_tests",
    "run_basic_test",
    "run_minimal_test",
    "run_final_test"
]

def run_minimal_test():
    """最小限のxTBテストを実行"""
    import subprocess
    import sys
    from pathlib import Path
    
    test_file = Path(__file__).parent / "minimal_xtb_test.py"
    result = subprocess.run([sys.executable, str(test_file)], capture_output=True, text=True)
    return result.returncode == 0, result.stdout, result.stderr

def run_final_test():
    """包括的なxTBテストを実行"""
    import subprocess
    import sys
    from pathlib import Path
    
    test_file = Path(__file__).parent / "final_xtb_test.py"
    result = subprocess.run([sys.executable, str(test_file)], capture_output=True, text=True)
    return result.returncode == 0, result.stdout, result.stderr

def run_basic_test():
    """基本的なxTBテストを実行 (minimal_testのエイリアス)"""
    return run_minimal_test()

def run_all_tests():
    """すべてのxTBテストを実行"""
    import subprocess
    import sys
    from pathlib import Path
    
    test_files = [
        "minimal_xtb_test.py",
        "safe_xtb_test.py", 
        "final_xtb_test.py"
    ]
    
    results = {}
    
    for test_file in test_files:
        test_path = Path(__file__).parent / test_file
        if test_path.exists():
            result = subprocess.run([sys.executable, str(test_path)], capture_output=True, text=True)
            results[test_file] = {
                "success": result.returncode == 0,
                "stdout": result.stdout,
                "stderr": result.stderr
            }
        else:
            results[test_file] = {
                "success": False,
                "stdout": "",
                "stderr": f"Test file not found: {test_path}"
            }
    
    return results
