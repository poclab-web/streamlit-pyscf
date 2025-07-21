#!/usr/bin/env python3
"""
xTB テストランナー

すべてのxTBテストを実行し、結果をまとめて表示します。

使用方法:
    python run_tests.py [オプション]

オプション:
    --minimal    最小限のテストのみ実行
    --final      包括的なテストのみ実行  
    --all        すべてのテストを実行 (デフォルト)
    --verbose    詳細な出力を表示
"""

import sys
import argparse
from pathlib import Path
import subprocess
from datetime import datetime

def run_single_test(test_file, verbose=False):
    """単一のテストファイルを実行"""
    test_path = Path(__file__).parent / test_file
    
    if not test_path.exists():
        return {
            "success": False,
            "stdout": "",
            "stderr": f"Test file not found: {test_path}",
            "runtime": 0
        }
    
    print(f"実行中: {test_file}...")
    start_time = datetime.now()
    
    try:
        result = subprocess.run(
            [sys.executable, str(test_path)],
            capture_output=True,
            text=True,
            timeout=300  # 5分でタイムアウト
        )
        
        end_time = datetime.now()
        runtime = (end_time - start_time).total_seconds()
        
        success = result.returncode == 0
        status = "✓ 成功" if success else "✗ 失敗"
        print(f"  {status} ({runtime:.1f}秒)")
        
        if verbose:
            if result.stdout:
                print(f"  出力: {result.stdout[:200]}...")
            if result.stderr and not success:
                print(f"  エラー: {result.stderr[:200]}...")
        
        return {
            "success": success,
            "stdout": result.stdout,
            "stderr": result.stderr,
            "runtime": runtime
        }
        
    except subprocess.TimeoutExpired:
        print(f"  ✗ タイムアウト (300秒)")
        return {
            "success": False,
            "stdout": "",
            "stderr": "Test timed out after 300 seconds",
            "runtime": 300
        }
    except Exception as e:
        print(f"  ✗ 実行エラー: {e}")
        return {
            "success": False,
            "stdout": "",
            "stderr": str(e),
            "runtime": 0
        }

def main():
    """メイン実行関数"""
    parser = argparse.ArgumentParser(description="xTB テストランナー")
    parser.add_argument("--minimal", action="store_true", help="最小限のテストのみ実行")
    parser.add_argument("--final", action="store_true", help="包括的なテストのみ実行")
    parser.add_argument("--all", action="store_true", help="すべてのテストを実行")
    parser.add_argument("--verbose", "-v", action="store_true", help="詳細な出力を表示")
    
    args = parser.parse_args()
    
    # デフォルトは--all
    if not (args.minimal or args.final):
        args.all = True
    
    print("=" * 60)
    print("xTB テストランナー")
    print("=" * 60)
    print(f"開始時刻: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # 実行するテストを決定
    test_files = []
    
    if args.minimal:
        test_files = ["minimal_xtb_test.py"]
    elif args.final:
        test_files = ["final_xtb_test.py"]
    elif args.all:
        test_files = [
            "minimal_xtb_test.py",
            "safe_xtb_test.py",
            "final_xtb_test.py"
        ]
    
    print(f"\n実行予定テスト: {len(test_files)}個")
    for test_file in test_files:
        print(f"  - {test_file}")
    
    # テスト実行
    print(f"\n{'-' * 40}")
    print("テスト実行")
    print(f"{'-' * 40}")
    
    results = {}
    total_start = datetime.now()
    
    for test_file in test_files:
        result = run_single_test(test_file, args.verbose)
        results[test_file] = result
    
    total_end = datetime.now()
    total_runtime = (total_end - total_start).total_seconds()
    
    # 結果サマリー
    print(f"\n{'-' * 40}")
    print("テスト結果サマリー")
    print(f"{'-' * 40}")
    
    successful_tests = sum(1 for result in results.values() if result["success"])
    total_tests = len(results)
    
    print(f"成功: {successful_tests}/{total_tests}")
    print(f"成功率: {successful_tests/total_tests*100:.1f}%")
    print(f"総実行時間: {total_runtime:.1f}秒")
    
    # 詳細結果
    if args.verbose or successful_tests < total_tests:
        print(f"\n詳細結果:")
        for test_file, result in results.items():
            status = "✓" if result["success"] else "✗"
            print(f"  {status} {test_file} ({result['runtime']:.1f}秒)")
            
            if not result["success"] and result["stderr"]:
                print(f"    エラー: {result['stderr'][:100]}...")
    
    # 失敗したテストの詳細（verbose時のみ）
    if args.verbose:
        failed_tests = [name for name, result in results.items() if not result["success"]]
        if failed_tests:
            print(f"\n失敗したテストの詳細:")
            for test_name in failed_tests:
                result = results[test_name]
                print(f"\n--- {test_name} ---")
                if result["stderr"]:
                    print(f"エラー出力:\n{result['stderr']}")
                if result["stdout"]:
                    print(f"標準出力:\n{result['stdout'][:500]}...")
    
    # 終了
    print(f"\n終了時刻: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    if successful_tests == total_tests:
        print("🎉 すべてのテストが成功しました！")
        exit_code = 0
    elif successful_tests > 0:
        print("⚠️ 一部のテストが失敗しました。")
        exit_code = 1
    else:
        print("❌ すべてのテストが失敗しました。")
        exit_code = 2
    
    return exit_code

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
