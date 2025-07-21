#!/usr/bin/env python3
"""
xTB セーフモードテスト

セグメンテーションフォルトを回避するための安全な設定でテスト
"""

import subprocess
import tempfile
import os
from pathlib import Path

def create_safe_water_xyz():
    """より安全な水分子のXYZ座標"""
    # 文献値を使用した安定な水分子構造
    xyz_content = """3
Water molecule - safe geometry
O    0.000000    0.000000    0.117176
H    0.000000    0.757200   -0.468706
H    0.000000   -0.757200   -0.468706
"""
    return xyz_content.strip()

def run_safe_xtb_test():
    """安全な設定でxTB計算を実行"""
    print("=== xTB セーフモードテスト ===")
    
    with tempfile.TemporaryDirectory(prefix="xtb_safe_") as temp_dir:
        temp_path = Path(temp_dir)
        xyz_file = temp_path / "water_safe.xyz"
        
        # XYZファイル作成
        with open(xyz_file, 'w') as f:
            f.write(create_safe_water_xyz())
        
        print(f"入力ファイル: {xyz_file}")
        print("分子: 水 (H2O) - 安定構造")
        
        # より安全な設定
        test_configs = [
            {
                "name": "GFN0 (最軽量)",
                "gfn": "0",
                "extra": ["--acc", "2.0"]  # 精度を下げて安定性向上
            },
            {
                "name": "GFN1 (標準)",
                "gfn": "1", 
                "extra": ["--acc", "1.0"]
            },
            {
                "name": "GFN2 (高精度)",
                "gfn": "2",
                "extra": ["--acc", "0.5"]
            }
        ]
        
        for config in test_configs:
            print(f"\n--- {config['name']} テスト ---")
            
            cmd = [
                "xtb", 
                str(xyz_file),
                "--gfn", config["gfn"],
                "--chrg", "0",
                "--uhf", "0",
                "--parallel", "1"  # シングルスレッドで実行
            ]
            cmd.extend(config["extra"])
            
            print(f"コマンド: {' '.join(cmd)}")
            
            try:
                result = subprocess.run(
                    cmd,
                    cwd=temp_dir,
                    capture_output=True,
                    text=True,
                    timeout=30  # 短時間でタイムアウト
                )
                
                if result.returncode == 0:
                    print("✓ 計算成功")
                    
                    # エネルギー抽出
                    energy = None
                    for line in result.stdout.splitlines():
                        if "TOTAL ENERGY" in line:
                            try:
                                parts = line.split()
                                for i, part in enumerate(parts):
                                    if "ENERGY" in part and i + 1 < len(parts):
                                        energy = float(parts[i + 1])
                                        break
                            except (ValueError, IndexError):
                                continue
                    
                    if energy is not None:
                        print(f"  エネルギー: {energy:.6f} Hartree")
                        print(f"  エネルギー: {energy * 27.2114:.3f} eV")
                    else:
                        print("  ⚠️ エネルギー抽出失敗")
                        
                else:
                    print("✗ 計算失敗")
                    print(f"  return code: {result.returncode}")
                    
                    # エラー分析
                    stderr = result.stderr.lower()
                    if "sigsegv" in stderr or "segmentation" in stderr:
                        print("  エラータイプ: セグメンテーションフォルト")
                        print("  → xTB内部エラーまたはメモリ問題")
                    elif "memory" in stderr:
                        print("  エラータイプ: メモリ不足")
                    elif "convergence" in stderr:
                        print("  エラータイプ: 収束失敗")
                    else:
                        print(f"  stderr: {result.stderr[:200]}...")
                
            except subprocess.TimeoutExpired:
                print("✗ タイムアウト (30秒)")
            except Exception as e:
                print(f"✗ 予期しないエラー: {e}")
        
        print(f"\n生成ファイル:")
        for file in temp_path.glob("*"):
            print(f"  {file.name}")

def main():
    """メイン実行"""
    
    # xTB確認
    try:
        result = subprocess.run(["xtb", "--version"], capture_output=True, text=True, timeout=5)
        if result.returncode != 0:
            print("xTBが利用できません")
            return
    except:
        print("xTBが見つかりません")
        return
    
    # セーフテスト実行
    run_safe_xtb_test()
    
    print("\n=== テスト完了 ===")

if __name__ == "__main__":
    main()
