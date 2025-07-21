#!/usr/bin/env python3
"""
最も単純なxTBテスト - 外部依存性なし

このスクリプトは：
1. xTBの基本動作確認
2. 最小限の水分子計算
3. 結果の検証

実行方法:
python minimal_xtb_test.py
"""

import subprocess
import tempfile
import os
from pathlib import Path

def create_water_xyz():
    """水分子のXYZファイル内容を作成"""
    xyz_content = """3
Water molecule
O    0.0000000    0.0000000    0.1173000
H    0.0000000    0.7572000   -0.4692000
H    0.0000000   -0.7572000   -0.4692000
"""
    return xyz_content.strip()

def check_xtb():
    """xTBの動作確認"""
    print("=== 最小限 xTB テスト ===")
    print("\n1. xTB実行可能性チェック...")
    
    try:
        result = subprocess.run(
            ["xtb", "--version"], 
            capture_output=True, 
            text=True, 
            timeout=10
        )
        
        if result.returncode == 0:
            version_info = result.stdout.strip().split('\n')[0]
            print(f"   ✓ xTB見つかりました: {version_info}")
            return True
        else:
            print(f"   ✗ xTBエラー (return code: {result.returncode})")
            print(f"   stderr: {result.stderr}")
            return False
            
    except FileNotFoundError:
        print("   ✗ xTB実行ファイルが見つかりません")
        print("   インストール方法:")
        print("     conda install -c conda-forge xtb")
        print("     または brew install xtb (macOS)")
        return False
    except subprocess.TimeoutExpired:
        print("   ✗ xTBバージョンチェックがタイムアウト")
        return False
    except Exception as e:
        print(f"   ✗ 予期しないエラー: {e}")
        return False

def run_minimal_calculation():
    """最小限のxTB計算を実行"""
    print("\n2. 最小計算実行...")
    
    # 一時ディレクトリで計算
    with tempfile.TemporaryDirectory(prefix="xtb_test_") as temp_dir:
        temp_path = Path(temp_dir)
        xyz_file = temp_path / "water.xyz"
        
        # XYZファイル作成
        print("   分子: 水 (H2O)")
        with open(xyz_file, 'w') as f:
            f.write(create_water_xyz())
        print(f"   入力ファイル: {xyz_file.name}")
        
        # xTB実行（シングルポイント、GFN2、気相）
        print("   計算設定:")
        print("     - シングルポイント計算")
        print("     - GFN2-xTB")
        print("     - 気相（溶媒なし）")
        print("     - 電荷: 0")
        print("     - 多重度: 1")
        
        cmd = [
            "xtb", 
            str(xyz_file),
            "--gfn", "2",
            "--chrg", "0",
            "--uhf", "0"
        ]
        
        try:
            print("   実行中...")
            result = subprocess.run(
                cmd,
                cwd=temp_dir,
                capture_output=True,
                text=True,
                timeout=60  # 1分でタイムアウト
            )
            
            print("\n3. 結果確認...")
            
            if result.returncode == 0:
                print("   ✓ 計算成功！")
                
                # エネルギー抽出
                energy = extract_energy(result.stdout)
                if energy is not None:
                    print(f"   全エネルギー: {energy:.8f} Hartree")
                    print(f"   全エネルギー: {energy * 27.2114:.4f} eV")
                    print(f"   全エネルギー: {energy * 627.509:.4f} kcal/mol")
                    
                    # 期待値との比較（水分子のGFN2エネルギーの概算値）
                    expected_range = (-5.5, -5.0)  # Hartree
                    if expected_range[0] <= energy <= expected_range[1]:
                        print("   ✓ エネルギー値は妥当な範囲内です")
                    else:
                        print(f"   ⚠️  エネルギー値が期待範囲外 ({expected_range[0]} ~ {expected_range[1]} Hartree)")
                else:
                    print("   ⚠️  エネルギー値を抽出できませんでした")
                
                # 追加情報
                dipole = extract_dipole(result.stdout)
                if dipole:
                    print(f"   双極子モーメント: {dipole:.3f} Debye")
                
                # 出力ファイル確認
                output_files = list(temp_path.glob("*"))
                print(f"   生成ファイル数: {len(output_files)}")
                
                return True
                
            else:
                print("   ✗ 計算失敗")
                print(f"   return code: {result.returncode}")
                print(f"   stdout: {result.stdout[:300]}...")
                print(f"   stderr: {result.stderr[:300]}...")
                return False
                
        except subprocess.TimeoutExpired:
            print("   ✗ 計算がタイムアウトしました（60秒）")
            return False
        except Exception as e:
            print(f"   ✗ 計算中にエラー: {e}")
            return False

def extract_energy(stdout):
    """xTB出力からエネルギーを抽出"""
    for line in stdout.splitlines():
        if "TOTAL ENERGY" in line:
            try:
                parts = line.split()
                # "TOTAL ENERGY" の後の数値を取得
                for i, part in enumerate(parts):
                    if "ENERGY" in part and i + 1 < len(parts):
                        return float(parts[i + 1])
            except (ValueError, IndexError):
                continue
    return None

def extract_dipole(stdout):
    """xTB出力から双極子モーメントを抽出"""
    for line in stdout.splitlines():
        if "molecular dipole:" in line.lower():
            try:
                import re
                match = re.search(r'(\d+\.?\d*)', line)
                if match:
                    return float(match.group(1))
            except (ValueError, AttributeError):
                continue
    return None

def main():
    """メイン実行"""
    
    # xTBチェック
    if not check_xtb():
        print("\n❌ xTBが利用できません")
        return False
    
    # 計算実行
    if run_minimal_calculation():
        print("\n🎉 テスト完了: 全て成功！")
        print("\nxTBは正常に動作しています。")
        return True
    else:
        print("\n❌ テスト失敗")
        return False

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)
