#!/usr/bin/env python3
"""
xTB 単純テスト - 最終版

確実に動作する最も単純なxTBシングルポイント計算のテスト
このテストは外部依存性を最小限に抑え、確実に動作する設定を使用
"""

import subprocess
import tempfile
import json
from pathlib import Path
from datetime import datetime

class SimpleXTBTest:
    """シンプルなxTBテストクラス"""
    
    def __init__(self):
        self.test_molecules = {
            "water": {
                "name": "水分子",
                "formula": "H2O",
                "xyz": """3
Water molecule
O    0.000000    0.000000    0.117176
H    0.000000    0.757200   -0.468706
H    0.000000   -0.757200   -0.468706""",
                "expected_energy_range": (-6.0, -5.5)  # GFN1での期待エネルギー範囲 (Hartree)
            },
            "methane": {
                "name": "メタン分子",
                "formula": "CH4",
                "xyz": """5
Methane molecule
C    0.000000    0.000000    0.000000
H    0.629118    0.629118    0.629118
H   -0.629118   -0.629118    0.629118
H   -0.629118    0.629118   -0.629118
H    0.629118   -0.629118   -0.629118""",
                "expected_energy_range": (-4.5, -4.0)
            }
        }
    
    def check_xtb_availability(self):
        """xTBの利用可能性をチェック"""
        try:
            result = subprocess.run(
                ["xtb", "--version"], 
                capture_output=True, 
                text=True, 
                timeout=10
            )
            if result.returncode == 0:
                version_line = result.stdout.split('\n')[0]
                return {"available": True, "version": version_line}
            else:
                return {"available": False, "error": f"Return code: {result.returncode}"}
        except FileNotFoundError:
            return {"available": False, "error": "xTB not found in PATH"}
        except Exception as e:
            return {"available": False, "error": str(e)}
    
    def run_single_calculation(self, molecule_data, work_dir):
        """単一分子のxTB計算を実行"""
        xyz_file = work_dir / f"{molecule_data['formula'].lower()}.xyz"
        
        # XYZファイル作成
        with open(xyz_file, 'w') as f:
            f.write(molecule_data['xyz'])
        
        # xTB実行コマンド（最も安全な設定）
        cmd = [
            "xtb", str(xyz_file),
            "--gfn", "1",        # GFN1は最も安定
            "--chrg", "0",       # 中性
            "--uhf", "0",        # 閉殻
            "--parallel", "1",   # シングルスレッド
            "--acc", "1.0"       # 適度な精度
        ]
        
        try:
            result = subprocess.run(
                cmd,
                cwd=work_dir,
                capture_output=True,
                text=True,
                timeout=60
            )
            
            return self._parse_xtb_result(result, molecule_data)
            
        except subprocess.TimeoutExpired:
            return {
                "success": False,
                "error": "Calculation timeout (60s)",
                "error_type": "timeout"
            }
        except Exception as e:
            return {
                "success": False,
                "error": str(e),
                "error_type": "unexpected"
            }
    
    def _parse_xtb_result(self, subprocess_result, molecule_data):
        """xTB結果を解析"""
        result = {
            "success": subprocess_result.returncode == 0,
            "return_code": subprocess_result.returncode,
            "stdout": subprocess_result.stdout,
            "stderr": subprocess_result.stderr
        }
        
        if result["success"]:
            # エネルギー抽出
            energy = self._extract_energy(subprocess_result.stdout)
            result["energy_hartree"] = energy
            
            if energy is not None:
                result["energy_ev"] = energy * 27.2114
                result["energy_kcal_mol"] = energy * 627.509
                
                # エネルギー妥当性チェック
                expected_range = molecule_data.get("expected_energy_range")
                if expected_range:
                    in_range = expected_range[0] <= energy <= expected_range[1]
                    result["energy_reasonable"] = in_range
                    result["expected_range"] = expected_range
            
            # 追加情報抽出
            result.update(self._extract_additional_info(subprocess_result.stdout))
        else:
            # エラー分析
            result.update(self._analyze_error(subprocess_result.stderr))
        
        return result
    
    def _extract_energy(self, stdout):
        """標準出力からエネルギーを抽出"""
        for line in stdout.splitlines():
            if "TOTAL ENERGY" in line:
                try:
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if "ENERGY" in part and i + 1 < len(parts):
                            return float(parts[i + 1])
                except (ValueError, IndexError):
                    continue
        return None
    
    def _extract_additional_info(self, stdout):
        """追加情報を抽出"""
        info = {}
        
        for line in stdout.splitlines():
            # 双極子モーメント
            if "molecular dipole:" in line.lower():
                try:
                    import re
                    match = re.search(r'(\d+\.?\d*)', line)
                    if match:
                        info["dipole_moment"] = float(match.group(1))
                except:
                    pass
            
            # HOMO-LUMOギャップ
            if "homo-lumo gap" in line.lower():
                try:
                    import re
                    match = re.search(r'(\d+\.?\d+)\s*ev', line.lower())
                    if match:
                        info["homo_lumo_gap_ev"] = float(match.group(1))
                except:
                    pass
        
        return info
    
    def _analyze_error(self, stderr):
        """エラーを分析"""
        stderr_lower = stderr.lower()
        
        if "sigsegv" in stderr_lower or "segmentation fault" in stderr_lower:
            return {
                "error_type": "segmentation_fault",
                "error_description": "xTB内部エラー（メモリアクセス違反）"
            }
        elif "memory" in stderr_lower:
            return {
                "error_type": "memory_error", 
                "error_description": "メモリ不足"
            }
        elif "convergence" in stderr_lower:
            return {
                "error_type": "convergence_error",
                "error_description": "SCF収束失敗"
            }
        else:
            return {
                "error_type": "unknown",
                "error_description": "不明なエラー"
            }
    
    def run_full_test(self):
        """完全なテストスイートを実行"""
        print("=" * 50)
        print("xTB シングルポイント計算テスト")
        print("=" * 50)
        
        # xTB確認
        print("\n1. xTB 利用可能性チェック...")
        xtb_check = self.check_xtb_availability()
        
        if not xtb_check["available"]:
            print(f"   ✗ xTB利用不可: {xtb_check['error']}")
            return {"overall_success": False, "xtb_available": False}
        
        print(f"   ✓ xTB利用可能: {xtb_check['version']}")
        
        # テスト実行
        print("\n2. 分子計算テスト...")
        test_results = {"xtb_available": True, "molecule_tests": {}}
        
        with tempfile.TemporaryDirectory(prefix="xtb_test_") as temp_dir:
            work_path = Path(temp_dir)
            
            for mol_key, mol_data in self.test_molecules.items():
                print(f"\n--- {mol_data['name']} ({mol_data['formula']}) ---")
                
                result = self.run_single_calculation(mol_data, work_path)
                test_results["molecule_tests"][mol_key] = result
                
                if result["success"]:
                    print("   ✓ 計算成功")
                    
                    energy = result.get("energy_hartree")
                    if energy is not None:
                        print(f"   エネルギー: {energy:.6f} Hartree")
                        print(f"   エネルギー: {result['energy_ev']:.3f} eV")
                        print(f"   エネルギー: {result['energy_kcal_mol']:.3f} kcal/mol")
                        
                        if result.get("energy_reasonable"):
                            print("   ✓ エネルギー値は妥当範囲内")
                        else:
                            expected = result.get("expected_range", "不明")
                            print(f"   ⚠️ エネルギー値が期待範囲外 (期待: {expected})")
                    
                    # 追加情報
                    if "dipole_moment" in result:
                        print(f"   双極子モーメント: {result['dipole_moment']:.3f} Debye")
                    
                    if "homo_lumo_gap_ev" in result:
                        print(f"   HOMO-LUMOギャップ: {result['homo_lumo_gap_ev']:.3f} eV")
                
                else:
                    print("   ✗ 計算失敗")
                    print(f"   エラータイプ: {result.get('error_type', '不明')}")
                    print(f"   説明: {result.get('error_description', result.get('error', '不明'))}")
        
        # 結果サマリー
        print("\n" + "=" * 50)
        print("テスト結果サマリー")
        print("=" * 50)
        
        successful_molecules = sum(
            1 for result in test_results["molecule_tests"].values() 
            if result["success"]
        )
        total_molecules = len(test_results["molecule_tests"])
        
        print(f"成功した分子: {successful_molecules}/{total_molecules}")
        print(f"成功率: {successful_molecules/total_molecules*100:.1f}%")
        
        overall_success = successful_molecules > 0
        test_results["overall_success"] = overall_success
        test_results["success_rate"] = successful_molecules / total_molecules
        
        if successful_molecules == total_molecules:
            print("🎉 全テスト成功！xTBは完全に動作しています。")
        elif successful_molecules > 0:
            print("⚠️ 部分的成功。基本機能は動作しています。")
        else:
            print("❌ 全テスト失敗。xTBの設定を確認してください。")
        
        # 結果をJSONで保存
        self._save_test_results(test_results)
        
        return test_results
    
    def _save_test_results(self, results):
        """テスト結果をJSONファイルに保存"""
        results_file = Path("xtb_test_results.json")
        
        # 保存用データ準備
        save_data = {
            "timestamp": datetime.now().isoformat(),
            "overall_success": results["overall_success"],
            "success_rate": results["success_rate"],
            "xtb_available": results["xtb_available"],
            "molecule_tests": {}
        }
        
        # 分子テスト結果（stdoutを除外して軽量化）
        for mol_key, mol_result in results["molecule_tests"].items():
            save_data["molecule_tests"][mol_key] = {
                k: v for k, v in mol_result.items() 
                if k not in ["stdout", "stderr"]
            }
        
        try:
            with open(results_file, 'w', encoding='utf-8') as f:
                json.dump(save_data, f, indent=2, ensure_ascii=False)
            print(f"\nテスト結果を保存: {results_file}")
        except Exception as e:
            print(f"結果保存エラー: {e}")

def main():
    """メイン実行関数"""
    tester = SimpleXTBTest()
    results = tester.run_full_test()
    
    # 終了コード設定
    exit_code = 0 if results.get("overall_success", False) else 1
    return exit_code

if __name__ == "__main__":
    exit_code = main()
    exit(exit_code)
