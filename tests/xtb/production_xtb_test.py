#!/usr/bin/env python3
"""
プロダクション対応 xTB テストスイート

実際のstreamlit-pyscfプロジェクトで使用するためのテスト
MoleculeHandlerクラスと連携したテスト
"""

import sys
import os
from pathlib import Path

# プロジェクトのパスを追加
current_dir = Path(__file__).parent
project_root = current_dir.parent.parent  # tests/xtb -> tests -> project_root
logic_dir = project_root / "logic"
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(logic_dir))

def test_xtb_with_molecule_handler():
    """MoleculeHandlerと連携したxTBテスト"""
    
    print("=== プロダクション xTB テスト ===")
    
    # 必要なモジュールをインポート
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        print("✓ RDKit インポート成功")
    except ImportError as e:
        print(f"✗ RDKit インポートエラー: {e}")
        print("RDKitをインストールしてください: conda install -c conda-forge rdkit")
        return False
    
    try:
        from xtb_calculation import check_xtb_installation, XTBCalculator
        print("✓ xtb_calculation モジュール インポート成功")
    except ImportError as e:
        print(f"✗ xtb_calculation インポートエラー: {e}")
        return False
    
    # xTBインストール確認
    print("\n1. xTB インストール確認...")
    xtb_status = check_xtb_installation()
    if not xtb_status['installed']:
        print(f"✗ xTB が利用できません: {xtb_status['error']}")
        return False
    print("✓ xTB 利用可能")
    
    # テスト用分子の作成
    print("\n2. テスト分子作成...")
    test_molecules = [
        {"smiles": "O", "name": "水", "expected_atoms": 3},
        {"smiles": "CC", "name": "エタン", "expected_atoms": 8},
        {"smiles": "c1ccccc1", "name": "ベンゼン", "expected_atoms": 12}
    ]
    
    # 実際のMoleculeHandlerクラスの模擬版
    class TestMoleculeHandler:
        def __init__(self, mol):
            self.mol = mol
        
        def get_xyz_coordinates(self):
            if not self.mol or self.mol.GetNumConformers() == 0:
                return []
            
            conf = self.mol.GetConformer()
            coords = []
            for i, atom in enumerate(self.mol.GetAtoms()):
                pos = conf.GetAtomPosition(i)
                symbol = atom.GetSymbol()
                coords.append((symbol, pos.x, pos.y, pos.z))
            return coords
    
    successful_tests = 0
    total_tests = 0
    
    for mol_data in test_molecules:
        print(f"\n--- {mol_data['name']} テスト ---")
        total_tests += 1
        
        try:
            # 分子作成
            mol = Chem.MolFromSmiles(mol_data["smiles"])
            if mol is None:
                print(f"✗ 分子作成失敗: {mol_data['smiles']}")
                continue
            
            # 水素追加と3D構造生成
            mol = Chem.AddHs(mol)
            if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
                print("✗ 3D埋め込み失敗")
                continue
            
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
            
            print(f"  分子: {mol_data['name']} ({mol_data['smiles']})")
            print(f"  原子数: {mol.GetNumAtoms()}")
            
            # MoleculeHandler作成
            molecule_handler = TestMoleculeHandler(mol)
            
            # XTBCalculator作成（安全な設定）
            work_dir = Path("test_results") / f"{mol_data['name']}_xtb"
            calculator = XTBCalculator(molecule_handler, work_dir=work_dir)
            
            # GFN1でシングルポイント計算（最も安定）
            print("  計算実行: GFN1 シングルポイント...")
            result = calculator.single_point_energy(
                gfn=1,  # GFN1は最も安定
                charge=0,
                uhf=0,
                solvent=None
            )
            
            if result['success']:
                print("  ✓ 計算成功")
                energy = result.get('energy')
                if energy is not None:
                    print(f"    エネルギー: {energy:.6f} Hartree")
                    print(f"    エネルギー: {energy * 27.2114:.3f} eV")
                
                # 追加情報
                summary = result.get('summary', {})
                key_results = summary.get('key_results', {})
                for key, value in key_results.items():
                    if key != 'total_energy':  # エネルギーは既に表示済み
                        print(f"    {key}: {value}")
                
                successful_tests += 1
                
            else:
                print("  ✗ 計算失敗")
                error = result.get('error', '不明なエラー')
                print(f"    エラー: {error}")
                
                # エラー分析
                if 'sigsegv' in error.lower() or 'segmentation' in error.lower():
                    print("    → セグメンテーションフォルト（xTB内部エラー）")
                elif 'convergence' in error.lower():
                    print("    → SCF収束失敗")
        
        except Exception as e:
            print(f"  ✗ 予期しないエラー: {e}")
            import traceback
            traceback.print_exc()
    
    # 結果サマリー
    print(f"\n=== テスト結果サマリー ===")
    print(f"成功: {successful_tests}/{total_tests}")
    print(f"成功率: {successful_tests/total_tests*100:.1f}%")
    
    if successful_tests == total_tests:
        print("🎉 全テスト成功！xTBは正常に動作しています。")
        return True
    elif successful_tests > 0:
        print("⚠️ 一部テスト成功。基本機能は動作しています。")
        return True
    else:
        print("❌ 全テスト失敗。設定を確認してください。")
        return False

def main():
    """メイン実行"""
    success = test_xtb_with_molecule_handler()
    
    # クリーンアップ（オプション）
    print("\nクリーンアップ...")
    test_dir = Path("test_results")
    if test_dir.exists():
        import shutil
        try:
            shutil.rmtree(test_dir)
            print(f"✓ {test_dir} を削除しました")
        except Exception as e:
            print(f"⚠️ {test_dir} の削除に失敗: {e}")
    
    return success

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
