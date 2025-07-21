#!/usr/bin/env python3
"""
xTB計算の最も単純なテストスクリプト

このスクリプトは：
1. xTBのインストールを確認
2. 水分子のシングルポイント計算を実行
3. 結果を表示

使用方法:
python test_xtb_simple.py
"""

import sys
import os
from pathlib import Path

# プロジェクトのlogicディレクトリをパスに追加
current_dir = Path(__file__).parent
project_root = current_dir.parent.parent  # tests/xtb -> tests -> project_root
logic_dir = project_root / "logic"
sys.path.insert(0, str(logic_dir))

def test_xtb_simple():
    """最も基本的なxTBテスト"""
    
    try:
        from xtb_calculation import check_xtb_installation, XTBCalculator
    except ImportError as e:
        print(f"インポートエラー: {e}")
        print("必要なモジュールが見つかりません")
        return False
    
    print("=== xTB 最小限テスト ===")
    
    # xTBインストール確認
    print("\n1. xTBインストール確認...")
    xtb_status = check_xtb_installation()
    print(f"   インストール状況: {'✓ OK' if xtb_status['installed'] else '✗ NG'}")
    
    if not xtb_status['installed']:
        print(f"   エラー: {xtb_status['error']}")
        print("\n   xTBをインストールしてください:")
        print("   conda install -c conda-forge xtb")
        print("   または")
        print("   brew install xtb  # macOS")
        return False
    
    print(f"   バージョン: {xtb_status['version'].split()[0] if xtb_status['version'] else '不明'}")
    
    # 簡単な分子ハンドラーの模擬
    class SimpleMoleculeHandler:
        def __init__(self):
            # 水分子の座標（Å）
            self.coords = [
                ("O", 0.0000, 0.0000, 0.1173),
                ("H", 0.0000, 0.7572, -0.4692),
                ("H", 0.0000, -0.7572, -0.4692)
            ]
            self.mol = "H2O"  # 文字列表現
        
        def get_xyz_coordinates(self):
            return self.coords
    
    # 分子ハンドラーを作成
    print("\n2. テスト分子準備...")
    molecule_handler = SimpleMoleculeHandler()
    print("   分子: 水 (H2O)")
    print(f"   原子数: {len(molecule_handler.coords)}")
    
    # XTBCalculatorでテスト
    print("\n3. xTB計算実行...")
    try:
        # 一時的な作業ディレクトリを指定
        work_dir = Path("test_xtb_temp")
        calculator = XTBCalculator(molecule_handler, work_dir=work_dir)
        print(f"   作業ディレクトリ: {calculator.work_dir}")
        
        # シングルポイント計算（最も単純な設定）
        print("   計算タイプ: シングルポイント")
        print("   理論レベル: GFN2-xTB")
        print("   電荷: 0")
        print("   多重度: 1")
        print("   溶媒: なし（気相）")
        
        result = calculator.single_point_energy(
            gfn=2,           # GFN2-xTB
            charge=0,        # 中性
            uhf=0,           # 一重項
            solvent=None     # 気相
        )
        
        # 結果表示
        print("\n4. 計算結果...")
        if result['success']:
            print("   ✓ 計算成功！")
            
            energy = result.get('energy')
            if energy is not None:
                print(f"   全エネルギー: {energy:.8f} Hartree")
                print(f"   全エネルギー: {energy * 27.2114:.4f} eV")
                print(f"   全エネルギー: {energy * 627.509:.4f} kcal/mol")
            
            # 追加情報
            dipole = result.get('dipole_moment')
            if dipole:
                print(f"   双極子モーメント: {dipole:.3f} Debye")
            
            gap = result.get('homo_lumo_gap')
            if gap:
                print(f"   HOMO-LUMOギャップ: {gap:.3f} eV")
            
            print(f"   入力ファイル: {Path(result.get('input_file', '')).name}")
            
        else:
            print("   ✗ 計算失敗")
            print(f"   エラー: {result.get('error', '不明')}")
            
            # デバッグ情報
            if 'stderr' in result and result['stderr']:
                print(f"   詳細エラー: {result['stderr'][:200]}...")
            
            return False
        
        # ファイル確認
        print("\n5. 生成ファイル確認...")
        files = calculator.list_calculation_files()
        total_files = sum(len(file_list) for file_list in files.values())
        print(f"   生成ファイル数: {total_files}")
        
        for file_type, file_list in files.items():
            if file_list:
                print(f"   {file_type}: {len(file_list)}")
        
        print("\n=== テスト完了: 成功 ===")
        return True
        
    except Exception as e:
        print(f"\n   ✗ 予期しないエラー: {e}")
        import traceback
        print("   詳細:")
        traceback.print_exc()
        return False

def main():
    """メイン関数"""
    success = test_xtb_simple()
    
    if success:
        print("\n🎉 xTBは正常に動作しています！")
        exit_code = 0
    else:
        print("\n❌ xTBテストが失敗しました")
        exit_code = 1
    
    print("\nクリーンアップ...")
    # テスト用ディレクトリの削除
    test_dir = Path("test_xtb_temp")
    if test_dir.exists():
        import shutil
        try:
            shutil.rmtree(test_dir)
            print(f"   ✓ {test_dir} を削除しました")
        except Exception as e:
            print(f"   ⚠️  {test_dir} の削除に失敗: {e}")
    
    sys.exit(exit_code)

if __name__ == "__main__":
    main()
