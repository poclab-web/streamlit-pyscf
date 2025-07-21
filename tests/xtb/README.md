# xTB テストスイート

このディレクトリには、xTB (extended tight-binding) 計算の動作確認のための複数のテストスクリプトが含まれています。

## ディレクトリ構造

```
tests/xtb/
├── README.md                    # このファイル
├── __init__.py                  # パッケージ初期化ファイル
├── run_tests.py                 # テストランナー（推奨実行方法）
├── minimal_xtb_test.py          # 最小限テスト
├── safe_xtb_test.py             # セーフモードテスト
├── final_xtb_test.py            # 包括的テスト（推奨）
├── test_xtb_simple.py           # プロジェクト統合テスト
├── production_xtb_test.py       # 本格的統合テスト
└── xtb_test_results.json       # テスト結果（自動生成）
```

## 実行方法

### 推奨: テストランナーを使用

```bash
# 基本的なテスト実行
cd tests/xtb
python run_tests.py

# 最小限のテストのみ
python run_tests.py --minimal

# 包括的なテストのみ
python run_tests.py --final

# 詳細出力付きですべてのテスト
python run_tests.py --all --verbose
```

### 個別テスト実行

```bash
cd tests/xtb

# 最小限テスト（xTBの基本動作確認）
python minimal_xtb_test.py

# 包括的テスト（複数分子でのベンチマーク）
python final_xtb_test.py

# セーフモードテスト（異なるGFNレベル）
python safe_xtb_test.py
```

### Python から実行

```python
# プロジェクトルートから
from tests.xtb import run_all_tests, run_basic_test

# 基本テスト
success, stdout, stderr = run_basic_test()

# 全テスト
results = run_all_tests()
```

## テスト結果の解釈

### 成功の場合
- ✓ 計算成功
- エネルギー値が妥当範囲内
- HOMO-LUMOギャップが出力される

### 失敗の場合
- セグメンテーションフォルト: xTB内部エラー（分子構造が原因の可能性）
- メモリエラー: システムリソース不足
- 収束エラー: SCF計算の収束失敗

## xTBインストール

テストが失敗する場合、xTBをインストールしてください：

```bash
# Conda環境の場合
conda install -c conda-forge xtb

# macOS (Homebrew)
brew install xtb

# その他のプラットフォーム
# https://github.com/grimme-lab/xtb からインストール
```

## 最も単純なテスト設定

確実に動作する最小設定：
- 理論レベル: GFN1
- 電荷: 0 (中性)
- 多重度: 1 (一重項)
- 溶媒: なし (気相)
- 精度: 標準 (--acc 1.0)

## 出力ファイル

- `xtb_test_results.json`: 詳細なテスト結果（final_xtb_test.pyで生成）

## トラブルシューティング

1. **xTB not found**: PATHにxTBが含まれていない
2. **Segmentation fault**: 分子構造が不安定、GFN1を使用
3. **Timeout**: 計算が複雑すぎる、より簡単な分子を使用
4. **Import errors**: 必要なPythonパッケージがインストールされていない
