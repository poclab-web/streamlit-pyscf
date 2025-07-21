# Troubleshooting

## Common Issues

### Installation Problems

#### PySCF Installation Error
```bash
# Ubuntu/Debian
sudo apt-get install libblas-dev liblapack-dev

# macOS
brew install openblas

# Then reinstall PySCF
pip install pyscf
```

#### RDKit Import Error
```bash
# Install via conda
conda install -c conda-forge rdkit

# Or via pip
pip install rdkit
```

### Runtime Errors

#### Memory Error
大きな分子系で計算する場合、メモリ不足が発生することがあります：

```python
# メモリ制限を調整
export OMP_NUM_THREADS=1
ulimit -v 8000000  # 8GB制限
```

#### Convergence Issues
SCF計算が収束しない場合：

1. より小さな基底関数を使用
2. 初期軌道を変更
3. レベルシフトを適用

### Performance Issues

#### Slow Calculations
- より小さな基底関数セット（STO-3G, 3-21G）を使用
- 対称性を利用（C1以外）
- 並列計算の設定確認

#### Large Output Files
- 出力レベルを調整（verbose設定）
- 不要なファイルの削除
- ディスク容量の確認

## Debug Mode

デバッグモードでの実行：

```bash
export PYSCF_DEBUG=1
streamlit run ComputationalChemistryTool.py
```

## Contact

問題が解決しない場合は、以下の情報を含めてお問い合わせください：

- エラーメッセージの全文
- 使用環境（OS、Python版、依存関係の版）
- 入力データ（SMILES、XYZ座標）
- 計算設定
