# Computational Chemistry Tool (streamlit-pyscf)

![Python](https://img.shields.io/badge/python-3.7+-blue.svg)
![PySCF](https://img.shields.io/badge/PySCF-2.x-green.svg)
![Gaussian](https://img.shields.io/badge/Gaussian-16/09-blue.svg)
![OpenMM](https://img.shields.io/badge/OpenMM-8.x-orange.svg)
![xTB](https://img.shields.io/badge/xTB-6.x-purple.svg)
![MOPAC](https://img.shields.io/badge/MOPAC-2016-red.svg)
![Streamlit](https://img.shields.io/badge/Streamlit-1.x-red.svg)

**多様な計算化学プログラムを統合した、包括的な計算化学プラットフォーム**

このプロジェクトは、PySCF、Gaussian、OpenMM、MOPAC、xTBなど、様々な計算化学プログラムを統合し、分子力場から量子化学計算まで幅広い計算手法を直感的なWebインターフェースで利用できるようにします。

## 🌟 特徴

- **統合プラットフォーム**: 複数の計算化学プログラムを一つのインターフェースで利用
- **GUIインストール不要**: ブラウザで動作する直感的なインターフェース
- **多階層計算手法**: 分子力場 → 半経験的 → 量子化学計算の階層的アプローチ
- **包括的な計算機能**: 構造最適化、振動解析、励起状態、NMR予測、分子動力学など
- **可視化サポート**: 分子構造と計算結果の3D表示
- **データベース管理**: 計算結果の自動保存と検索機能
- **プログラミング不要**: 化学の知識だけで使用可能
- **カスタマイズ可能**: ページ表示設定とユーザー設定

## 🚀 デモ

インストール不要のオンライン版（ローカル環境推奨）:

**[https://conputational-chemistry-tool.streamlit.app/](https://conputational-chemistry-tool.streamlit.app/)**

## 📦 機能一覧

### 🔧 前処理・ユーティリティ
- **構造チェック**: 分子構造の確認・バリデーション
- **フォーマット変換**: 様々な分子ファイル形式の変換

### 🔀 統合・ワークフロー
- **総合計算**: 複数の計算を統合した自動化ワークフロー

### 🧲 分子力場計算
- **配座探索**: RDKitを使用した配座解析
- **OpenMM分子動力学**: 高速分子動力学シミュレーション
  - 古典力場計算（AMBER、CHARMM、OPLSなど）
  - 溶媒効果シミュレーション
  - エネルギー最小化

### 🧬 半経験的手法
- **MOPAC計算**: PM6/PM7を使用した高速計算
  - 一点エネルギー計算
  - 構造最適化
  - 振動解析
  - 結果可視化
- **xTB計算**: 拡張タイトバインディング法による高速計算
  - GFN1-xTB、GFN2-xTB
  - 溶媒効果（GBSA、ALPB）
  - 配座解析

### 🧪 量子化学計算
- **PySCF計算**: オープンソース量子化学パッケージ
  - 一点エネルギー計算: エネルギー、軌道エネルギー、双極子モーメント
  - 構造最適化: 分子幾何の最適化
  - 振動解析: 振動数、熱力学量、IR/ラマンスペクトル
- **Gaussian計算**: 商用量子化学パッケージ（要ライセンス）
  - 高精度DFT計算
  - 相関エネルギー法（MP2、CCSD(T)など）
  - 大規模分子系の計算

### 🔍 可視化と解析
- **分子軌道可視化**: 3D分子構造と軌道の表示
- **エネルギー分解解析（EDA）**: 分子間相互作用の詳細解析
- **フラグメント分解**: xyz座標からのフラグメント解析
- **配座エネルギー分解**: 配座異性体のエネルギー解析

### ⚡ 物性計算
- **イオン化ポテンシャル**: 電子親和力の理論計算
- **溶媒効果**: PCM/ddCOSMOモデルを使用した溶媒計算
- **結合解離エネルギー**: 化学結合の強さの評価
- **pKa計算**: プロトン親和力の理論予測

### 📊 スペクトル計算
- **IR/ラマンスペクトル**: 振動スペクトルの予測
- **NMR予測**: 化学シフトの理論計算
- **分極率計算**: 分子の分極率
- **UV-Vis予測**: 時間依存DFT計算による励起状態

### 🔄 遷移状態計算
- **NEB計算**: 複雑な反応経路の探索
- **遷移状態探索**: 化学反応の活性化エネルギー
- **IRC計算**: 反応経路の追跡

### ⚙️ システム・設定
- **データベース**: 計算結果の管理と検索
- **結果集計**: 計算データの統計と可視化
- **システム設定**: アプリケーション設定
- **ページ設定**: 表示するページのカスタマイズ

## 📋 システム要件

### 必須要件
- Python 3.7以上
- 推奨メモリ: 8GB以上
- 推奨CPU: マルチコア（計算高速化のため）

### 計算エンジンの要件
- **PySCF**: Python環境に自動インストール
- **OpenMM**: Python環境に自動インストール
- **xTB**: 実行可能ファイルのパス設定が必要
- **MOPAC**: 実行可能ファイルのパス設定が必要
- **Gaussian**: 商用ライセンスと実行環境の設定が必要（オプション）

## 🛠️ インストール

### 1. リポジトリのクローン
```bash
git clone https://github.com/poclab-web/streamlit-pyscf.git
cd streamlit-pyscf
```

### 2. 仮想環境の作成（推奨）
```bash
# Condaを使用
conda create -n streamlit-pyscf python=3.11
conda activate streamlit-pyscf

# venvを使用
python -m venv streamlit-pyscf
source streamlit-pyscf/bin/activate  # Linux/Mac
# streamlit-pyscf\Scripts\activate  # Windows
```

### 3. 依存関係のインストール
```bash
pip install -r requirements.txt
```

**主要な依存関係:**
- `pyscf`: 量子化学計算エンジン（オープンソース）
- `streamlit`: Webアプリケーションフレームワーク
- `rdkit`: ケモインフォマティクスライブラリ
- `matplotlib`, `plotly`: データ可視化
- `py3Dmol`, `stmol`: 3D分子表示
- `openmm`: 分子動力学シミュレーション
- `geometric`, `pyberny`: 構造最適化
- `ase`: 原子シミュレーション環境

### 4. 外部プログラムの設定（オプション）
以下のプログラムは必要に応じて別途インストール・設定してください：

```bash
# xTB（半経験的計算）
# https://github.com/grimme-lab/xtb からダウンロード・インストール

# MOPAC（半経験的計算）  
# http://openmopac.net/ からダウンロード・インストール

# Gaussian（量子化学計算、要ライセンス）
# ライセンス購入後、公式サイトからインストール
```

設定ファイル `config/external_software_config.py` で各プログラムのパスを設定できます。

## 🚀 使用方法

### 基本的な起動
```bash
streamlit run ComputationalChemistryTool.py
```

ブラウザが自動的に開き、`http://localhost:8501`でアプリケーションにアクセスできます。

### 基本的な使用手順
1. **分子入力**: SMILES記法またはXYZ座標で分子を入力
2. **計算設定**: 理論レベル（HF、DFTなど）と基底関数を選択
3. **計算実行**: 計算ボタンをクリックして実行
4. **結果確認**: 可視化された結果を確認・ダウンロード

### サポートされている計算手法と理論レベル

#### 分子力場
- **OpenMM力場**: AMBER、CHARMM、OPLS-AA、GAFF など

#### 半経験的手法
- **MOPAC**: PM6、PM7、PM3、AM1、MNDO
- **xTB**: GFN1-xTB、GFN2-xTB、GFN-FF

#### 量子化学計算
- **Hartree-Fock（HF）**
- **密度汎関数理論（DFT）**: B3LYP、PBE、M06-2X、ωB97X-D など
- **相関理論**: MP2、CCSD、CCSD(T)（Gaussianで利用可能）
- **分散力補正**: D3、D4補正

### サポートされている基底関数
- STO-3G、3-21G、6-31G(d,p)
- cc-pVDZ、cc-pVTZ、aug-cc-pVDZ
- def2-SVP、def2-TZVP

## 📁 プロジェクト構造

```
streamlit-pyscf/
├── ComputationalChemistryTool.py   # メインアプリケーション
├── pages/                          # 機能ページ
│   ├── 001_StructureCheck.py       # 構造チェック
│   ├── 101_GeneralCalculation.py   # 総合計算
│   ├── 201_RDKit_ConfSearch.py     # RDKit配座探索
│   ├── 211_openmm_singlepoint.py   # OpenMM計算
│   ├── 301-304_mopac_*.py          # MOPAC計算群
│   ├── 311_xtb_calculation_suite.py# xTB計算
│   ├── 401-403_*.py                # 量子化学計算
│   ├── 501-504_*.py                # 可視化・解析
│   ├── 601-604_*.py                # 物性計算
│   ├── 701-704_*.py                # スペクトル計算
│   ├── 801-803_*.py                # 遷移状態計算
│   └── 901-904_*.py                # システム設定
├── logic/                          # 計算ロジック
│   ├── calculation.py              # 量子化学計算
│   ├── molecule_handler.py         # 分子データ処理
│   ├── visualization.py            # 結果可視化
│   └── database.py                 # データベース管理
├── config/                         # 設定ファイル
│   ├── page_visibility.json        # ページ表示設定
│   ├── user_preferences.py         # ユーザー設定
│   └── external_software_config.py # 外部ソフトウェア設定
├── data/                          # 計算データ保存
├── controllers/                   # 計算制御
├── utils/                         # ユーティリティ
└── tests/                         # テストコード
```

## 💡 使用例

### 例1: 階層的計算アプローチ
1. **分子力場での初期構造最適化**
   - 「分子力場 > RDKit配座探索」で初期配座を生成
   - 「分子力場 > OpenMM計算」で粗い最適化
2. **半経験的手法での精密化**
   - 「半経験的 > xTB計算」で中程度の精度で最適化
3. **量子化学計算での高精度計算**
   - 「量子化学計算 > PySCF構造最適化」で最終的な高精度計算

### 例2: 水分子の構造最適化（PySCF）
1. 「量子化学計算 > 構造最適化」ページを選択
2. SMILES入力: `O`
3. 理論レベル: `B3LYP`、基底関数: `6-31G(d,p)`を選択
4. 「最適化実行」をクリック

### 例3: ベンゼンの分子動力学シミュレーション（OpenMM）
1. 「分子力場 > OpenMM計算」ページを選択
2. SMILES入力: `c1ccccc1`
3. 力場: `GAFF`、温度: `300K`を設定
4. 「MD実行」をクリック

### 例4: 有機分子の半経験的計算（xTB）
1. 「半経験的 > xTB計算」ページを選択
2. SMILES入力: `CCO`（エタノール）
3. 手法: `GFN2-xTB`、溶媒: `水`を選択
4. 「計算実行」をクリック

### 例5: 分子間相互作用のEDA解析
1. 「可視化と解析 > エネルギー分解解析」ページを選択
2. 複合体分子を入力
3. フラグメント分割を設定
4. EDA解析を実行

## 🔬 計算例とベンチマーク

### 推定計算時間（タンパク質アミノ酸残基1個程度の分子）

| 計算手法 | 理論レベル | 基底関数/設定 | 一点計算 | 構造最適化 |
|---------|-----------|-------------|---------|-----------|
| **分子力場** | OpenMM | GAFF | ~0.1秒 | ~1秒 |
| **半経験的** | xTB | GFN2-xTB | ~0.5秒 | ~5秒 |
| **半経験的** | MOPAC | PM7 | ~1秒 | ~10秒 |
| **量子化学** | PySCF | B3LYP/6-31G(d,p) | ~10秒 | ~5分 |
| **量子化学** | Gaussian | B3LYP/6-31G(d,p) | ~5秒 | ~3分 |

### 大きな分子系での計算時間（ベンゼン = 12原子）

| 計算手法 | 理論レベル | 設定 | 一点計算 | 構造最適化 |
|---------|-----------|-----|---------|-----------|
| **分子力場** | OpenMM | GAFF | ~0.1秒 | ~2秒 |
| **半経験的** | xTB | GFN2-xTB | ~1秒 | ~10秒 |
| **半経験的** | MOPAC | PM7 | ~2秒 | ~30秒 |
| **量子化学** | PySCF | B3LYP/6-31G(d,p) | ~30秒 | ~15分 |

*計算時間は環境により大きく変動します

### 精度とコストのトレードオフ
| 計算手法 | 計算コスト | 構造精度 | スペクトル精度 | 適用範囲 |
|---------|-----------|---------|-------------|---------|
| **分子力場** | 非常に低 | 中 | × | 大規模系、MD |
| **半経験的** | 低 | 中-高 | 中 | 中規模系、スクリーニング |
| **DFT** | 中 | 高 | 高 | 一般的な化学系 |
| **ab initio** | 高 | 非常に高 | 非常に高 | 高精度が必要な系 |

## 🐛 トラブルシューティング

### よくある問題

**Q: 計算が失敗する**
- 分子構造が正しいか確認
- より高速な計算手法から試行（分子力場 → 半経験的 → 量子化学）
- より小さな基底関数を試行（量子化学計算の場合）
- 電荷とスピン多重度の設定を確認

**Q: メモリ不足エラー**
- より高速な計算手法を使用（OpenMM、xTB、MOPAC）
- より小さな基底関数を使用（量子化学計算の場合）
- 分子サイズを小さくする
- 利用可能メモリを増やす

**Q: 外部プログラムが見つからない**
- `config/external_software_config.py`でパス設定を確認
- 該当プログラムが正しくインストールされているか確認
- 環境変数の設定を確認

**Q: 収束しない**
- より粗い計算手法で初期構造を最適化
- 初期構造を変更
- 収束条件を緩和
- 異なる理論レベルを試行

**Q: Gaussianライセンスエラー**
- 有効なライセンスがあることを確認
- Gaussianの環境設定を確認
- PySCFでの代替計算を検討

## 🎨 ページカスタマイズ

アプリケーションでは表示するページをカスタマイズできます：

1. 「システム・設定 > ページ設定」にアクセス
2. カテゴリ別にページの表示/非表示を設定
3. ユーザーごとに異なる設定が可能

利用可能なカテゴリ：
- 🔧 前処理・ユーティリティ
- 🔀 統合・ワークフロー  
- 🧲 分子力場
- 🧬 半経験的
- 🧪 量子化学計算
- 🔍 可視化と解析
- ⚡ 物性計算
- 📊 スペクトル計算
- 🔄 遷移状態計算
- ⚙️ システム・設定

## 🤝 貢献

プロジェクトへの貢献を歓迎します！

### 開発環境のセットアップ
```bash
# 開発用依存関係のインストール
pip install -r requirements-dev.txt

# pre-commitフックのインストール
pre-commit install
```

### 貢献ガイドライン
1. このリポジトリをフォーク
2. 機能ブランチを作成: `git checkout -b feature/new-feature`
3. 変更をコミット: `git commit -m "新機能を追加"`
4. ブランチをプッシュ: `git push origin feature/new-feature`
5. プルリクエストを作成

### バグ報告・機能リクエスト
Issuesページで報告してください：
**[https://github.com/poclab-web/streamlit-pyscf/issues](https://github.com/poclab-web/streamlit-pyscf/issues)**

## 📖 ドキュメント

### 詳細ドキュメント
- [アーキテクチャ概要](docs/architecture.md)
- [API リファレンス](docs/api.md)
- [開発者ガイド](docs/development.md)

### 学習リソース
- [PySCF公式ドキュメント](https://pyscf.org/)
- [量子化学計算入門](https://example.com/quantum-chemistry)
- [DFT理論の基礎](https://example.com/dft-basics)

## 🏆 引用

この ソフトウェアを研究で使用する場合は、以下のように引用してください：

```bibtex
@software{computational_chemistry_tool,
  title = {Computational Chemistry Tool: Integrated Multi-Engine Platform},
  author = {PocLab},
  url = {https://github.com/poclab-web/streamlit-pyscf},
  year = {2024}
}
```

### 依存ライブラリ・プログラムの引用
- **RDKit**: Landrum, G. RDKit: Open-source cheminformatics.
- **OpenMM**: Eastman, P. et al. *PLoS Comput Biol* **2017**, *13*, e1005659.
- **xTB**: Bannwarth, C. et al. *J Chem Theory Comput* **2019**, *15*, 1652-1671.
- **MOPAC**: Stewart, J. J. P. *J Mol Model* **2013**, *19*, 1-32.
- **PySCF**: Sun, Q. et al. *WIREs Comput Mol Sci* **2018**, *8*, e1340.
- **Gaussian**: Frisch, M. J. et al. Gaussian 16, Revision C.01, Gaussian, Inc., Wallingford CT, 2016.

## 🛡️ セキュリティ

セキュリティ問題を発見した場合は、公開Issuesではなく開発者に直接メールしてください。

## 👥 開発チーム

- **PocLab** - プロジェクトリード
- 全ての貢献者は[Contributors](https://github.com/poclab-web/streamlit-pyscf/contributors)をご覧ください

## 🔗 関連リンク

### 統合プログラム
- [PySCF公式サイト](https://pyscf.org/)
- [OpenMM公式サイト](http://openmm.org/)
- [xTB公式サイト](https://xtb-docs.readthedocs.io/)
- [MOPAC公式サイト](http://openmopac.net/)
- [Gaussian公式サイト](https://gaussian.com/)

### フレームワーク・ライブラリ
- [Streamlit公式サイト](https://streamlit.io/)
- [RDKit公式サイト](https://www.rdkit.org/)
- [ASE公式サイト](https://wiki.fysik.dtu.dk/ase/)

### 学習リソース
- [量子化学計算ベストプラクティス](https://example.com/qc-best-practices)
- [分子動力学入門](https://example.com/md-basics)
- [計算化学手法比較](https://example.com/computational-methods)

## 📊 統計

![GitHub stars](https://img.shields.io/github/stars/poclab-web/streamlit-pyscf?style=social)
![GitHub forks](https://img.shields.io/github/forks/poclab-web/streamlit-pyscf?style=social)
![GitHub issues](https://img.shields.io/github/issues/poclab-web/streamlit-pyscf)
![GitHub last commit](https://img.shields.io/github/last-commit/poclab-web/streamlit-pyscf)

---

**🔬 分子力場から量子化学まで、包括的な計算化学プラットフォーム。**
