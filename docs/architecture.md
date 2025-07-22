# Architecture

## システム概要

streamlit-pyscfは、PySCFとRDKitを基盤とした量子化学計算のWebアプリケーションです。

### 主要コンポーネント

1. **Streamlitフロントエンド**: ユーザーインターフェース
2. **Logic Layer**: 計算ロジックとデータ処理
3. **Controller Layer**: 高度な計算の制御
4. **Configuration**: 設定管理
5. **Data Storage**: 計算結果の永続化

## ディレクトリ構成

```
streamlit-pyscf/
├── .github/                            # GitHub設定
│   ├── ISSUE_TEMPLATE/                 # Issue テンプレート
│   │   ├── bug_report.md               # バグレポート
│   │   ├── calculation_issue.md        # 計算関連の問題
│   │   ├── config.yml                  # Issue設定
│   │   ├── documentation.md            # ドキュメント改善
│   │   ├── feature_request.md          # 機能要求
│   │   └── improvement.md              # 改善提案
│   └── workflows/                      # CI/CDワークフロー
│
├── config/                             # 計算設定と定数を管理
│   ├── __init__.py                     # パッケージ初期化
│   ├── config.py                       # 計算パラメータの設定
│   ├── external_software_config.py     # 外部ソフトウェア設定
│   ├── page_visibility.json            # ページ表示設定
│   ├── user_page_visibility.json       # ユーザー別ページ表示設定
│   ├── user_preferences.py             # ユーザー設定管理
│   ├── pKa_reference_list.csv           # pKa計算の参考データ
│   ├── solvents_epsilon.csv             # 溶媒の誘電率データ
│   ├── styles.css                       # streamlitアプリのスタイル設定
│   └── README.md                        # 設定ファイルの説明
│
├── controllers/                        # 高度な計算処理を制御
│   ├── __init__.py                     # パッケージ初期化
│   ├── energydecompositionanalysis.py  # エネルギー分解解析の制御
│   ├── excited_state_calculation.py    # 励起状態計算の制御
│   ├── fragment_calculation.py         # フラグメント計算（BDE、pKaなど）の制御
│   └── xyz_fragment_decomposition.py   # XYZ座標からの分子分解処理
│
├── data/                               # 計算結果データが保存される場所
│   ├── energy_db.sqlite                # エネルギー計算結果のデータベース
│   └── [分子別フォルダ]/                # 各分子の計算データフォルダ
│
├── docs/                               # ドキュメント
│   ├── _static/                        # 静的ファイル（画像など）
│   ├── _templates/                     # Sphinxテンプレート
│   ├── api_reference.md                # API参考文書
│   ├── architecture.md                 # この設計書（アーキテクチャ概要）
│   ├── conf.py                         # Sphinx設定
│   ├── contributing.md                 # 貢献ガイド
│   ├── examples.md                     # 使用例
│   ├── index.md                        # ドキュメント索引（MkDocs）
│   ├── index.rst                       # ドキュメント索引（Sphinx）
│   ├── installation.md                 # インストールガイド
│   ├── requirements.txt                # ドキュメント生成用依存関係
│   ├── troubleshooting.md              # トラブルシューティング
│   ├── usage.md                        # 使用方法
│   └── user_guide.md                   # ユーザーガイド
│
├── jobs/                               # 未計算ジョブ管理（現在未使用）
│
├── logic/                              # 計算の核となるロジック
│   ├── calculation.py                  # PySCFでの量子化学計算の主要ロジック
│   ├── data_loader.py                  # データの読み込み用ユーティリティ
│   ├── database.py                     # 計算結果のSQLデータベース操作
│   ├── logger.py                       # ログ出力・管理ユーティリティ
│   ├── molecule_handler.py             # RDKitを使った分子情報の変換・入力ファイル作成
│   ├── output_handler.py               # PySCFの出力結果の抽出・整形
│   ├── pubchem_api.py                  # PubChem APIを使った文献情報アクセス
│   └── visualization.py                # 計算結果の可視化・グラフ作成
│
├── pages/                              # streamlitアプリの各ページ（新しい番号体系）
│   # === 前処理・ユーティリティ (000-099) ===
│   ├── 001_StructureCheck.py           # 分子構造の確認・変換
│   
│   # === 統合・ワークフロー (100-199) ===
│   ├── 101_GeneralCalculation.py       # 一般的な量子化学計算
│   
│   # === 分子力場 (200-299) ===
│   ├── 201_RDKit_ConfSearch.py         # RDKitコンフォマー探索
│   ├── 211_openmm_singlepoint.py       # OpenMM分子力場計算
│   
│   # === 半経験的 (300-399) ===
│   ├── 301_mopac_singlepoint.py        # MOPAC一点計算
│   ├── 302_mopac_Opt.py                # MOPAC構造最適化
│   ├── 303_mopac_OptandFreq.py         # MOPAC最適化+振動解析
│   ├── 304_mopac_Visualization.py      # MOPAC結果可視化
│   ├── 311_xtb_calculation_suite.py    # xTB計算スイート
│   
│   # === 量子化学計算 (400-499) ===
│   ├── 401_Singlepointcalculation.py   # 一点エネルギー計算
│   ├── 402_Optimization.py             # 構造最適化
│   ├── 403_OPTandFreq.py               # 最適化と振動解析
│   
│   # === 可視化と解析 (500-599) ===
│   ├── 501_Visualization.py            # 計算結果の可視化
│   ├── 502_EDA.py                      # エネルギー分解解析
│   ├── 503_xyz_fragment_decomposition.py # XYZ分子分解
│   ├── 504_ConformationalEnergyDecomposition.py # コンフォマーEDA
│   
│   # === 物性計算 (600-699) ===
│   ├── 601_IP.py                       # イオン化ポテンシャル計算
│   ├── 602_Solvation.py                # 溶媒効果計算
│   ├── 603_BDE.py                      # 結合解離エネルギー計算
│   ├── 604_Pka.py                      # pKa計算
│   
│   # === スペクトル計算 (700-799) ===
│   ├── 701_IR.py                       # IR スペクトル予測
│   ├── 702_NMR.py                      # NMR スペクトル予測
│   ├── 703_Polar.py                    # 分極率計算
│   ├── 704_UV.py                       # UV-Vis スペクトル予測
│   
│   # === 遷移状態計算 (800-899) ===
│   ├── 801_NEB.py                      # 反応経路計算（NEB）
│   ├── 802_TS.py                       # 遷移状態探索
│   ├── 803_IRC.py                      # 内在反応座標計算
│   
│   # === システム・設定 (900-999) ===
│   ├── 901_DataBase.py                 # データベース管理
│   ├── 902_Summarization.py            # 結果要約
│   ├── 903_SystemSettings.py           # システム設定
│   ├── 904_PageSettings.py             # ページ設定管理
│   └── __init__.py                     # パッケージ初期化
│
├── site/                               # MkDocsビルド出力ディレクトリ
│
├── tests/                              # テストコード
│
├── utils/                              # 共通ユーティリティ
│   └── module.py                       # streamlit表示や共通処理のユーティリティ
│
├── .gitignore                          # Git管理対象外ファイル設定
├── .readthedocs.yaml                   # Read the Docs設定
├── ComputationalChemistryTool.py       # streamlitアプリのエントリポイント
├── main.py                             # コマンドライン・バッチ処理用エントリポイント
├── mkdocs.yml                          # MkDocs設定ファイル
├── packages.txt                        # システムパッケージ一覧
├── README.md                           # プロジェクト概要（英語）
├── README_jp.md                        # プロジェクト概要（日本語）
└── requirements.txt                    # Python依存パッケージ
```

## 主なモジュールと役割
RDKitやPySCFを用いて、プログラミング技術がなくても計算化学を行うために以下の構成で作成。
logic部分で基本的な計算処理、controllers部分で複雑な計算の制御、pages部分でユーザーインターフェースを提供。

## ページ体系と番号システム

アプリケーションのページは、機能別に3桁の番号体系で整理されています：

### 番号体系
- **000-099**: 🔧 前処理・ユーティリティ - 分子構造の確認・変換、入力ファイル形式チェック、フォーマット変換
- **100-199**: 🔀 統合・ワークフロー - 複数プログラムの統合実行、自動化ワークフロー、総合計算
- **200-299**: 🧲 分子力場 - 分子力場を用いた分子動力学計算、OpenMM、AMBER、CHARMM
- **300-399**: 🧬 半経験的 - 半経験的手法による高速計算、PM6/PM7、AM1、MNDO、DFTB、xTB
- **400-499**: 🧪 量子化学計算 - 分子構造最適化、一点エネルギー計算、配座解析
- **500-599**: 🔍 可視化と解析 - 分子軌道可視化、エネルギー分解解析、フラグメント解析
- **600-699**: ⚡ 物性計算 - イオン化ポテンシャル、溶媒効果、結合解離エネルギー、pKa計算
- **700-799**: 📊 スペクトル計算 - IR、NMR、UV-Visスペクトル予測、分極率計算
- **800-899**: 🔄 遷移状態計算 - 遷移状態探索、反応経路計算、IRC解析
- **900-999**: ⚙️ システム・設定 - 設定管理、データベース、結果集計

この番号体系により、機能の分類が明確になり、新機能の追加も体系的に行えます。

### controllers（高度計算制御）
#### controllers/energydecompositionanalysis.py
- エネルギー分解解析（EDA）の実行を制御
- 複数分子間の相互作用エネルギーを成分ごとに分析

#### controllers/fragment_calculation.py
- BDE（結合解離エネルギー）やpKaなど、フラグメントを用いた計算の制御

#### controllers/excited_state_calculation.py  
- TD-DFT計算による励起状態計算の制御
- UV/Visスペクトル予測

#### controllers/xyz_fragment_decomposition.py
- XYZ座標ファイルから2分子への自動分解処理
- 距離ベースでの分子分離機能


### logic（計算処理の核）
#### logic/calculation.py
- 量子化学計算の中心的な実装を担当
- 主な関数
    - `run_quantum_calculation` : 一点計算などのエントリポイント
    - `run_geometry_optimization` : 構造最適化
    - `calculate_vibrational_frequencies` : 振動数・熱力学量計算
    - `calc_nmr_and_shift` : NMRシールド値計算
    - `compute_electric_properties` : 双極子モーメント・分極率など
    - `calculate_ionization_potential` : VIP（イオン化ポテンシャル）計算
- 各関数はPySCFを用いて分子の計算を行い、必要に応じてファイル保存やログ出力も行う

#### logic/database.py
- 計算結果や分子情報をSQLデータベースに保存・取得するための関数群

#### logic/molecule_handler.py
- RDKitを用いて分子構造の変換や入力ファイルの作成を行う
- SMILES、XYZ座標の処理と、2次元・3次元表示、Z-matrix、PySCF入力形式への変換

#### logic/output_handler.py
- PySCFの計算結果から必要な情報を抽出・整形する

#### logic/visualization.py
- 計算結果を可視化し、streamlit上で表示するための図を作成
- エネルギー収束グラフ、振動スペクトル、EDA解析、UVスペクトル、NMRスペクトルなどの可視化

#### logic/logger.py
- 計算ログやエラー情報の記録・管理

#### logic/data_loader.py
- データの読み込みや前処理用のユーティリティ

#### logic/pubchem_api.py
- PubChem APIを使った文献情報や分子データへのアクセス

### config/（設定管理）
- 計算に使用する定数やパラメータを管理
- ページ表示設定、ユーザー設定、溶媒データ、pKa参考データなど

### utils/（共通ユーティリティ）
- streamlit表示や共通処理のためのユーティリティ関数

### pages/（ユーザーインターフェース）
- streamlitアプリの各ページ（UI部分）
- 3桁の番号体系で機能別に整理された計算ページ群
- 各カテゴリで特定の計算機能を提供（前処理、分子力場、半経験的、量子化学など）

### 補助的なユーティリティ
- ファイル名生成、ログ保存、分子情報の抽出などの補助関数も含まれる
- configディレクトリにユーザー設定やページ表示制御機能
- dataディレクトリにSQLiteデータベースと分子別の計算結果フォルダが自動生成される

## 主要な機能カテゴリ

### 前処理・ユーティリティ機能
- 分子構造の確認・変換、フォーマット変換
- 入力ファイル形式のチェック

### 統合・ワークフロー機能
- 複数プログラムの統合実行
- 自動化された計算ワークフロー

### 分子力場・半経験的計算
- OpenMM分子力場計算、RDKitコンフォマー探索
- MOPAC、xTB を用いた高速半経験的計算

### 量子化学計算機能
- 一点エネルギー計算、構造最適化、振動解析
- PySCFを用いた高精度ab initio計算

### 可視化と解析機能
- 分子軌道可視化、エネルギー分解解析（EDA）
- フラグメント解析、コンフォマー解析

### 物性・スペクトル計算機能  
- イオン化ポテンシャル、溶媒効果、結合解離エネルギー、pKa計算
- IR、NMR、UV-Visスペクトル予測、分極率計算

### 遷移状態・反応経路解析
- NEB反応経路計算、遷移状態探索
- IRC（内在反応座標）解析

### システム管理・データベース
- 計算結果のデータベース管理
- 結果の集計・要約、システム設定

---
## ドキュメント体系

プロジェクトでは複数のドキュメント形式をサポートしています：

### MkDocs（推奨）
- `mkdocs.yml`: MkDocs設定ファイル
- `docs/index.md`: メインドキュメント
- Material Designテーマによるモダンなドキュメントサイト

### Sphinx
- `docs/conf.py`: Sphinx設定ファイル  
- `docs/index.rst`: Sphinxメインドキュメント
- Python開発者向けの標準的なドキュメント形式

### CI/CD
- `.github/workflows/`: GitHub Actions設定
- `.readthedocs.yaml`: Read the Docs自動ビルド設定

## 今後の改善ポイント
- calculation.pyの機能分割・モジュール化の継続検討
- controllersディレクトリの機能充実（より複雑な計算制御への対応）
- ユーザー設定機能の拡充とページカスタマイゼーションの向上
- データベーススキーマの最適化とパフォーマンス改善
- テストカバレッジの向上と自動テストの充実
- ドキュメントの多言語対応とAPI文書の自動生成

## 参考
- 詳細なAPIや各関数の仕様は各ファイルのdocstringやコメントを参照
