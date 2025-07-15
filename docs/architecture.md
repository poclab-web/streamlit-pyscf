# streamlit-pyscf アーキテクチャ概要

ここの部分は、開発者用に記載

## ディレクトリ構成

```
streamlit-pyscf/
├── config/                             # 計算設定と定数を管理
│   ├── config.py                       # 計算パラメータの設定
│   ├── page_visibility.json            # ページ表示設定
│   ├── user_page_visibility.json       # ユーザー別ページ表示設定
│   ├── user_preferences.py             # ユーザー設定管理
│   ├── pKa_reference_list.csv           # pKa計算の参考データ
│   ├── solvents_epsilon.csv             # 溶媒の誘電率データ
│   └── styles.css                       # streamlitアプリのスタイル設定
│
├── controllers/                        # 高度な計算処理を制御
│   ├── energydecompositionanalysis.py  # エネルギー分解解析の制御
│   ├── excited_state_calculation.py    # 励起状態計算の制御
│   ├── fragment_calculation.py         # フラグメント計算（BDE、pKaなど）の制御
│   └── xyz_fragment_decomposition.py   # XYZ座標からの分子分解処理
│
├── data/                               # 計算結果データが保存される場所
│   ├── energy_db.sqlite                # エネルギー計算結果のデータベース
│   └── [分子別フォルダ]/                # 各分子の計算データフォルダ
│
├── docs/
│   └── architecture.md                 # この設計書（アーキテクチャ概要）
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
├── pages/                              # streamlitアプリの各ページ
│   ├── 01_GeneralCalculation.py        # 一般的な量子化学計算
│   ├── 02_StructureCheck.py            # 分子構造の確認
│   ├── 03_Singlepointcalculation.py    # 一点エネルギー計算
│   ├── 04_Optimization.py              # 構造最適化
│   ├── 05_ConformationalSearch.py      # コンフォマー探索
│   ├── 06_OPTandFreq.py                # 最適化と振動解析
│   ├── 07_Visualization.py             # 計算結果の可視化
│   ├── 08_PropertyCalculation*.py      # 物性計算（IR、NMR、分極率）
│   ├── 09_IonizationPotencial.py       # イオン化ポテンシャル計算
│   ├── 09_Solvation.py                 # 溶媒効果計算
│   ├── 10_BondDissociationEnergy.py    # 結合解離エネルギー計算
│   ├── 11_Pka.py                       # pKa計算
│   ├── 12_UV_SpectrumPrediction.py     # UV吸収スペクトル予測
│   ├── 13_*EnergyDecomposition*.py     # エネルギー分解解析関連
│   ├── 14_*.py                         # 反応経路計算（NEB、遷移状態）
│   ├── 15_IRC.py                       # 内在反応座標計算
│   ├── 16_*.py                         # データベース管理・結果要約
│   ├── 99_SystemSettings.py            # システム設定
│   └── PageSettings.py                 # ページ設定管理
│
├── utils/
│   └── module.py                       # streamlit表示や共通処理のユーティリティ
│
├── .gitignore                          # Git管理対象外ファイル設定
├── ComputationalChemistryTool.py       # streamlitアプリのエントリポイント
├── main.py                             # コマンドライン・バッチ処理用エントリポイント
├── packages.txt                        # システムパッケージ一覧
├── README.md                           # プロジェクト概要
└── requirements.txt                    # Python依存パッケージ
```

## 主なモジュールと役割
RDKitやPySCFを用いて、プログラミング技術がなくても計算化学を行うために以下の構成で作成。
logic部分で基本的な計算処理、controllers部分で複雑な計算の制御、pages部分でユーザーインターフェースを提供。

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
- 01-16番台で機能別に整理された計算ページ群
- 各ページで特定の計算機能を提供

### 補助的なユーティリティ
- ファイル名生成、ログ保存、分子情報の抽出などの補助関数も含まれる
- configディレクトリにユーザー設定やページ表示制御機能
- dataディレクトリにSQLiteデータベースと分子別の計算結果フォルダが自動生成される

## 主要な機能カテゴリ

### 基本計算機能
- 一点エネルギー計算、構造最適化、振動解析
- コンフォマー探索、溶媒効果計算

### 物性計算機能  
- NMRシールド値、IR周波数、分極率計算
- イオン化ポテンシャル、UV/Visスペクトル予測

### 高度な解析機能
- エネルギー分解解析（EDA）
- 結合解離エネルギー（BDE）、pKa計算
- 反応経路解析（NEB、遷移状態、IRC）

### データ管理・可視化
- 計算結果のデータベース管理
- グラフ・スペクトルの可視化
- 分子構造の2D/3D表示

---
## 今後の改善ポイント
- calculation.pyの機能分割・モジュール化の継続検討
- controllersディレクトリの機能充実（より複雑な計算制御への対応）
- ユーザー設定機能の拡充とページカスタマイゼーションの向上
- データベーススキーマの最適化とパフォーマンス改善

## 参考
- 詳細なAPIや各関数の仕様は各ファイルのdocstringやコメントを参照
