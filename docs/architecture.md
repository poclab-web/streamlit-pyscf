# streamlit-pyscf アーキテクチャ概要

ここの部分は、開発者用に記載

## ディレクトリ構成（抜粋）

```
streamlit-pyscf/
├── config/                             # 計算に使用する定数を保存する場所
├── controllers/                        # 複数計算を制御する場所
│   └── excited_state_calculation.py    # BDEやpKaなどフラグメント計算用の関数
│   └── fragment_calculation.py         # 励起状態計算用の関数
├── data                                # 計算したデータが保存される場所
├── docs/
│   └── architecture.md                 # この設計書(ディレクトリ構成などを記載)
├── jobs/                               # 未計算ジョブ
│   └── pending_jobs.csv                # Excel等で編集可能なCSVで計算したい分子とパラメータの一覧が作成される
├── logic/
│   └── calculation.py                  # pyscfでの量子化学計算の主要なロジック
│   └── data_loader.py                  # データの読み込み用ユーティリティ
│   └── database.py                     # 計算結果などをSQLデータベースに保存する関数
│   └── logger.py                       # ログ出力・管理用のユーティリティ
│   └── molecular_handler.py            # RDKitを使った分子情報の変換・入力ファイル作成
│   └── output_handler.py               # pyscfの出力結果の抽出・整形
│   └── pubchem_api.py                  # pubchemのapiを使って、文献情報へアクセス
│   └── plot_handler.py                 # 結果をstreamlit上に表示するための図作成
│
├── pages/                              # streamlitアプリの各ページ Docstringに内容を記載
├── utils/                              # streamlit表示や共通処理のユーティリティ
├── .gitignore                          # Git管理対象外ファイル設定
├── ComputationalChemistryTool.py       # streamlitアプリのエントリポイント
├── main.py                             # コマンドラインやバッチ処理で分子計算を実行するためのエントリポイント
├── packages.txt                        # 使用パッケージ一覧
├── README.txt                          # プロジェクト概要
├── requirements.txt                    # Python依存パッケージ
```

## 主なモジュールと役割
RDKitやpyscfを用いて計算化学をプログラミング技術がなくても、行うために以下の構成で作成
logic部分は基本的な計算のためのもの、controllersの部分で計算部分をまとめて、pagesの部分でユーザーへ表示したり入力を受け取るようにしている。

### controllers
#### controllers/fragment_calculation.py
- BDE（結合解離エネルギー）やpKaなど、フラグメントを用いた計算の関数を提供


### logic
#### logic/calculation.py
- 量子化学計算の中心的な実装を担当
- 主な関数
    - `run_quantum_calculation` : 1点計算などのエントリポイント
    - `run_geometry_optimization` : 構造最適化
    - `calculate_vibrational_frequencies` : 振動数・熱力学量計算
    - `calc_nmr_and_shift` : NMRシールド値計算
    - `compute_electric_properties` : 双極子モーメント・分極率など
    - `calculate_ionization_potential` : VIP（イオン化ポテンシャル）計算

- 各関数はPySCFを用いて分子の計算を行い、必要に応じてファイル保存やログ出力も行う

#### logic/database.py
- 計算結果や分子情報をSQLデータベースに保存・取得するための関数群

#### logic/molecular_handler.py
- RDKitを用いて分子構造の変換や入力ファイルの作成を行う

#### logic/output_handler.py
- PySCFの計算結果から必要な情報を抽出・整形する

#### logic/plot_handler.py
- 計算結果を可視化し、streamlit上で表示するための図を作成

#### logic/logger.py
- 計算ログやエラー情報の記録・管理

#### logic/data_loadre.py
- データの読み込みや前処理用のユーティリティ

### config/
- 計算に使用する定数やパラメータを管理

### utils/
- streamlit表示や共通処理のためのユーティリティ関数

### pages/
- streamlitアプリの各ページ（UI部分）

### 補助的なユーティリティ
- ファイル名生成、ログ保存、分子情報の抽出などの補助関数も含まれる


---
## 今後の改善ポイント
- calculation.pyが肥大化しているため、機能ごとに分割・共通化を検討中
- 共通処理やユーティリティの整理、テスト容易性の向上も今後の課題

## 参考
- 詳細なAPIや各関数の仕様は各ファイルのdocstringやコメントを参照
