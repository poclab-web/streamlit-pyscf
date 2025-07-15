import pandas as pd


# 量子化学計算の理論手法と基底関数の設定値を定義
theory_info = {
    "HF": {
        "group": "Hartree-Fock",
        "desc": "電子相関なし、基準計算に用いられる波動関数法"
    },
    "MP2": {
        "group": "Post-HF",
        "desc": "2次摂動論による電子相関補正",
    },
    "B3LYP": {
        "group": "DFT",
        "desc": "ハイブリッド汎用汎関数",
        "xc": "b3lyp"
    },
    "B3LYPD3": {
        "group": "DFT + Dispersion",
        "desc": "B3LYPにD3分散補正を加えたモデル",
        "xc": "b3lyp"
    },
    "PBE": {
        "group": "DFT",
        "desc": "GGA汎関数、周期系や高速計算向き",
        "xc": "pbe"
    },
    "M06-2X": {
        "group": "DFT",
        "desc": "中長距離相互作用に強いmeta-GGA",
        "xc": "m06-2x"
    },
    "B97X-D": {
        "group": "DFT",
        "desc": "分散補正込みハイブリッド汎関数",
        "xc": "wb97x-d"  # PySCFでは "w" を含む綴り
    },
    "CAMB3LYP": {
        "group": "Advanced (Experimental)",
        "desc": "長距離補正付きB3LYP",
        "xc": "camb3lyp"
    },
    "PBE0": {
        "group": "Advanced",
        "desc": "PBEにHF交換を混ぜたハイブリッドGGA",
        "xc": "pbe0"
    },
    "LC-ωPBE": {
        "group": "Advanced (Experimental)",
        "desc": "長距離補正型PBE",
        "xc": "lc-wpbe"
    },
    "ωB97X-D": {
        "group": "Advanced",
        "desc": "レンジ分離＋分散補正型ハイブリッド",
        "xc": "wb97x-d"
    },
    "TPSS": {
        "group": "Advanced",
        "desc": "中庸なmeta-GGA",
        "xc": "tpss"
    },
    "SCAN": {
        "group": "Advanced",
        "desc": "近年注目される高精度meta-GGA",
        "xc": "scan"
    }
}

# 量子化学計算の基底関数の設定値を定義
basis_set_info = {
    "Pople系列": [
    {"name": "sto-3g",        "desc": "最小基底（STO-3G）：教育・初学者向けの簡易モデル"},
    {"name": "3-21g",         "desc": "3-21G：内側3つ＋外側2+1の縮約、軽量な分極なし基底"},
    {"name": "6-31g",         "desc": "6-31G：二重ζ基底、価電子の柔軟な記述に対応"},
    {"name": "6-31g*",        "desc": "6-31G(d)：重原子にd型分極関数を付加"},
    {"name": "6-31g**",       "desc": "6-31G(d,p)：全原子に分極関数（Hにもp軌道）を付加"},
    {"name": "6-31+g**",      "desc": "6-31+G(d,p)：陰イオン・励起状態向けの拡散関数つき"},
    {"name": "6-31++g**",     "desc": "6-31++G(d,p)：すべての原子に拡散関数を付加（H含む）"},
    ],
    "Dunning系列": [
        {"name": "cc-pVDZ",       "desc": "二重ζ（DFT・MP2向け）"},
        {"name": "cc-pVTZ",       "desc": "三重ζ（より精度重視）"},
        {"name": "cc-pVQZ",       "desc": "四重ζ（高精度、基準計算向け）"},
        {"name": "aug-cc-pVDZ",   "desc": "+拡散関数（アニオン・励起状態）"},
        {"name": "aug-cc-pVTZ",   "desc": "+拡散関数（精度高）"},
        {"name": "aug-cc-pVQZ",   "desc": "+拡散付き四重ζ（Rydberg状態、TDDFT）"},
    ],
    "def2系列（Karlsruhe）": [
        {"name": "def2-SVP",      "desc": "軽量＋分極（実用的な初期計算）"},
        {"name": "def2-SVPD",     "desc": "+拡散関数（アニオン向け）"},
        {"name": "def2-SV(P)",    "desc": "軽量版（H, C, N, O 向け）"},
        {"name": "def2-TZVP",     "desc": "triple-ζ 分極付き（標準精度）"},
        {"name": "def2-TZVPD",    "desc": "+拡散関数（高精度アニオン向け）"},
        {"name": "def2-QZVP",     "desc": "quadruple-ζ（基準計算向け）"},
    ],
    "Jensen系列（pcseg）": [
        {"name": "pcseg-1",       "desc": "分極一貫性あり（軽量）"},
        {"name": "pcseg-2",       "desc": "中等精度（triple-ζ相当）"},
        {"name": "pcseg-3",       "desc": "triple-ζ精度（高精度計算用）"},
    ],
    "その他": [
        {"name": "ma-def2-SVP",   "desc": "最小限拡張付きdef2-SVP（高速化重視）"},
        {"name": "cc-pwCVTZ",     "desc": "コアバレンス相関考慮（遷移金属・重原子）"},
    ]
}

# 量子化学計算の構造最適化のときの設定値を定義
conv_preset_values = {
    "Loose":  {"energy": 1.0e-4, "grms": 1.0e-3, "gmax": 3.0e-3, "drms": 4.0e-3, "dmax": 6.0e-3},
    "Normal": {"energy": 1.0e-5, "grms": 5.0e-4, "gmax": 1.5e-3, "drms": 2.0e-3, "dmax": 3.0e-3},
    "Tight":  {"energy": 1.0e-6, "grms": 3.0e-4, "gmax": 1.2e-3, "drms": 1.2e-3, "dmax": 1.8e-3},
}

# solvent_models
solvent_models = ["None", "PCM", "ddCOSMO"]
solvents_file = "config/solvents_epsilon.csv"
solvents_data = pd.read_csv(solvents_file)


# データベースのテーブル定義
# 各列の名前、型、説明を定義
columns_info = [
    ("id", "INTEGER PRIMARY KEY AUTOINCREMENT", "自動採番ID"),
    ("inchi", "TEXT NOT NULL", "InChI表記"),
    ("inchikey", "TEXT NOT NULL", "InChIKey"),
    ("mol", "TEXT", "MolBlock形式"),
    ("mw", "REAL", "分子量"),
    ("formula", "TEXT", "化学式"),
    ("charge", "INTEGER NOT NULL", "電荷"),
    ("spin", "INTEGER NOT NULL", "スピン多重度"),
    ("g_tot", "REAL", "全エネルギー"),
    ("zpe", "REAL", "ゼロ点振動エネルギー"),
    ("method", "TEXT NOT NULL", "計算手法"),
    ("basis", "TEXT NOT NULL", "基底関数"),
    ("solvent", "TEXT", "溶媒モデル"),
    ("dielectric", "REAL", "誘電率"),
    ("temperature", "REAL", "温度"),
    ("pressure", "REAL", "圧力"),
    ("frequencies", "TEXT", "振動数リスト（JSON）"),
    ("chk_file", "BLOB", "チェックポイントファイル"),
    ("num_imaginary", "INTEGER", "虚数振動数の数"),
    ("nstates", "INTEGER", "励起状態数"),
    ("excited_spin", "TEXT", "励起状態のスピン（singlet/triplet）"),
    ("tda", "INTEGER", "TDA法の有無（0/1）"),
    ("excited_energies", "TEXT", "励起状態エネルギーリスト（JSON）"),  
    ("oscillator_strengths", "TEXT", "振動子強度リスト（JSON）"),      
    ("timestamp", "TEXT DEFAULT CURRENT_TIMESTAMP", "登録日時"),
]