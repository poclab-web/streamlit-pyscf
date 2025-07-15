"""
セッティングのパラメーターのデフォルトを変更する。
"""

# 保存先


# エネルギーの変換


import streamlit as st
import platform
import psutil
import multiprocessing
import sys
import importlib
import os
import pandas as pd

# タイトル
st.title("PySCF および環境情報の表示")

# 環境情報
st.subheader("🔧 システム情報")
st.write(f"Python バージョン: {sys.version}")
st.write(f"プラットフォーム: {platform.platform()}")

# CPU情報
st.subheader("💡 CPU 情報")
st.write(f"物理コア数: {psutil.cpu_count(logical=False)}")
st.write(f"論理コア数: {psutil.cpu_count(logical=True)}")
st.write(f"最大周波数: {psutil.cpu_freq().max:.2f} MHz")
st.write(f"現在のCPU使用率: {psutil.cpu_percent()} %")

# 並列処理のサポート
st.subheader("🔄 並列計算のサポート")
cpu_count = multiprocessing.cpu_count()
st.write(f"multiprocessing による並列数: {cpu_count}")

# メモリ情報
st.subheader("🧠 メモリ情報")
virtual_mem = psutil.virtual_memory()
st.write(f"総メモリ: {virtual_mem.total / (1024 ** 3):.2f} GB")
st.write(f"使用中メモリ: {virtual_mem.used / (1024 ** 3):.2f} GB")
st.write(f"空きメモリ: {virtual_mem.available / (1024 ** 3):.2f} GB")
st.write(f"メモリ使用率: {virtual_mem.percent} %")

# OpenMPサポートの確認（PySCFなどのC backendに重要）
st.subheader("🧵 OpenMP サポート確認")
omp_threads = os.environ.get("OMP_NUM_THREADS", "未設定（自動）")
st.write(f"OMP_NUM_THREADS: {omp_threads}")

# チェックする主なモジュール
modules = [
    "pyscf",
    "streamlit",
    "numpy",
    "matplotlib",
    "rdkit",  # PySCFの内部Cライブラリ
]

st.subheader("📦 モジュールのバージョン")
for mod in modules:
    try:
        m = importlib.import_module(mod)
        st.write(f"{mod}: {m.__version__}")
    except ImportError:
        st.warning(f"{mod} はインポートできませんでした。")
    except AttributeError:
        st.write(f"{mod}: バージョン情報が取得できませんでした。")

# すべてのインストール済みパッケージを表示（必要であれば）
if st.checkbox("全パッケージ一覧も表示"):
    import pkg_resources
    all_packages = sorted([(dist.project_name, dist.version) for dist in pkg_resources.working_set])
    for name, version in all_packages:
        st.write(f"{name}: {version}")

st.subheader("🧪 溶媒誘電率データ（solvents_epsion.csv）")

# CSVファイルのパス
csv_path =  "config/solvents_epsilon.csv"

# CSV読み込み
if os.path.exists(csv_path):
    df = pd.read_csv(csv_path)
    st.dataframe(df)
else:
    st.warning(f"{csv_path} が見つかりません。")

# 新規追加フォーム
with st.expander("溶媒データを追加"):
    new_name = st.text_input("溶媒名")
    new_eps = st.number_input("誘電率 (ε)", min_value=0.0, step=0.01, format="%.2f")
    if st.button("追加"):
        if new_name and new_eps > 0:
            # 既存データに追加
            new_row = pd.DataFrame([[new_name, new_eps]], columns=df.columns)
            df = pd.concat([df, new_row], ignore_index=True)
            df.to_csv(csv_path, index=False)
            st.success(f"{new_name} (ε={new_eps}) を追加しました。")
            st.rerun()
        else:
            st.error("溶媒名と誘電率を正しく入力してください。")

st.subheader("📂 pka list ")

# pka listのCSVファイルのパス
pka_csv_path = "config/pKa_reference_list.csv"

# pka listのCSV読み込み
if os.path.exists(pka_csv_path):
    pka_df = pd.read_csv(pka_csv_path)
    st.dataframe(pka_df)
else:
    st.warning(f"{pka_csv_path} が見つかりません。")    
