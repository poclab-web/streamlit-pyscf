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
