"""
Controllersの部分にある計算を連続実行するためのメインプログラムです。
jobs/pending_jobs.csvに記載された分子の計算を自動で行います。
このプログラムは、コマンドラインやバッチ処理での実行を想定しています。
計算ジョブはCSV形式で記載され、各ジョブは独立したプロセスで計算されます。
PCのCPUコア数に応じて並列数を自動で設定し、効率的に計算を行います。
計算結果は標準出力に表示されますが、必要に応じてファイルやデータベースに保存することも可能です。    
"""

import os
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
from logic.calculation import compute_molecule_properties

# from logic.data_loader import ensure_all_columns_exist

def load_pending_jobs(filepath="jobs/pending_jobs.csv"):
    """CSVから未計算ジョブを読み込む（1行1ジョブ）"""
    with open(filepath, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        return list(reader)

def parse_bool(val):
    if isinstance(val, bool):
        return val
    if isinstance(val, str):
        return val.lower() in ("true", "1", "yes")
    return False

def get_max_workers():
    """PCのCPUコア数から並列数を自動決定（コア数の半分、最大8）"""
    cpu_count = os.cpu_count() or 1
    return min(max(1, cpu_count // 2), 8)

def run_job(job):
    """1ジョブ分の計算を実行"""
    name = job.get("name")
    smiles = job.get("smiles")
    xyz = job.get("xyz")
    charge = int(job.get("charge", 0))
    spin = int(job.get("spin", 0))
    # conv_paramsは辞書形式で記入（例: {"max_cycle": 100}）
    conv_params = eval(job.get("conv_params", "{}"))
    theory = job.get("theory", "B3LYP")
    basis_set = job.get("basis_set", "def2-SVP")
    opt_theory = job.get("opt_theory") or None
    opt_basis_set = job.get("opt_basis_set") or None
    solvent_model = job.get("solvent_model") or None
    eps = float(job["eps"]) if job.get("eps") else None
    maxsteps = int(job.get("maxsteps", 100))
    optimize_with_qc = parse_bool(job.get("optimize_with_qc", "True"))

    print(f"=== {name} の計算を開始 ===")
    result = compute_molecule_properties(
        name, smiles, xyz, charge, spin, conv_params,
        theory, basis_set, opt_theory, opt_basis_set,
        solvent_model, eps, maxsteps, optimize_with_qc
    )
    print(f"=== {name} の計算が完了 ===")
    print(result)
    return name, result

def main():
    jobs = load_pending_jobs()
    max_workers = get_max_workers()
    print(f"並列計算数: {max_workers}")
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(run_job, job) for job in jobs]
        for future in as_completed(futures):
            name, result = future.result()
            print(f"[{name}] 完了")

import sqlite3
from config.config import columns_info

def ensure_all_columns_exist(db_path="data/energy_db.sqlite"):
    """データベースのmoleculesテーブルに必要なカラムが存在するか確認し、存在しない場合は追加する。"""

    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute("PRAGMA table_info(molecules);")
    existing_columns = [row[1] for row in cur.fetchall()]

    added = False
    for col_name, col_type, _ in columns_info:
        if col_name not in existing_columns:
            # timestampのDEFAULTはALTER TABLEで追加できないため、型のみ追加
            if col_name == "timestamp":
                cur.execute(f"ALTER TABLE molecules ADD COLUMN {col_name} TEXT;")
            else:
                # AUTOINCREMENTやNOT NULL制約はALTER TABLEでは追加できないので型のみ
                base_type = col_type.split()[0]
                cur.execute(f"ALTER TABLE molecules ADD COLUMN {col_name} {base_type};")
            print(f"✅ カラム {col_name} を追加しました。")
            added = True

    if not added:
        print("✅ すべてのカラムが既に存在しています。追加はありません。")

    conn.commit()
    conn.close()

if __name__ == "__main__":
    print("=== 計算ジョブの実行を開始 ===")
    ensure_all_columns_exist()
    # main()