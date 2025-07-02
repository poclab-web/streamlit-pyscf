import sqlite3
import json
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import inchi
import csv

from config.config import columns_info

# データベースのカラムが存在するか確認し、必要なカラムを追加
# 既存のカラム情報を取得し、必要なカラムが存在しない場合は追加する

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

# データベース接続と統計取得
def get_summary_statistics(db_path="energy_db.sqlite"):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    cur.execute("SELECT COUNT(*) FROM molecules")
    total = cur.fetchone()[0]

    cur.execute("SELECT DISTINCT method FROM molecules")
    methods = [row[0] for row in cur.fetchall()]

    cur.execute("SELECT DISTINCT basis FROM molecules")
    bases = [row[0] for row in cur.fetchall()]

    cur.execute("SELECT DISTINCT solvent FROM molecules WHERE solvent IS NOT NULL")
    solvents = [row[0] for row in cur.fetchall()]

    conn.close()
    return total, methods, bases, solvents

# データベースへのデータの登録部分
def insert_molecule_with_frequencies(
    inchi,
    inchikey,
    g_tot=None,
    zpe=None,
    method="B3LYP",
    basis="def2-SVP",
    charge=0,
    spin=0,
    solvent=None,
    dielectric=None,
    temperature=298.15,
    pressure=1.0,
    frequencies=None,
    num_imaginary=None,
    chk_file_path=None,
    nstates=None,
    excited_spin=None,
    tda=None,
    excited_energies=None,         
    oscillator_strengths=None,     
    db_path="data/energy_db.sqlite"
):
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    import sqlite3
    import json
    from datetime import datetime
    import os

    # inchi, inchikey からmolを生成
    mol = Chem.MolFromInchi(inchi)
    molblock = Chem.MolToMolBlock(mol)
    mw = Descriptors.MolWt(mol)
    formula = rdMolDescriptors.CalcMolFormula(mol)

    # .chkファイルのバイナリデータを読み込み
    chk_file_data = None
    if chk_file_path and os.path.isfile(chk_file_path):
        with open(chk_file_path, "rb") as f:
            chk_file_data = f.read()

    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    # テーブル作成（存在しなければ）
    cur.execute("""
    CREATE TABLE IF NOT EXISTS molecules (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        inchi TEXT NOT NULL,
        inchikey TEXT NOT NULL,
        mol TEXT,
        mw REAL,
        formula TEXT,
        charge INTEGER NOT NULL,
        spin INTEGER NOT NULL,
        g_tot REAL,
        zpe REAL,
        method TEXT NOT NULL,
        basis TEXT NOT NULL,
        solvent TEXT,
        dielectric REAL,
        temperature REAL,
        pressure REAL,
        frequencies TEXT,
        chk_file BLOB,
        num_imaginary INTEGER,
        nstates INTEGER,                
        excited_spin TEXT,              
        tda INTEGER,                    
        timestamp TEXT DEFAULT CURRENT_TIMESTAMP
    );
    """)

    freq_json = json.dumps(frequencies) if frequencies is not None else None
    excited_energies_json = json.dumps(excited_energies) if excited_energies is not None else None
    oscillator_strengths_json = json.dumps(oscillator_strengths) if oscillator_strengths is not None else None

    # データ挿入
    cur.execute("""
    INSERT INTO molecules (
        inchi, inchikey, mol, mw, formula, charge, spin, g_tot, zpe,
        method, basis, solvent, dielectric, temperature, pressure,
        frequencies, chk_file, num_imaginary,
        nstates, excited_spin, tda,
        excited_energies, oscillator_strengths,
        timestamp
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, (
        inchi, inchikey, molblock, mw, formula, charge, spin, g_tot, zpe,
        method, basis, solvent, dielectric, temperature, pressure,
        freq_json, chk_file_data, num_imaginary,
        nstates, excited_spin, int(tda) if tda is not None else None,
        excited_energies_json, oscillator_strengths_json,   # 追加
        datetime.now().isoformat()
    ))

    conn.commit()
    conn.close()
    print(f"✅ Inserted molecule: {inchikey} ")


# データベースから分子データを取得する関数
def get_molecule_from_sqlite(
    inchikey, method, basis, spin, charge,
    solvent=None, dielectric=None,
    temperature=298.15, pressure=1.0,
    nstates=None, excited_spin=None, tda=None,
    db_path="energy_db.sqlite"
):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    # テーブルがなければ作成
    cur.execute("""
    CREATE TABLE IF NOT EXISTS molecules (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        inchi TEXT NOT NULL,
        inchikey TEXT NOT NULL,
        mol TEXT,
        mw REAL,
        formula TEXT,
        charge INTEGER NOT NULL,
        spin INTEGER NOT NULL,
        g_tot REAL,
        zpe REAL,
        method TEXT NOT NULL,
        basis TEXT NOT NULL,
        solvent TEXT,
        dielectric REAL,
        temperature REAL,
        pressure REAL,
        frequencies TEXT,
        chk_file BLOB,
        num_imaginary INTEGER,
        nstates INTEGER,
        excited_spin TEXT,
        tda INTEGER,
        timestamp TEXT DEFAULT CURRENT_TIMESTAMP
    );
    """)

    query = """
    SELECT * FROM molecules
    WHERE inchikey = ? AND method = ? AND basis = ?
      AND spin = ? AND charge = ? AND temperature = ? AND pressure = ?
    """
    params = [inchikey, method, basis, spin, charge, temperature, pressure]

    if solvent is None:
        query += " AND solvent IS NULL"
    else:
        query += " AND solvent = ?"
        params.append(solvent)

    if dielectric is None:
        query += " AND dielectric IS NULL"
    else:
        query += " AND dielectric = ?"
        params.append(dielectric)

    # 追加: nstates, excited_spin, tda
    if nstates is not None:
        query += " AND nstates = ?"
        params.append(nstates)
    if excited_spin is not None:
        query += " AND excited_spin = ?"
        params.append(excited_spin)
    if tda is not None:
        query += " AND tda = ?"
        params.append(int(tda))

    cur.execute(query, params)
    row = cur.fetchone()
    conn.close()
    if row:
        (molecule_id, g_tot, zpe, freq_json, chk_file, num_imaginary,
         excited_energies_json, oscillator_strengths_json) = row
        frequencies = json.loads(freq_json) if freq_json else None
        excited_energies = json.loads(excited_energies_json) if excited_energies_json else None
        oscillator_strengths = json.loads(oscillator_strengths_json) if oscillator_strengths_json else None
        return {
            "id": molecule_id,
            "g_tot": g_tot,
            "zpe": zpe,
            "frequencies": frequencies,
            "chk_file": chk_file,
            "num_imaginary": num_imaginary,
            "excited_energies": excited_energies,
            "oscillator_strengths": oscillator_strengths
        }
    else:
        return None

# データベースから安定な分子を取得する関数
def get_stable_molecule_from_sqlite(
    inchikey,
    method, basis, spin, charge,
    solvent=None, dielectric=None,
    temperature=298.15, pressure=1.0,
    nstates=None, excited_spin=None, tda=None,
    db_path="energy_db.sqlite"
):
    import sqlite3
    import json

    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    # テーブルがなければ作成
    cur.execute("""
    CREATE TABLE IF NOT EXISTS molecules (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        inchi TEXT NOT NULL,
        inchikey TEXT NOT NULL,
        mol TEXT,
        mw REAL,
        formula TEXT,
        charge INTEGER NOT NULL,
        spin INTEGER NOT NULL,
        g_tot REAL,
        zpe REAL,
        method TEXT NOT NULL,
        basis TEXT NOT NULL,
        solvent TEXT,
        dielectric REAL,
        temperature REAL,
        pressure REAL,
        frequencies TEXT,
        chk_file BLOB,
        num_imaginary INTEGER,
        nstates INTEGER,
        excited_spin TEXT,
        tda INTEGER,
        timestamp TEXT DEFAULT CURRENT_TIMESTAMP
    );
    """)

    query = """
    SELECT id, g_tot, zpe, frequencies, chk_file, num_imaginary
    FROM molecules
    WHERE inchikey = ? AND method = ? AND basis = ?
      AND spin = ? AND charge = ? AND temperature = ? AND pressure = ?
      AND num_imaginary = 0
    """
    params = [inchikey, method, basis, spin, charge, temperature, pressure]

    if solvent is None:
        query += " AND solvent IS NULL"
    else:
        query += " AND solvent = ?"
        params.append(solvent)

    if dielectric is None:
        query += " AND dielectric IS NULL"
    else:
        query += " AND dielectric = ?"
        params.append(dielectric)

    if nstates is None:
        query += " AND nstates IS NULL"
    else:
        query += " AND nstates = ?"
        params.append(nstates)

    if excited_spin is None:
        query += " AND excited_spin IS NULL"
    else:
        query += " AND excited_spin = ?"
        params.append(excited_spin)

    if tda is None:
        query += " AND tda IS NULL"
    else:
        query += " AND tda = ?"
        params.append(int(tda))

    query += " ORDER BY g_tot ASC LIMIT 1"

    cur.execute(query, params)
    row = cur.fetchone()
    conn.close()

    if row:
        molecule_id, g_tot, zpe, freq_json, chk_file, num_imaginary = row
        frequencies = json.loads(freq_json) if freq_json else None
        return {
            "id": molecule_id,
            "g_tot": g_tot,
            "zpe": zpe,
            "frequencies": frequencies,
            "chk_file": chk_file,
            "num_imaginary": num_imaginary
        }
    else:
        return None


# データベースから最新の分子データを取得する関数
def get_latest_molecules(n=5, db_path="data/energy_db.sqlite"):
    """
    データベースから最新n件の分子データを取得
    """
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()
    cur.execute(
        "SELECT * FROM molecules ORDER BY id DESC LIMIT ?", (n,)
    )
    rows = cur.fetchall()
    conn.close()
    return [dict(row) for row in rows]

# データベースから分子データをJSON形式でエクスポートする関数
def export_molecule_data_to_json(export_path="exported_molecules.json", db_path="energy_db.sqlite"):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    # 全レコードを取得（条件を絞ってもOK）
    cur.execute("SELECT * FROM molecules")
    rows = cur.fetchall()

    # カラム名を取得
    colnames = [desc[0] for desc in cur.description]

    # レコードを辞書に変換
    data = [dict(zip(colnames, row)) for row in rows]

    # バイナリや日時などを扱いやすくする
    for entry in data:
        if isinstance(entry["chk_file"], bytes):
            entry["chk_file"] = entry["chk_file"].hex()  # 16進文字列に変換
        if isinstance(entry["timestamp"], str):
            entry["timestamp"] = entry["timestamp"]

    with open(export_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

    conn.close()
    print(f"✅ Exported {len(data)} molecules to {export_path}")


# データベースから分子データをCSV形式でエクスポートする関数
def export_molecule_data_to_csv(export_path="exported_molecules.csv", db_path="energy_db.sqlite"):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute("SELECT * FROM molecules")
    rows = cur.fetchall()
    headers = [desc[0] for desc in cur.description]

    with open(export_path, mode="w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        for row in rows:
            # chk_fileなどバイナリは省略またはhex変換
            row = [col.hex() if isinstance(col, bytes) else col for col in row]
            writer.writerow(row)

    conn.close()
    print(f"✅ Exported {len(rows)} rows to {export_path}")

# データベースから特定の計算方法の分子データをJSON形式でエクスポートする関数
def export_filtered_data(method="B3LYP", export_path="filtered_b3lyp.json", db_path="energy_db.sqlite"):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute("SELECT * FROM molecules WHERE method = ?", (method,))
    rows = cur.fetchall()
    headers = [desc[0] for desc in cur.description]

    data = [dict(zip(headers, row)) for row in rows]
    for entry in data:
        if isinstance(entry["chk_file"], bytes):
            entry["chk_file"] = entry["chk_file"].hex()

    with open(export_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

    conn.close()
    print(f"✅ Exported {len(data)} entries with method={method} to {export_path}")

# データベースから分子データをJSON形式でインポートする関数
def import_molecules_from_json(json_path, db_path="energy_db.sqlite"):
    import sqlite3
    import json

    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    # テーブルが存在しない場合は作成
    cur.execute("""
    CREATE TABLE IF NOT EXISTS molecules (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        inchi TEXT NOT NULL,
        inchikey TEXT NOT NULL,
        mol TEXT,
        mw REAL,
        formula TEXT,
        charge INTEGER NOT NULL,
        spin INTEGER NOT NULL,
        g_tot REAL,
        zpe REAL,
        method TEXT NOT NULL,
        basis TEXT NOT NULL,
        solvent TEXT,
        dielectric REAL,
        temperature REAL,
        pressure REAL,
        frequencies TEXT,
        chk_file BLOB,
        num_imaginary INTEGER,
        nstates INTEGER,                -- 追加
        excited_spin TEXT,              -- 追加 (singlet/triplet)
        tda INTEGER,                    -- 追加 (0/1)
        timestamp TEXT DEFAULT CURRENT_TIMESTAMP
    );
    """)

    for entry in data:
        # hex文字列だったバイナリを復元
        chk_file_data = bytes.fromhex(entry["chk_file"]) if entry["chk_file"] else None

        cur.execute("""
        INSERT INTO molecules (
            inchi, inchikey, mol, mw, formula, charge, spin,
            g_tot, zpe, method, basis, solvent, dielectric,
            temperature, pressure, frequencies, chk_file,
            num_imaginary, nstates, excited_spin, tda,
            excited_energies, oscillator_strengths,   -- 追加
            timestamp
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            entry["inchi"], entry["inchikey"], entry["mol"], entry["mw"], entry["formula"],
            entry["charge"], entry["spin"], entry["g_tot"], entry["zpe"],
            entry["method"], entry["basis"], entry["solvent"], entry["dielectric"],
            entry["temperature"], entry["pressure"],
            json.dumps(entry["frequencies"]) if entry.get("frequencies") else None,
            chk_file_data, entry.get("num_imaginary"),
            entry.get("nstates"), entry.get("excited_spin"), entry.get("tda"),
            json.dumps(entry.get("excited_energies")) if entry.get("excited_energies") else None,
            json.dumps(entry.get("oscillator_strengths")) if entry.get("oscillator_strengths") else None,
            entry["timestamp"]
        ))

    conn.commit()
    conn.close()
    print(f"✅ Imported {len(data)} entries from {json_path}")

# データベースから分子データをCSV形式でインポートする関数
def import_molecules_from_csv(csv_path, db_path="energy_db.sqlite"):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    cur.execute("""
    CREATE TABLE IF NOT EXISTS molecules (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        inchi TEXT NOT NULL,
        inchikey TEXT NOT NULL,
        mol TEXT,
        mw REAL,
        formula TEXT,
        charge INTEGER NOT NULL,
        spin INTEGER NOT NULL,
        g_tot REAL,
        zpe REAL,
        method TEXT NOT NULL,
        basis TEXT NOT NULL,
        solvent TEXT,
        dielectric REAL,
        temperature REAL,
        pressure REAL,
        frequencies TEXT,
        chk_file BLOB,
        num_imaginary INTEGER,
        nstates INTEGER,
        excited_spin TEXT,
        tda INTEGER,
        excited_energies TEXT,
        oscillator_strengths TEXT,
        timestamp TEXT DEFAULT CURRENT_TIMESTAMP
    );
    """)

    with open(csv_path, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        entries = list(reader)

    for entry in entries:
        chk_file_data = bytes.fromhex(entry["chk_file"]) if entry["chk_file"] else None
        frequencies = entry["frequencies"]
        if frequencies:
            try:
                frequencies = json.dumps(json.loads(frequencies))
            except Exception:
                frequencies = None

        excited_energies = entry.get("excited_energies")
        oscillator_strengths = entry.get("oscillator_strengths")
        if excited_energies:
            try:
                excited_energies = json.dumps(json.loads(excited_energies))
            except Exception:
                excited_energies = None
        if oscillator_strengths:
            try:
                oscillator_strengths = json.dumps(json.loads(oscillator_strengths))
            except Exception:
                oscillator_strengths = None

        cur.execute("""
        INSERT INTO molecules (
            inchi, inchikey, mol, mw, formula, charge, spin,
            g_tot, zpe, method, basis, solvent, dielectric,
            temperature, pressure, frequencies, chk_file,
            num_imaginary, nstates, excited_spin, tda,
            excited_energies, oscillator_strengths, timestamp
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            entry["inchi"], entry["inchikey"], entry["mol"], float(entry["mw"]), entry["formula"],
            int(entry["charge"]), int(entry["spin"]), float_or_none(entry["g_tot"]),
            float_or_none(entry["zpe"]), entry["method"], entry["basis"], entry["solvent"],
            float_or_none(entry["dielectric"]), float(entry["temperature"]),
            float(entry["pressure"]), frequencies, chk_file_data,
            int(entry["num_imaginary"]) if entry["num_imaginary"] else None,
            int(entry["nstates"]) if entry.get("nstates") else None,
            entry.get("excited_spin"),
            int(entry["tda"]) if entry.get("tda") else None,
            excited_energies,
            oscillator_strengths,
            entry["timestamp"]
        ))

    conn.commit()
    conn.close()
    print(f"✅ Imported {len(entries)} entries from {csv_path}")

def float_or_none(val):
    try:
        return float(val)
    except (TypeError, ValueError):
        return None

# データベースの分子エネルギーを更新する関数
def update_molecule_energy(
    molecule_id,
    g_tot=None,
    zpe=None,
    db_path="energy_db.sqlite"
):

    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    # 更新するカラムと値を動的に構築
    fields = []
    values = []
    if g_tot is not None:
        fields.append("g_tot = ?")
        values.append(g_tot)
    if zpe is not None:
        fields.append("zpe = ?")
        values.append(zpe)

    if not fields:
        print("⚠️ 更新するフィールドが指定されていません。")
        conn.close()
        return

    values.append(molecule_id)
    query = f"UPDATE molecules SET {', '.join(fields)} WHERE id = ?"
    cur.execute(query, values)

    conn.commit()
    conn.close()
    print(f"✅ Molecule ID {molecule_id} を更新しました。")

# データベースから分子をIDで削除する関数
def delete_molecule_by_id(molecule_id, db_path="energy_db.sqlite"):
    import sqlite3

    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    cur.execute("DELETE FROM molecules WHERE id = ?", (molecule_id,))
    conn.commit()
    conn.close()

    print(f"🗑️ Molecule ID {molecule_id} を削除しました。")


def get_molecules_from_sqlite(
    inchikey, method, basis, spin, charge,
    solvent=None, dielectric=None,
    temperature=298.15, pressure=1.0,
    nstates=None, excited_spin=None, tda=None,
    db_path="energy_db.sqlite"
):
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row  # ← これを追加
    cur = conn.cursor()
    # テーブルがなければ作成
    cur.execute("""
    CREATE TABLE IF NOT EXISTS molecules (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        inchi TEXT NOT NULL,
        inchikey TEXT NOT NULL,
        mol TEXT,
        mw REAL,
        formula TEXT,
        charge INTEGER NOT NULL,
        spin INTEGER NOT NULL,
        g_tot REAL,
        zpe REAL,
        method TEXT NOT NULL,
        basis TEXT NOT NULL,
        solvent TEXT,
        dielectric REAL,
        temperature REAL,
        pressure REAL,
        frequencies TEXT,
        chk_file BLOB,
        num_imaginary INTEGER,
        nstates INTEGER,
        excited_spin TEXT,
        tda INTEGER,
        timestamp TEXT DEFAULT CURRENT_TIMESTAMP
    );
    """)

    query = """
    SELECT id, g_tot, zpe, frequencies, chk_file, num_imaginary,
           excited_energies, oscillator_strengths, timestamp
    FROM molecules
    WHERE inchikey = ? AND method = ? AND basis = ?
      AND spin = ? AND charge = ? AND temperature = ? AND pressure = ?
    """
    params = [inchikey, method, basis, spin, charge, temperature, pressure]

    if solvent is None:
        query += " AND solvent IS NULL"
    else:
        query += " AND solvent = ?"
        params.append(solvent)

    if dielectric is None:
        query += " AND dielectric IS NULL"
    else:
        query += " AND dielectric = ?"
        params.append(dielectric)

    # nstates, excited_spin, tda の条件をNoneのとき無視
    if nstates is not None:
        query += " AND nstates = ?"
        params.append(nstates)
    if excited_spin is not None:
        query += " AND excited_spin = ?"
        params.append(excited_spin)
    if tda is not None:
        query += " AND tda = ?"
        params.append(int(tda))

    query += " ORDER BY id DESC"
    cur.execute(query, params)
    rows = cur.fetchall()

    # DBの全データを確認
    cur.execute("SELECT id, inchikey, method, basis, spin, charge, solvent, dielectric, temperature, pressure, nstates, excited_spin, tda FROM molecules")
    all_rows = cur.fetchall()

    # 検索実行
    cur.execute(query, params)
    rows = cur.fetchall()
    print("DEBUG: rows =", rows)
    conn.close()
    results = [dict(row) for row in rows]  # ここでdict化
    return results


# 使用例
if __name__ == "__main__":
    from rdkit import Chem
    from rdkit.Chem import inchi

    smiles = "CCO"
    mol = Chem.MolFromSmiles(smiles)
    inchi_str = inchi.MolToInchi(mol)
    inchikey_str = inchi.InchiToInchiKey(inchi_str)

    insert_molecule_with_frequencies(
        inchi=inchi_str,
        inchikey=inchikey_str,
        g_tot=-154.1234,
        zpe=0.0456,
        frequencies=[1234.5, 1500.0, 2980.2],
        method="B3LYP",
        basis="def2-SVP",
        charge=0,
        spin=0,
        solvent="PCM",
        dielectric=78.4,
        temperature=298.15,
        pressure=1.0
    )