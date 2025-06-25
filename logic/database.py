import sqlite3
import json
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import inchi
import csv

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
    db_path="data/energy_db.sqlite"
):
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    import sqlite3
    import json
    from datetime import datetime

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
        timestamp TEXT DEFAULT CURRENT_TIMESTAMP
    );
    """)

    freq_json = json.dumps(frequencies) if frequencies is not None else None

    # データ挿入
    cur.execute("""
    INSERT INTO molecules (
        inchi, inchikey, mol, mw, formula, charge, spin, g_tot, zpe,
        method, basis, solvent, dielectric, temperature, pressure,
        frequencies, chk_file, num_imaginary, timestamp
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, (
        inchi, inchikey, molblock, mw, formula, charge, spin, g_tot, zpe,
        method, basis, solvent, dielectric, temperature, pressure,
        freq_json, chk_file_data, num_imaginary, datetime.now().isoformat()
    ))

    conn.commit()
    conn.close()
    print(f"✅ Inserted molecule: {inchikey} ")

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

# データベースから分子データを取得する関数
def get_molecule_from_sqlite(inchikey, method, basis, spin, charge,
                              solvent=None, dielectric=None,
                              temperature=298.15, pressure=1.0,
                              db_path="energy_db.sqlite"):
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
        timestamp TEXT DEFAULT CURRENT_TIMESTAMP
    );
    """)
    query = """
    SELECT id, g_tot, zpe, frequencies, chk_file, num_imaginary
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

# データベースから安定な分子を取得する関数
def get_stable_molecule_from_sqlite(
    inchikey,  # ←ここをsmilesからinchikeyに変更
    method, basis, spin, charge,
    solvent=None, dielectric=None,
    temperature=298.15, pressure=1.0,
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
            num_imaginary, timestamp
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            entry["inchi"], entry["inchikey"], entry["mol"], entry["mw"], entry["formula"],
            entry["charge"], entry["spin"], entry["g_tot"], entry["zpe"],
            entry["method"], entry["basis"], entry["solvent"], entry["dielectric"],
            entry["temperature"], entry["pressure"],
            json.dumps(entry["frequencies"]) if entry["frequencies"] else None,
            chk_file_data, entry["num_imaginary"], entry["timestamp"]
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
                frequencies = json.dumps(json.loads(frequencies))  # json文字列として整形
            except Exception:
                frequencies = None

        cur.execute("""
        INSERT INTO molecules (
            inchi, inchikey, mol, mw, formula, charge, spin,
            g_tot, zpe, method, basis, solvent, dielectric,
            temperature, pressure, frequencies, chk_file,
            num_imaginary, timestamp
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            entry["inchi"], entry["inchikey"], entry["mol"], float(entry["mw"]), entry["formula"],
            int(entry["charge"]), int(entry["spin"]), float_or_none(entry["g_tot"]),
            float_or_none(entry["zpe"]), entry["method"], entry["basis"], entry["solvent"],
            float_or_none(entry["dielectric"]), float(entry["temperature"]),
            float(entry["pressure"]), frequencies, chk_file_data,
            int(entry["num_imaginary"]) if entry["num_imaginary"] else None,
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
    import sqlite3

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
