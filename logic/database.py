import sqlite3
import json
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

def insert_molecule_with_frequencies(
    smiles,
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
    db_path="data/energy_db.sqlite"
):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    molblock = Chem.MolToMolBlock(mol)
    mw = Descriptors.MolWt(mol)
    formula = rdMolDescriptors.CalcMolFormula(mol)

    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    # テーブル作成（存在しなければ）
    cur.execute("""
    CREATE TABLE IF NOT EXISTS molecules (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        smiles TEXT NOT NULL,
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
        timestamp TEXT DEFAULT CURRENT_TIMESTAMP
    );
    """)

    # JSON文字列に変換
    freq_json = json.dumps(frequencies) if frequencies else None

    # データ挿入
    cur.execute("""
    INSERT INTO molecules (
        smiles, mol, mw, formula, charge, spin, g_tot, zpe,
        method, basis, solvent, dielectric, temperature, pressure,
        frequencies, timestamp
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, (
        smiles, molblock, mw, formula, charge, spin, g_tot, zpe,
        method, basis, solvent, dielectric, temperature, pressure,
        freq_json, datetime.now().isoformat()
    ))

    conn.commit()
    conn.close()
    print(f"✅ Inserted molecule: {smiles} ")

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

def get_molecule_from_sqlite(smiles, method, basis, spin, charge,
                              solvent=None, dielectric=None,
                              temperature=298.15, pressure=1.0,
                              db_path="energy_db.sqlite"):
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()

    query = """
    SELECT id, g_tot, zpe, frequencies
    FROM molecules
    WHERE smiles = ? AND method = ? AND basis = ?
      AND spin = ? AND charge = ? AND temperature = ? AND pressure = ?
    """
    params = [smiles, method, basis, spin, charge, temperature, pressure]

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
        molecule_id, g_tot, zpe, freq_json = row
        frequencies = json.loads(freq_json) if freq_json else None
        return {
            "id": molecule_id,
            "g_tot": g_tot,
            "zpe": zpe,
            "frequencies": frequencies
        }
    else:
        return None


# 使用例
if __name__ == "__main__":
    insert_molecule_with_frequencies(
        smiles="CCO",
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
