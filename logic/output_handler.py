"""
pySCFの計算を実行するための出力フォーマットに変更する部分

1. 結果の整形
PySCFの計算結果は多くの場合、数値や配列形式で出力される。それらを分かりやすいキーを持つ辞書形式やJSON形式に変換します。
例:
エネルギー値の整形
構造最適化後の座標情報の整形
フォノン周波数や振動モードの整理

"""

import os

def format_energy_output(energy):
    """
    シングルポイントエネルギーの出力を整形。
    """
    return {"energy (Hartree)": energy}

def convert_energy_units(energy_hartree, unit="eV"):
    """
    エネルギーを指定された単位に変換。
    """
    conversion_factors = {"eV": 27.2114, "kcal/mol": 627.509}
    if unit in conversion_factors:
        return energy_hartree * conversion_factors[unit]
    else:
        return energy_hartree

def extract_orbital_energies(molden_file_path):
    """
    MOLDENファイルから軌道エネルギーを抽出し、HOMOとLUMOを特定する。

    Parameters:
        molden_file_path (str): MOLDENファイルのパス。

    Returns:
        dict: 軌道エネルギー、HOMO、LUMOを含む情報。
    """
    if not os.path.exists(molden_file_path):
        raise FileNotFoundError(f"MOLDENファイルが見つかりません: {molden_file_path}")

    orbital_energies = []
    homo_index = -1

    with open(molden_file_path, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if line.startswith(" Ene="):
            # 軌道エネルギーを抽出
            energy = float(line.split("=")[1].strip())
            orbital_energies.append(energy)

            # Occupancy（占有数）を確認してHOMOを特定
            occupancy_line = lines[i + 2]
            occupancy = float(occupancy_line.split("=")[1].strip())
            if occupancy > 0:
                homo_index = len(orbital_energies) - 1

    if homo_index == -1:
        raise ValueError("HOMOが見つかりませんでした。")

    lumo_index = homo_index + 1 if homo_index + 1 < len(orbital_energies) else None

    homo_index = max(i for i, e in enumerate(orbital_energies) if e < 0)  # HOMOのインデックス
    lumo_index = homo_index + 1 if homo_index + 1 < len(orbital_energies) else None  # LUMOのインデックス
    return {
        "orbital_energies": orbital_energies,
        "homo_index": homo_index,
        "lumo_index": lumo_index,
    }
