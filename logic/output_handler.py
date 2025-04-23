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
import re
import streamlit as st

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

    return {
        "orbital_energies": orbital_energies,
        "homo_index": homo_index,
        "lumo_index": lumo_index,
    }

# フォルダパスを解析して inchikey とファイル名を抽出する関数
def parse_folder_and_filename(file_path):
    try:
        # フォルダを "/" で分割
        parts = file_path.split(os.sep)
        inchikey = parts[-2]  # data の次のフォルダ名が inchikey
        filename = parts[-1]  # ファイル名
        return inchikey, filename
    except Exception as e:
        st.write(f"Error parsing path {file_path}: {e}")
        return None, None

# MOLDENファイル名を解析して情報を抽出する関数
def parse_filename(filename):
    try:
        # ファイル名を "_" と "." で分割
        parts = filename.split("_")
        theory = parts[0]
        basis_set = parts[1]
        index, extension = parts[2].split(".")
        return theory, basis_set, index, extension
    except Exception as e:
        st.write(f"Error parsing filename {filename}: {e}")
        return None, None, None, None

# .outファイルから最後に出現するHOMOおよびLUMOエネルギーを抽出する関数
def extract_homo_lumo_from_out(out_file):
    try:
        with open(out_file, "r") as f:
            lines = f.readlines()
        
        # HOMOとLUMOのエネルギーを格納する変数
        homo_energy = None
        lumo_energy = None

        # 正規表現でHOMOとLUMOのエネルギーを検索
        homo_pattern = r"HOMO.*?=\s*(-?\d+\.\d+)"
        lumo_pattern = r"LUMO.*?=\s*(-?\d+\.\d+)"

        for line in lines:
            homo_match = re.search(homo_pattern, line)
            lumo_match = re.search(lumo_pattern, line)
            if homo_match:
                homo_energy = float(homo_match.group(1))
            if lumo_match:
                lumo_energy = float(lumo_match.group(1))
        
        return homo_energy, lumo_energy
    except Exception as e:
        st.write(f"Error processing {out_file}: {e}")
        return None, None
