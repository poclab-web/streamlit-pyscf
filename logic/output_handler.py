"""
pySCFの計算を実行するための出力フォーマットに変更する部分

1. 結果の整形
PySCFの計算結果は多くの場合、数値や配列形式で出力されます。それらを分かりやすいキーを持つ辞書形式やJSON形式に変換します。
例:
エネルギー値の整形
構造最適化後の座標情報の整形
フォノン周波数や振動モードの整理

"""

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
