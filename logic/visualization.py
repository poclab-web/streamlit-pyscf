"""
pySCFの計算を実行した結果をわかりやすく表示する部分

1. エネルギー関連の可視化
エネルギー収束グラフ:
PySCFの収束過程をプロットし、計算の収束性を確認。
横軸に繰り返し回数（iteration）、縦軸にエネルギーをプロット。

2. 周波数スペクトルの可視化
周波数解析結果をスペクトル形式でプロット。
横軸を振動数（cm⁻¹）、縦軸を強度（Intensity）としてプロット。

3. エネルギー分解解析 (EDA) の可視化
エネルギー分解結果を棒グラフで表示し、どの項が支配的かを可視化。

4. UVスペクトル予測
TD-DFT計算の結果を可視化し、吸収スペクトルをプロット。
横軸を波長（nm）、縦軸を吸収強度（a.u.）として描画。

5. NMRスペクトルの可視化
NMR計算結果をスペクトル形式でプロット。
横軸を化学シフト（ppm）、縦軸を強度（Intensity）としてプロット。

"""

import matplotlib.pyplot as plt
import json
import io
import numpy as np
import pandas as pd



def plot_energy_convergence(energies):
    """
    # エネルギー関連の可視化
    エネルギー収束のグラフをプロット。
    """
    plt.figure(figsize=(8, 5))
    plt.plot(range(len(energies)), energies, marker="o", linestyle="-")
    plt.xlabel("Iteration")
    plt.ylabel("Energy (Hartree)")
    plt.title("Energy Convergence")
    plt.grid()
    plt.show()

def plot_frequency_spectrum(frequencies, intensities):
    """
    周波数スペクトルの可視化
    周波数スペクトルをプロット。
    """
    plt.figure(figsize=(10, 6))
    plt.bar(frequencies, intensities, width=10, color="blue", alpha=0.7)
    plt.xlabel("Frequency (cm^-1)")
    plt.ylabel("Intensity")
    plt.title("Vibrational Frequency Spectrum")
    plt.grid(axis="y")
    plt.show()

def plot_energy_decomposition(terms, values):
    """
    エネルギー分解解析 (EDA) の可視化
    エネルギー分解解析の結果をプロット。
    """
    plt.figure(figsize=(8, 5))
    plt.bar(terms, values, color="green", alpha=0.8)
    plt.xlabel("Energy Term")
    plt.ylabel("Energy Contribution (Hartree)")
    plt.title("Energy Decomposition Analysis")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()


def plot_uv_spectrum(wavelengths, intensities):
    """
    UVスペクトルをプロット。
    """
    plt.figure(figsize=(10, 6))
    plt.plot(wavelengths, intensities, color="purple", marker="o", linestyle="-")
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("Absorption Intensity (a.u.)")
    plt.title("UV Spectrum")
    plt.grid()
    plt.show()

def plot_nmr_spectrum(chemical_shifts, intensities):
    """
    5. NMRスペクトルの可視化
    """
    plt.figure(figsize=(10, 6))
    plt.stem(chemical_shifts, intensities, linefmt="C1-", markerfmt="C1o", basefmt="C0-")
    plt.xlabel("Chemical Shift (ppm)")
    plt.ylabel("Intensity")
    plt.title("NMR Spectrum")
    plt.grid()
    plt.show()

def simulate_nmr_spectrum(shifts_ppm, j_values_hz, intensities=None, lw=1.0, points=8000, width=0.02):
    if intensities is None:
        intensities = [1.0] * len(shifts_ppm)
    x = np.linspace(0, 10, points)
    y = np.zeros_like(x)
    for shift, j, inten in zip(shifts_ppm, j_values_hz, intensities):
        delta = j / 2 / 400  # 400 MHz NMR換算でHz→ppmに
        positions = [shift - delta, shift + delta]
        for pos in positions:
            line = inten * np.exp(-((x - pos) ** 2) / (2 * width ** 2))
            y += line
    fig, ax = plt.subplots(figsize=(10, 3))
    ax.plot(x, y, color='black')
    ax.invert_xaxis()
    ax.set_xlabel("Chemical Shift (ppm)")
    ax.set_ylabel("Intensity")
    ax.set_title("Simulated 1H-NMR Spectrum")
    plt.tight_layout()
    return fig

def simulate_nmr_spectrum_from_csv(
    csv_path,
    target_element=None,
    TMS_reference_csv_path=None,
    # chcl3_csv_path=None,
    # chcl3_exp_shift=7.26,
    j_values_hz=None,
    intensities=None,
    lw=1.0,
    points=8000,
    width=0.02
):
    """
    NMRシールド値CSVファイルを読み込み、指定元素のみでNMRスペクトルをシミュレートして表示する
    chcl3_csv_pathを指定した場合、CHCl3の計算値を使ってスケーリング補正を行う
    """
    df = pd.read_csv(csv_path)
    if target_element is not None:
        df = df[df["Element"] == target_element]

    # TMS基準値の取得
    if TMS_reference_csv_path is not None:
        tms_df = pd.read_csv(TMS_reference_csv_path)
        sigma_tms = tms_df["NMR Shielding"].mean()
    else:
        sigma_tms = df["NMR Shielding"].mean()

    # # CHCl3基準値の取得とスケーリングファクター
    # scaling_factor = 1.0
    # if chcl3_csv_path is not None:
    #     chcl3_df = pd.read_csv(chcl3_csv_path)
    #     sigma_chcl3 = chcl3_df["NMR Shielding"].mean()
    #     # 実験値7.26ppm, TMSは0ppm
    #     scaling_factor = chcl3_exp_shift / (sigma_tms - sigma_chcl3)

    # ケミカルシフト計算＋スケーリング
    df["Chemical Shift (ppm)"] = sigma_tms - df["NMR Shielding"]
    shifts_ppm = df["Chemical Shift (ppm)"]

    if j_values_hz is None:
        j_values_hz = [0.0] * len(shifts_ppm)
    if intensities is None:
        intensities = [1.0] * len(shifts_ppm)

    x = np.linspace(0, 10, points)
    y = np.zeros_like(x)
    for shift, j, inten in zip(shifts_ppm, j_values_hz, intensities):
        delta = j / 2 / 400
        positions = [shift - delta, shift + delta]
        for pos in positions:
            line = inten * np.exp(-((x - pos) ** 2) / (2 * width ** 2))
            y += line
    fig, ax = plt.subplots(figsize=(10, 3))
    ax.plot(x, y, color='black')
    ax.invert_xaxis()
    ax.set_xlabel("Chemical Shift (ppm)")
    ax.set_ylabel("Intensity")
    ax.set_title("Simulated 1H-NMR Spectrum (scaled)")
    plt.tight_layout()

    table_df = df[["Atom Index", "Element", "Chemical Shift (ppm)"]].reset_index(drop=True)
    return fig, table_df

def generate_cjson(mol, freqs, modes):
    BOHR_TO_ANGSTROM = 0.529177
    n_atoms = mol.natm
    atomic_numbers = [mol.atom_charge(i) for i in range(n_atoms)]
    coords = (mol.atom_coords() * BOHR_TO_ANGSTROM).flatten().tolist()

    modes = np.array(modes)
    n_modes = modes.shape[0]
    n_freqs = len(freqs)
    n = min(n_modes, n_freqs)

    vibrational_modes = []
    for i in range(n):
        # shapeが(n_modes, n_atoms, 3)の場合
        if modes.shape == (n_modes, n_atoms, 3):
            vec = modes[i]  # shape: (n_atoms, 3)
        # shapeが(n_atoms*3, n_modes)の場合
        elif modes.shape == (n_atoms * 3, n_modes):
            vec = modes[:, i].reshape(n_atoms, 3)
        # shapeが(n_modes, n_atoms*3)の場合
        elif modes.shape == (n_modes, n_atoms * 3):
            vec = modes[i, :].reshape(n_atoms, 3)
        else:
            raise ValueError(f"modes shape {modes.shape} is not compatible with n_atoms={n_atoms}")
        
        # 🔧 Å単位に変換
        vec = vec * BOHR_TO_ANGSTROM
        vec_list = vec.tolist()  # shape: (n_atoms, 3)

        vibrational_modes.append({
            "frequency": freqs[i],
            "eigenVectors": vec_list
        })

    data = {
        "chemical json": 0,
        "atoms": {
            "elements": {"number": atomic_numbers},
            "coords": {"3d": coords}
        },
        "vibrationalModes": vibrational_modes
    }

    return json.dumps(data, indent=2)


# TODO : gaussianで読み取れない。
def write_gaussian_log(mol, freqs, modes, filename="vibrations.log"):
    BOHR_TO_ANGSTROM = 0.529177
    n_atoms = mol.natm
    atom_symbols = [mol.atom_symbol(i) for i in range(n_atoms)]
    atom_numbers = [mol.atom_charge(i) for i in range(n_atoms)]

    # 原子座標取得（省略可能）
    coords = mol.atom_coords() * BOHR_TO_ANGSTROM

    modes = np.array(modes)
    n_modes = len(freqs)
    assert modes.shape[0] in [n_modes, n_atoms * 3]

    with open(filename, 'w') as f:
        f.write(" Entering Gaussian System, Link 801\n")
        f.write("\n Frequencies -- " + "  ".join(f"{freqs[i]:10.4f}" for i in range(min(3, n_modes))) + "\n")
        f.write(" Red. masses --" + "  ".join(f"{1.0:10.4f}" for _ in range(min(3, n_modes))) + "\n")
        f.write(" Frc consts  --" + "  ".join(f"{0.5:10.4f}" for _ in range(min(3, n_modes))) + "\n")
        f.write(" IR Inten    --" + "  ".join(f"{10.0:10.4f}" for _ in range(min(3, n_modes))) + "\n\n")

        f.write(" Atom  AN      X      Y      Z\n")

        # 最初の3モードだけ表示（Avogadroでは1モードずつでOK）
        for i_atom in range(n_atoms):
            f.write(f" {atom_symbols[i_atom]:>4} {int(atom_numbers[i_atom]):>3}")
            for i_mode in range(min(3, n_modes)):
                if modes.shape[0] == n_modes:
                    vec = modes[i_mode].reshape(n_atoms, 3)[i_atom]
                else:
                    vec = modes[:, i_mode].reshape(n_atoms, 3)[i_atom]
                vec = vec * BOHR_TO_ANGSTROM
                f.write(f" {vec[0]:>8.4f} {vec[1]:>8.4f} {vec[2]:>8.4f}")
            f.write("\n")

        f.write("\n Normal termination of Gaussian 09\n")

def read_cube_values(cube_path):
    with open(cube_path) as f:
        lines = f.readlines()
    n_atoms = int(lines[2].split()[0])
    data_lines = lines[6 + n_atoms:]
    values = np.array([float(v) for line in data_lines for v in line.strip().split()])
    return values

def save_cube_from_data(filename, cube_data):
    try:
        with open(filename, "w") as f:
            f.write(cube_data)
        st.success(f"✅ {filename} を保存しました。")
    except Exception as e:
        st.error(f"❌ 保存に失敗しました: {e}")