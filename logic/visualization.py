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
from scipy.stats import norm
import streamlit as st


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


# streamlitの警告メッセージを表示する関数
def show_imaginary_frequency_warning(freq_data):
    ""
    " 振動数データから虚振動の有無をチェックし、警告メッセージを表示する。"
    freq_info = freq_data.get("frequencies", {})

    # 文字列なら辞書に変換
    if isinstance(freq_info, str):
        try:
            freq_info = eval(freq_info)
        except Exception:
            st.warning("⚠️ 振動数データの解析に失敗しました。")
            return

    freq_wavenumbers = freq_info.get("freq_wavenumber", None)

    if freq_wavenumbers is None:
        st.warning("⚠️ 振動数データが存在しません。")
        return

    # 虚数成分の有無をチェック
    try:
        imag_freqs = [f for f in freq_wavenumbers if isinstance(f, complex) and f.imag != 0]
        if imag_freqs:
            st.error(f"❌ 虚振動あり: {len(imag_freqs)} 個（例: {imag_freqs[0]:.2f}）")
        else:
            st.success("✅ 虚振動なし")
    except Exception:
        st.warning("⚠️ 虚振動の判定中にエラーが発生しました。")

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

# UVスペクトルの可視化
def plot_uv_spectrum(
    wavelengths,
    wavelength_min,
    wavelength_max,
    width=10,
    n_points=2000,
    ax=None
):

    # Noneを除外
    wavelengths = [(w, o) for w, o in wavelengths if w is not None and o is not None]
    wavelengths = sorted(wavelengths, key=lambda x: x[0])
    filtered = [(w, o) for w, o in wavelengths if wavelength_min <= w <= wavelength_max]

    x = np.linspace(wavelength_min, wavelength_max, n_points)
    spectrum = np.zeros_like(x)

    for w, o in filtered:
        spectrum += o.real * norm.pdf(x, w, width) * width * np.sqrt(2 * np.pi)

    stick_wavelengths = [w for w, o in filtered]
    stick_intensities = [o for w, o in filtered]

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 4))
    else:
        fig = ax.figure

    ax.plot(x, spectrum, color='blue', alpha=0.9, linewidth=2, label='Gaussian Broadening')
    ax.vlines(stick_wavelengths, 0, stick_intensities, color='red', alpha=0.3, linewidth=1, label='Transitions')
    ax.set_title("UV-Vis Spectrum (Gaussian + Sticks)")
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Intensity (a.u.)")
    ax.set_xlim(wavelength_min, wavelength_max)
    ax.grid(True)
    ax.legend()
    return fig, ax

def find_allowed_excited_state_from_lists(excitation_energies, oscillator_strengths, threshold=0.01):
    """
    エネルギーとオシレーター強度のリストから、最初に有意な遷移を返す。

    Parameters:
        excitation_energies: List[float]  # eV単位
        oscillator_strengths: List[float]
        threshold: float = 0.01

    Returns:
        (state_index, energy_in_eV, oscillator_strength) or (None, None, None)
    """
    for i, (e, f) in enumerate(zip(excitation_energies, oscillator_strengths)):
        f = float(f)  # オシレーター強度をfloatに変換
        e = float(e)  # エネルギーをfloatに変換
        if f >= threshold:
            return i + 1, e, f  # 1-indexed
    return None, None, None

def prepare_excited_states_table(excitation_energies, oscillator_strengths, threshold=0.01):
    """
    励起状態のリストをDataFrameに整形し、有意な遷移状態を選定する。
    """
    df = pd.DataFrame({
        "State": list(range(1, len(excitation_energies) + 1)),
        "Energy (eV)": excitation_energies,
        "Oscillator Strength": oscillator_strengths
    })

    for i, (e, f) in enumerate(zip(excitation_energies, oscillator_strengths)):
        try:
            f_value = float(f)
            e_value = float(e)
        except Exception:
            continue
        if f_value >= threshold:
            selected_state = {"index": i + 1, "energy": e_value, "strength": f_value}
            return df, selected_state

    return df, None

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

import numpy as np
import json

def make_json_safe(obj):
    """再帰的に complex128 や numpy 型を JSON で扱える形式に変換"""
    if isinstance(obj, (np.complexfloating, complex)):
        return obj.real  # または str(obj) にしてもOK
    elif isinstance(obj, (np.integer, np.floating)):
        return obj.item()
    elif isinstance(obj, np.ndarray):
        return [make_json_safe(x) for x in obj.tolist()]
    elif isinstance(obj, (list, tuple)):
        return [make_json_safe(x) for x in obj]
    elif isinstance(obj, dict):
        return {k: make_json_safe(v) for k, v in obj.items()}
    return obj

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
        if modes.shape == (n_modes, n_atoms, 3):
            vec = modes[i]
        elif modes.shape == (n_atoms * 3, n_modes):
            vec = modes[:, i].reshape(n_atoms, 3)
        elif modes.shape == (n_modes, n_atoms * 3):
            vec = modes[i, :].reshape(n_atoms, 3)
        else:
            raise ValueError(f"modes shape {modes.shape} is not compatible with n_atoms={n_atoms}")

        vec = vec * BOHR_TO_ANGSTROM
        vec_list = vec.tolist()

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

    # ✅ JSONセーフ化してから出力
    return json.dumps(make_json_safe(data), indent=2)


# def generate_cjson(mol, freqs, modes):
#     BOHR_TO_ANGSTROM = 0.529177
#     n_atoms = mol.natm
#     atomic_numbers = [mol.atom_charge(i) for i in range(n_atoms)]
#     coords = (mol.atom_coords() * BOHR_TO_ANGSTROM).flatten().tolist()

#     modes = np.array(modes)
#     n_modes = modes.shape[0]
#     n_freqs = len(freqs)
#     n = min(n_modes, n_freqs)

#     vibrational_modes = []
#     for i in range(n):
#         # shapeが(n_modes, n_atoms, 3)の場合
#         if modes.shape == (n_modes, n_atoms, 3):
#             vec = modes[i]  # shape: (n_atoms, 3)
#         # shapeが(n_atoms*3, n_modes)の場合
#         elif modes.shape == (n_atoms * 3, n_modes):
#             vec = modes[:, i].reshape(n_atoms, 3)
#         # shapeが(n_modes, n_atoms*3)の場合
#         elif modes.shape == (n_modes, n_atoms * 3):
#             vec = modes[i, :].reshape(n_atoms, 3)
#         else:
#             raise ValueError(f"modes shape {modes.shape} is not compatible with n_atoms={n_atoms}")
        
#         # 🔧 Å単位に変換
#         vec = vec * BOHR_TO_ANGSTROM
#         vec_list = vec.tolist()  # shape: (n_atoms, 3)

#         vibrational_modes.append({
#             "frequency": freqs[i],
#             "eigenVectors": vec_list
#         })

#     data = {
#         "chemical json": 0,
#         "atoms": {
#             "elements": {"number": atomic_numbers},
#             "coords": {"3d": coords}
#         },
#         "vibrationalModes": vibrational_modes
#     }

#     return json.dumps(data, indent=2)


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