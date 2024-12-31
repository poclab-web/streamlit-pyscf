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