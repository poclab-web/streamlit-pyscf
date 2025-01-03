"""
pySCFの計算を実行する部分
1. 1点計算
2. 構造最適化計算
3. 振動計算
4. 励起状態計算（時間依存密度汎関数理論: TD-DFT）
5. 分子の応答性（極性率や双極子モーメント）
"""


from pyscf import gto, scf, dft, hessian
import numpy as np


# 理論と基底関数の選択肢
# TODO: 今後のことを考えて、自動取得もしくはファイルを分割したい。
theory_options = ["HF", "B3LYP", "PBE", "M06-2X", "B97X-D"]
basis_set_options = [
    "sto-3g", "6-31g", "6-31g*", "6-31+g(d,p)", "cc-pVDZ",
    "cc-pVTZ", "aug-cc-pVDZ", "aug-cc-pVTZ", "def2-SVP", "def2-TZVP"
]
hartree_to_cm1 = 219474.63  # 1 Hartree = 219474.63 cm^-1

# 1. 1点計算
def run_quantum_calculation(atom_input, basis_set, theory):
    """
    Execute a quantum chemistry calculation.

    Parameters:
        atom_input (str): Atomic coordinates in XYZ format.
        basis_set (str): Basis set for the calculation.
        theory (str): Theory to be used for the calculation (HF, DFT, MP2, etc.).

    Returns:
        float: The calculated energy.
    """
    try:
        # 分子のセットアップ
        mol = gto.M(atom=atom_input, basis=basis_set)
        
        # 理論に応じて計算を実行
        if theory == "HF":
            mf = scf.RHF(mol)
        elif theory == "MP2":
            mf = scf.RHF(mol).run()
            from pyscf import mp
            mf = mp.MP2(mf).run()
        elif theory in ["B3LYP", "PBE", "M06-2X", "B97X-D"]:
            mf = dft.RKS(mol)
            mf.xc = theory
        else:
            raise ValueError(f"Unsupported theory: {theory}")
        
        energy = mf.kernel()
        return energy
    except Exception as e:
        raise RuntimeError(f"Calculation failed: {e}")


# 2. 構造最適化計算
def run_geometry_optimization(atom_input, basis_set, theory, conv_params):
    """
    Perform geometry optimization for a molecule using PySCF and geometric.

    Parameters:
        atom_input (str): Atomic coordinates in XYZ format.
        basis_set (str): Basis set for the calculation.
        theory (str): Theory to be used for the calculation (HF, DFT, etc.).
        conv_params (dict): Convergence parameters for optimization.

    Returns:
        list: Energies during optimization.
        list: Geometries during optimization (list of numpy arrays).
    """
    try:
        # 分子のセットアップ
        mol = gto.M(atom=atom_input, basis=basis_set)

        # 理論の設定
        if theory == "HF":
            mf = scf.RHF(mol)
        elif theory in ["B3LYP", "PBE", "M06-2X", "B97X-D"]:
            mf = dft.RKS(mol)
            mf.xc = theory
        else:
            raise ValueError(f"Unsupported theory: {theory}")

        # エネルギーとジオメトリを記録するためのリスト
        energies = []
        geometries = []

        # コールバック関数
        def callback(update):
            """
            Callback function to capture optimization steps.

            Parameters:
                update (dict): A dictionary containing geometry and energy information.
            """
            if "coords" in update and "energy" in update:
                geometries.append(update["coords"])
                energies.append(update["energy"])

        # 構造最適化を実行
        geometric_solver.optimize(mf, callback=callback, **conv_params)

        return energies, geometries
    
    except Exception as e:
        raise RuntimeError(f"Optimization failed: {e}")


# 3. 振動計算
def calculate_vibrational_frequencies(molecule: str, basis: str = 'cc-pVDZ'):
    """
    Calculate vibrational frequencies for a given molecule.
    
    Parameters:
        molecule (str): Atomic coordinates in PySCF format.
        basis (str): Basis set to use for the calculation.
        
    Returns:
        frequencies (list): List of vibrational frequencies (cm^-1).
    """
    # 定義する分子
    mol = gto.M(atom=molecule, basis=basis, symmetry=True)

    # SCF計算
    mf = scf.RHF(mol)
    mf.kernel()

    # ヘシアン行列の計算
    hess = hessian.RHF(mf).kernel()

    # 質量重み付きヘシアン行列の作成
    natoms = mol.natm
    mass = np.repeat(mol.atom_mass_list(), 3)
    sqrt_mass = np.sqrt(mass)
    mass_weighted_hess = hess / np.outer(sqrt_mass, sqrt_mass)

    # 固有値と固有ベクトルを計算
    eigenvalues, _ = np.linalg.eigh(mass_weighted_hess)

    # 周波数（cm^-1）に変換
    frequencies = np.sqrt(np.abs(eigenvalues)) * np.sign(eigenvalues) * hartree_to_cm1
    frequencies = np.real_if_close(frequencies)  # 虚数成分を除去

    return frequencies.tolist()


# 4. 励起状態計算（時間依存密度汎関数理論：TD-DFT）
from pyscf import tdscf

def calculate_excited_states(molecule: str, basis: str = 'cc-pVDZ', theory: str = 'B3LYP', num_states: int = 5):
    """
    Calculate excited states using TD-DFT.

    Parameters:
        molecule (str): Atomic coordinates in PySCF format.
        basis (str): Basis set to use for the calculation.
        theory (str): DFT functional to use for the calculation.
        num_states (int): Number of excited states to compute.

    Returns:
        list: Excitation energies (eV).
    """
    mol = gto.M(atom=molecule, basis=basis)
    mf = dft.RKS(mol)
    mf.xc = theory
    mf.kernel()

    # TD-DFT計算
    td = tdscf.TDDFT(mf)
    td.nstates = num_states
    excitation_energies = td.kernel()

    # ハートリー単位から電子ボルト（eV）に変換
    hartree_to_ev = 27.2114
    excitation_energies_ev = [exc[0] * hartree_to_ev for exc in excitation_energies]

    return excitation_energies_ev

# 5. 分子の応答性（極性率や双極子モーメント）
def calculate_polarizability(molecule: str, basis: str = 'cc-pVDZ'):
    """
    Calculate molecular polarizability using SCF.

    Parameters:
        molecule (str): Atomic coordinates in PySCF format.
        basis (str): Basis set to use for the calculation.

    Returns:
        numpy.ndarray: Polarizability tensor (au).
    """
    mol = gto.M(atom=molecule, basis=basis)
    mf = scf.RHF(mol)
    mf.kernel()

    # 応答計算
    polarizability = mf.Polarizability().kernel()
    return polarizability

# 6. 溶媒効果を考慮した計算

from pyscf import solvent

def calculate_in_solvent(molecule: str, basis: str = 'cc-pVDZ', theory: str = 'B3LYP', solvent_model: str = 'pcm'):
    """
    Perform a calculation considering solvent effects using PCM.

    Parameters:
        molecule (str): Atomic coordinates in PySCF format.
        basis (str): Basis set to use for the calculation.
        theory (str): DFT functional to use for the calculation.
        solvent_model (str): Solvent model to use (e.g., 'pcm').

    Returns:
        float: Solvated energy.
    """
    mol = gto.M(atom=molecule, basis=basis)
    mf = dft.RKS(mol)
    mf.xc = theory

    # 溶媒モデルの設定
    solvent.set_pcm(mf)
    energy = mf.kernel()

    return energy
