"""
pySCFの計算を実行して、実行した結果を出力させる部分
ファイルの入力については、molecule_handler.pyに行わせている。
1. 1点計算
2. 構造最適化計算
3. 振動計算
4. 励起状態計算（時間依存密度汎関数理論: TD-DFT）
5. 分子の応答性（極性率や双極子モーメント）
出てきたファイルの解析については、output_handler.pyで数値の変換し、visualization.pyでグラフなどに変換する。
"""

import os
from pyscf import gto, scf, dft, solvent
import json
from datetime import datetime
import tempfile
import sys

from pyscf.geomopt.berny_solver import optimize


# 理論と基底関数の選択肢
theory_options = ["HF", "B3LYP", "PBE", "M06-2X", "B97X-D"]
basis_set_options = [
    "sto-3g", "6-31g", "6-31g*", "6-31+g(d,p)", "cc-pVDZ",
    "cc-pVTZ", "aug-cc-pVDZ", "aug-cc-pVTZ", "def2-SVP", "def2-TZVP"
]
hartree_to_cm1 = 219474.63  # 1 Hartree = 219474.63 cm^-1


# ログ保存関数（化合物ごとのフォルダに保存）
def save_log(compound_name, data):
    """
    Save calculation data to a JSON log file in the compound-specific folder.

    Parameters:
        compound_name (str): Name of the compound.
        data (dict): A dictionary containing calculation parameters and results.
    """
    # ディレクトリ作成
    directory = os.path.join("data", compound_name)
    os.makedirs(directory, exist_ok=True)

    # ログファイルのパス
    log_file = os.path.join(directory, "calculation_log.json")

    # JSONデータを保存
    with open(log_file, 'a') as f:
        f.write(json.dumps(data, indent=4) + '\n')
    
    return directory  # ディレクトリパスを返す

def setup_molecule(atom_input, basis_set, charge=0, spin=0, solvent_model=None, eps=None):
    """
    Set up a PySCF molecule with specified parameters.

    Parameters:
        atom_input (str): Atomic coordinates in XYZ format.
        basis_set (str): Basis set for the calculation.
        charge (int): Molecular charge.
        spin (int): Spin multiplicity (2S + 1).
        solvent_model (str): Solvent model to be used (e.g., PCM).
        eps (float): Dielectric constant for solvent.

    Returns:
        gto.Mole: A PySCF molecule object.
    """
    mol = gto.M(atom=atom_input, basis=basis_set, charge=charge, spin=spin)
    
    if solvent_model == "PCM":
        pcm = solvent.PCM(mol)
        if eps is not None:
            pcm.eps = eps
        mol = pcm

    return mol

# 1. 1点計算
def run_quantum_calculation(compound_name, smiles, atom_input, basis_set, theory, charge=0, spin=0, solvent_model=None, eps=None):
    """
    Execute a quantum chemistry calculation.

    Parameters:
        compound_name (str): Name of the compound.
        atom_input (str): Atomic coordinates in XYZ format.
        basis_set (str): Basis set for the calculation.
        theory (str): Theory to be used for the calculation (HF, DFT, MP2, etc.).
        charge (int): Molecular charge.
        spin (int): Spin multiplicity (2S + 1).
        solvent_model (str): Solvent model to be used (e.g., PCM).
        eps (float): Dielectric constant for solvent.

    Returns:
        float: The calculated energy.
    """
    log_entry = {
        "time": "Start_Time",
        "timestamp": datetime.now().isoformat(),
        "compound": compound_name,
        "smiles": smiles,
        "calculation_type": "single_point",
        "parameters": {
            "atom_input": atom_input,
            "basis_set": basis_set,
            "theory": theory,
            "charge": charge,
            "spin": spin,
            "solvent_model": solvent_model,
            "dielectric": eps
        }
    }
    try:
        directory = save_log(compound_name, log_entry)
        mol = setup_molecule(atom_input, basis_set, charge, spin, solvent_model, eps)

        chkfile_name = os.path.join(directory, f"{compound_name}_{theory}.chk")

        with open(os.path.join(directory, f"output_{theory}.out"), 'w') as f:
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

            mf.chkfile = chkfile_name
            mf.verbose = 4
            mf.stdout = f

            energy = mf.kernel()
            log_entry["time"] = "End_Time"
            log_entry["result"] = {
                "energy": energy,
                "converged": mf.converged,
                "chkfile": os.path.abspath(chkfile_name)
            }
            save_log(compound_name, log_entry)

        return energy
    except Exception as e:
        log_entry["time"] = "Error_End"
        log_entry["error"] = str(e)
        save_log(compound_name, log_entry)
        raise RuntimeError(f"Calculation failed: {e}")


# 2. 構造最適化計算

def run_geometry_optimization(compound_name, smiles, atom_input, basis_set, theory, charge=0, spin=0, solvent_model=None, eps=None, conv_params=None):
    """
    Perform geometry optimization for a molecule using PySCF's native optimizer.

    Parameters:
        compound_name (str): Name of the compound.
        atom_input (str): Atomic coordinates in XYZ format.
        basis_set (str): Basis set for the calculation.
        theory (str): Theory to be used for the calculation (HF, DFT, etc.).
        charge (int): Molecular charge.
        spin (int): Spin multiplicity (2S + 1).
        solvent_model (str): Solvent model to be used (e.g., PCM).
        eps (float): Dielectric constant for solvent.
        conv_params (dict, optional): Convergence parameters for optimization.

    Returns:
        list: Final optimized geometry (list of numpy arrays).
    """
    log_entry = {
        "time": "Start_Time",
        "timestamp": datetime.now().isoformat(),
        "compound": compound_name,
        "smiles": smiles,
        "calculation_type": "geometry_optimization",
        "parameters": {
            "atom_input": atom_input,
            "basis_set": basis_set,
            "theory": theory,
            "charge": charge,
            "spin": spin,
            "solvent_model": solvent_model,
            "dielectric": eps,
            "conv_params": conv_params
        }
    }
    try:
        directory = save_log(compound_name, log_entry)
        mol = setup_molecule(atom_input, basis_set, charge, spin, solvent_model, eps)

        chkfile_name = os.path.join(directory, f"{compound_name}_{theory}.chk")

        with open(os.path.join(directory, f"output_{theory}.txt"), 'w') as f:
            if theory == "HF":
                mf = scf.RHF(mol)
            elif theory in ["B3LYP", "PBE", "M06-2X", "B97X-D"]:
                mf = dft.RKS(mol)
                mf.xc = theory
            else:
                raise ValueError(f"Unsupported theory: {theory}")

            mf.chkfile = chkfile_name
            mf.kernel()

            options = conv_params if conv_params else {}
            optimized_mol = optimize(mf, solver='geomeTRIC', options=options)

            final_geometry = optimized_mol.atom_coords()
            log_entry["time"] = "End_Time"
            log_entry["result"] = {"final_geometry": final_geometry.tolist()}
            save_log(compound_name, log_entry)

        return final_geometry

    except Exception as e:
        log_entry["time"] = "Error_End"
        log_entry["error"] = str(e)
        save_log(compound_name, log_entry)
        raise RuntimeError(f"Optimization failed: {e}")


# def run_geometry_optimization(compound_name, smiles, atom_input, basis_set, theory, conv_params=None):
#     """
#     Perform geometry optimization for a molecule using PySCF's native optimizer.

#     Parameters:
#         compound_name (str): Name of the compound.
#         atom_input (str): Atomic coordinates in XYZ format.
#         basis_set (str): Basis set for the calculation.
#         theory (str): Theory to be used for the calculation (HF, DFT, etc.).
#         conv_params (dict, optional): Convergence parameters for optimization.

#     Returns:
#         list: Final optimized geometry (list of numpy arrays).
#     """
#     log_entry = {
#         "time": "Start_Time",
#         "timestamp": datetime.now().isoformat(),
#         "compound": compound_name,
#         "smiles": smiles,
#         "calculation_type": "geometry_optimization",
#         "parameters": {
#             "atom_input": atom_input,
#             "basis_set": basis_set,
#             "theory": theory,
#             "conv_params": conv_params
#         }
#     }
#     try:
#         # ログファイルとチェックポイント保存用ディレクトリを作成
#         directory = save_log(compound_name, log_entry)

#         with open(os.path.join(directory, f"output_{theory}.txt"), 'w') as f:

#             # 分子のセットアップ
#             mol = gto.M(atom=atom_input, basis=basis_set)

#             # チェックポイントファイルのパスを作成
#             chkfile_name = os.path.join(directory, f"{compound_name}_{theory}.chk")

#             # 理論の設定
#             if theory == "HF":
#                 mf = scf.RHF(mol)
#                 mf.chkfile = chkfile_name  # チェックポイントファイルを設定
#             elif theory in ["B3LYP", "PBE", "M06-2X", "B97X-D"]:
#                 mf = dft.RKS(mol)
#                 mf.xc = theory
#                 mf.chkfile = chkfile_name  # チェックポイントファイルを設定
#             else:
#                 raise ValueError(f"Unsupported theory: {theory}")

#             # SCFを実行して初期エネルギーを計算
#             mf.kernel()

#             # 収束条件をgeomeTRIC用に設定
#             options = conv_params if conv_params else {}

#             # 構造最適化の実行（geomeTRICを使用）
#             optimized_mol = optimize(mf, solver='geomeTRIC', options=options)

#             # 最適化後の座標
#             final_geometry = optimized_mol.atom_coords()
#             log_entry["time"] = "End_Time"
#             log_entry["result"] = {"final_geometry": final_geometry.tolist()}
#             save_log(compound_name, log_entry)

#         return final_geometry

#     except Exception as e:
#         log_entry["time"] = "Error_End"
#         log_entry["error"] = str(e)
#         save_log(compound_name, log_entry)
#         raise RuntimeError(f"Optimization failed: {e}")


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
