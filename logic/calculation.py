"""
pySCFの計算を実行して、実行した結果を出力させる部分
ファイルの入力については、molecule_handler.pyに行わせている。
1. 1点計算
2. 構造最適化計算
3. 振動計算
4. NMR計算
5. 励起状態計算（時間依存密度汎関数理論: TD-DFT）
6. 分子の応答性（極性率や双極子モーメント）
出てきたファイルの解析については、output_handler.pyで数値の変換し、visualization.pyでグラフなどに変換する。
"""

import os
import json
from datetime import datetime
import tempfile
import sys
import re
import csv

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

from pyscf import gto, scf, dft, solvent, tools, hessian
from pyscf.geomopt.berny_solver import optimize
from pyscf.hessian import thermo, rhf, rks
from pyscf.prop import nmr
from pyscf.prop.polarizability import rhf
from pyscf.prop.polarizability.rhf import polarizability as get_polarizability

from logic.molecule_handler import MoleculeHandler
from logic.visualization import generate_cjson


# 現在の日時を取得してフォーマット
current_time = datetime.now().strftime("%Y%m%d_%H%M%S")

# 理論と基底関数の選択肢
theory_options = ["HF", "B3LYP", "PBE", "M06-2X", "B97X-D"]
basis_set_options = [
    "sto-3g", "6-31g", "6-31g*", "6-31+g(d,p)", "cc-pVDZ",
    "cc-pVTZ", "aug-cc-pVDZ", "aug-cc-pVTZ", "def2-SVP", "def2-TZVP"
]
hartree_to_cm1 = 219474.63  # 1 Hartree = 219474.63 cm^-1

# 定数を自前で定義
HARTREE2CM = 219474.6313705
HARTREE2KCAL = 627.5094740631

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

    # JSONデータを配列形式で保存
    if os.path.exists(log_file):
        # 既存のデータを読み込む
        with open(log_file, 'r') as f:
            try:
                existing_data = json.load(f)
            except json.JSONDecodeError:
                existing_data = []
    else:
        existing_data = []

    # 新しいデータを追加
    existing_data.append(data)

    # ファイルに書き込む
    with open(log_file, 'w') as f:
        json.dump(existing_data, f, indent=4)

    return directory  # ディレクトリパスを返す


def get_sequential_filename(directory, theory, basis_set, opt_theory="none", opt_basis_set="none", extension="molden", purpose=None):
    """
    Generate a sequential filename to avoid overwriting existing files.
    Optionally include a purpose (e.g., 'nmr', 'vib') in the filename.
    """
    opt_theory = opt_theory or "none"
    opt_basis_set = opt_basis_set or "none"
    # 用途（purpose）があればファイル名に追加
    if purpose:
        filename_base = f"{theory}_{basis_set}__{opt_theory}_{opt_basis_set}_{purpose}"
    else:
        filename_base = f"{theory}_{basis_set}__{opt_theory}_{opt_basis_set}"

    pattern = re.compile(rf"{re.escape(filename_base)}_(\d+)\.{re.escape(extension)}$")
    existing_files = os.listdir(directory)
    max_index = 0
    for file in existing_files:
        match = pattern.match(file)
        if match:
            max_index = max(max_index, int(match.group(1)))
    new_index = max_index + 1
    new_filename = f"{filename_base}_{new_index}.{extension}"
    return os.path.join(directory, new_filename)

# 分子の設定
def setup_molecule(atom_input, basis_set, charge=0, spin=0, solvent_model=None, eps=None, symmetry=False):
    """
    Set up a PySCF molecule with specified parameters.

    Parameters:
        atom_input (str): Atomic coordinates in XYZ format.
        basis_set (str): Basis set for the calculation.
        charge (int): Molecular charge.
        spin (int): Spin multiplicity (2S + 1).
        solvent_model (str): Solvent model to be used (e.g., PCM).
        eps (float): Dielectric constant for solvent.
        symmetry (bool): Whether to consider molecular symmetry.

    Returns:
        gto.Mole: A PySCF molecule object.
    """
    mol = gto.M(atom=atom_input, basis=basis_set, charge=charge, spin=spin, symmetry=symmetry)
    mol.unit = 'Angstrom'
    
    if solvent_model == "PCM":
        pcm = solvent.PCM(mol)
        if eps is not None:
            pcm.eps = eps
        mol = pcm

    return mol

def extract_mol_mf_params(mol, mf):
    """
    mol, mfオブジェクトからログ用のパラメータを抽出するユーティリティ関数
    """
    params = {
        "atom_input": mol.atom,  # 原子座標情報
        "basis_set": mol.basis if isinstance(mol.basis, str) else ",".join(sorted(mol.basis.keys())),
        "theory": getattr(mf, "xc", "HF") if hasattr(mf, "xc") else "HF",
        "charge": mol.charge,
        "spin": mol.spin,
        "solvent_model": getattr(mol, "solvent", None),
        "dielectric": getattr(mol, "eps", None),
        "symmetry": mol.symmetry
    }
    return params

# 1. 1点計算
def run_quantum_calculation(compound_name, smiles, atom_input, basis_set, theory, opt_theory=None, opt_basis_set=None, charge=0, spin=0, solvent_model=None, eps=None, symmetry=False):
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
        symmetry (bool): Whether to consider molecular symmetry.

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
            "dielectric": eps,
            "symmetry": symmetry
        }
    }
    try:
        directory = save_log(compound_name, log_entry)
        mol = setup_molecule(atom_input, basis_set, charge, spin, solvent_model, eps, symmetry)

        # ファイル名を get_sequential_filename を使って生成
        chkfile_name = get_sequential_filename(directory, theory, basis_set, opt_theory, opt_basis_set, extension="chk")
        molden_file = get_sequential_filename(directory, theory, basis_set, opt_theory, opt_basis_set, extension="molden")
        output_file = get_sequential_filename(directory, theory, basis_set, opt_theory, opt_basis_set, extension="out")

        # ファイル名情報をログに追加
        log_entry["file_info"] = {
            "chkfile": os.path.basename(chkfile_name),
            "molden_file": os.path.basename(molden_file),
            "output_file": os.path.basename(output_file)
        }

        with open(output_file, 'w') as f:
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
            log_entry["time"] = "Normal_End_Time"
            log_entry["result"] = {
                "energy": energy,
                "converged": mf.converged,
                "chkfile": os.path.abspath(chkfile_name)
            }
            
            # molden.dump_scfを使用してエネルギー情報を含む出力を保存
            tools.molden.dump_scf(mf, molden_file)

            save_log(compound_name, log_entry)

        return energy, molden_file
    except Exception as e:
        log_entry["time"] = "Error_End_Time"
        log_entry["error"] = str(e)
        save_log(compound_name, log_entry)
        raise RuntimeError(f"Calculation failed: {e}")


# 2. 構造最適化計算

def run_geometry_optimization(compound_name, smiles, atom_input, basis_set, theory, charge=0, spin=0, solvent_model=None, eps=None, symmetry=False, conv_params=None, maxsteps=100):
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
            "symmetry": symmetry,
            "conv_params": conv_params
        }
    }

    try:
        directory = save_log(compound_name, log_entry)
        mol = setup_molecule(atom_input, basis_set, charge, spin, solvent_model, eps, symmetry)

        if theory == "HF":
            mf = scf.RHF(mol)
        elif theory in ["B3LYP", "PBE", "M06-2X", "B97X-D"]:
            mf = dft.RKS(mol)
            mf.xc = theory
        else:
            raise ValueError(f"Unsupported theory: {theory}")

        # mf.chkfile = chkfile_name
        mf.kernel()

        conv_params = conv_params if conv_params else {}
        optimized_mol = optimize(mf, conv_params=conv_params, maxsteps=maxsteps)        
        # optimized_mol = optimize(mf, solver='geomeTRIC', options=options)

        bohr_to_angstrom = 0.52917721092

        final_geometry = "\n".join(
            f"{atom[0]} {coord[0] * bohr_to_angstrom:.6f} {coord[1] * bohr_to_angstrom:.6f} {coord[2] * bohr_to_angstrom:.6f}"
            for atom, coord in zip(optimized_mol.atom, optimized_mol.atom_coords())
            )

        log_entry["time"] = "End_Time"
        log_entry["result"] = {"final_geometry": final_geometry}
        save_log(compound_name, log_entry)

        return final_geometry

    except Exception as e:
        log_entry["time"] = "Error_End"
        log_entry["error"] = str(e)
        save_log(compound_name, log_entry)
        raise RuntimeError(f"Optimization failed: {e}")


# 3. 振動計算

def calculate_vibrational_frequencies(atom_input: str, theory: str = 'B3LYP', basis_set: str = 'cc-pVDZ',opt_theory=None, opt_basis_set=None, 
                                       charge=0, spin=0, solvent_model=None, eps=None, symmetry=False,
                                       temperature=298.15, compound_name=None, smiles=None):
    """
    Calculate vibrational frequencies and thermodynamic properties using PySCF.
    """
    log_entry = {
        "time": "Start_Time",
        "timestamp": datetime.now().isoformat(),
        "compound": compound_name,
        "smiles": smiles,
        "calculation_type": "vibrational_frequency",
        "parameters": {
            "atom_input": atom_input,
            "basis_set": basis_set,
            "theory": theory,
            "charge": charge,
            "spin": spin,
            "solvent_model": solvent_model,
            "dielectric": eps,
            "symmetry": symmetry,
            "temperature": temperature
        }
    }

    try:
        directory = None
        if compound_name:
            directory = save_log(compound_name, log_entry)

        chkfile_name = get_sequential_filename(directory, theory, basis_set, opt_theory, opt_basis_set, extension="chk")
        molden_file = get_sequential_filename(directory, theory, basis_set, opt_theory, opt_basis_set, extension="molden")
        output_file = get_sequential_filename(directory, theory, basis_set, opt_theory, opt_basis_set, extension="out")
        cjson_file = get_sequential_filename(directory, theory, basis_set, opt_theory, opt_basis_set, extension="cjson")

        mol = setup_molecule(atom_input, basis_set, charge, spin, solvent_model, eps, symmetry)

        theory = theory.upper()
        if theory == 'HF':
            mf = scf.RHF(mol)
            hess_class = rhf.Hessian
        elif theory in ['B3LYP', 'PBE', 'M06-2X', 'B97X-D']:
            mf = dft.RKS(mol)
            mf.xc = theory.lower()
            hess_class = rks.Hessian
        else:
            raise ValueError(f"Unsupported theory: {theory}")

        mf.chkfile = chkfile_name
        mf.kernel()
        
        # ヘシアン計算と振動解析
        hess_calc = hess_class(mf)
        hess = hess_calc.kernel()
        vib_info = thermo.harmonic_analysis(mf.mol, hess)  # ← 修正
        # norm_mode（固有ベクトル）があれば modes として追加
        modes = vib_info.get('norm_mode', None)

        # Thermochemistry analysis at 298.15 K and 1 atmospheric pressure
        freq_info = thermo.harmonic_analysis(mf.mol, hess)  # ← 修正
        thermo_info = thermo.thermo(mf, freq_info['freq_au'], temperature=temperature)

        # --- ここからfreq情報をoutfileに出力 ---
        with open(output_file, 'w') as f:
            f.write("Vibrational Frequency Analysis\n")
            f.write("====================================\n")
            if "freq_wavenumber" in freq_info:
                f.write("Frequencies (cm^-1):\n")
                for i, freq in enumerate(freq_info["freq_wavenumber"]):
                    f.write(f"  Mode {i+1}: {freq:.2f} cm^-1\n")
            else:
                f.write("No frequency information found.\n")
            f.write("\n")
            if "IR_intensity" in freq_info:
                f.write("IR Intensities:\n")
                for i, inten in enumerate(freq_info["IR_intensity"]):
                    f.write(f"  Mode {i+1}: {inten:.4f}\n")
            f.write("\n")
            # --- thermo_infoの出力 ---
            f.write("Thermodynamic Properties (thermo_info):\n")
            for key, value in thermo_info.items():
                if isinstance(value, (float, int)):
                    f.write(f"{key}: {value}\n")
                elif isinstance(value, (list, tuple, np.ndarray)):
                    # 値と単位のペアの場合は値だけ出力
                    if len(value) == 2 and isinstance(value[0], (float, int)):
                        f.write(f"{key}: {value[0]} {value[1]}\n")
                    else:
                        # それ以外は文字列化して出力
                        f.write(f"{key}: {str(value)}\n")
                else:
                    f.write(f"{key}: {value}\n")
            f.write("\n")


        log_entry["time"] = "End_Time"
        if compound_name:
            save_log(compound_name, log_entry)

        # molden.dump_scfを使用してエネルギー情報を含む出力を保存
        tools.molden.dump_scf(mf, molden_file)

        # generate_cjsonの戻り値をファイルに保存
        cjson_data = generate_cjson(mol, freq_info["freq_wavenumber"], modes)
        with open(cjson_file, 'w') as f:
            f.write(cjson_data)

        # ここでmolとmodesを返す
        return {
            'frequencies': freq_info,
            'thermo_info': thermo_info,
            'mol': mol,
            'modes': modes
        }

    except Exception as e:
        log_entry["time"] = "Error_End"
        log_entry["error"] = str(e)
        if compound_name:
            save_log(compound_name, log_entry)
        raise RuntimeError(f"Vibrational frequency calculation failed: {e}")

# 4. NMR計算（化学シフト）
def calc_nmr_and_shift(mf, mol, target_element=None, basis_set=None, directory=None, theory="HF"):
    """
    NMRシールド値のみ計算（Jカップリングなし, 指定元素のみ）
    basis_setがmolの基底関数と異なる場合は再生成して再計算
    """

    # PySCFのmolからxyz座標を取得
    xyz_str = mol.atom  
    # MoleculeHandlerでRDKit Molに変換
    handler = MoleculeHandler(xyz_str, input_type="xyz")
    compound_name = Chem.MolToInchiKey(handler.mol)
    smiles = Chem.MolToSmiles(handler.mol)

    # mol, mfからパラメータを抽出
    log_params = extract_mol_mf_params(mol, mf)
    # theory, basis_set, target_elementなど引数で上書きしたい場合はここで上書き
    if basis_set is not None:
        log_params["basis_set"] = basis_set

    log_entry = {
        "time": "Start_Time",
        "timestamp": datetime.now().isoformat(),
        "compound": compound_name,
        "smiles": smiles,
        "calculation_type": "NMR",
        "parameters": log_params
    }

    directory = save_log(compound_name, log_entry)

    # molの基底関数と引数のbasis_setが異なる場合は再生成
    if basis_set is not None:
        theory="HF"  # NMR計算はHFで行う
        mol_basis = mol.basis
        # mol.basisはdict型なので、比較用にstr化
        if isinstance(mol_basis, dict):
            mol_basis_str = ",".join(sorted(mol_basis.keys()))
        else:
            mol_basis_str = str(mol_basis)
        if basis_set != mol_basis_str:
            # 元の座標・電荷・スピン・対称性を使って新しいmolを作成
            mol = gto.M(
                atom=mol.atom,
                charge=mol.charge,
                spin=mol.spin,
                symmetry=mol.symmetry,
                basis=basis_set
            )
            mf = scf.RHF(mol)
            mf.kernel()

    # NMRシールド値保存用ファイル名を生成
    nmr_calc_result_name = None

    if directory and theory and basis_set:
        nmr_calc_result_name = get_sequential_filename(directory, theory, basis_set, extension="csv", purpose="nmr")

    nmr_calculator = nmr.RHF(mf)
    nmr_tensors = nmr_calculator.kernel()

    results = []
    for i, tensor in enumerate(nmr_tensors):
        symbol = mol.atom_symbol(i)
        if target_element is not None and symbol != target_element:
            continue
        shield = tensor[0, 0]
        results.append(
            {
                "Atom Index": i,
                "Element": symbol,
                "Theory": theory if theory else "",
                "Basis Set": basis_set if basis_set else "",
                "NMR Shielding": shield
            }
        )

    # --- NMRシールド値をCSVで保存 ---
    print("nmr_calc_result_name:", nmr_calc_result_name)
    if nmr_calc_result_name is not None:
        print("directory exists:", os.path.exists(os.path.dirname(nmr_calc_result_name)))
        with open(nmr_calc_result_name, "w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=["Atom Index", "Element", "Theory", "Basis Set", "NMR Shielding"])
            writer.writeheader()
            for row in results:
                writer.writerow(row)

    log_entry["time"] = "End_Time"
    log_entry["result"] = {
        "nmr_shielding": results,
    }

    if compound_name:
        save_log(compound_name, log_entry)

    return results, nmr_calc_result_name

# 5. 分子の応答性（極性率や双極子モーメント）
def compute_electric_properties(mol, basis_set='6-31g', density_g_cm3=1.0, slope = 278.26, intercept = -299.10):

    # PySCFのmolからxyz座標を取得
    xyz_str = mol.atom  
    # MoleculeHandlerでRDKit Molに変換
    handler = MoleculeHandler(xyz_str, input_type="xyz")
    compound_name = Chem.MolToInchiKey(handler.mol)
    smiles = Chem.MolToSmiles(handler.mol)

    # 分子量（g/mol）をRDKitから取得
    mol_weight = rdMolDescriptors.CalcExactMolWt(handler.mol)

    # PySCFで分子構築
    mol = gto.Mole(
                atom=mol.atom,
                charge=mol.charge,
                spin=mol.spin,
                symmetry=mol.symmetry,
                basis=basis_set
                )
    # 単位をAngstromに設定
    mol.unit = 'Angstrom'
    mol.build()

    # SCF計算
    mf = scf.RHF(mol).run()

    # 双極子モーメント（Debye）
    dipole = mf.dip_moment()
    dipole_norm = np.linalg.norm(dipole)

    # 分極率テンソル（a.u.）
    pol = rhf.Polarizability(mf)
    alpha_tensor = get_polarizability(pol)
    alpha_iso = np.trace(alpha_tensor) / 3

    # 誘電率計算
    NA = 6.022e23
    N_cm3 = (density_g_cm3 / mol_weight) * NA
    alpha_cm3 = alpha_iso * 0.1482e-24
    eps_calc = 1 + (4 * np.pi / 3) * N_cm3 * alpha_cm3

    # ✅ 補正式（回帰式）による実験値推定

    eps_pred = slope * eps_calc + intercept

    return {
        "smiles": smiles,
        "dipole_moment": dipole_norm,
        "polarizability": alpha_iso,
        "dielectric_constant_calc": eps_calc,
        "dielectric_constant_pred": eps_pred,
    }


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


# 6. 励起状態計算（時間依存密度汎関数理論：TD-DFT）
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


# 7. 溶媒効果を考慮した計算

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


