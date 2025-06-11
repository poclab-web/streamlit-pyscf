"""
pySCFの計算を実行して、実行した結果を出力させる部分
ファイルの入力については、molecule_handler.pyに行わせている。
1. 1点計算
2. 構造最適化計算
3. 振動計算
4. NMR計算
5. 励起状態計算（時間依存密度汚関数理論: TD-DFT）
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

from pyscf import gto, scf, dft, solvent, mp, tools, hessian, geomopt
from pyscf.geomopt.berny_solver import optimize
from pyscf.geomopt import geometric_solver
from pyscf.hessian import thermo, rhf, rks, uhf  # これだけでOK


from logic.molecule_handler import MoleculeHandler
from logic.visualization import generate_cjson
import dftd3.pyscf as d3  # 追加

# 現在の日時を取得してフォーマット
current_time = datetime.now().strftime("%Y%m%d_%H%M%S")

# 理論と基底関数の選択肢
theory_options = ["HF", "MP2", "B3LYP", "B3LYPD3", "PBE", "M06-2X", "B97X-D"]
basis_set_options = [
    "sto-3g",          # 最小基底
    "6-31g",           # 二重ζ基底（DZ）
    "6-31g*",          # +分極関数（重原子にd軌道）
    "6-31g**",         # +分極関数（全原子にd,p軌道）
    "6-31+g**",        # +分散関数（陰イオン・励起状態向け）
    "6-31++g**",       # +分散関数×2（水素にも分散関数）
    "cc-pVDZ",         # correlation-consistent (Dunning系)
    "cc-pVTZ",
    "aug-cc-pVDZ",     # 分散関数付き（augmented）
    "aug-cc-pVTZ",
    "def2-SVP",        # def2 系（Karlsruhe基底）
    "def2-TZVP"
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

    Returns:
        gto.Mole: A PySCF molecule object.
        dict: 溶媒モデル情報（solvent_model, eps）を返す（必要なら）
    """
    mol = gto.M(atom=atom_input, basis=basis_set, charge=charge, spin=spin, symmetry=symmetry)
    mol.unit = 'Angstrom'
    # PCMなどの溶媒モデルはmolではなくSCF計算時にセットする
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

# --- ZPEとGibbsの取得を柔軟に ---
def extract_float_from_dict(d, keys):
    for k in keys:
        v = d.get(k)
        if v is not None:
            if isinstance(v, (list, tuple, np.ndarray)):
                v = v[0]
            try:
                return float(v)
            except Exception:
                continue
    return None

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
        solvent_model (str): Solvent model to be used (e.g., PCM, DDCOSMO).
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

        chkfile_name = get_sequential_filename(directory, theory, basis_set, opt_theory, opt_basis_set, extension="chk")
        molden_file = get_sequential_filename(directory, theory, basis_set, opt_theory, opt_basis_set, extension="molden")
        output_file = get_sequential_filename(directory, theory, basis_set, opt_theory, opt_basis_set, extension="out")

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
            elif theory == "B3LYPD3":
                mf = dft.RKS(mol)
                mf.xc = "B3LYP"
                mf = d3.energy(mf)  # D3分散補正を付加
            else:
                raise ValueError(f"Unsupported theory: {theory}")

            # 溶媒モデルをセット
            if solvent_model is not None:
                if solvent_model.upper() == "PCM":
                    mf = solvent.PCM(mf)
                    if eps is not None:
                        mf.with_solvent.eps = eps
                elif solvent_model.upper() == "DDCOSMO":
                    mf = solvent.ddCOSMO(mf)
                    if eps is not None:
                        mf.with_solvent.eps = eps

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
        elif theory == "B3LYPD3":
            mf = dft.RKS(mol)
            mf.xc = "B3LYP"
            mf = d3.energy(mf)  # D3分散補正を付加
        else:
            raise ValueError(f"Unsupported theory: {theory}")
        
        # 溶媒モデルをセット
        if solvent_model is not None:
            if solvent_model.upper() == "PCM":
                mf = solvent.PCM(mf)
                if eps is not None:
                    mf.with_solvent.eps = eps
            elif solvent_model.upper() == "DDCOSMO":
                mf = solvent.ddCOSMO(mf)
                if eps is not None:
                    mf.with_solvent.eps = eps

        mf.kernel()

        conv_params = conv_params if conv_params else {}
        optimized_mol = optimize(mf, conv_params=conv_params, maxsteps=maxsteps)        

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

        theory_upper = theory.upper()
        if theory_upper == 'HF':
            mf = scf.RHF(mol)
            hess_class = rhf.Hessian
        elif theory_upper in ['B3LYP', 'PBE', 'M06-2X', 'B97X-D']:
            mf = dft.RKS(mol)
            mf.xc = theory.lower()
            hess_class = rks.Hessian
        elif theory_upper == "B3LYPD3":
            mf = dft.RKS(mol)
            mf.xc = "b3lyp"
            mf = d3.energy(mf)  # D3分散補正を付加
            hess_class = rks.Hessian
        else:
            raise ValueError(f"Unsupported theory: {theory}")

        # 溶媒モデルをセット
        if solvent_model is not None:
            if solvent_model.upper() == "PCM":
                mf = solvent.PCM(mf)
                if eps is not None:
                    mf.with_solvent.eps = eps
            elif solvent_model.upper() == "DDCOSMO":
                mf = solvent.ddCOSMO(mf)
                if eps is not None:
                    mf.with_solvent.eps = eps

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


# IP計算
def calculate_ionization_potential(
    mol, theory, basis, optimize_neutral, optimize_radical_cation,
    solvent_model=None, eps=None, compound_name=None, smiles=None, opt_theory=None, opt_basis_set=None
):
    """
    分子のVIP（Vertical/Adiabatic Ionization Potential）を計算する。
    """
    hartree_to_ev = 27.2114
    temperature = 298.15

    # ログ保存用ディレクトリ
    directory = None
    log_entry = None
    if compound_name:
        log_entry = {
            "time": "Start_Time",
            "timestamp": datetime.now().isoformat(),
            "compound": compound_name,
            "smiles": smiles,
            "calculation_type": "ionization_potential",
            "parameters": {
                "basis_set": basis,
                "theory": theory,
                "solvent_model": solvent_model,
                "dielectric": eps,
                "optimize_neutral": optimize_neutral,
                "optimize_radical_cation": optimize_radical_cation
            }
        }
        directory = save_log(compound_name, log_entry)

    # chkファイル名生成はdirectoryがある場合のみ
    if directory is not None:
        chkfile_neutral = get_sequential_filename(
            directory, theory, basis, opt_theory, opt_basis_set, extension="chk"
        )
        chkfile_cation = get_sequential_filename(
            directory, theory, basis, opt_theory, opt_basis_set, extension="chk",
            purpose=f"{solvent_model}_cation" if solvent_model else "cation"
        )
    else:
        chkfile_neutral = None
        chkfile_cation = None

    try:
        # --- 理論ごとのSCFオブジェクト生成 ---
        def make_mf(mol_obj, charge=0, spin=0, is_cation=False):
            theory_upper = theory.upper()
            if theory_upper == "HF":
                mf = scf.RHF(mol_obj) if not is_cation else scf.UHF(mol_obj)
                hess_class = rhf.Hessian if not is_cation else uhf.Hessian
            elif theory_upper == "MP2":
                mf = scf.RHF(mol_obj).run() if not is_cation else scf.UHF(mol_obj).run()
                mf = mp.MP2(mf).run()
                hess_class = rhf.Hessian if not is_cation else uhf.Hessian
            elif theory_upper in ["B3LYP", "PBE", "M06-2X", "B97X-D"]:
                if not is_cation:
                    mf = dft.RKS(mol_obj)
                    mf.xc = theory
                    hess_class = rks.Hessian
                else:
                    mf = dft.UKS(mol_obj)
                    mf.xc = theory
                    hess_class = uhf.Hessian
            elif theory_upper == "B3LYPD3":
                if not is_cation:
                    mf = dft.RKS(mol_obj)
                    mf.xc = "B3LYP"
                    mf = d3.energy(mf)
                    hess_class = rks.Hessian
                else:
                    mf = dft.UKS(mol_obj)
                    mf.xc = "B3LYP"
                    mf = d3.energy(mf)
                    hess_class = uhf.Hessian
            else:
                raise ValueError(f"Unsupported theory: {theory}")
            # 溶媒モデル
            if solvent_model is not None:
                if solvent_model.upper() == "PCM":
                    mf = solvent.PCM(mf)
                    if eps is not None:
                        mf.with_solvent.eps = eps
                elif solvent_model.upper() == "DDCOSMO":
                    mf = solvent.ddCOSMO(mf)
                    if eps is not None:
                        mf.with_solvent.eps = eps
            return mf, hess_class

        # 1. 中性分子の最適化
        mf_neutral, hess_class = make_mf(mol, charge=mol.charge, spin=mol.spin, is_cation=False)
        if chkfile_neutral:
            mf_neutral.chkfile = chkfile_neutral

        if optimize_neutral:
            mol_opt = geometric_solver.optimize(mf_neutral)
            mf_neutral, hess_class = make_mf(mol_opt, charge=mol_opt.charge, spin=mol_opt.spin, is_cation=False)
            if hasattr(mf_neutral, "xc"):
                mf_neutral.xc = theory if theory.upper() != "B3LYPD3" else "B3LYP"
            if chkfile_neutral:
                mf_neutral.chkfile = chkfile_neutral
            mf_neutral.kernel()
            mol_for_cation = mol_opt
        else:
            mf_neutral.kernel()
            mol_for_cation = mol

        e_neutral = mf_neutral.e_tot

        # 振動解析・ZPE・Gibbs
        hess_neutral = hess_class(mf_neutral).kernel()
        vib_info_neutral = thermo.harmonic_analysis(mol_for_cation, hess_neutral)
        freqs_neutral = vib_info_neutral.get("freq_au", [])
        thermo_info_neutral = thermo.thermo(mf_neutral, freqs_neutral, temperature=temperature)

        zpe_neutral = extract_float_from_dict(
            thermo_info_neutral,
            ['Zero Point Energy', 'ZPE', 'ZeroPointEnergy']
        )
        gibbs_neutral = extract_float_from_dict(
            thermo_info_neutral,
            ['Gibbs Free Energy', 'GibbsFreeEnergy', 'Gibbs']
        )

        # 2. カチオン分子の準備
        mol_cation = mol_for_cation.copy()
        mol_cation.charge += 1
        mol_cation.spin = 1
        mol_cation.build()

        # カチオンの構造最適化 or 中性構造流用
        mf_cation, hess_class_cation = make_mf(mol_cation, charge=mol_cation.charge, spin=mol_cation.spin, is_cation=True)
        if optimize_radical_cation:
            mol_cation_opt = geometric_solver.optimize(mf_cation)
            mf_cation, hess_class_cation = make_mf(mol_cation_opt, charge=mol_cation_opt.charge, spin=mol_cation_opt.spin, is_cation=True)
            if hasattr(mf_cation, "xc"):
                mf_cation.xc = theory if theory.upper() != "B3LYPD3" else "B3LYP"
            if chkfile_cation:
                mf_cation.chkfile = chkfile_cation
            mf_cation.kernel()
            mol_cation = mol_cation_opt
        else:
            if chkfile_cation:
                mf_cation.chkfile = chkfile_cation
            mf_cation.kernel()

        e_cation = mf_cation.e_tot

        # 振動解析・ZPE・Gibbs（カチオン）
        hess_cation = hess_class_cation(mf_cation).kernel()
        vib_info_cation = thermo.harmonic_analysis(mol_cation, hess_cation)
        freqs_cation = vib_info_cation.get("freq_au", [])
        thermo_info_cation = thermo.thermo(mf_cation, freqs_cation, temperature=temperature)
        zpe_cation = extract_float_from_dict(
            thermo_info_cation,
            ['Zero Point Energy', 'ZPE', 'ZeroPointEnergy']
        )
        gibbs_cation = extract_float_from_dict(
            thermo_info_cation,
            ['Gibbs Free Energy', 'GibbsFreeEnergy', 'Gibbs']
        )

        # --- 型チェックを追加 ---
        def safe_float(val):
            try:
                return float(val)
            except Exception:
                return None

        zpe_neutral_f = safe_float(zpe_neutral)
        zpe_cation_f = safe_float(zpe_cation)
        gibbs_neutral_f = safe_float(gibbs_neutral)
        gibbs_cation_f = safe_float(gibbs_cation)

        # VIP計算
        vip_e = (e_cation - e_neutral) * hartree_to_ev
        vip_zpe = None
        vip_gibbs = None
        if zpe_neutral_f is not None and zpe_cation_f is not None:
            vip_zpe = ((e_cation + zpe_cation_f) - (e_neutral + zpe_neutral_f)) * hartree_to_ev
        if gibbs_neutral_f is not None and gibbs_cation_f is not None:
            vip_gibbs = (gibbs_cation_f - gibbs_neutral_f) * hartree_to_ev

        result = {
            "VIP_electronic": vip_e,
            "VIP_ZPE": vip_zpe,
            "VIP_Gibbs": vip_gibbs,
            "neutral": {
                "energy": e_neutral,
                "zpe": zpe_neutral_f,
                "gibbs": gibbs_neutral_f,
                "freqs": freqs_neutral,
                "thermo": thermo_info_neutral,
                "chkfile": chkfile_neutral
            },
            "cation": {
                "energy": e_cation,
                "zpe": zpe_cation_f,
                "gibbs": gibbs_cation_f,
                "freqs": freqs_cation,
                "thermo": thermo_info_cation,
                "chkfile": chkfile_cation
            }
        }

        if compound_name:
            log_entry["time"] = "End_Time"
            log_entry["result"] = result
            save_log(compound_name, log_entry)

        return result

    except Exception as e:
        if compound_name and log_entry is not None:
            log_entry["time"] = "Error_End"
            log_entry["error"] = str(e)
            save_log(compound_name, log_entry)
        raise RuntimeError(f"Ionization potential calculation failed: {e}")
# ...existing code...
def calculate_solvation_energy(
    mol, theory, basis,
    optimize_gas=False, optimize_solvent=False,
    solvent_model="PCM", eps=78.39,
    compound_name=None, smiles=None, opt_theory=None, opt_basis_set=None,
    solvent_name=None  
):
    """
    溶媒和エネルギー（ΔG_solv）を計算する。
    """

    HARTREE2KCAL = 627.5094740631
    temperature = 298.15

    # 入力がxyz文字列の場合はmolを作成
    if isinstance(mol, str):
        mol = gto.M(atom=mol, basis=basis)
        mol.unit = 'Angstrom'

    # ログ保存用ディレクトリ
    directory = None
    if compound_name:
        log_entry = {
            "time": "Start_Time",
            "timestamp": datetime.now().isoformat(),
            "compound": compound_name,
            "smiles": smiles,
            "calculation_type": "solvation_energy",
            "parameters": {
                "basis_set": basis,
                "theory": theory,
                "solvent_model": solvent_model,
                "dielectric": eps,
                "optimize_gas": optimize_gas,
                "optimize_solvent": optimize_solvent
            }
        }
        directory = save_log(compound_name, log_entry)

    # chkファイル名生成はdirectoryがある場合のみ
    if directory is not None:
        chkfile_gas = get_sequential_filename(
            directory, theory, basis, opt_theory, opt_basis_set, extension="chk"
        )
        chkfile_solvent = get_sequential_filename(
            directory, theory, basis, opt_theory, opt_basis_set, extension="chk",
            purpose=solvent_name if solvent_name else solvent_model
        )
    else:
        chkfile_gas = None
        chkfile_solvent = None

    # --- 理論ごとのSCFオブジェクト生成 ---
    def make_mf(mol_obj, is_solvent=False):
        theory_upper = theory.upper()
        if theory_upper == "HF":
            mf = scf.RHF(mol_obj)
            hess_class = rhf.Hessian
        elif theory_upper == "MP2":
            mf = scf.RHF(mol_obj).run()
            mf = mp.MP2(mf).run()
            hess_class = rhf.Hessian
        elif theory_upper in ["B3LYP", "PBE", "M06-2X", "B97X-D"]:
            mf = dft.RKS(mol_obj)
            mf.xc = theory
            hess_class = rks.Hessian
        elif theory_upper == "B3LYPD3":
            mf = dft.RKS(mol_obj)
            mf.xc = "B3LYP"
            import dftd3.pyscf as d3
            mf = d3.energy(mf)
            hess_class = rks.Hessian
        else:
            raise ValueError(f"Unsupported theory: {theory}")
        if is_solvent:
            if solvent_model.upper() == "PCM":
                mf = solvent.PCM(mf)
                if eps is not None:
                    mf.with_solvent.eps = eps
            elif solvent_model.upper() == "DDCOSMO":
                mf = solvent.ddCOSMO(mf)
                if eps is not None:
                    mf.with_solvent.eps = eps
        return mf, hess_class

    # 気相構造最適化
    if optimize_gas:
        mf_gas, hess_class = make_mf(mol, is_solvent=False)
        mol_gas_opt = geometric_solver.optimize(mf_gas)
        mf_gas, hess_class = make_mf(mol_gas_opt, is_solvent=False)
        if hasattr(mf_gas, "xc"):
            mf_gas.xc = theory if theory.upper() != "B3LYPD3" else "B3LYP"
        if chkfile_gas:
            mf_gas.chkfile = chkfile_gas
        mf_gas.kernel()
        mol_gas = mol_gas_opt
    else:
        mf_gas, hess_class = make_mf(mol, is_solvent=False)
        if chkfile_gas:
            mf_gas.chkfile = chkfile_gas
        mf_gas.kernel()
        mol_gas = mol

    gas_phase_energy = mf_gas.e_tot

    # 気相ZPE・Gibbs
    hess_gas = hess_class(mf_gas).kernel()
    vib_info_gas = thermo.harmonic_analysis(mol_gas, hess_gas)
    freqs_gas = vib_info_gas.get("freq_au", [])
    thermo_info_gas = thermo.thermo(mf_gas, freqs_gas, temperature=temperature)
    zpe_gas = extract_float_from_dict(thermo_info_gas, ['Zero Point Energy', 'ZPE', 'ZeroPointEnergy'])
    gibbs_gas = extract_float_from_dict(thermo_info_gas, ['Gibbs Free Energy', 'GibbsFreeEnergy', 'Gibbs'])

    # 溶媒中構造最適化
    if optimize_solvent:
        mf_solvent, hess_class_sol = make_mf(mol_gas, is_solvent=True)
        mol_solvent_opt = geometric_solver.optimize(mf_solvent)
        mf_solvent, hess_class_sol = make_mf(mol_solvent_opt, is_solvent=True)
        if hasattr(mf_solvent, "xc"):
            mf_solvent.xc = theory if theory.upper() != "B3LYPD3" else "B3LYP"
        if chkfile_solvent:
            mf_solvent.chkfile = chkfile_solvent
        mf_solvent.kernel()
        mol_solvent = mol_solvent_opt
    else:
        mf_solvent, hess_class_sol = make_mf(mol_gas, is_solvent=True)
        if chkfile_solvent:
            mf_solvent.chkfile = chkfile_solvent
        mf_solvent.kernel()
        mol_solvent = mol_gas

    solvent_phase_energy = mf_solvent.e_tot

    # 溶媒中ZPE・Gibbs
    hess_solvent = hess_class_sol(mf_solvent).kernel()
    vib_info_solvent = thermo.harmonic_analysis(mol_solvent, hess_solvent)
    freqs_solvent = vib_info_solvent.get("freq_au", [])
    thermo_info_solvent = thermo.thermo(mf_solvent, freqs_solvent, temperature=temperature)
    zpe_solvent = extract_float_from_dict(thermo_info_solvent, ['Zero Point Energy', 'ZPE', 'ZeroPointEnergy'])
    gibbs_solvent = extract_float_from_dict(thermo_info_solvent, ['Gibbs Free Energy', 'GibbsFreeEnergy', 'Gibbs'])

    # ΔG_solv計算（Gibbs Free Energyベース、なければエネルギー差）
    if gibbs_gas is not None and gibbs_solvent is not None:
        delta_g_solv = (gibbs_solvent - gibbs_gas) * HARTREE2KCAL
    else:
        delta_g_solv = (solvent_phase_energy - gas_phase_energy) * HARTREE2KCAL

    # ログ保存
    if compound_name:
        log_entry["time"] = "End_Time"
        log_entry["result"] = {
            "solvation_energy": delta_g_solv,
            "gas_phase_energy": gas_phase_energy,
            "solvent_phase_energy": solvent_phase_energy,
            "zpe_gas": zpe_gas,
            "zpe_solvent": zpe_solvent,
            "gibbs_gas": gibbs_gas,
            "gibbs_solvent": gibbs_solvent,
            "chkfile_gas": os.path.basename(chkfile_gas) if chkfile_gas else None,
            "chkfile_solvent": os.path.basename(chkfile_solvent) if chkfile_solvent else None,
        }
        save_log(compound_name, log_entry)

    return {
        "solvation_energy": delta_g_solv,
        "gas_phase_energy": gas_phase_energy,
        "solvent_phase_energy": solvent_phase_energy,
        "zpe_gas": zpe_gas,
        "zpe_solvent": zpe_solvent,
        "gibbs_gas": gibbs_gas,
        "gibbs_solvent": gibbs_solvent,
        "chkfile_gas": chkfile_gas,
        "chkfile_solvent": chkfile_solvent,
    }



