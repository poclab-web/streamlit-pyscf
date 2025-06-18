import pandas as pd
import rdkit
import rdkit.Chem
import rdkit.Chem.AllChem
from collections import Counter

from logic.calculation import (
    run_geometry_optimization,
    calculate_vibrational_frequencies
)
import logging


def fragment_iterator(smiles, skip_warnings=False, skip_rings=False):

    mol_stereo = enumerate_stereocenters(smiles)
    if (mol_stereo["atom_unassigned"] != 0) or (mol_stereo["bond_unassigned"] != 0):
        logging.warning(f"Molecule {smiles} has undefined stereochemistry")
        if skip_warnings:
            return

    mol = rdkit.Chem.MolFromSmiles(smiles)
    mol = rdkit.Chem.rdmolops.AddHs(mol)
    rdkit.Chem.Kekulize(mol, clearAromaticFlags=True)

    for bond in mol.GetBonds():

        if skip_rings and bond.IsInRing():
            continue

        if bond.GetBondTypeAsDouble() > 1.9999:
            continue

        try:

            # Use RDkit to break the given bond
            mh = rdkit.Chem.RWMol(mol)
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            mh.RemoveBond(a1, a2)

            mh.GetAtomWithIdx(a1).SetNoImplicit(True)
            mh.GetAtomWithIdx(a2).SetNoImplicit(True)

            # Call SanitizeMol to update radicals
            rdkit.Chem.SanitizeMol(mh)

            # Convert the two molecules into a SMILES string
            fragmented_smiles = rdkit.Chem.MolToSmiles(mh)

            # Split fragment and canonicalize
            split_smiles = fragmented_smiles.split(".")
            if len(split_smiles) == 2:
                frag1, frag2 = sorted(split_smiles)
            elif len(split_smiles) == 1:
                frag1, frag2 = split_smiles[0], ""
            else:
                raise ValueError("Too many fragments")

            frag1 = canonicalize_smiles(frag1)
            frag2 = canonicalize_smiles(frag2)

            # Stoichiometry check
            assert (
                count_atom_types(frag1) + count_atom_types(frag2)
            ) == count_atom_types(smiles), "Error with {}; {}; {}".format(
                frag1, frag2, smiles
            )

            # Check introduction of new stereocenters
            is_valid_stereo = check_stereocenters(frag1) and check_stereocenters(frag2)

            yield pd.Series(
                {
                    "molecule": smiles,
                    "bond_index": bond.GetIdx(),
                    "bond_type": get_bond_type(bond),
                    "fragment1": frag1,
                    "fragment2": frag2,
                    "is_valid_stereo": is_valid_stereo,
                }
            )

        except ValueError:
            logging.error(
                "Fragmentation error with {}, bond {}".format(smiles, bond.GetIdx())
            )
            continue

def enumerate_stereocenters(smiles):
    """ Returns a count of both assigned and unassigned stereocenters in the
    given molecule """

    mol = rdkit.Chem.MolFromSmiles(smiles)
    rdkit.Chem.FindPotentialStereoBonds(mol)

    stereocenters = rdkit.Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    stereobonds = [
        bond
        for bond in mol.GetBonds()
        if bond.GetStereo() is not rdkit.Chem.rdchem.BondStereo.STEREONONE
    ]

    atom_assigned = len([center for center in stereocenters if center[1] != "?"])
    atom_unassigned = len([center for center in stereocenters if center[1] == "?"])

    bond_assigned = len(
        [
            bond
            for bond in stereobonds
            if bond.GetStereo() is not rdkit.Chem.rdchem.BondStereo.STEREOANY
        ]
    )
    bond_unassigned = len(
        [
            bond
            for bond in stereobonds
            if bond.GetStereo() is rdkit.Chem.rdchem.BondStereo.STEREOANY
        ]
    )

    return pd.Series(
        {
            "atom_assigned": atom_assigned,
            "atom_unassigned": atom_unassigned,
            "bond_assigned": bond_assigned,
            "bond_unassigned": bond_unassigned,
        }
    )

def canonicalize_smiles(smiles):
    """ Return a consistent SMILES representation for the given molecule """
    mol = rdkit.Chem.MolFromSmiles(smiles)
    return rdkit.Chem.MolToSmiles(mol)

def count_atom_types(smiles):
    """ Return a dictionary of each atom type in the given fragment or molecule
    """
    mol = rdkit.Chem.MolFromSmiles(smiles, sanitize=True)
    mol = rdkit.Chem.rdmolops.AddHs(mol)
    return Counter([atom.GetSymbol() for atom in mol.GetAtoms()])

def check_stereocenters(smiles):
    """Check the given SMILES string to determine whether accurate
    enthalpies can be calculated with the given stereochem information
    """
    stereocenters = enumerate_stereocenters(smiles)
    if stereocenters["bond_unassigned"] > 0:
        return False

    max_unassigned = 1 if stereocenters["atom_assigned"] == 0 else 1
    if stereocenters["atom_unassigned"] <= max_unassigned:
        return True
    else:
        return False

def get_bond_type(bond):
    return "{}-{}".format(
        *tuple(sorted((bond.GetBeginAtom().GetSymbol(), bond.GetEndAtom().GetSymbol())))
    )

def get_fragment_dataframe(smiles_list, **kwargs):
    rows = []
    for smiles in smiles_list:
        for frag in fragment_iterator(smiles, **kwargs):
            rows.append(frag)
    return pd.DataFrame(rows)


from pyscf import gto
from logic.calculation import run_geometry_optimization, calculate_vibrational_frequencies, run_quantum_calculation

def is_single_atom(xyz: str, charge: int, spin: int) -> bool:
    """
    与えられたxyz形式の構造が1原子のみかを判定する。
    PySCFの分子オブジェクトを構築し、原子数を確認する。

    Parameters:
        xyz (str): XYZ形式の原子座標。
        charge (int): 分子全体の電荷。
        spin (int): Nalpha - Nbeta （2S ではない）で指定するスピン。

    Returns:
        bool: 原子数が1つならTrue、それ以外はFalse。
    """
    try:
        mol = gto.M(atom=xyz, basis='sto-3g', charge=charge, spin=spin)
        return mol.natm == 1
    except Exception as e:
        print(f"[ERROR] Failed to parse molecule in is_single_atom(): {e}")
        return False

def compute_neutral_molecule_properties(name, smiles, xyz,conv_params,
                                        theory="B3LYP", basis_set="def2-SVP",
                                        solvent_model=None, eps=None,  maxsteps=100):
    """
    中性分子に対して、構造最適化＋振動計算を行い、Gibbs自由エネルギーなどの情報を辞書形式で返す。
    ただし、1原子系の場合はSCFエネルギーのみを返す。
    """
    if is_single_atom(xyz, charge=0, spin=0):
        print(f"[INFO] {name} is a single-atom molecule. Skipping optimization and frequency.")
        energy, _ = run_quantum_calculation(
            compound_name=name,
            smiles=smiles,
            atom_input=xyz,
            basis_set=basis_set,
            theory=theory,
            charge=0,
            spin=0,
            solvent_model=solvent_model,
            eps=eps,
            symmetry=False
        )
        return {
            'G_tot': energy,
            'frequencies': {},
            'thermo_info': {'E_scf_only': energy}
        }

    # 通常の分子（2原子以上）に対する処理
    xyz_opt = run_geometry_optimization(
        compound_name=name,
        smiles=smiles,
        atom_input=xyz,
        theory=theory,
        basis_set=basis_set,
        charge=0,
        spin=0,
        solvent_model=solvent_model,
        eps=eps,
        conv_params=conv_params, 
        maxsteps=maxsteps
    )

    vib_result = calculate_vibrational_frequencies(
        atom_input=xyz_opt,
        theory=theory,
        basis_set=basis_set,
        charge=0,
        spin=0,
        solvent_model=solvent_model,
        eps=eps,
        compound_name=name,
        smiles=smiles
    )

    thermo = vib_result['thermo_info']
    G_tot = thermo.get('G_tot', None)

    return {
        'G_tot': G_tot,
        'frequencies': vib_result.get('frequencies', {}),
        'thermo_info': thermo
    }


def compute_radical_fragment_properties(name, smiles, xyz, theory, basis_set, charge, spin, conv_params, maxsteps=100, solvent_model=None, eps=None):
    """
    ラジカルフラグメントに対して、構造最適化+振動計算を行い、Gibbs自由エネルギーを含む辞書を返す。
    """

    if is_single_atom(xyz, charge, spin):
        print(f"[INFO] {name} is a single-atom molecule. Skipping optimization and frequency.")
        energy, _ = run_quantum_calculation(
            compound_name=name,
            smiles=smiles,
            atom_input=xyz,
            basis_set=basis_set,
            theory=theory,
            charge=charge,
            spin=spin,
            solvent_model=solvent_model,
            eps=eps,
            symmetry=False
        )
        return {
            'G_tot': energy,
            'frequencies': {},
            'thermo_info': {'E_scf_only': energy}
        }

    print(f"[INFO] Running geometry optimization for {name} with SMILES: {smiles}")
    xyz_opt = run_geometry_optimization(
        compound_name=name,
        smiles=smiles,
        atom_input=xyz,
        theory=theory,
        basis_set=basis_set,
        charge=charge,
        spin=spin, 
        solvent_model=solvent_model,
        eps=eps,
        conv_params=conv_params, 
        maxsteps=maxsteps
    )

    print(f"[INFO] Running vibrational frequency calculation for {name} with SMILES: {smiles}")
    vib_result = calculate_vibrational_frequencies(
        atom_input=xyz_opt,
        theory=theory,
        basis_set=basis_set,
        charge=charge,
        spin=spin,
        solvent_model=solvent_model,
        eps=eps,
        compound_name=name,
        smiles=smiles
    )

    thermo = vib_result['thermo_info']
    G_tot = thermo.get('G_tot', None)

    return {
        'G_tot': G_tot,
        'frequencies': vib_result.get('frequencies', {}),
        'thermo_info': thermo
    }

def calculate_bde_from_gibbs(G_molecule, G_frag1, G_frag2):
    """
    Gibbs自由エネルギーを用いてBDEを計算（単位：kcal/mol）
    """
    bde_au = (G_frag1 + G_frag2) - G_molecule
    return bde_au * 627.509  # au → kcal/mol


if __name__ == "__main__":
    smiles_list = ['CCO', 'c1ccccc1O']  # 例：エタノールとフェノール
    df = get_fragment_dataframe(smiles_list, skip_rings=False, skip_warnings=True)
    print(df)