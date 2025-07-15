from logic.calculation import (
    run_geometry_optimization,
    calculate_vibrational_frequencies
)
from pyscf import tdscf  # 追加

def calculate_excited_state(name, smiles, xyz, conv_params,
                            theory="B3LYP", basis_set="def2-SVP",
                            charge=0, spin=0,
                            solvent_model=None, eps=None, maxsteps=100,
                            nstates=10, excited_spin="singlet", tda=False):

    # 通常の分子（2原子以上）に対する処理
    xyz_opt, mf = run_geometry_optimization(
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


    # UVスペクトル（遷移エネルギー）の計算（TDDFT）
    if tda:
        td = tdscf.TDA(mf)
    else:
        td = tdscf.TDHF(mf)
    td.nstates = nstates

    # singlet/tripletの分岐
    if excited_spin == "triplet":
        td.singlet = False
    else:
        td.singlet = True

    excitation_energies = td.kernel()

    result = {
        'td': td,
        'excitation_energies': excitation_energies,
        'xyz_opt': xyz_opt,
        'mf': mf,
        'G_tot': G_tot,
        'frequencies': vib_result.get('frequencies', {}),
        'thermo_info': thermo
    }
    return result



if __name__== "__main__":
    # テスト用のコード
    name = "test_molecule"
    smiles = "CCO"
    xyz = "C 0.000000 0.000000 0.000000\nC 1.529000 0.000000 0.000000\nO 2.329000 0.000000 0.000000"
    conv_params = {"max_cycle": 50, "tol": 1e-6}
    
    td, energies, xyz_opt, mf = calculate_excited_state(
        name, smiles, xyz, conv_params,
        theory="HF", basis_set="sto-3g",
        charge=0, spin=0,
        solvent_model=None, eps=None, maxsteps=100,
        nstates=10, excited_spin="singlet", tda=True
    )
    
    print("Excitation Energies:", energies)
    print("Optimized Geometry:", xyz_opt)
    print("Molecular Orbital Energies:", mf.mo_energy)
    print("Molecular Orbital Coefficients:", mf.mo_coeff)