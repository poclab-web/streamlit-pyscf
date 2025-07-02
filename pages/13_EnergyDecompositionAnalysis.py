"""
個々のフラグメント（分子 A, 分子 B）のエネルギーを計算。
フラグメントが結合した複合体（A + B）のエネルギーを計算。
相互作用エネルギーを計算：
"""

import streamlit as st
# from utils.module import load_css

# # カスタムCSSを適用
# load_css("config/styles.css")



from pyscf import gto, scf
import numpy as np

def get_hf_energy_decomposition(atom_str, basis="sto-3g"):
    mol = gto.Mole()
    mol.atom = atom_str
    mol.basis = basis
    mol.build()
    
    mf = scf.RHF(mol)
    e_total = mf.kernel()
    
    dm = mf.make_rdm1()
    hcore = mf.get_hcore()
    vj, vk = mf.get_jk()
    
    e_nuc = mol.energy_nuc()
    e_core = np.einsum('ij,ji', dm, hcore)
    e_J = np.einsum('ij,ji', dm, vj) * 0.5
    e_K = np.einsum('ij,ji', dm, vk) * 0.5
    e_elec = e_core + e_J - e_K
    
    return {
        "E_total": e_total,
        "E_nuc": e_nuc,
        "E_core": e_core,
        "E_J": e_J,
        "E_K": e_K,
        "E_elec": e_elec
    }

# HF
hf = "H 0 0 0; F 0 0 0.9"

# NH3
nh3 = "N 0 0 0; H 0 0.9 0; H 0.9 0 0; H 0 0 -0.9"

# HF + NH3 複合体
complex_ab = """
N 0.0000 0.0000 0.0000
H 0.9000 0.0000 0.0000
H 0.0000 0.9000 0.0000
H 0.0000 0.0000 -0.9000
F 0.0000 0.0000 2.5000
H 0.0000 0.0000 3.4000
"""

# 各分子について計算
energy_a = get_hf_energy_decomposition(hf)
energy_b = get_hf_energy_decomposition(nh3)
energy_ab = get_hf_energy_decomposition(complex_ab)

# 差分計算
def print_energy_diff(energy_a, energy_b, energy_ab):
    print("===== 相互作用エネルギーの分解 =====")
    for key in ["E_total", "E_nuc", "E_core", "E_J", "E_K", "E_elec"]:
        delta = energy_ab[key] - (energy_a[key] + energy_b[key])
        kcal = delta * 627.509
        print(f"{key:10s} : ΔE = {delta:+.6f} Ha = {kcal:+.2f} kcal/mol")

print_energy_diff(energy_a, energy_b, energy_ab)
