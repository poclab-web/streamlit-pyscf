
import os
import numpy as np
import streamlit as st
from pyscf import gto, scf, lib



def list_chk_files(data_dir="data"):
    chk_paths = []
    for root, _, files in os.walk(data_dir):
        for f in files:
            if f.endswith(".chk"):
                chk_paths.append(os.path.join(root, f))
    return chk_paths

def load_mo_info(chk_path):
    try:
        mol = lib.chkfile.load_mol(chk_path)
        mf = scf.RHF(mol)
        mf.__dict__.update(lib.chkfile.load(chk_path, 'scf'))
        mo_occ = np.array(mf.mo_occ)
        homo_idx = np.where(mo_occ > 0)[0][-1]
        lumo_idx = homo_idx + 1
        return mol, mf.mo_coeff, mo_occ, homo_idx, lumo_idx, mf
    except Exception as e:
        st.error(f"❌ 読み込み失敗: {e}")
        return None, None, None, None, None