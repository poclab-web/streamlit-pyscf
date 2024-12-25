import streamlit as st
from pyscf import gto, scf

st.title("PySCF + Streamlit: Quantum Chemistry on the Cloud")

# ユーザー入力
st.sidebar.header("Molecular Input")
atom_input = st.sidebar.text_area("Atoms (XYZ format)", "H 0 0 0\nH 0 0 0.74")
basis_set = st.sidebar.selectbox("Basis Set", ["sto-3g", "cc-pVDZ", "cc-pVTZ"])

# 計算を実行
if st.sidebar.button("Run Calculation"):
    try:
        st.write("Setting up molecule...")
        mol = gto.M(atom=atom_input, basis=basis_set)
        
        st.write("Running Hartree-Fock calculation...")
        mf = scf.RHF(mol)
        energy = mf.kernel()
        
        st.success(f"Calculated HF Energy: {energy} Hartree")
    except Exception as e:
        st.error(f"Error: {e}")
