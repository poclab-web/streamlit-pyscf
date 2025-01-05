"""
分子力場計算を用いた構造探索
TODO: 改修中 molecule_handler.pyに組み入れる。
"""

import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

# Function to generate and optimize conformers
def generate_conformers(mol, num_conformers, force_field, prune_rms_thresh):
    conformers = []
    try:
        # Generate multiple conformers with pruning
        params = AllChem.ETKDGv3()
        params.numThreads = 0  # Use all available threads
        params.pruneRmsThresh = prune_rms_thresh
        conformer_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, params=params)

        # Optimize each conformer
        for conf_id in conformer_ids:
            if force_field == "UFF":
                energy = AllChem.UFFOptimizeMolecule(mol, confId=conf_id)
            elif force_field == "MMFF":
                mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94')
                if mmff_props is not None:
                    energy = AllChem.MMFFOptimizeMolecule(mol, mmff_props, confId=conf_id)
                else:
                    raise ValueError("MMFF properties could not be calculated for the molecule.")
            else:
                raise ValueError("Unsupported force field")
            conformers.append((conf_id, energy))
    except Exception as e:
        st.error(f"Error generating conformers: {e}")
    return conformers


# Streamlit interface
st.title("Molecule Conformer Generator")

# Input molecule in SMILES format
smiles = st.text_input("Enter SMILES string:", value="CCO")

# Number of conformers to generate
num_conformers = st.slider("Number of Conformers", min_value=1, max_value=100000, value=1000)

# Choose force field
force_field = st.selectbox("Select Force Field", ["UFF", "MMFF"])

# RMS threshold for pruning similar structures
prune_rms_thresh = st.slider("Prune RMS Threshold", min_value=0.1, max_value=2.0, value=0.5, step=0.1)

# Energy difference threshold for filtering conformers
energy_diff_thresh = st.slider("Energy Difference Threshold (kcal/mol)", min_value=-.5, max_value=50.0, value=5.0, step=1.0)

if st.button("Generate Conformers"):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        mol = Chem.AddHs(mol)
        conformers = generate_conformers(mol, num_conformers, force_field, prune_rms_thresh)

        # Sort conformers by energy
        conformers = sorted(conformers, key=lambda x: x[1])

        # Filter conformers by energy difference
        min_energy = conformers[0][1] if conformers else None
        filtered_conformers = [
            (conf_id, energy)
            for conf_id, energy in conformers
            if energy - min_energy <= energy_diff_thresh
        ]

        st.write(f"Generated {len(conformers)} conformers. {len(filtered_conformers)} conformers remain after filtering:")
        for conf_id, energy in filtered_conformers:
            st.write(f"Conformer ID: {conf_id}, Energy: {energy:.2f} kcal/mol")
    else:
        st.error("Invalid SMILES string.")
