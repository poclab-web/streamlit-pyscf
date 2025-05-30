"""
分子力場計算を用いた構造探索
"""

import streamlit as st
import py3Dmol  # 3D可視化用ライブラリ
import stmol
import streamlit.components.v1 as components  # StreamlitでHTML埋め込み用
from logic.molecule_handler import MoleculeHandler  # MoleculeHandlerクラスをインポート


# Function to display 3D structure using py3Dmol
def show_3d_structure(mol_block):
    try:        
        viewer = py3Dmol.view(width=400, height=400)
        viewer.addModel(mol_block, "mol")  
        viewer.setStyle({"stick": {}})
        viewer.zoomTo()
        stmol.showmol(viewer, height=400)
    except Exception as e:
        st.warning(f"Unable to generate 3D structure: {e}")

                
# Streamlit interface
st.title("Molecule Conformer Generator")

# Input molecule in SMILES format
smiles = st.text_input("Enter SMILES string:", value="CCO")

# Choose force field
force_field = st.selectbox("Select Force Field", ["UFF", "MMFF"])

# Number of conformers to generate
num_conformers = st.number_input("Number of Conformers", value=100)

# RMS threshold for pruning similar structures
prune_rms_thresh = st.slider("Prune RMS Threshold", min_value=0.1, max_value=2.0, value=0.5, step=0.1)

# Energy difference threshold for filtering conformers
energy_diff_thresh = st.slider("Energy Difference Threshold (kcal/mol)", min_value=0.0, max_value=50.0, value=5.0, step=1.0)

if st.button("Generate Conformers"):
    try:
        # Initialize MoleculeHandler
        handler = MoleculeHandler(smiles, input_type="smiles")

        # Generate conformers
        conformers = handler.generate_conformers(
            num_conformers=num_conformers,
            prune_rms_threshold=prune_rms_thresh,
            forcefield=force_field
        )

        # Sort conformers by energy
        conformers_with_energy = []
        for conf in conformers:
            energy = handler.mol.GetProp(f"Energy_{conf.GetId()}")
            conformers_with_energy.append((conf.GetId(), float(energy)))

        conformers_with_energy = sorted(conformers_with_energy, key=lambda x: x[1])

        # Filter conformers by energy difference
        min_energy = conformers_with_energy[0][1] if conformers_with_energy else None
        filtered_conformers = [
            (conf_id, energy)
            for conf_id, energy in conformers_with_energy
            if energy - min_energy <= energy_diff_thresh
        ]

        st.write(f"Generated {len(conformers_with_energy)} conformers. {len(filtered_conformers)} conformers remain after filtering:")
        for conf_id, energy in filtered_conformers:
            energy_diff = energy - min_energy
            st.write(f"Conformer ID: {conf_id}, Energy: {energy:.2f} kcal/mol, ΔE: {energy_diff:.2f} kcal/mol")
            mol_block = handler.generate_3d_molblock(conf_id=conf_id)
            show_3d_structure(mol_block)
            
            # XYZ座標を取得してフォーマットして表示
            xyz_coordinates = handler.get_xyz_coordinates(conf_id=conf_id)
            st.write("XYZ Coordinates:")
            xyz_text = "\n".join([f"{atom:<2} {x:>10.6f} {y:>10.6f} {z:>10.6f}" for atom, x, y, z in xyz_coordinates])
            st.text(xyz_text)

    except Exception as e:
        st.error(f"Error: {e}")
