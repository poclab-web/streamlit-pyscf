# Examples

## Basic Usage Examples

### Single Point Energy Calculation

```python
from logic.molecule_handler import MoleculeHandler
from logic.calculation import run_quantum_calculation

# 分子の準備
handler = MoleculeHandler("CCO", input_type="smiles")
smiles = "CCO"
xyz = handler.to_pyscf_input()

# 計算実行
result = run_quantum_calculation(
    compound_name="ethanol",
    smiles=smiles,
    atom_input=xyz,
    basis_set="6-31G",
    theory="B3LYP",
    charge=0,
    spin=0
)

print(f"Energy: {result['energy']} Hartree")
```

### Geometry Optimization

```python
from logic.calculation import run_geometry_optimization

result = run_geometry_optimization(
    compound_name="ethanol",
    smiles="CCO",
    atom_input=xyz,
    basis_set="6-31G",
    theory="B3LYP",
    charge=0,
    spin=0,
    maxsteps=50
)

print(f"Optimized structure saved to: {result['xyz_file']}")
```

### Frequency Analysis

```python
from logic.calculation import calculate_vibrational_frequencies

freq_result = calculate_vibrational_frequencies(
    atom_input=optimized_xyz,
    theory="B3LYP",
    basis_set="6-31G",
    charge=0,
    spin=0,
    compound_name="ethanol"
)

print(f"Frequencies: {freq_result['frequencies']}")
print(f"ZPE: {freq_result['zpe']} Hartree")
```

## Advanced Examples

### Energy Decomposition Analysis

```python
from controllers.energydecompositionanalysis import run_eda_calculation

# 分子複合体の準備
complex_smiles = "CCO.O"  # エタノール-水複合体
handler = MoleculeHandler(complex_smiles, input_type="smiles")

# EDA計算
eda_result = run_eda_calculation(
    molecule_handler=handler,
    theory="B3LYP",
    basis_set="6-31G*",
    charge=0,
    spin=0
)

print(f"Interaction Energy: {eda_result['interaction_energy']} kcal/mol")
```

### Excited State Calculation

```python
from logic.calculation import run_tddft_calculation

td_result = run_tddft_calculation(
    compound_name="benzene",
    smiles="c1ccccc1",
    atom_input=xyz,
    basis_set="6-31G*",
    theory="B3LYP",
    charge=0,
    spin=0,
    nstates=10
)

print(f"Excitation energies: {td_result['excitation_energies']}")
```

## MOPAC Examples

### MOPAC Optimization

```python
from logic.mopac_calculation import run_mopac_calculation

handler = MoleculeHandler("CCO", input_type="smiles")

result = run_mopac_calculation(
    molecule_handler=handler,
    calculation_type="optimize",
    theory="PM7",
    charge=0,
    multiplicity=1
)

print(f"Optimization successful: {result['success']}")
```

## Database Examples

### Save Calculation Results

```python
from logic.database import insert_molecule_with_frequencies
from rdkit import Chem
from rdkit.Chem import inchi

mol = Chem.MolFromSmiles("CCO")
inchi_str = inchi.MolToInchi(mol)
inchikey_str = inchi.InchiToInchiKey(inchi_str)

insert_molecule_with_frequencies(
    inchi=inchi_str,
    inchikey=inchikey_str,
    g_tot=-154.1234,
    zpe=0.0456,
    frequencies=[1234.5, 1500.0, 2980.2],
    method="B3LYP",
    basis="6-31G",
    charge=0,
    spin=0
)
```

### Retrieve from Database

```python
from logic.database import get_molecule_from_sqlite

result = get_molecule_from_sqlite(
    inchikey="LFQSCWFLJHTTHZ-UHFFFAOYSA-N",  # エタノール
    method="B3LYP",
    basis="6-31G",
    spin=0,
    charge=0
)

if result:
    print(f"Found molecule with energy: {result['g_tot']} Hartree")
```
