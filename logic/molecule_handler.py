"""
pySCFの計算を実行するための入力フォーマットに変更する部分
"""


from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.rdmolfiles import MolFromXYZFile
from rdkit.Geometry import Point3D

class MoleculeHandler:
    def __init__(self, input_data, input_type="smiles"):
        self.mol = None
        if input_type.lower() == "smiles":
            self.mol = Chem.MolFromSmiles(input_data)
            if self.mol:
                AllChem.Compute2DCoords(self.mol)
                self.mol = Chem.AddHs(self.mol)  # Add hydrogens to the molecule
        elif input_type.lower() == "xyz":
            self.mol = self._load_from_xyz(input_data)
        elif input_type.lower() == "z-matrix":
            raise NotImplementedError("Z-Matrix parsing is not supported natively by RDKit.")
        else:
            raise ValueError("Unsupported input type. Use 'smiles', 'xyz', or 'z-matrix'.")

        if self.mol:
            AllChem.EmbedMolecule(self.mol)
            AllChem.UFFOptimizeMolecule(self.mol)

    def _load_from_xyz(self, xyz):
        try:
            # Parse the XYZ data manually
            lines = xyz.strip().split("\n")
            atoms = []
            coords = []
            
            for line in lines:
                parts = line.split()
                if len(parts) != 4:
                    raise ValueError(f"Invalid XYZ format: {line}")
                atom = parts[0]
                x, y, z = map(float, parts[1:])
                atoms.append(atom)
                coords.append(Point3D(x, y, z))

            # Create an editable molecule
            mol = Chem.RWMol()

            # Add atoms
            atom_indices = []
            for atom in atoms:
                a = Chem.Atom(atom)
                atom_indices.append(mol.AddAtom(a))

            # Set 3D coordinates
            conf = Chem.Conformer(len(atoms))
            for i, coord in enumerate(coords):
                conf.SetAtomPosition(i, coord)
            mol.AddConformer(conf)

            # Attempt to determine bonds
            Chem.rdDetermineBonds.DetermineBonds(mol)

            # Add hydrogens
            mol = Chem.AddHs(mol)

            return mol
        except Exception as e:
            print(f"Error loading XYZ data: {e}")
            return None

    def generate_2d_image(self, output_path="molecule_2d.png"):
        if not self.mol:
            raise ValueError("Molecule is not initialized.")
        smiles = Chem.MolToSmiles(self.mol)
        mol = Chem.MolFromSmiles(smiles)
        img = Draw.MolToImage(mol, size=(300, 300))
        img.save(output_path)

    def generate_3d_molblock(self):
        if not self.mol:
            raise ValueError("Molecule is not initialized.")
        return Chem.MolToXYZBlock(self.mol)

    def save_to_sdf(self, output_path="molecule.sdf"):
        if not self.mol:
            raise ValueError("Molecule is not initialized.")
        self.mol.SetProp("_Name", Chem.MolToSmiles(self.mol))
        writer = Chem.SDWriter(output_path)
        writer.write(self.mol)
        writer.close()

    def save_to_xyz(self, output_path="molecule.xyz"):
        if not self.mol:
            raise ValueError("Molecule is not initialized.")
        conf = self.mol.GetConformer()
        with open(output_path, "w") as f:
            f.write(f"{self.mol.GetNumAtoms()}\n\n")
            for atom in self.mol.GetAtoms():
                pos = conf.GetAtomPosition(atom.GetIdx())
                f.write(f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n")

    def save_to_zmatrix(self, output_path="molecule.zmatrix"):
        if not self.mol:
            raise ValueError("Molecule is not initialized.")
        conf = self.mol.GetConformer()
        atoms = [atom.GetSymbol() for atom in self.mol.GetAtoms()]
        with open(output_path, "w") as f:
            f.write(f"{len(atoms)}\nZ-Matrix Format\n")
            for i, atom in enumerate(atoms):
                pos = conf.GetAtomPosition(i)
                if i == 0:
                    f.write(f"{atom}\n")
                elif i == 1:
                    f.write(f"{atom} 1 {pos.x:.4f}\n")
                elif i == 2:
                    f.write(f"{atom} 1 {pos.x:.4f} 2 {pos.y:.4f}\n")
                else:
                    f.write(f"{atom} 1 {pos.x:.4f} 2 {pos.y:.4f} 3 {pos.z:.4f}\n")

    def to_pyscf_input(self):
        """Convert RDKit molecule to PySCF gto.M(atom=...) format."""
        if not self.mol:
            raise ValueError("Molecule is not initialized.")

        conf = self.mol.GetConformer()
        pyscf_atoms = []

        for atom in self.mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            pyscf_atoms.append(f"{atom.GetSymbol()} {pos.x:.8f} {pos.y:.8f} {pos.z:.8f}")

        return "\n".join(pyscf_atoms)