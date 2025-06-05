"""
pySCFの計算を実行するための入力フォーマットに変更する部分
入力は、smilesとXYZに対応
出力は、2次元画像、3次元、zmatrixやpyscfのinput形式などに変換する。
"""

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdDetermineBonds
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
                AllChem.EmbedMolecule(self.mol)
                AllChem.UFFOptimizeMolecule(self.mol)  # Optimize the molecule
        elif input_type.lower() == "xyz":
            self.mol = self._load_from_xyz(input_data)
            # No embedding or optimization for XYZ input
        elif input_type.lower() == "z-matrix":
            raise NotImplementedError("Z-Matrix parsing is not supported natively by RDKit.")
        else:
            raise ValueError("Unsupported input type. Use 'smiles', 'xyz', or 'z-matrix'.")

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
            for atom in atoms:
                mol.AddAtom(Chem.Atom(atom))

            # Set 3D coordinates
            conf = Chem.Conformer(len(atoms))
            for i, coord in enumerate(coords):
                conf.SetAtomPosition(i, coord)
            mol.AddConformer(conf)

            # Use DetermineBonds to infer bonds
            rdDetermineBonds.DetermineBonds(mol)

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

    def generate_2d_image_with_atom_index(self, out_path):
        from rdkit.Chem import rdDepictor
        from rdkit.Chem.Draw import rdMolDraw2D

        # molをコピーして編集
        mol = Chem.Mol(self.mol)
        for atom in mol.GetAtoms():
            # molAtomMapNumberをセット（これで「C:0」などと表示される）
            atom.SetProp('molAtomMapNumber', str(atom.GetIdx()))

        rdDepictor.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DCairo(300, 300)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        with open(out_path, "wb") as f:
            f.write(drawer.GetDrawingText())

    def generate_3d_molblock(self, conf_id=0):
        """
        指定された配座の3D構造をMolBlock形式で返す。

        Args:
            conf_id (int): 配座ID（デフォルトは0）。

        Returns:
            str: MolBlock形式の3D構造。
        """
        if not self.mol:
            raise ValueError("Molecule is not initialized.")
        
        try:
            # 指定された配座IDのMolBlockを生成
            mol_block = Chem.MolToMolBlock(self.mol, confId=conf_id)
            if not mol_block:
                raise ValueError("MolBlock generation failed.")
            return mol_block
        except Exception as e:
            raise ValueError(f"Failed to generate 3D MolBlock: {e}")

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

    def get_xyz_coordinates(self, conf_id=0):
        """
        指定した配座IDのXYZ座標を取得する。

        Args:
            conf_id (int): 配座ID（デフォルトは0）。

        Returns:
            List[Tuple[str, float, float, float]]: 各原子のシンボルと座標のリスト。
        """
        if not self.mol:
            raise ValueError("Molecule is not initialized.")
        
        try:
            conf = self.mol.GetConformer(conf_id)
        except ValueError:
            raise ValueError(f"Conformer with ID {conf_id} does not exist.")
        
        xyz_coordinates = []
        for atom in self.mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            xyz_coordinates.append((atom.GetSymbol(), pos.x, pos.y, pos.z))
        
        return xyz_coordinates

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

    def generate_conformers(self, num_conformers=10, max_iterations=200, prune_rms_threshold=0.5, forcefield="UFF"):
        """
        配座探索を実行し、複数の配座を生成する。
        
        Args:
            num_conformers (int): 生成する配座の数。
            max_iterations (int): 配座生成の最大反復回数。
            prune_rms_threshold (float): 配座間のRMSDのしきい値（類似配座を削除するため）。
            forcefield (str): エネルギー最適化に使用する力場 ("UFF" または "MMFF")。
        
        Returns:
            List[Chem.Conformer]: 生成された配座のリスト。
        """
        if not self.mol:
            raise ValueError("Molecule is not initialized.")
        
        # 配座生成の設定
        params = AllChem.ETKDGv3()
        params.numThreads = 0  # 自動でスレッド数を設定
        params.maxIterations = max_iterations
        params.pruneRmsThresh = prune_rms_threshold
        
        # 配座を生成
        conformer_ids = AllChem.EmbedMultipleConfs(self.mol, numConfs=num_conformers, params=params)
        
        # 配座のエネルギーを最適化し、エネルギーを計算
        for conf_id in conformer_ids:
            if forcefield.upper() == "UFF":
                success = AllChem.UFFOptimizeMolecule(self.mol, confId=conf_id)
                if success == 0:  # 最適化成功
                    ff = AllChem.UFFGetMoleculeForceField(self.mol, confId=conf_id)
                    energy = ff.CalcEnergy()
                    self.mol.SetProp(f"Energy_{conf_id}", str(energy))
            elif forcefield.upper() == "MMFF":
                mmff_props = AllChem.MMFFGetMoleculeProperties(self.mol, mmffVariant="MMFF94")
                if mmff_props is None:
                    raise ValueError("MMFF parameters could not be assigned to the molecule.")
                success = AllChem.MMFFOptimizeMolecule(self.mol, mmff_props, confId=conf_id)
                if success == 0:  # 最適化成功
                    ff = AllChem.MMFFGetMoleculeForceField(self.mol, mmff_props, confId=conf_id)
                    energy = ff.CalcEnergy()
                    self.mol.SetProp(f"Energy_{conf_id}", str(energy))
            else:
                raise ValueError("Unsupported forcefield. Use 'UFF' or 'MMFF'.")
        
        # 配座をリストとして返す
        return [self.mol.GetConformer(conf_id) for conf_id in conformer_ids]


