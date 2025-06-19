"""
pySCFの計算を実行するための入力フォーマットに変更する部分
入力は、smilesとXYZに対応
出力は、2次元画像、3次元、zmatrixやpyscfのinput形式などに変換する。
"""
from io import BytesIO

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdDetermineBonds
from rdkit.Chem.rdmolfiles import MolFromXYZFile
from rdkit.Geometry import Point3D, Point2D
from rdkit.Chem import rdchem

from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

class MoleculeHandler:
    def __init__(self, input_data, input_type="smiles"):
        self.mol = None

        if input_type.lower() == "smiles":
            raw_smiles = input_data
            base_mol = Chem.MolFromSmiles(raw_smiles)

            if base_mol:
                AllChem.Compute2DCoords(base_mol)
                bracket_indices = self._get_bracket_atom_indices(raw_smiles, base_mol)

                # Hを追加
                mol_with_H = Chem.AddHs(base_mol)
                mol_with_H = Chem.RWMol(mol_with_H)

                # 幾何最適化
                AllChem.EmbedMolecule(mol_with_H)
                AllChem.UFFOptimizeMolecule(mol_with_H)
                self.mol = mol_with_H

        elif input_type.lower() == "xyz":
            self.mol = self._load_from_xyz(input_data)

        elif input_type.lower() == "z-matrix":
            raise NotImplementedError("Z-Matrix parsing is not supported.")

        else:
            raise ValueError("Unsupported input type. Use 'smiles', 'xyz', or 'z-matrix'.")

    def _get_bracket_atom_indices(self, smiles, mol):
        """
        SMILES文字列中の[]に囲まれた原子のインデックスを返す
        """
        bracket_atoms = []
        import re
        bracket_elements = re.findall(r'\[.*?\]', smiles)
        for be in bracket_elements:
            temp = Chem.MolFromSmiles(be)
            if temp:
                match_idx = self._find_substructure_index(mol, temp)
                if match_idx is not None:
                    bracket_atoms.extend(match_idx)
        return list(set(bracket_atoms))

    def _find_substructure_index(self, mol, query):
        match = mol.GetSubstructMatch(query)
        return list(match) if match else None

    def _contains_radical_atoms(self, smiles):
        return "[" in smiles and "]" in smiles

    def _load_from_xyz(self, xyz):
        try:
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

            mol = Chem.RWMol()
            for atom in atoms:
                mol.AddAtom(Chem.Atom(atom))

            conf = Chem.Conformer(len(atoms))
            for i, coord in enumerate(coords):
                conf.SetAtomPosition(i, coord)
            mol.AddConformer(conf)

            rdDetermineBonds.DetermineBonds(mol)
            mol = Chem.AddHs(mol)
            return mol
        except Exception as e:
            print(f"Error loading XYZ data: {e}")
            return None

    def is_radical(self):
        if not self.mol:
            raise ValueError("Molecule is not initialized.")
        return any(atom.GetNumRadicalElectrons() > 0 for atom in self.mol.GetAtoms())




    def generate_conformers(self, num_conformers=10, max_iterations=200, prune_rms_threshold=0.5, forcefield="MMFF"):
        if not self.mol:
            raise ValueError("Molecule is not initialized.")
        
        # [H] のような単原子分子はスキップ
        if self.mol.GetNumAtoms() <= 1:
            print("[INFO] Single-atom molecule detected. Skipping conformer generation.")
            return []

        params = AllChem.ETKDGv3()
        params.numThreads = 0
        params.maxIterations = max_iterations
        params.pruneRmsThresh = prune_rms_threshold

        conformer_ids = AllChem.EmbedMultipleConfs(self.mol, numConfs=num_conformers, params=params)

        for conf_id in conformer_ids:
            if forcefield.upper() == "UFF":
                success = AllChem.UFFOptimizeMolecule(self.mol, confId=conf_id)
                if success == 0:
                    ff = AllChem.UFFGetMoleculeForceField(self.mol, confId=conf_id)
                    energy = ff.CalcEnergy()
                    self.mol.SetProp(f"Energy_{conf_id}", str(energy))
                else:
                    print(f"[WARNING] UFF optimization failed for conformer {conf_id}.")

            elif forcefield.upper() == "MMFF":
                mmff_props = AllChem.MMFFGetMoleculeProperties(self.mol, mmffVariant="MMFF94")
                if mmff_props is None:
                    print(f"[WARNING] MMFF parameters could not be assigned to the molecule. Conformer {conf_id} kept without optimization.")
                    continue  

                success = AllChem.MMFFOptimizeMolecule(self.mol, confId=conf_id)
                if success == 0:
                    ff = AllChem.MMFFGetMoleculeForceField(self.mol, mmff_props, confId=conf_id)
                    energy = ff.CalcEnergy()
                    self.mol.SetProp(f"Energy_{conf_id}", str(energy))
                else:
                    print(f"[WARNING] MMFF optimization failed for conformer {conf_id}.")
                    self.mol.SetProp(f"Energy_{conf_id}", "NaN")

            else:
                raise ValueError("Unsupported forcefield. Use 'UFF' or 'MMFF'.")

        return [self.mol.GetConformer(conf_id) for conf_id in conformer_ids]


    def keep_lowest_energy_conformer(self):
        """
        self.mol 内の conformer のうち、最もエネルギーが低いもの 1 つだけ残す。
        エネルギーは self.mol.GetProp(f"Energy_{conf_id}") に保存されていることを前提とする。
        """
        if not self.mol:
            raise ValueError("Molecule is not initialized.")

        energies = []
        for conf in self.mol.GetConformers():
            conf_id = conf.GetId()
            if self.mol.HasProp(f"Energy_{conf_id}"):
                energy = float(self.mol.GetProp(f"Energy_{conf_id}"))
                energies.append((conf_id, energy))

        if not energies:
            raise RuntimeError("No conformer energies found.")

        # 最もエネルギーが低い conformer を取得
        min_conf_id, _ = min(energies, key=lambda x: x[1])
        best_conf = rdchem.Conformer(self.mol.GetConformer(min_conf_id))  # コピー
        self.mol.RemoveAllConformers()
        self.mol.AddConformer(best_conf, assignId=True)  # confId = 0 に再設定される


    def get_xyz_coordinates(self, conf_id=0):
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
    
    def generate_2d_image(self, output_path="molecule_2d.png"):
        if not self.mol:
            raise ValueError("Molecule is not initialized.")
        smiles = Chem.MolToSmiles(self.mol)
        mol = Chem.MolFromSmiles(smiles)
        img = Draw.MolToImage(mol, size=(300, 300))
        img.save(output_path)


    def generate_2d_image_with_atom_index(self, out_path):
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


    def generate_2d_image_with_bond_index(self):
        if self.mol is None:
            raise ValueError("Molecule not initialized.")

        mol = Chem.Mol(self.mol)
        rdDepictor.Compute2DCoords(mol)

        drawer = rdMolDraw2D.MolDraw2DCairo(500, 500)
        drawer.drawOptions().addAtomIndices = False
        drawer.drawOptions().fontSize = 0.1
        drawer.DrawMolecule(mol)

        conf = mol.GetConformer()
        for bond in mol.GetBonds():
            idx = bond.GetIdx()
            begin = conf.GetAtomPosition(bond.GetBeginAtomIdx())
            end = conf.GetAtomPosition(bond.GetEndAtomIdx())

            x = (begin.x + end.x) / 2
            y = (begin.y + end.y) / 2

            offset = 0.3
            dx = end.x - begin.x
            dy = end.y - begin.y
            length = (dx ** 2 + dy ** 2) ** 0.5
            perp_x = -dy / length * offset
            perp_y = dx / length * offset
            x_shifted = x + perp_x
            y_shifted = y + perp_y

            drawer.DrawString(str(idx), Point2D(x_shifted, y_shifted))

        drawer.FinishDrawing()

        png = drawer.GetDrawingText()
        return BytesIO(png)


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
        if not self.mol:
            raise ValueError("Molecule is not initialized.")
        conf = self.mol.GetConformer()
        pyscf_atoms = []
        for atom in self.mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            pyscf_atoms.append(f"{atom.GetSymbol()} {pos.x:.8f} {pos.y:.8f} {pos.z:.8f}")
        return "\n".join(pyscf_atoms)

