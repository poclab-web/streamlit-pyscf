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
from rdkit.Chem import rdmolops

from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from scipy.spatial.distance import cdist
import numpy as np

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

        elif input_type.lower() == "rdkit":
            self.mol = input_data

        else:
            raise ValueError("Unsupported input type. Use 'smiles', 'xyz', or 'rdkit'.")

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
            # 改行と空白で分割し、1行に複数原子があっても対応
            lines = [l for l in xyz.strip().split("\n") if l.strip()]
            # 1行目が原子数ならスキップ
            try:
                natoms = int(lines[0])
                if len(lines) == natoms + 2:
                    lines = lines[2:]  # コメント行もスキップ
                elif len(lines) == natoms:
                    lines = lines[1:]  # コメントなし
                else:
                    lines = lines[2:]  # フォーマット不一致時もスキップ
            except ValueError:
                # 1行目が原子数でなければそのまま
                pass

            atoms = []
            coords = []
            for line in lines:
                parts = line.split()
                # 1行に複数原子が並ぶ場合も対応
                if len(parts) % 4 != 0:
                    raise ValueError(f"Invalid XYZ format: {line}")
                for i in range(0, len(parts), 4):
                    atom = parts[i]
                    x, y, z = map(float, parts[i+1:i+4])
                    atoms.append(atom)
                    coords.append((x, y, z))

            # RDKit Mol作成
            mol = Chem.RWMol()
            for atom in atoms:
                mol.AddAtom(Chem.Atom(atom))

            conf = Chem.Conformer(len(atoms))
            for i, (x, y, z) in enumerate(coords):
                conf.SetAtomPosition(i, Point3D(x, y, z))
            mol.AddConformer(conf)

            # 必要なら結合推定
            rdDetermineBonds.DetermineBonds(mol)
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

    def get_fragments(self):
        """
        分子をフラグメントごとに分割してMolオブジェクトのリストで返す
        """
        if not self.mol:
            raise ValueError("Molecule is not initialized.")
        return list(Chem.GetMolFrags(self.mol, asMols=True, sanitizeFrags=True))

    def get_fragments_smiles(self):
        """
        分子をフラグメントごとに分割し、それぞれのSMILESをリストで返す
        """
        fragments = self.get_fragments()
        return [Chem.MolToSmiles(frag) for frag in fragments]

    @staticmethod
    def combine_mols_with_3d(smiles_list, offset=5.0, forcefield="UFF"):
        """
        複数のSMILESから個別に3D構造を生成し、重ならないように配置して
        1つの複合体のMolオブジェクトとして返す。
        """
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mols = []
        for smi in smiles_list:
            mol = Chem.AddHs(Chem.MolFromSmiles(smi))
            AllChem.EmbedMolecule(mol)
            if forcefield.upper() == "UFF":
                AllChem.UFFOptimizeMolecule(mol)
            elif forcefield.upper() == "MMFF":
                mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
                if mmff_props is not None:
                    AllChem.MMFFOptimizeMolecule(mol, mmff_props)
            mols.append(mol)
        # 3D座標をずらして重ならないように配置
        for i, mol in enumerate(mols):
            conf = mol.GetConformer()
            for j in range(mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(j)
                conf.SetAtomPosition(j, (pos.x + i * offset, pos.y, pos.z))
        # 複合体Molを作成
        combo = mols[0]
        for m in mols[1:]:
            combo = Chem.CombineMols(combo, m)
        # Molオブジェクトとして返す
        return combo

    @staticmethod
    def place_mol_by_closest_distance(mol1, mol2, target_distance=2.5):
        """
        mol2をmol1の近くに、最近接原子間距離がtarget_distanceになるように平行移動して配置する。
        mol1, mol2: RDKit Molオブジェクト（各1配座のみ）
        target_distance: 配置したい最小距離（Å）
        戻り値: 配置後のmol2（mol1は変更しません）
        """
        conf1 = mol1.GetConformer()
        conf2 = mol2.GetConformer()
        coords1 = np.array([list(conf1.GetAtomPosition(i)) for i in range(mol1.GetNumAtoms())])
        coords2 = np.array([list(conf2.GetAtomPosition(i)) for i in range(mol2.GetNumAtoms())])

        # 最近接原子ペアを探す
        dists = cdist(coords1, coords2)
        idx1, idx2 = np.unravel_index(np.argmin(dists), dists.shape)
        min_dist = dists[idx1, idx2]

        # 移動ベクトルを計算
        vec = coords2[idx2] - coords1[idx1]
        norm = np.linalg.norm(vec)
        if norm == 0:
            # 万一同じ座標なら適当にずらす
            vec = np.array([1.0, 0.0, 0.0])
            norm = 1.0
        vec = vec / norm

        # 必要な移動量
        move = (target_distance - min_dist) * vec

        # mol2全体を平行移動
        for i in range(mol2.GetNumAtoms()):
            pos = conf2.GetAtomPosition(i)
            new_pos = np.array([pos.x, pos.y, pos.z]) + move
            conf2.SetAtomPosition(i, new_pos.tolist())

        return mol2

    @staticmethod
    def place_mol_for_ch_pi_interaction(mol1, mol2, target_distance=2.5, approach_angle=0.0, rotation_angle=0.0):
        """
        C-H π相互作用用の分子配置関数
        mol1: ベンゼン環などの π系分子
        mol2: メタンなどのC-H結合を持つ分子
        target_distance: ベンゼン環の重心から配置分子までの距離（Å）
        approach_angle: 接近角度（度）- 0°は垂直（π面に垂直）、90°は平行
        rotation_angle: 回転角度（度）- ベンゼン環の周りでの回転
        """
        from rdkit.Chem import Descriptors
        import math
        
        conf1 = mol1.GetConformer()
        conf2 = mol2.GetConformer()
        
        # ベンゼン環の重心を計算
        coords1 = np.array([list(conf1.GetAtomPosition(i)) for i in range(mol1.GetNumAtoms())])
        centroid1 = np.mean(coords1, axis=0)
        
        # ベンゼン環の法線ベクトルを計算（最初の3つの原子から平面を定義）
        if mol1.GetNumAtoms() >= 3:
            v1 = coords1[1] - coords1[0]
            v2 = coords1[2] - coords1[0]
            normal = np.cross(v1, v2)
            normal = normal / np.linalg.norm(normal)
        else:
            normal = np.array([0.0, 0.0, 1.0])  # デフォルト
        
        # 接近角度を適用
        approach_rad = math.radians(approach_angle)
        
        # 回転角度を適用
        rotation_rad = math.radians(rotation_angle)
        
        # 配置ベクトルを計算
        # 垂直成分
        vertical_component = normal * math.cos(approach_rad)
        
        # 水平成分（任意の方向）
        # ベンゼン環の平面内での方向を決定
        tangent = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(normal, tangent)) > 0.9:
            tangent = np.array([0.0, 1.0, 0.0])
        
        # 法線ベクトルに垂直な方向を計算
        tangent = tangent - np.dot(tangent, normal) * normal
        tangent = tangent / np.linalg.norm(tangent)
        
        # 回転を適用
        cos_rot = math.cos(rotation_rad)
        sin_rot = math.sin(rotation_rad)
        
        # 法線ベクトルの周りで回転
        bitangent = np.cross(normal, tangent)
        rotated_tangent = cos_rot * tangent + sin_rot * bitangent
        
        horizontal_component = rotated_tangent * math.sin(approach_rad)
        
        # 最終的な配置方向
        placement_direction = vertical_component + horizontal_component
        placement_direction = placement_direction / np.linalg.norm(placement_direction)
        
        # mol2の重心を計算
        coords2 = np.array([list(conf2.GetAtomPosition(i)) for i in range(mol2.GetNumAtoms())])
        centroid2 = np.mean(coords2, axis=0)
        
        # 目標位置を計算
        target_position = centroid1 + placement_direction * target_distance
        
        # 移動ベクトルを計算
        move_vector = target_position - centroid2
        
        # mol2全体を移動
        for i in range(mol2.GetNumAtoms()):
            pos = conf2.GetAtomPosition(i)
            new_pos = np.array([pos.x, pos.y, pos.z]) + move_vector
            conf2.SetAtomPosition(i, new_pos.tolist())
        
        return mol2

    @staticmethod
    def place_mol_by_specific_atoms(mol1, mol2, atom_idx1, atom_idx2, target_distance=2.5):
        """
        特定の原子間距離で分子を配置する
        mol1, mol2: RDKit Molオブジェクト
        atom_idx1: mol1の原子インデックス
        atom_idx2: mol2の原子インデックス
        target_distance: 目標距離（Å）
        """
        conf1 = mol1.GetConformer()
        conf2 = mol2.GetConformer()
        
        # 指定された原子の座標を取得
        pos1 = np.array(list(conf1.GetAtomPosition(atom_idx1)))
        pos2 = np.array(list(conf2.GetAtomPosition(atom_idx2)))
        
        # 現在の距離
        current_distance = np.linalg.norm(pos2 - pos1)
        
        # 移動ベクトルを計算
        if current_distance > 0:
            direction = (pos2 - pos1) / current_distance
            required_move = (target_distance - current_distance) * direction
        else:
            # 同じ位置にある場合は適当にずらす
            required_move = np.array([target_distance, 0.0, 0.0])
        
        # mol2全体を移動
        for i in range(mol2.GetNumAtoms()):
            pos = conf2.GetAtomPosition(i)
            new_pos = np.array([pos.x, pos.y, pos.z]) + required_move
            conf2.SetAtomPosition(i, new_pos.tolist())
        
        return mol2

