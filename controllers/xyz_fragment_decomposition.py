"""
XYZ座標から2分子を自動分解するためのコントローラー
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.distance import cdist
import networkx as nx
from collections import defaultdict


def separate_molecules_by_distance(mol, cutoff_distance=2.5):
    """
    距離ベースで分子を分離する
    
    Parameters:
    - mol: RDKit Mol object
    - cutoff_distance: 分子間距離の閾値 (Å)
    
    Returns:
    - fragments: 分離された分子のリスト
    """
    if mol.GetNumConformers() == 0:
        # 3D座標が無い場合は結合情報から分離
        return list(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True))
    
    # 原子座標を取得
    conf = mol.GetConformer()
    coords = []
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        coords.append([pos.x, pos.y, pos.z])
    coords = np.array(coords)
    
    # 距離行列を計算
    dist_matrix = cdist(coords, coords)
    
    # 距離閾値以下の原子間にエッジを作成
    graph = nx.Graph()
    for i in range(mol.GetNumAtoms()):
        graph.add_node(i)
    
    for i in range(mol.GetNumAtoms()):
        for j in range(i + 1, mol.GetNumAtoms()):
            if dist_matrix[i][j] <= cutoff_distance:
                graph.add_edge(i, j)
    
    # 連結成分を取得（これが各分子に対応）
    components = list(nx.connected_components(graph))
    
    # 各成分から分子を作成
    fragments = []
    for component in components:
        if len(component) == 0:
            continue
        
        # 新しい分子を作成
        fragment_mol = Chem.RWMol()
        atom_mapping = {}
        
        # 原子を追加
        for old_idx in component:
            atom = mol.GetAtomWithIdx(old_idx)
            new_idx = fragment_mol.AddAtom(Chem.Atom(atom.GetSymbol()))
            atom_mapping[old_idx] = new_idx
        
        # 結合を追加
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if begin_idx in component and end_idx in component:
                fragment_mol.AddBond(
                    atom_mapping[begin_idx],
                    atom_mapping[end_idx],
                    bond.GetBondType()
                )
        
        # 配座を追加
        if mol.GetNumConformers() > 0:
            fragment_conf = Chem.Conformer(len(component))
            for i, old_idx in enumerate(component):
                pos = conf.GetAtomPosition(old_idx)
                fragment_conf.SetAtomPosition(i, pos)
            fragment_mol.AddConformer(fragment_conf)
        
        # 分子を仕上げ
        try:
            fragment_mol = fragment_mol.GetMol()
            Chem.SanitizeMol(fragment_mol)
            fragments.append(fragment_mol)
        except:
            # サニタイズに失敗した場合はスキップ
            continue
    
    return fragments


def separate_molecules_by_clustering(mol, n_clusters=2, method="simple"):
    """
    シンプルな距離ベースクラスタリングによる分子分離
    
    Parameters:
    - mol: RDKit Mol object
    - n_clusters: クラスター数
    - method: クラスタリング手法 ("simple")
    
    Returns:
    - fragments: 分離された分子のリスト
    """
    if mol.GetNumConformers() == 0:
        # 3D座標が無い場合は結合情報から分離
        return list(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True))
    
    # 原子座標を取得
    conf = mol.GetConformer()
    coords = []
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        coords.append([pos.x, pos.y, pos.z])
    coords = np.array(coords)
    
    # シンプルなクラスタリング：重心からの距離でグループ化
    if method == "simple":
        # 全体の重心を計算
        centroid = np.mean(coords, axis=0)
        
        # 各原子から重心への距離を計算
        distances = np.linalg.norm(coords - centroid, axis=1)
        
        # 距離で2グループに分ける（単純な閾値ベース）
        threshold = np.median(distances)
        cluster_labels = (distances > threshold).astype(int)
        
        # 少なくとも各グループに1つの原子があることを確認
        if np.all(cluster_labels == 0) or np.all(cluster_labels == 1):
            # 距離でソートして半分に分ける
            sorted_indices = np.argsort(distances)
            cluster_labels = np.zeros(len(coords), dtype=int)
            cluster_labels[sorted_indices[len(sorted_indices)//2:]] = 1
    else:
        raise ValueError(f"Unknown clustering method: {method}")
    
    # 各クラスターから分子を作成
    fragments = []
    for cluster_id in range(n_clusters):
        component = [i for i, label in enumerate(cluster_labels) if label == cluster_id]
        
        if len(component) == 0:
            continue
        
        # 新しい分子を作成
        fragment_mol = Chem.RWMol()
        atom_mapping = {}
        
        # 原子を追加
        for old_idx in component:
            atom = mol.GetAtomWithIdx(old_idx)
            new_idx = fragment_mol.AddAtom(Chem.Atom(atom.GetSymbol()))
            atom_mapping[old_idx] = new_idx
        
        # 結合を追加
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if begin_idx in component and end_idx in component:
                fragment_mol.AddBond(
                    atom_mapping[begin_idx],
                    atom_mapping[end_idx],
                    bond.GetBondType()
                )
        
        # 配座を追加
        if mol.GetNumConformers() > 0:
            fragment_conf = Chem.Conformer(len(component))
            for i, old_idx in enumerate(component):
                pos = conf.GetAtomPosition(old_idx)
                fragment_conf.SetAtomPosition(i, pos)
            fragment_mol.AddConformer(fragment_conf)
        
        # 分子を仕上げ
        try:
            fragment_mol = fragment_mol.GetMol()
            Chem.SanitizeMol(fragment_mol)
            fragments.append(fragment_mol)
        except:
            # サニタイズに失敗した場合はスキップ
            continue
    
    return fragments


def analyze_fragment_separation(fragments):
    """
    分子分離の結果を分析する
    
    Parameters:
    - fragments: 分離された分子のリスト
    
    Returns:
    - analysis: 分析結果の辞書
    """
    analysis = {
        "num_fragments": len(fragments),
        "fragment_sizes": [mol.GetNumAtoms() for mol in fragments],
        "fragment_formulas": [],
        "fragment_charges": [],
        "separation_quality": "unknown"
    }
    
    # 各フラグメントの分子式を計算
    for mol in fragments:
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
        analysis["fragment_formulas"].append(formula)
        
        # 形式電荷を計算
        total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        analysis["fragment_charges"].append(total_charge)
    
    # 分離品質の評価
    if len(fragments) == 2:
        size_diff = abs(analysis["fragment_sizes"][0] - analysis["fragment_sizes"][1])
        total_atoms = sum(analysis["fragment_sizes"])
        
        if size_diff / total_atoms < 0.1:
            analysis["separation_quality"] = "excellent"
        elif size_diff / total_atoms < 0.3:
            analysis["separation_quality"] = "good"
        else:
            analysis["separation_quality"] = "poor"
    elif len(fragments) == 1:
        analysis["separation_quality"] = "failed"
    else:
        analysis["separation_quality"] = "too_many_fragments"
    
    return analysis


def get_fragment_interaction_distance(mol1, mol2):
    """
    2つの分子間の最短距離を計算
    
    Parameters:
    - mol1, mol2: RDKit Mol objects
    
    Returns:
    - min_distance: 最短距離 (Å)
    """
    if mol1.GetNumConformers() == 0 or mol2.GetNumConformers() == 0:
        return None
    
    # 座標を取得
    conf1 = mol1.GetConformer()
    conf2 = mol2.GetConformer()
    
    coords1 = []
    for i in range(mol1.GetNumAtoms()):
        pos = conf1.GetAtomPosition(i)
        coords1.append([pos.x, pos.y, pos.z])
    
    coords2 = []
    for i in range(mol2.GetNumAtoms()):
        pos = conf2.GetAtomPosition(i)
        coords2.append([pos.x, pos.y, pos.z])
    
    coords1 = np.array(coords1)
    coords2 = np.array(coords2)
    
    # 距離行列を計算
    dist_matrix = cdist(coords1, coords2)
    
    return np.min(dist_matrix)
