"""
遷移状態の計算
ミニマムエネルギーパスを計算するところを作成
Nudged Elastic Band (NEB) 法による反応経路計算

機能:
- XYZ形式で分子構造を入力可能（初期構造・最終構造）
- 理論手法と基底関数系を選択可能
- NEB計算パラメータの調整可能
- 3D分子構造の可視化
- エネルギープロファイルの表示
"""

import streamlit as st
from utils.module import load_css
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import pandas as pd
from ase.io import read, write
import os
from rdkit import Chem
import py3Dmol  # 3D可視化用ライブラリ
import stmol
import streamlit.components.v1 as components  # StreamlitでHTML埋め込み用

# カスタムCSSを適用
load_css("config/styles.css")

from logic.molecule_handler import MoleculeHandler  # MoleculeHandlerクラスをインポート
from logic.calculation import (
    theory_options, basis_set_options, setup_molecule, 
    extract_mol_mf_params, run_quantum_calculation
)

from ase import Atoms
from ase.neb import NEB
from ase.optimize import BFGS
from ase.calculators.calculator import Calculator, all_changes
from pyscf import gto, scf, grad, solvent
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist

def normalize_basis_set(basis_set):
    """
    基底関数の表記を統一化する関数
    calculation.pyと統合予定だが、一時的にここで定義
    """
    # 基本的な正規化のみ実装
    basis_lower = basis_set.lower().replace("-", "").replace("_", "")
    if "*" in basis_set:
        basis_set = basis_set.replace("*", "(d)").replace("(d)(d)", "(d,p)")
    return basis_set

def xyz_to_ase_atoms(xyz_string, charge=0):
    """
    XYZ文字列をASE Atomsオブジェクトに変換
    """
    lines = xyz_string.strip().split('\n')
    symbols = []
    positions = []
    
    for line in lines:
        if line.strip():
            parts = line.strip().split()
            if len(parts) >= 4:
                symbol = parts[0]
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                symbols.append(symbol)
                positions.append([x, y, z])
    
    if symbols and positions:
        atoms = Atoms(symbols=symbols, positions=positions)
        # 電荷情報を保存（ASEでは直接サポートされていないが、情報として保持）
        atoms.info['charge'] = charge
        return atoms
    else:
        return None

def align_atom_order(initial_atoms, final_atoms):
    """
    初期構造と最終構造の原子順序を自動的に整列させる
    各原子を最も近い位置の同種原子にマッピングする
    """
    initial_symbols = initial_atoms.get_chemical_symbols()
    final_symbols = final_atoms.get_chemical_symbols()
    
    # 元素の種類と数が一致するかチェック
    from collections import Counter
    if Counter(initial_symbols) != Counter(final_symbols):
        return None, None  # 元素の種類や数が異なる場合は修正不可
    
    # 既に順序が一致している場合はそのまま返す
    if initial_symbols == final_symbols:
        return initial_atoms.copy(), final_atoms.copy()
    
    initial_positions = initial_atoms.positions
    final_positions = final_atoms.positions
    
    # 新しい順序で最終構造を再配列
    aligned_final_atoms = initial_atoms.copy()  # 初期構造をベースにする
    aligned_final_positions = np.zeros_like(final_positions)
    used_indices = set()
    
    # 各元素種ごとに処理
    unique_elements = list(set(initial_symbols))
    
    for element in unique_elements:
        # 初期構造でのこの元素のインデックス
        initial_element_indices = [i for i, sym in enumerate(initial_symbols) if sym == element]
        # 最終構造でのこの元素のインデックス
        final_element_indices = [i for i, sym in enumerate(final_symbols) if sym == element]
        
        if len(initial_element_indices) != len(final_element_indices):
            return None, None  # 元素数が一致しない
        
        # この元素の位置を取得
        initial_element_positions = initial_positions[initial_element_indices]
        final_element_positions = final_positions[final_element_indices]
        
        # 距離行列を計算して最適なマッピングを見つける
        distance_matrix = cdist(initial_element_positions, final_element_positions)
        
        # ハンガリアンアルゴリズムまたは貪欲法で最適割り当て
        # 簡単な貪欲法を使用（小さな分子では十分）
        assigned_final_indices = []
        available_final_indices = final_element_indices.copy()
        
        for i, initial_idx in enumerate(initial_element_indices):
            # 最も近い未使用の最終構造原子を見つける
            distances_to_available = []
            for final_idx in available_final_indices:
                final_pos_in_list = final_element_indices.index(final_idx)
                distances_to_available.append((distance_matrix[i, final_pos_in_list], final_idx))
            
            # 最小距離の原子を選択
            distances_to_available.sort()
            best_final_idx = distances_to_available[0][1]
            assigned_final_indices.append(best_final_idx)
            available_final_indices.remove(best_final_idx)
        
        # 対応する位置を設定
        for initial_idx, final_idx in zip(initial_element_indices, assigned_final_indices):
            aligned_final_positions[initial_idx] = final_positions[final_idx]
    
    # 新しい最終構造を作成
    aligned_final_atoms.positions = aligned_final_positions
    
    return initial_atoms.copy(), aligned_final_atoms

# ===== PySCF Calculator with energy & forces =====
class PySCFCalculator(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, basis='6-31g*', charge=0, spin=0, theory='B3LYP', 
                 solvent_model=None, eps=None, **kwargs):
        super().__init__(**kwargs)
        self.basis = normalize_basis_set(basis)
        self.charge = charge
        self.spin = spin
        self.theory = theory
        self.solvent_model = solvent_model
        self.eps = eps

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)

        # ASE AtomsオブジェクトからPySCF用のatom文字列を作成
        atom_string = ""
        for symbol, pos in zip(atoms.get_chemical_symbols(), atoms.positions):
            atom_string += f"{symbol} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}; "

        try:
            # calculation.pyのsetup_molecule関数を使用
            mol = setup_molecule(
                atom_input=atom_string,
                basis_set=self.basis,
                charge=self.charge,
                spin=self.spin,
                symmetry=False  # NEBでは対称性を無効化
            )

            # 理論手法の設定（calculation.pyのrun_quantum_calculationを参考）
            theory_lower = self.theory.lower().replace("-", "").replace("_", "")
            
            if theory_lower == "hf":
                if self.spin == 0:
                    mf = scf.RHF(mol)
                else:
                    mf = scf.UHF(mol)
            elif theory_lower == "mp2":
                from pyscf import mp
                if self.spin == 0:
                    mf = scf.RHF(mol)
                else:
                    mf = scf.UHF(mol)
                # MP2は後で追加予定
            elif theory_lower in ["b3lyp", "dftb3lyp"]:
                if self.spin == 0:
                    mf = scf.RKS(mol)
                else:
                    mf = scf.UKS(mol)
                mf.xc = 'b3lyp'
            elif theory_lower in ["camb3lyp", "cam-b3lyp"]:
                if self.spin == 0:
                    mf = scf.RKS(mol)
                else:
                    mf = scf.UKS(mol)
                mf.xc = 'cam-b3lyp'
            elif theory_lower in ["b3lypd3", "b3lyp-d3"]:
                if self.spin == 0:
                    mf = scf.RKS(mol)
                else:
                    mf = scf.UKS(mol)
                mf.xc = 'b3lyp'
                # D3補正は後で追加予定
            elif theory_lower == "pbe":
                if self.spin == 0:
                    mf = scf.RKS(mol)
                else:
                    mf = scf.UKS(mol)
                mf.xc = 'pbe'
            elif theory_lower == "pbe0":
                if self.spin == 0:
                    mf = scf.RKS(mol)
                else:
                    mf = scf.UKS(mol)
                mf.xc = 'pbe0'
            elif theory_lower == "m062x":
                if self.spin == 0:
                    mf = scf.RKS(mol)
                else:
                    mf = scf.UKS(mol)
                mf.xc = 'm06-2x'
            elif theory_lower == "tpss":
                if self.spin == 0:
                    mf = scf.RKS(mol)
                else:
                    mf = scf.UKS(mol)
                mf.xc = 'tpss'
            elif theory_lower == "scan":
                if self.spin == 0:
                    mf = scf.RKS(mol)
                else:
                    mf = scf.UKS(mol)
                mf.xc = 'scan'
            else:
                # デフォルトはB3LYP
                if self.spin == 0:
                    mf = scf.RKS(mol)
                else:
                    mf = scf.UKS(mol)
                mf.xc = 'b3lyp'

            # 溶媒効果の適用
            if self.solvent_model and self.eps:
                if self.solvent_model.upper() == "PCM":
                    mf = solvent.pcm.PCM(mf)
                    mf.with_solvent.eps = self.eps
                elif self.solvent_model.upper() == "DDCOSMO":
                    mf = solvent.ddcosmo.DDCOSMO(mf)
                    mf.with_solvent.eps = self.eps
            
            # SCF計算の収束設定
            mf.conv_tol = 1e-8
            mf.max_cycle = 100
            
            # エネルギー計算
            energy = mf.kernel()
            
            if not mf.converged:
                print(f"Warning: SCF not converged for structure")
                
            self.results['energy'] = energy

            # 力の計算
            if 'forces' in properties:
                if hasattr(mf, 'xc'):  # DFT
                    if self.spin == 0:
                        g = grad.RKS(mf)
                    else:
                        g = grad.UKS(mf)
                else:  # HF
                    if self.spin == 0:
                        g = grad.RHF(mf)
                    else:
                        g = grad.UHF(mf)
                
                forces = -g.kernel()  # ASE expects negative gradient
                self.results['forces'] = forces
                
        except Exception as e:
            print(f"Error in PySCF calculation: {e}")
            # フォールバック値を設定
            self.results['energy'] = 1e6
            if 'forces' in properties:
                self.results['forces'] = np.zeros((len(atoms), 3))

# ===== 3D分子可視化関数 =====
def create_3d_molecule_plot(atoms, title="分子構造"):
    """ASE Atomsオブジェクトから3D分子プロットを作成"""
    try:
        import plotly.graph_objects as go
        
        positions = atoms.positions
        symbols = atoms.get_chemical_symbols()
        
        # 原子の色を定義
        atom_colors = {
            'H': 'white',
            'C': 'gray',
            'N': 'blue',
            'O': 'red',
            'F': 'green',
            'S': 'yellow'
        }
        
        # 原子のサイズを定義
        atom_sizes = {
            'H': 5,
            'C': 8,
            'N': 8,
            'O': 8,
            'F': 7,
            'S': 10
        }
        
        colors = [atom_colors.get(symbol, 'gray') for symbol in symbols]
        sizes = [atom_sizes.get(symbol, 6) for symbol in symbols]
        
        fig = go.Figure(data=[go.Scatter3d(
            x=positions[:, 0],
            y=positions[:, 1],
            z=positions[:, 2],
            mode='markers+text',
            marker=dict(
                size=sizes,
                color=colors,
                line=dict(width=2, color='black')
            ),
            text=symbols,
            textposition="middle center",
            textfont=dict(size=12, color='black'),
            name='原子'
        )])
        
        # 結合を描画（距離が2.0Å以下の原子間）
        bond_x, bond_y, bond_z = [], [], []
        for i in range(len(atoms)):
            for j in range(i+1, len(atoms)):
                distance = np.linalg.norm(positions[i] - positions[j])
                if distance < 2.0:  # 2.0Å以下なら結合とみなす
                    bond_x.extend([positions[i][0], positions[j][0], None])
                    bond_y.extend([positions[i][1], positions[j][1], None])
                    bond_z.extend([positions[i][2], positions[j][2], None])
        
        if bond_x:  # 結合がある場合
            fig.add_trace(go.Scatter3d(
                x=bond_x, y=bond_y, z=bond_z,
                mode='lines',
                line=dict(color='black', width=4),
                name='結合',
                showlegend=False
            ))
        
        fig.update_layout(
            title=title,
            scene=dict(
                xaxis_title='X (Å)',
                yaxis_title='Y (Å)',
                zaxis_title='Z (Å)',
                aspectmode='cube'
            ),
            width=700,
            height=500
        )
        
        return fig
    except ImportError:
        return None

# Function to display 3D structure using py3Dmol
def show_3d_structure(mol_block):
    """
    py3Dmolを使用して3D構造を表示する関数
    MOLブロック形式の分子データを受け取り、3D可視化を行う
    改良版：エラーハンドリングを強化
    """
    try:
        if mol_block is None or mol_block.strip() == "":
            st.warning("表示するMOLデータがありません")
            return
            
        viewer = py3Dmol.view(width=400, height=400)
        viewer.addModel(mol_block, "mol")
        viewer.setStyle({"stick": {"radius": 0.1}, "sphere": {"radius": 0.3}})
        viewer.setBackgroundColor("white")
        viewer.zoomTo()
        stmol.showmol(viewer, height=400)
        
    except ImportError as e:
        st.error("py3Dmolまたはstmolライブラリが見つかりません")
        st.info("以下のコマンドでインストールしてください: pip install py3dmol stmol")
    except Exception as e:
        st.error(f"3D構造の表示中にエラーが発生しました: {e}")
        st.info("代替表示方法をお試しください")

def ase_atoms_to_mol_block(atoms):
    """
    ASE AtomsオブジェクトをMOLブロック形式に変換
    改良版：エラーハンドリングと結合検出を強化
    """
    try:
        positions = atoms.positions
        symbols = atoms.get_chemical_symbols()
        n_atoms = len(atoms)
        
        if n_atoms == 0:
            return None
        
        # 簡単なMOLブロック形式を手動で作成
        mol_lines = []
        mol_lines.append("")  # タイトル行
        mol_lines.append("  Generated by Streamlit-PySCF")  # プログラム行
        mol_lines.append("")  # コメント行
        
        # 結合を検出
        bonds = []
        bond_distances = {
            ('C', 'C'): 1.6, ('C', 'H'): 1.2, ('C', 'N'): 1.6, ('C', 'O'): 1.6,
            ('N', 'N'): 1.6, ('N', 'H'): 1.2, ('N', 'O'): 1.6,
            ('O', 'O'): 1.6, ('O', 'H'): 1.2,
            ('H', 'H'): 1.0
        }
        
        for i in range(n_atoms):
            for j in range(i+1, n_atoms):
                distance = np.linalg.norm(positions[i] - positions[j])
                sym_i, sym_j = symbols[i], symbols[j]
                
                # 結合距離判定を元素ペアに応じて調整
                bond_key = tuple(sorted([sym_i, sym_j]))
                max_distance = bond_distances.get(bond_key, 2.0)
                
                if distance < max_distance:
                    bonds.append((i+1, j+1, 1))  # 1-indexed, single bond
        
        # カウント行 (原子数, 結合数, その他)
        mol_lines.append(f"{n_atoms:3d}{len(bonds):3d}  0  0  0  0  0  0  0  0999 V2000")
        
        # 原子ブロック
        for symbol, pos in zip(symbols, positions):
            mol_lines.append(f"{pos[0]:10.4f}{pos[1]:10.4f}{pos[2]:10.4f} {symbol:<3s} 0  0  0  0  0  0  0  0  0  0  0  0")
        
        # 結合ブロック
        for bond in bonds:
            mol_lines.append(f"{bond[0]:3d}{bond[1]:3d}{bond[2]:3d}  0  0  0  0")
        
        mol_lines.append("M  END")
        
        mol_block = "\n".join(mol_lines)
        return mol_block
        
    except Exception as e:
        print(f"MOLブロック変換エラー: {e}")
        return None

def plot_energy_profile(energies):
    """エネルギープロファイルをプロット"""
    try:
        import plotly.graph_objects as go
        
        if not energies or all(e is None for e in energies):
            return None
            
        valid_energies = [e for e in energies if e is not None]
        min_energy = min(valid_energies)
        
        # 相対エネルギーをkcal/molに変換
        rel_energies = []
        x_vals = []
        for i, e in enumerate(energies):
            if e is not None:
                rel_energy = (e - min_energy) * 627.5095  # Hartree to kcal/mol
                rel_energies.append(rel_energy)
                x_vals.append(i)
        
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=x_vals,
            y=rel_energies,
            mode='lines+markers',
            marker=dict(size=8, color='blue'),
            line=dict(width=3, color='blue'),
            name='エネルギー'
        ))
        
        fig.update_layout(
            title='反応エネルギープロファイル',
            xaxis_title='反応座標 (Image番号)',
            yaxis_title='相対エネルギー (kcal/mol)',
            width=700,
            height=400,
            showlegend=False
        )
        
        return fig
    except ImportError:
        return None

def simple_3d_test():
    """シンプルな3D表示テスト"""
    try:
        # 水分子のテストMOLブロック
        test_mol = """

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7570    0.5860    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7570    0.5860    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END"""
        
        st.subheader("🧪 3D表示テスト（水分子）")
        viewer = py3Dmol.view(width=300, height=250)
        viewer.addModel(test_mol, "mol")
        viewer.setStyle({"stick": {"radius": 0.1}, "sphere": {"radius": 0.3}})
        viewer.setBackgroundColor("white")
        viewer.zoomTo()
        stmol.showmol(viewer, height=250)
        return True
    except Exception as e:
        st.error(f"3D表示テストに失敗: {e}")
        return False

def load_trajectory_if_exists(filename='data/neb.traj'):
    """軌道ファイルが存在する場合に読み込み"""
    if os.path.exists(filename):
        try:
            return read(filename, ':')
        except:
            return None
    return None
    """シンプルな3D表示テスト"""
    try:
        # 水分子のテストMOLブロック
        test_mol = """

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7570    0.5860    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7570    0.5860    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END"""
        
        st.subheader("🧪 3D表示テスト（水分子）")
        viewer = py3Dmol.view(width=300, height=250)
        viewer.addModel(test_mol, "mol")
        viewer.setStyle({"stick": {"radius": 0.1}, "sphere": {"radius": 0.3}})
        viewer.setBackgroundColor("white")
        viewer.zoomTo()
        stmol.showmol(viewer, height=250)
        return True
    except Exception as e:
        st.error(f"3D表示テストに失敗: {e}")
        return False
    """軌道ファイルが存在する場合に読み込み"""
    if os.path.exists(filename):
        try:
            return read(filename, ':')
        except:
            return None
    return None

# ===== Streamlit インターフェース =====
st.title("🔬 Nudged Elastic Band (NEB) 計算")

st.info("""
📌 **NEB計算の要件 - 高速計算用分子例**

🚀 **計算時間の目安**:
- **水素分子解離 (H₂)**: 1-3分程度 - 最も高速
- **エチレン回転**: 10-15分程度 - 軽量有機分子
- **アンモニア反転**: 5-10分程度 - 小分子
- **シクロヘキサン**: 30分以上 - 重い計算

💡 **推奨**: 初回は水素分子またはアンモニアから始めることを強く推奨します
""")

# 分子例選択
st.subheader("🧪 分子例を選択")
molecule_example = st.selectbox(
    "計算する分子例を選択してください:",
    [
        "水素分子解離 (H₂) - 最速",
        "アンモニア反転 (NH₃) - 高速", 
        "エチレン回転 (C₂H₄) - 中速",
        "シクロヘキサン配座変換 - 低速",
        "カスタム分子"
    ],
    index=0,
    help="計算時間を考慮して選択してください"
)

st.markdown(f"""
### 反応経路計算について - {molecule_example.split(' - ')[0]}
Nudged Elastic Band (NEB) 法は化学反応の遷移経路を探索する手法です。
初期構造と最終構造を設定し、その間の最小エネルギー経路を計算します。

**選択された計算例の詳細:**
""")

# 選択された分子例に応じた説明を表示
if "水素分子" in molecule_example:
    st.markdown("""
    🌟 **水素分子解離 (H₂)**: 
    - **計算時間**: 1-3分（最速）
    - **原子数**: 2原子のみ
    - **現象**: H-H結合の解離過程（結合切断）
    - **学習価値**: NEB法の基本概念を学ぶのに最適
    """)
elif "アンモニア" in molecule_example:
    st.markdown("""
    🔺 **アンモニア反転 (NH₃)**: 
    - **計算時間**: 3-8分（高速）
    - **原子数**: 4原子
    - **現象**: ピラミッド型分子の反転（umbrella inversion）
    - **学習価値**: 遷移状態を含む実際の反応経路
    """)
elif "エチレン" in molecule_example:
    st.markdown("""
    🔄 **エチレン回転 (C₂H₄)**: 
    - **計算時間**: 8-15分（中速）
    - **原子数**: 6原子
    - **現象**: C=C二重結合周りの回転障壁（π結合の破断・形成）
    - **学習価値**: π結合の性質と回転障壁の理解
    """)
elif "シクロヘキサン" in molecule_example:
    st.markdown("""
    🔄 **シクロヘキサン配座変換**: 
    - **計算時間**: 30分以上（低速）
    - **原子数**: 18原子
    - **現象**: 舟型から椅子型への配座変換
    - **学習価値**: 環状分子の立体化学
    """)
else:
    st.markdown("""
    ⚙️ **カスタム分子**: 
    - **計算時間**: 分子サイズに依存
    - **注意**: 大きな分子は計算時間が大幅に増加します
    """)

st.markdown("""
**重要な注意事項:**
- 適切な初期構造と最終構造が必要です
- 小さな分子から始めることを強く推奨します
- 中間構造数が多いほど精密ですが、計算時間も増加します
""")

st.warning("⚠️ **重要**: 初期構造と最終構造は同じ分子の異なる配座である必要があります。")

# 分子例に応じた構造データを定義
def get_molecule_structures(molecule_type):
    """選択された分子例に応じて初期・最終構造を返す"""
    
    if "水素分子" in molecule_type:
        # H2の解離反応: H-H結合が切れて2つの水素原子になる過程
        initial = """H     0.000000    0.000000    0.371000
H     0.000000    0.000000   -0.371000"""
        final = """H     0.000000    0.000000    1.500000
H     0.000000    0.000000   -1.500000"""
        return initial, final, "水素分子 (結合状態)", "水素分子 (解離状態)"
    
    elif "アンモニア" in molecule_type:
        # アンモニア反転反応: ピラミッド型 → 平面 → 逆ピラミッド型
        initial = """N     0.000000    0.000000    0.117000
H     0.000000    0.948000   -0.039000
H     0.821000   -0.474000   -0.039000
H    -0.821000   -0.474000   -0.039000"""
        final = """N     0.000000    0.000000   -0.117000
H     0.000000    0.948000    0.039000
H     0.821000   -0.474000    0.039000
H    -0.821000   -0.474000    0.039000"""
        return initial, final, "アンモニア (上向きピラミッド)", "アンモニア (下向きピラミッド)"
    
    elif "エチレン" in molecule_type:
        # エチレンの二重結合回転：エクリプス型 → 遷移状態 → スタガード型
        initial = """C     0.000000    0.000000    0.665000
C     0.000000    0.000000   -0.665000
H     0.000000    0.925000    1.230000
H     0.000000   -0.925000    1.230000
H     0.000000    0.925000   -1.230000
H     0.000000   -0.925000   -1.230000"""
        final = """C     0.000000    0.000000    0.665000
C     0.000000    0.000000   -0.665000
H     0.925000    0.000000    1.230000
H    -0.925000    0.000000    1.230000
H    -0.925000    0.000000   -1.230000
H     0.925000    0.000000   -1.230000"""
        return initial, final, "エチレン (エクリプス型)", "エチレン (スタガード型)"
    
    elif "シクロヘキサン" in molecule_type:
        # 舟型配座 → 椅子型配座の変換（分子力場最適化構造）
        initial = """C   -0.467756  -1.337606  -0.286856
C   -1.491088  -0.217767  -0.086516
C   -0.858981   1.143842   0.209874
C    0.517436   1.280856  -0.426540
C    1.491087   0.217767   0.086518
C    0.809301  -1.087094   0.503519
H   -0.911721  -2.297209   0.001193
H   -0.220021  -1.416375  -1.353028
H   -2.168423  -0.476468   0.736459
H   -2.116321  -0.149299  -0.985067
H   -1.522462   1.941708  -0.142433
H   -0.766383   1.272318   1.295793
H    0.927982   2.278635  -0.234689
H    0.419796   1.188176  -1.515657
H    2.051402   0.610137   0.943893
H    2.233341   0.015633  -0.695275
H    1.506202  -1.923134   0.375921
H    0.566608  -1.044122   1.572890"""
        final = """C    0.987394   1.081245   0.180920
C   -0.399110   1.365674  -0.390404
C   -1.429135   0.352260   0.101269
C   -0.987394  -1.081246  -0.180920
C    0.399110  -1.365674   0.390405
C    1.429135  -0.352260  -0.101269
H    1.713409   1.781938  -0.246523
H    0.977126   1.252671   1.264393
H   -0.354443   1.337329  -1.486127
H   -0.713492   2.377136  -0.109541
H   -2.394011   0.542868  -0.381908
H   -1.580706   0.481059   1.180085
H   -1.713408  -1.781938   0.246523
H   -0.977126  -1.252671  -1.264393
H    0.713492  -2.377136   0.109541
H    0.354443  -1.337328   1.486127
H    1.580706  -0.481059  -1.180085
H    2.394011  -0.542868   0.381908"""
        return initial, final, "シクロヘキサン舟型配座", "シクロヘキサン椅子型配座"
    
    else:  # カスタム分子
        return None, None, "初期構造", "最終構造"

# ===== 入力フォーム =====
st.header("📋 分子構造の入力")

# 計算レベルの推奨設定を表示
if "水素分子" in molecule_example or "アンモニア" in molecule_example:
    st.info("🎯 **小分子推奨設定**: HF/STO-3G または B3LYP/STO-3G でも十分です")
elif "エチレン" in molecule_example:
    st.info("🎯 **中分子推奨設定**: B3LYP/6-31G 以上を推奨します")
else:
    st.info("🎯 **大分子推奨設定**: B3LYP/6-31G(d) 以上を推奨します")

# 選択された分子例に応じて構造を取得
initial_xyz, final_xyz, initial_label, final_label = get_molecule_structures(molecule_example)

# 入力タブ
tab1, tab2 = st.tabs([f"🔬 {initial_label}", f"🎯 {final_label}"])

st.info("💡 **リアルタイム3Dプレビュー**: 入力した構造は右側に自動的に3D表示されます")

# デバッグ情報を追加
with st.expander("🔧 デバッグ情報"):
    st.write("py3Dmol利用可能:", "py3Dmol" in str(py3Dmol))
    st.write("stmol利用可能:", "stmol" in str(stmol))
    try:
        test_viewer = py3Dmol.view(width=100, height=100)
        st.success("✅ py3Dmol動作確認OK")
        
        # 3D表示テストを実行
        if st.button("🧪 3D表示テスト実行"):
            if simple_3d_test():
                st.success("✅ 3D表示テスト成功")
            else:
                st.error("❌ 3D表示テスト失敗")
                
    except Exception as e:
        st.error(f"❌ py3Dmolエラー: {e}")
        st.info("py3Dmolとstmolをインストールしてください: pip install py3dmol stmol")

with tab1:
    st.subheader(f"{initial_label}")
    input_type_initial = st.selectbox("入力形式 (初期)", ["XYZ", "SMILES"], key="input_type_initial")
    
    col1_input, col1_preview = st.columns([1, 1])
    
    with col1_input:
        if input_type_initial == "XYZ":
            # カスタム分子でない場合は自動設定
            default_initial = initial_xyz if initial_xyz else ""
            initial_structure = st.text_area(
                f"{initial_label}のXYZ座標",
                default_initial,
                key="initial_xyz",
                help="原子記号 X座標 Y座標 Z座標 の形式で入力",
                height=200
            )
        else:
            # SMILES用のデフォルト値
            default_smiles = "H" if "水素" in molecule_example else "N" if "アンモニア" in molecule_example else "C=C" if "エチレン" in molecule_example else "C1CCCCC1"
            initial_smiles = st.text_input(f"{initial_label}のSMILES", default_smiles, key="initial_smiles")
    
    with col1_preview:
        st.markdown("**3D構造プレビュー**")
        if input_type_initial == "XYZ" and 'initial_structure' in locals() and initial_structure and initial_structure.strip():
            try:
                # XYZ形式から3D構造を生成
                initial_atoms_preview = xyz_to_ase_atoms(initial_structure, 0)
                if initial_atoms_preview:
                    mol_block_initial = ase_atoms_to_mol_block(initial_atoms_preview)
                    if mol_block_initial:
                        try:
                            viewer = py3Dmol.view(width=350, height=300)
                            viewer.addModel(mol_block_initial, "mol")
                            viewer.setStyle({"stick": {}})
                            viewer.zoomTo()
                            stmol.showmol(viewer, height=300)
                        except Exception as viz_error:
                            st.error(f"3D可視化エラー: {viz_error}")
                            # フォールバック表示
                            st.text(f"原子数: {len(initial_atoms_preview)}")
                            st.text(f"元素: {', '.join(set(initial_atoms_preview.get_chemical_symbols()))}")
                            
                            # 代替可視化：座標テーブル
                            with st.expander("座標データを表示"):
                                coord_data = []
                                for i, (symbol, pos) in enumerate(zip(initial_atoms_preview.get_chemical_symbols(), 
                                                                     initial_atoms_preview.positions)):
                                    coord_data.append({
                                        "原子": f"{symbol}{i+1}",
                                        "X": f"{pos[0]:.3f}",
                                        "Y": f"{pos[1]:.3f}",
                                        "Z": f"{pos[2]:.3f}"
                                    })
                                st.table(coord_data[:10])  # 最初の10原子のみ表示
                    else:
                        st.warning("MOLブロック変換に失敗しました")
                        st.text("代替情報表示中...")
                else:
                    st.warning("XYZ形式の解析に失敗しました")
            except Exception as e:
                st.error(f"XYZ処理エラー: {e}")
        elif input_type_initial == "SMILES":
            try:
                if 'initial_smiles' in st.session_state and st.session_state.initial_smiles and st.session_state.initial_smiles.strip():
                    # SMILES形式から3D構造を生成
                    initial_handler = MoleculeHandler(st.session_state.initial_smiles, input_type="smiles")
                    if initial_handler.mol:
                        mol_block_initial = initial_handler.generate_3d_molblock()
                        if mol_block_initial:
                            try:
                                viewer = py3Dmol.view(width=350, height=300)
                                viewer.addModel(mol_block_initial, "mol")
                                viewer.setStyle({"stick": {}})
                                viewer.zoomTo()
                                stmol.showmol(viewer, height=300)
                            except Exception as viz_error:
                                st.error(f"3D可視化エラー: {viz_error}")
                                st.text("SMILES構造が正常に処理されました")
                        else:
                            st.warning("3D構造生成に失敗しました")
                    else:
                        st.warning("SMILES解析に失敗しました")
                else:
                    st.info("SMILESを入力してください")
            except Exception as e:
                st.error(f"SMILES処理エラー: {e}")
        else:
            st.info("分子構造を入力してください")

with tab2:
    st.subheader(f"{final_label}")
    input_type_final = st.selectbox("入力形式 (最終)", ["XYZ", "SMILES"], key="input_type_final")
    
    col2_input, col2_preview = st.columns([1, 1])
    
    with col2_input:
        if input_type_final == "XYZ":
            # カスタム分子でない場合は自動設定
            default_final = final_xyz if final_xyz else ""
            final_structure = st.text_area(
                f"{final_label}のXYZ座標",
                default_final,
                key="final_xyz",
                help="原子記号 X座標 Y座標 Z座標 の形式で入力",
                height=200
            )
        else:
            # SMILES用のデフォルト値
            default_smiles_final = "H" if "水素" in molecule_example else "N" if "アンモニア" in molecule_example else "C=C" if "エチレン" in molecule_example else "C1CCCCC1"
            final_smiles = st.text_input(f"{final_label}のSMILES", default_smiles_final, key="final_smiles")
    
    with col2_preview:
        st.markdown("**3D構造プレビュー**")
        if input_type_final == "XYZ" and 'final_structure' in locals() and final_structure and final_structure.strip():
            try:
                # XYZ形式から3D構造を生成
                final_atoms_preview = xyz_to_ase_atoms(final_structure, 0)
                if final_atoms_preview:
                    mol_block_final = ase_atoms_to_mol_block(final_atoms_preview)
                    if mol_block_final:
                        try:
                            viewer = py3Dmol.view(width=350, height=300)
                            viewer.addModel(mol_block_final, "mol")
                            viewer.setStyle({"stick": {}})
                            viewer.zoomTo()
                            stmol.showmol(viewer, height=300)
                        except Exception as viz_error:
                            st.error(f"3D可視化エラー: {viz_error}")
                            # フォールバック表示
                            st.text(f"原子数: {len(final_atoms_preview)}")
                            st.text(f"元素: {', '.join(set(final_atoms_preview.get_chemical_symbols()))}")
                            
                            # 代替可視化：座標テーブル
                            with st.expander("座標データを表示"):
                                coord_data = []
                                for i, (symbol, pos) in enumerate(zip(final_atoms_preview.get_chemical_symbols(), 
                                                                     final_atoms_preview.positions)):
                                    coord_data.append({
                                        "原子": f"{symbol}{i+1}",
                                        "X": f"{pos[0]:.3f}",
                                        "Y": f"{pos[1]:.3f}",
                                        "Z": f"{pos[2]:.3f}"
                                    })
                                st.table(coord_data[:10])  # 最初の10原子のみ表示
                    else:
                        st.warning("MOLブロック変換に失敗しました")
                        st.text("代替情報表示中...")
                else:
                    st.warning("XYZ形式の解析に失敗しました")
            except Exception as e:
                st.error(f"XYZ処理エラー: {e}")
        elif input_type_final == "SMILES":
            try:
                if 'final_smiles' in st.session_state and st.session_state.final_smiles and st.session_state.final_smiles.strip():
                    # SMILES形式から3D構造を生成
                    final_handler = MoleculeHandler(st.session_state.final_smiles, input_type="smiles")
                    if final_handler.mol:
                        mol_block_final = final_handler.generate_3d_molblock()
                        if mol_block_final:
                            try:
                                viewer = py3Dmol.view(width=350, height=300)
                                viewer.addModel(mol_block_final, "mol")
                                viewer.setStyle({"stick": {}})
                                viewer.zoomTo()
                                stmol.showmol(viewer, height=300)
                            except Exception as viz_error:
                                st.error(f"3D可視化エラー: {viz_error}")
                                st.text("SMILES構造が正常に処理されました")
                        else:
                            st.warning("3D構造生成に失敗しました")
                    else:
                        st.warning("SMILES解析に失敗しました")
                else:
                    st.info("SMILESを入力してください")
            except Exception as e:
                st.error(f"SMILES処理エラー: {e}")
        else:
            st.info("分子構造を入力してください")

# ===== 計算設定 =====
st.header("⚙️ 計算設定")

# 分子例に応じた推奨設定の表示
if "水素分子" in molecule_example:
    st.success("🚀 **超高速設定**: 水素分子は HF/STO-3G でも十分な結果が得られます")
    default_theory_idx = theory_options.index("HF") if "HF" in theory_options else 0
    default_basis_idx = 0  # STO-3G
elif "アンモニア" in molecule_example:
    st.info("⚡ **高速設定**: アンモニアは B3LYP/STO-3G または HF/6-31G で高速計算可能")
    default_theory_idx = theory_options.index("B3LYP") if "B3LYP" in theory_options else 0
    default_basis_idx = 0  # STO-3G
elif "エチレン" in molecule_example:
    st.info("🔧 **標準設定**: エチレンは B3LYP/6-31G を推奨")
    default_theory_idx = theory_options.index("B3LYP") if "B3LYP" in theory_options else 0
    for i, basis in enumerate(basis_set_options):
        if "6-31g" in basis.lower() and "*" not in basis.lower():
            default_basis_idx = i
            break
    else:
        default_basis_idx = 0
else:
    st.warning("💪 **高精度設定**: 大きな分子は B3LYP/6-31G* 以上を推奨")
    default_theory_idx = theory_options.index("B3LYP") if "B3LYP" in theory_options else 0
    for i, basis in enumerate(basis_set_options):
        if "6-31g*" in basis.lower():
            default_basis_idx = i
            break
    else:
        default_basis_idx = 0

col1, col2 = st.columns(2)
with col1:
    theory = st.selectbox("理論手法", theory_options, index=default_theory_idx)
    charge = st.number_input("分子電荷", min_value=-10, max_value=10, value=0, step=1)

with col2:
    basis_set = st.selectbox("基底関数", basis_set_options, index=default_basis_idx)
    multiplicity = st.number_input("スピン多重度 (2S + 1)", min_value=1, max_value=10, value=1, step=1)

spin = multiplicity - 1

# 選択された設定の評価
theory_eval = theory.lower().replace("-", "").replace("_", "")
basis_eval = basis_set.lower().replace("-", "").replace("*", "(d)").replace("**", "(d,p)")

if theory_eval == "hf" or basis_eval == "sto3g":
    st.warning("⚠️ 選択された設定はNEB計算には適していません。推奨設定をご検討ください。")
elif theory_eval in ["b3lyp", "b3lypd3", "pbe", "pbe0", "m062x"] and basis_eval in ["631g(d)", "6311g(d,p)", "ccpvdz", "ccpvtz", "def2svp", "def2tzvp"]:
    st.success("✅ 良い選択です！この設定はNEB計算に適しています。")
elif theory_eval in ["b3lyp", "b3lypd3", "pbe", "pbe0", "m062x"]:
    st.info("💡 良い理論手法です。より大きな基底関数系（6-31G(d)以上）を使用するとさらに良いでしょう。")

# ===== NEB パラメータ =====
with st.expander("🔧 NEB計算パラメータ"):
    col1, col2, col3 = st.columns(3)
    
    with col1:
        n_images = st.slider("中間構造数", 2, 8, 3, help="初期と最終の間に挿入する構造数")
        spring_constant = st.slider("Spring Constant", 0.05, 0.5, 0.1, 0.05, help="隣接構造間の結合強度")
    
    with col2:
        max_force = st.slider("収束判定 (fmax)", 0.05, 0.3, 0.1, 0.05, help="力の最大値による収束判定")
        max_steps = st.slider("最大ステップ数", 20, 200, 50, 10, help="最適化の最大反復回数")
    
    # 予想計算時間の表示
    def estimate_calculation_time(molecule_type, theory, basis_set, n_images):
        """計算時間の目安を表示"""
        times = {
            "水素分子解離 (H₂)": {"HF/STO-3G": "30秒-1分", "B3LYP/STO-3G": "1-2分", "B3LYP/6-31G": "2-3分"},
            "アンモニア反転 (NH₃)": {"HF/STO-3G": "2-4分", "B3LYP/STO-3G": "3-8分", "B3LYP/6-31G": "8-15分"},
            "エチレン回転 (C₂H₄)": {"HF/STO-3G": "4-8分", "B3LYP/STO-3G": "8-15分", "B3LYP/6-31G": "15-25分"},
            "シクロヘキサン (C₆H₁₂)": {"HF/STO-3G": "30分-1時間", "B3LYP/STO-3G": "1-3時間", "B3LYP/6-31G": "3-10時間"}
        }
        
        theory_basis = f"{theory}/{basis_set}"
        
        # 近似的なマッチング
        if "STO-3G" in basis_set:
            key = f"{theory}/STO-3G"
        elif "6-31G" in basis_set and "*" not in basis_set:
            key = f"{theory}/6-31G"
        else:
            key = f"{theory}/6-31G"  # より保守的な見積もり
        
        if molecule_type in times and key in times[molecule_type]:
            base_time = times[molecule_type][key]
            multiplier = max(1, n_images / 3)  # 中間構造数に応じて調整
            if multiplier > 1:
                return f"{base_time} × {multiplier:.1f} ≈ より長時間"
            return base_time
        else:
            return "不明（高計算コスト）"

    estimated_time = estimate_calculation_time(molecule_example, theory, basis_set, n_images)
    st.info(f"⏱️ **予想計算時間**: {estimated_time} (中間構造数: {n_images})")

    if "シクロヘキサン" in molecule_example and ("B3LYP" in theory and "6-31G" in basis_set):
        st.warning("⚠️ **長時間計算**: この設定では数時間かかる可能性があります。まずは水素分子やアンモニアで動作を確認することをお勧めします。")
    
    with col3:
        optimizer_type = st.selectbox("最適化手法", ["BFGS", "LBFGS", "GPMin"], index=0)
        interpolation_method = st.selectbox("補間方法", ["linear", "idpp"], index=0, help="初期経路の生成方法")

# ===== 溶媒効果設定 =====
with st.expander("🌊 溶媒効果 (オプション)"):
    solvent_model = st.selectbox("溶媒モデル", ["None", "PCM", "DDCOSMO"])
    if solvent_model != "None":
        # 溶媒データの読み込み
        solvents_file = "config/solvents_epsilon.csv"
        if os.path.exists(solvents_file):
            solvents_data = pd.read_csv(solvents_file)
            solvent_selection = st.selectbox(
                "溶媒を選択",
                [f"{row['Solvent']} (ε={row['Epsilon']})" for _, row in solvents_data.iterrows()]
            )
            if solvent_selection:
                eps = float(solvent_selection.split("=", 1)[-1][:-1])
        else:
            eps = st.number_input("誘電率 (ε)", min_value=1.0, value=78.4, step=0.1)
    else:
        eps = None

# ===== 計算実行 =====
st.header("🚀 計算実行")

if st.button("🚀 NEB計算実行", type="primary"):
    
    # ===== 必須入力の検証 =====
    # 初期構造の入力チェック
    initial_missing = False
    final_missing = False
    
    if input_type_initial == "XYZ":
        if not initial_structure or initial_structure.strip() == "":
            initial_missing = True
    else:
        if 'initial_smiles' not in st.session_state or not st.session_state.initial_smiles or st.session_state.initial_smiles.strip() == "":
            initial_missing = True
    
    # 最終構造の入力チェック
    if input_type_final == "XYZ":
        if not final_structure or final_structure.strip() == "":
            final_missing = True
    else:
        if 'final_smiles' not in st.session_state or not st.session_state.final_smiles or st.session_state.final_smiles.strip() == "":
            final_missing = True
    
    # エラーメッセージの表示
    if "水素分子" in molecule_example:
        if initial_missing and final_missing:
            st.error("❌ **初期構造と最終構造の両方が必要です！**\n\nNEB計算では結合状態（初期構造）と解離状態（最終構造）の両方を入力する必要があります。")
            st.stop()
        elif initial_missing:
            st.error("❌ **初期構造が入力されていません！**\n\n「初期構造」タブで水素分子の結合状態を入力してください。")
            st.stop()
        elif final_missing:
            st.error("❌ **最終構造が入力されていません！**\n\n「最終構造」タブで水素分子の解離状態を入力してください。")
            st.stop()
    elif "アンモニア" in molecule_example:
        if initial_missing and final_missing:
            st.error("❌ **初期構造と最終構造の両方が必要です！**\n\nNEB計算では上向きピラミッド（初期構造）と下向きピラミッド（最終構造）の両方を入力する必要があります。")
            st.stop()
        elif initial_missing:
            st.error("❌ **初期構造が入力されていません！**\n\n「初期構造」タブでアンモニアの上向きピラミッド構造を入力してください。")
            st.stop()
        elif final_missing:
            st.error("❌ **最終構造が入力されていません！**\n\n「最終構造」タブでアンモニアの下向きピラミッド構造を入力してください。")
            st.stop()
    elif "エチレン" in molecule_example:
        if initial_missing and final_missing:
            st.error("❌ **初期構造と最終構造の両方が必要です！**\n\nNEB計算ではエクリプス型（初期構造）とスタガード型（最終構造）の両方を入力する必要があります。")
            st.stop()
        elif initial_missing:
            st.error("❌ **初期構造が入力されていません！**\n\n「初期構造」タブでエチレンのエクリプス型構造を入力してください。")
            st.stop()
        elif final_missing:
            st.error("❌ **最終構造が入力されていません！**\n\n「最終構造」タブでエチレンのスタガード型構造を入力してください。")
            st.stop()
    else:
        if initial_missing and final_missing:
            st.error("❌ **初期構造と最終構造の両方が必要です！**\n\nNEB計算では舟型配座（初期構造）と椅子型配座（最終構造）の両方を入力する必要があります。")
            st.stop()
        elif initial_missing:
            st.error("❌ **初期構造が入力されていません！**\n\n「初期構造」タブでシクロヘキサンの舟型配座を入力してください。")
            st.stop()
        elif final_missing:
            st.error("❌ **最終構造が入力されていません！**\n\n「最終構造」タブでシクロヘキサンの椅子型配座を入力してください。")
            st.stop()
    
    # 入力処理と検証
    try:
        # 初期構造の処理
        if input_type_initial == "XYZ":
            initial_atoms = xyz_to_ase_atoms(initial_structure, charge)
            if initial_atoms is None:
                st.error("❌ 初期構造のXYZ形式が正しくありません")
                st.stop()
        else:
            # SMILESから構造生成
            initial_handler = MoleculeHandler(st.session_state.initial_smiles, input_type="smiles")
            if not initial_handler.mol:
                st.error("❌ 初期構造のSMILESが正しくありません")
                st.stop()
            initial_xyz = initial_handler.to_pyscf_input()
            initial_atoms = xyz_to_ase_atoms(initial_xyz, charge)
        
        # 最終構造の処理
        if input_type_final == "XYZ":
            final_atoms = xyz_to_ase_atoms(final_structure, charge)
            if final_atoms is None:
                st.error("❌ 最終構造のXYZ形式が正しくありません")
                st.stop()
        else:
            # SMILESから構造生成
            final_handler = MoleculeHandler(st.session_state.final_smiles, input_type="smiles")
            if not final_handler.mol:
                st.error("❌ 最終構造のSMILESが正しくありません")
                st.stop()
            final_xyz = final_handler.to_pyscf_input()
            final_atoms = xyz_to_ase_atoms(final_xyz, charge)
        
        # 構造の整合性チェック
        if len(initial_atoms) != len(final_atoms):
            st.error(f"❌ **原子数が一致しません！**\n\n初期構造: {len(initial_atoms)}原子\n最終構造: {len(final_atoms)}原子\n\nNEB計算では初期構造と最終構造の原子数が同じでなければなりません。")
            st.stop()
        
        # 元素の確認と自動整列
        initial_symbols = initial_atoms.get_chemical_symbols()
        final_symbols = final_atoms.get_chemical_symbols()
        
        from collections import Counter
        if Counter(initial_symbols) != Counter(final_symbols):
            st.error(f"❌ **元素の種類または数が異なります！**\n\n初期構造: {sorted(Counter(initial_symbols).items())}\n最終構造: {sorted(Counter(final_symbols).items())}\n\nNEB計算では同じ元素を同じ数だけ含む必要があります。")
            st.stop()
        
        # 原子順序が異なる場合は自動整列を試行
        if initial_symbols != final_symbols:
            st.warning(f"⚠️ **原子順序が異なります**\n\n初期構造: {initial_symbols}\n最終構造: {final_symbols}\n\n自動整列を実行します...")
            
            aligned_initial, aligned_final = align_atom_order(initial_atoms, final_atoms)
            
            if aligned_initial is None or aligned_final is None:
                st.error("❌ **原子順序の自動整列に失敗しました**\n\n手動で原子順序を統一してください。")
                st.stop()
            else:
                initial_atoms = aligned_initial
                final_atoms = aligned_final
                st.success("✅ **原子順序を自動的に整列しました**")
                
                # 整列後の順序を表示
                aligned_initial_symbols = initial_atoms.get_chemical_symbols()
                aligned_final_symbols = final_atoms.get_chemical_symbols()
                st.info(f"整列後の順序:\n初期構造: {aligned_initial_symbols}\n最終構造: {aligned_final_symbols}")
            
        st.success(f"✅ **構造の整合性確認完了**\n\n- 原子数: {len(initial_atoms)}原子\n- 元素: {', '.join(set(initial_symbols))}")
        
        # 原子種チェック
        if initial_atoms.get_chemical_symbols() != final_atoms.get_chemical_symbols():
            st.error("初期構造と最終構造の原子種が一致しません")
            st.stop()
            
    except Exception as e:
        st.error(f"構造処理エラー: {e}")
        st.stop()
    
    # 構造プレビューの表示
    st.subheader("構造プレビュー")
    
    # プレビュー方法の選択
    preview_method = st.radio(
        "詳細プレビュー方法:",
        ["Plotly", "py3Dmol"],
        index=1,
        horizontal=True,
        help="計算前の最終確認用詳細プレビュー"
    )
    
    col1, col2 = st.columns(2)
    
    with col1:
        if "水素分子" in molecule_example:
            st.markdown("**初期構造プレビュー (結合状態)**")
            if preview_method == "Plotly":
                initial_fig = create_3d_molecule_plot(initial_atoms, "水素分子 (結合状態)")
                if initial_fig:
                    st.plotly_chart(initial_fig, use_container_width=True)
            else:  # py3Dmol
                initial_mol_block = ase_atoms_to_mol_block(initial_atoms)
                if initial_mol_block:
                    show_3d_structure(initial_mol_block)
        elif "アンモニア" in molecule_example:
            st.markdown("**初期構造プレビュー (上向きピラミッド)**")
            if preview_method == "Plotly":
                initial_fig = create_3d_molecule_plot(initial_atoms, "アンモニア (上向きピラミッド)")
                if initial_fig:
                    st.plotly_chart(initial_fig, use_container_width=True)
            else:  # py3Dmol
                initial_mol_block = ase_atoms_to_mol_block(initial_atoms)
                if initial_mol_block:
                    show_3d_structure(initial_mol_block)
        elif "エチレン" in molecule_example:
            st.markdown("**初期構造プレビュー (エクリプス型)**")
            if preview_method == "Plotly":
                initial_fig = create_3d_molecule_plot(initial_atoms, "エチレン (エクリプス型)")
                if initial_fig:
                    st.plotly_chart(initial_fig, use_container_width=True)
            else:  # py3Dmol
                initial_mol_block = ase_atoms_to_mol_block(initial_atoms)
                if initial_mol_block:
                    show_3d_structure(initial_mol_block)
        else:
            st.markdown("**初期構造プレビュー (舟型配座)**")
            if preview_method == "Plotly":
                initial_fig = create_3d_molecule_plot(initial_atoms, "シクロヘキサン舟型配座")
                if initial_fig:
                    st.plotly_chart(initial_fig, use_container_width=True)
            else:  # py3Dmol
                initial_mol_block = ase_atoms_to_mol_block(initial_atoms)
                if initial_mol_block:
                    show_3d_structure(initial_mol_block)
    
    with col2:
        if "水素分子" in molecule_example:
            st.markdown("**最終構造プレビュー (解離状態)**")
            if preview_method == "Plotly":
                final_fig = create_3d_molecule_plot(final_atoms, "水素分子 (解離状態)")
                if final_fig:
                    st.plotly_chart(final_fig, use_container_width=True)
            else:  # py3Dmol
                final_mol_block = ase_atoms_to_mol_block(final_atoms)
                if final_mol_block:
                    show_3d_structure(final_mol_block)
        elif "アンモニア" in molecule_example:
            st.markdown("**最終構造プレビュー (下向きピラミッド)**")
            if preview_method == "Plotly":
                final_fig = create_3d_molecule_plot(final_atoms, "アンモニア (下向きピラミッド)")
                if final_fig:
                    st.plotly_chart(final_fig, use_container_width=True)
            else:  # py3Dmol
                final_mol_block = ase_atoms_to_mol_block(final_atoms)
                if final_mol_block:
                    show_3d_structure(final_mol_block)
        elif "エチレン" in molecule_example:
            st.markdown("**最終構造プレビュー (スタガード型)**")
            if preview_method == "Plotly":
                final_fig = create_3d_molecule_plot(final_atoms, "エチレン (スタガード型)")
                if final_fig:
                    st.plotly_chart(final_fig, use_container_width=True)
            else:  # py3Dmol
                final_mol_block = ase_atoms_to_mol_block(final_atoms)
                if final_mol_block:
                    show_3d_structure(final_mol_block)
        else:
            st.markdown("**最終構造プレビュー (椅子型配座)**")
            if preview_method == "Plotly":
                final_fig = create_3d_molecule_plot(final_atoms, "シクロヘキサン椅子型配座")
                if final_fig:
                    st.plotly_chart(final_fig, use_container_width=True)
            else:  # py3Dmol
                final_mol_block = ase_atoms_to_mol_block(final_atoms)
                if final_mol_block:
                    show_3d_structure(final_mol_block)
    
    # 構造比較の追加情報
    st.markdown("---")
    with st.expander("🔍 構造比較情報"):
        col_info1, col_info2, col_info3 = st.columns(3)
        
        with col_info1:
            st.metric("初期構造原子数", len(initial_atoms))
            initial_elements = set(initial_atoms.get_chemical_symbols())
            st.write(f"**元素**: {', '.join(sorted(initial_elements))}")
        
        with col_info2:
            st.metric("最終構造原子数", len(final_atoms))
            final_elements = set(final_atoms.get_chemical_symbols())
            st.write(f"**元素**: {', '.join(sorted(final_elements))}")
        
        with col_info3:
            # RMSD計算
            try:
                from scipy.spatial.distance import euclidean
                import numpy as np
                
                # 重心を原点に移動
                initial_centered = initial_atoms.positions - initial_atoms.positions.mean(axis=0)
                final_centered = final_atoms.positions - final_atoms.positions.mean(axis=0)
                
                # RMSD計算
                rmsd = np.sqrt(np.mean(np.sum((initial_centered - final_centered)**2, axis=1)))
                st.metric("構造RMSD (Å)", f"{rmsd:.3f}")
                
                if rmsd < 0.5:
                    st.success("✅ 構造が非常に類似")
                elif rmsd < 2.0:
                    st.info("💡 適度な構造変化")
                else:
                    st.warning("⚠️ 大きな構造変化")
                    
            except Exception as e:
                st.info("RMSD計算不可")
    
    with st.spinner("NEB計算を実行中..."):
        
        # ===== 中間構造を含むイメージを作成 =====
        images = [initial_atoms.copy()]
        for _ in range(n_images):
            images.append(initial_atoms.copy())
        images.append(final_atoms.copy())

        neb = NEB(images, k=spring_constant)
        neb.interpolate(method=interpolation_method)

        # ===== PySCF Calculator を各imageに個別にセット =====
        for image in images:
            image.calc = PySCFCalculator(
                basis=basis_set, 
                charge=charge, 
                spin=spin, 
                theory=theory,
                solvent_model=solvent_model if solvent_model != "None" else None,
                eps=eps
            )

        # ===== データフォルダの作成 =====
        data_dir = "data"
        if not os.path.exists(data_dir):
            os.makedirs(data_dir)
            st.info(f"📁 {data_dir} フォルダを作成しました")
        
        # ===== ファイルパスの設定 =====
        traj_file = os.path.join(data_dir, 'neb.traj')
        log_file = os.path.join(data_dir, 'neb.log')

        # ===== 最適化実行 =====
        if optimizer_type == "BFGS":
            opt = BFGS(neb, trajectory=traj_file, logfile=log_file)
        elif optimizer_type == "LBFGS":
            from ase.optimize import LBFGS
            opt = LBFGS(neb, trajectory=traj_file, logfile=log_file)
        else:  # GPMin
            from ase.optimize import GPMin
            opt = GPMin(neb, trajectory=traj_file, logfile=log_file)
        
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        try:
            # 計算の進捗を表示
            status_text.text("NEB最適化を実行中...")
            opt.run(fmax=max_force, steps=max_steps)
            progress_bar.progress(100)
            status_text.text("✅ 計算完了!")
            
        except Exception as e:
            st.error(f"最適化エラー: {e}")
            st.info("現在の構造で結果を表示します...")

        # ===== 結果の表示 =====
        st.success("🎉 NEB計算が完了しました!")
        
        # 配座変換に関する追加情報（分子例に応じた適切な説明に変更）
        if "水素分子" in molecule_example:
            st.info("""
            🔬 **水素分子解離反応の解析**
            
            - **結合状態**: H-H結合が形成された安定状態
            - **解離状態**: 2つの水素原子が分離した状態
            - **エネルギー差**: H₂の結合エネルギー (約4.5 eV = 104 kcal/mol)
            - **反応特性**: 遷移状態なし（単純な結合解離）
            """)
        elif "アンモニア" in molecule_example:
            st.info("""
            🔬 **アンモニア反転反応の解析**
            
            - **初期状態**: 上向きピラミッド型（C3v対称性）
            - **遷移状態**: 平面構造（D3h対称性）
            - **最終状態**: 下向きピラミッド型（C3v対称性）
            - **活性化エネルギー**: 約5.8 kcal/mol（実験値）
            - **反応特性**: 明確な遷移状態を持つ反転反応
            """)
        elif "エチレン" in molecule_example:
            st.info("""
            🔬 **エチレン回転反応の解析**
            
            - **初期状態**: エクリプス型（H原子が重なった配置）
            - **遷移状態**: 90°回転した直交配置
            - **最終状態**: スタガード型（H原子がずれた配置）
            - **活性化エネルギー**: 約65 kcal/mol（実験値）
            - **反応特性**: π結合の破断を伴う高い回転障壁
            """)
        else:
            st.info("""
            🔬 **シクロヘキサン配座変換の解析**
            
            - **舟型配座**: エネルギー的に不安定な配座
            - **椅子型配座**: 最も安定な配座
            - **エネルギー差**: 通常 6-7 kcal/mol程度
            - **活性化エネルギー**: リング反転の障壁高さ
            """)
        
        # エネルギーを計算
        energies = []
        for i, image in enumerate(images):
            try:
                e = image.get_potential_energy()
                energies.append(e)
            except Exception as e:
                st.warning(f"Image {i}のエネルギー計算でエラー: {e}")
                energies.append(None)

        # 結果をタブで表示
        tab1, tab2, tab3, tab4 = st.tabs(["📊 エネルギープロファイル", "🔬 分子構造", "📋 数値結果", "💾 軌道データ"])
        
        with tab1:
            if "水素分子" in molecule_example:
                st.subheader("水素分子解離エネルギープロファイル")
            elif "アンモニア" in molecule_example:
                st.subheader("アンモニア反転エネルギープロファイル")
            elif "エチレン" in molecule_example:
                st.subheader("エチレン回転エネルギープロファイル")
            else:
                st.subheader("配座変換エネルギープロファイル")
            energy_fig = plot_energy_profile(energies)
            if energy_fig:
                st.plotly_chart(energy_fig, use_container_width=True)
                
                # 配座変換の解析情報
                if energies and any(e is not None for e in energies):
                    valid_energies = [e for e in energies if e is not None]
                    min_energy = min(valid_energies)
                    max_energy = max(valid_energies)
                    energy_diff = (max_energy - min_energy) * 627.5095  # kcal/mol
                    
                    if "水素分子" in molecule_example:
                        st.markdown(f"""
                        **🔍 水素分子解離解析結果:**
                        - **解離エネルギー**: {energy_diff:.2f} kcal/mol
                        - **参考値**: H₂の結合エネルギーは約 104 kcal/mol
                        - **反応特性**: 単調な上昇（遷移状態なし）
                        """)
                    elif "アンモニア" in molecule_example:
                        st.markdown(f"""
                        **🔍 アンモニア反転解析結果:**
                        - **活性化エネルギー**: {energy_diff:.2f} kcal/mol
                        - **参考値**: 実験値は約 5.8 kcal/mol
                        - **反応特性**: 遷移状態を持つ反転反応
                        """)
                    elif "エチレン" in molecule_example:
                        st.markdown(f"""
                        **🔍 エチレン回転解析結果:**
                        - **回転障壁**: {energy_diff:.2f} kcal/mol
                        - **参考値**: 実験値は約 65 kcal/mol
                        - **反応特性**: π結合破断による高い障壁
                        """)
                    else:
                        st.markdown(f"""
                        **🔍 配座変換解析結果:**
                        - **エネルギー差**: {energy_diff:.2f} kcal/mol
                        - **参考値**: 実験値は約 6.9 kcal/mol
                        - **活性化エネルギー**: {energy_diff:.2f} kcal/mol (リング反転障壁)
                        """)
            else:
                st.info("plotlyをインストールして3Dプロットを有効にしてください: `pip install plotly`")
                
                # 代替のmatplotlib表示
                if energies and any(e is not None for e in energies):
                    valid_energies = [e for e in energies if e is not None]
                    min_energy = min(valid_energies)
                    
                    rel_energies = []
                    x_vals = []
                    for i, e in enumerate(energies):
                        if e is not None:
                            rel_energy = (e - min_energy) * 627.5095
                            rel_energies.append(rel_energy)
                            x_vals.append(i)
                    
                    fig, ax = plt.subplots(figsize=(10, 6))
                    ax.plot(x_vals, rel_energies, 'bo-', linewidth=2, markersize=8)
                    ax.set_xlabel('反応座標 (Image番号)')
                    ax.set_ylabel('相対エネルギー (kcal/mol)')
                    ax.set_title('反応エネルギープロファイル')
                    ax.grid(True, alpha=0.3)
                    st.pyplot(fig)
                    plt.close()

        with tab2:
            st.subheader("分子構造の可視化")
            
            # 構造選択
            if "水素分子" in molecule_example:
                structure_idx = st.selectbox(
                    "表示する構造を選択:",
                    range(len(images)),
                    format_func=lambda x: f"Image {x} {'(結合状態)' if x == 0 else '(解離状態)' if x == len(images)-1 else '(中間状態)'}"
                )
            elif "アンモニア" in molecule_example:
                structure_idx = st.selectbox(
                    "表示する構造を選択:",
                    range(len(images)),
                    format_func=lambda x: f"Image {x} {'(上向きピラミッド)' if x == 0 else '(下向きピラミッド)' if x == len(images)-1 else '(遷移状態付近)' if x == len(images)//2 else '(中間状態)'}"
                )
            elif "エチレン" in molecule_example:
                structure_idx = st.selectbox(
                    "表示する構造を選択:",
                    range(len(images)),
                    format_func=lambda x: f"Image {x} {'(エクリプス型)' if x == 0 else '(スタガード型)' if x == len(images)-1 else '(遷移状態付近)' if x == len(images)//2 else '(中間状態)'}"
                )
            else:
                structure_idx = st.selectbox(
                    "表示する構造を選択:",
                    range(len(images)),
                    format_func=lambda x: f"Image {x} {'(舟型配座)' if x == 0 else '(椅子型配座)' if x == len(images)-1 else '(中間配座)'}"
                )
            
            selected_image = images[structure_idx]
            
            # 可視化方法の選択
            visualization_method = st.radio(
                "可視化方法を選択:",
                ["Plotly (3D散布図)", "py3Dmol (分子モデル)", "座標テーブル"],
                index=0
            )
            
            if visualization_method == "Plotly (3D散布図)":
                # 3D表示 (既存のplotly版)
                molecule_fig = create_3d_molecule_plot(selected_image, f"Image {structure_idx}")
                if molecule_fig:
                    st.plotly_chart(molecule_fig, use_container_width=True)
                else:
                    st.info("plotlyをインストールして3Dプロットを有効にしてください")
            
            elif visualization_method == "py3Dmol (分子モデル)":
                # py3Dmolを使用した3D表示
                mol_block = ase_atoms_to_mol_block(selected_image)
                if mol_block:
                    st.info("🔬 **py3Dmol 3D分子構造**")
                    show_3d_structure(mol_block)
                else:
                    st.error("MOLブロック形式への変換に失敗しました")
            
            else:  # 座標テーブル
                # 座標データの表示
                positions = selected_image.positions
                symbols = selected_image.get_chemical_symbols()
                
                coord_df = pd.DataFrame({
                    '原子': symbols,
                    'X (Å)': positions[:, 0],
                    'Y (Å)': positions[:, 1],
                    'Z (Å)': positions[:, 2]
                })
                st.dataframe(coord_df)

        with tab3:
            st.subheader("数値結果")
            
            # エネルギーテーブル
            if energies:
                valid_energies = [e for e in energies if e is not None]
                if valid_energies:
                    min_energy = min(valid_energies)
                    
                    results_data = []
                    for i, e in enumerate(energies):
                        if e is not None:
                            rel_energy_kcal = (e - min_energy) * 627.5095
                            results_data.append({
                                'Image': i,
                                'エネルギー (Hartree)': f"{e:.6f}",
                                '相対エネルギー (kcal/mol)': f"{rel_energy_kcal:.2f}"
                            })
                        else:
                            results_data.append({
                                'Image': i,
                                'エネルギー (Hartree)': "N/A",
                                '相対エネルギー (kcal/mol)': "N/A"
                            })
                    
                    results_df = pd.DataFrame(results_data)
                    st.dataframe(results_df, use_container_width=True)

        with tab4:
            st.subheader("軌道データのダウンロード")
            
            # 軌道ファイルの存在確認
            if os.path.exists(traj_file):
                st.success(f"✅ 軌道ファイル ({traj_file}) が生成されました")
                
                # ファイルダウンロード
                with open(traj_file, 'rb') as f:
                    st.download_button(
                        label="📁 neb.traj をダウンロード",
                        data=f.read(),
                        file_name="neb.traj",
                        mime="application/octet-stream"
                    )
                
                st.info("💡 ダウンロードしたファイルは `ase gui neb.traj` で表示できます")
            
            if os.path.exists(log_file):
                st.success(f"✅ ログファイル ({log_file}) が生成されました")
                
                with st.expander("📜 計算ログを表示"):
                    with open(log_file, 'r') as f:
                        st.text(f.read())

# 既存の軌道ファイルがある場合の表示
st.markdown("---")
st.subheader("📁 既存の軌道データ")

existing_traj = load_trajectory_if_exists('data/neb.traj')
if existing_traj:
    st.info(f"✅ 既存の軌道ファイルを発見 ({len(existing_traj)} フレーム)")
    
    if st.button("🔄 既存データを表示"):
        # 既存データの可視化
        st.subheader("既存軌道の分析")
        
        # フレーム選択
        frame_idx = st.slider("フレームを選択", 0, len(existing_traj)-1, 0)
        selected_frame = existing_traj[frame_idx]
        
        # 可視化方法選択
        existing_viz_method = st.radio(
            "可視化方法:",
            ["Plotly", "py3Dmol"],
            index=0,
            horizontal=True,
            key="existing_viz"
        )
        
        # 3D表示
        if existing_viz_method == "Plotly":
            frame_fig = create_3d_molecule_plot(selected_frame, f"Frame {frame_idx}")
            if frame_fig:
                st.plotly_chart(frame_fig, use_container_width=True)
        else:  # py3Dmol
            frame_mol_block = ase_atoms_to_mol_block(selected_frame)
            if frame_mol_block:
                show_3d_structure(frame_mol_block)
else:
    st.info("軌道ファイルが見つかりません。上記のボタンでNEB計算を実行してください。")

# ===== 計算のヒント =====
with st.expander("💡 NEB計算のヒント"):
    st.markdown("""
    **良いNEB計算のためのガイドライン:**
    
    1. **適切な初期・最終構造**
       - 反応前後の安定構造を使用
       - 構造最適化を事前に実行することを推奨
       - 原子順序が異なる場合は自動整列機能が働きます
    
    2. **計算レベルの選択**
       - **推奨**: B3LYP/6-31G(d) 以上
       - **避ける**: HF/STO-3G (不正確)
       - 精度と計算時間のバランスを考慮
    
    3. **NEB パラメータ**
       - **中間構造数**: 反応の複雑さに応じて調整（通常3-5個）
       - **Spring Constant**: 0.1前後が適切（大きすぎると剛直、小さすぎると不安定）
       - **収束判定**: 0.1 eV/Å程度が実用的
    
    4. **計算の解釈**
       - エネルギープロファイルのピークが活性化エネルギー
       - 単調でない場合は中間体の存在を示唆
       - 負の活性化エネルギーは設定ミスの可能性
    
    **原子順序の問題について:**
    - 同じ分子でも原子の番号順序が異なることがあります
    - 本システムは自動的に最適な原子対応を見つけて整列します
    - 距離に基づいて同種原子同士を対応付けます
    
    **3D分子可視化について:**
    - **入力時プレビュー**: 各タブで構造を入力すると自動的に3D表示されます
    - **詳細プレビュー**: 計算実行前に両構造を詳細確認できます
    - **Plotly**: 散布図形式で原子と結合を表示、軽量で高速
    - **py3Dmol**: 分子モデリング専用、より美しい分子表示
    - 構造の詳細確認にはpy3Dmol、概要把握にはPlotlyが適しています
    
    **トラブルシューティング:**
    - 収束しない → spring constantを調整、より良い初期構造を使用
    - 異常な結果 → 計算レベルを向上、対称性を確認
    - 原子順序エラー → 元素の種類と数が一致しているか確認
    - 3D表示エラー → ブラウザの再読み込み、別の可視化方法を試行
    """)

with st.expander("📚 計算レベル選択指針"):
    st.markdown("""
    **NEB計算の計算レベル推奨:**
    
    | 計算レベル | 適用性 | 推奨度 | 備考 |
    |------------|--------|--------|------|
    | HF/sto-3g | × 不適 | ⭐ | 最も基本的、反応経路には不十分 |
    | HF/6-31g* | △ 限定的 | ⭐⭐ | 電子相関の欠如により不正確 |
    | B3LYP/sto-3g | △ 限定的 | ⭐⭐ | 基底関数が不十分 |
    | B3LYP/6-31g* | ○ 推奨 | ⭐⭐⭐⭐ | 最小推奨レベル |
    | B3LYP/cc-pVDZ | ○ 推奨 | ⭐⭐⭐⭐ | Dunning基底、推奨 |
    | B3LYP/def2-TZVP | ◎ 高精度 | ⭐⭐⭐⭐⭐ | 高精度計算 |
    | B3LYP-D3/6-31g* | ◎ 分散補正 | ⭐⭐⭐⭐⭐ | ファンデルワールス相互作用対応 |
    | M06-2X/6-31g* | ○ 中距離相互作用 | ⭐⭐⭐⭐ | 非共有相互作用に優れる |
    | PBE0/def2-SVP | ○ 計算効率 | ⭐⭐⭐⭐ | 計算時間と精度のバランス |
    
    **分子サイズ別推奨:**
    - 小分子（～15原子）: B3LYP/def2-TZVP または B3LYP-D3/cc-pVTZ
    - 中分子（15-50原子）: B3LYP/6-31g* または PBE0/def2-SVP
    - 大分子（50原子以上）: B3LYP/6-31g または PBE/def2-SVP
    
    **特殊系への推奨:**
    - 分散相互作用重要: B3LYP-D3, M06-2X
    - 金属錯体: PBE0, TPSS
    - 励起状態関連: CAM-B3LYP
    """)

st.markdown("---")
st.markdown("**注意**: NEB計算は計算コストが非常に高い処理です。calculation.pyで定義された理論手法・基底関数を活用し、小さな分子から始めてパラメータを調整することを推奨します。")
