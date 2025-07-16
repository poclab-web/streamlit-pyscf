"""
MOPAC計算結果の分子軌道可視化用ロジック

機能:
- arcファイルまたはoutファイルから分子軌道データを読み取り
- HOMO/LUMO軌道の3D可視化
- 軌道エネルギー図の表示
- 電子密度表面の可視化
- 分子軌道係数の解析
"""

import streamlit as st
import pandas as pd
import numpy as np
import os
import tempfile
import matplotlib.pyplot as plt
import py3Dmol
import stmol
from pathlib import Path
import re
import io

from logic.mopac_calculation import MopacCalculator
from logic.molecule_handler import MoleculeHandler


@st.cache_resource
def get_mopac_calculator(xyz_data):
    """
    MopacCalculatorインスタンスを作成（キャッシュ付き）
    """
    molecule_handler = MoleculeHandler(xyz_data, input_type="xyz")
    return MopacCalculator(molecule_handler)


def list_out_files(data_dir="data"):
    """
    dataフォルダ内のMOPACファイル（mopac_workディレクトリ内のみ）を再帰的に探索
    
    Args:
        data_dir: 検索するディレクトリ
        
    Returns:
        dict: ファイル種別をキーとするファイルパスリスト
    """
    files = {"out": []}
    for root, _, file_list in os.walk(data_dir):
        # mopac_workディレクトリ内のファイルのみを対象とする
        if "mopac_work" in root:
            for f in file_list:
                file_path = os.path.join(root, f)
                if f.endswith(".out"):
                    files["out"].append(file_path)
    return files


def display_file_info(file_path):
    """
    ファイルの基本情報を表示
    
    Args:
        file_path: ファイルのパス
    """
    file_size = os.path.getsize(file_path)
    file_dir = os.path.dirname(file_path)
    file_name = os.path.basename(file_path)
    
    st.write(f"**ファイル名**: {file_name}")
    st.write(f"**ディレクトリ**: {file_dir}")
    st.write(f"**ファイルサイズ**: {file_size:,} bytes")
    
    # ファイルの最初の数行を表示（プレビュー）
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            preview_lines = []
            for i, line in enumerate(f):
                if i >= 10:  # 最初の10行のみ
                    break
                preview_lines.append(line.rstrip())
        
        if preview_lines:
            st.text_area(
                "ファイルプレビュー（最初の10行）",
                "\n".join(preview_lines),
                height=200,
                disabled=True
            )
    except Exception as e:
        st.warning(f"ファイルプレビューの読み取りに失敗しました: {e}")


def parse_mopac_out_file(out_file_path):
    """
    MOPACのoutファイルから詳細な分子軌道情報を抽出
    
    Args:
        out_file_path: outファイルのパス
        
    Returns:
        dict: 分子軌道データ
    """
    data = {
        'coordinates': [],
        'atoms': [],
        'orbital_energies': [],
        'orbital_coefficients': [],
        'homo_index': None,
        'total_energy': None,
        'dipole_moment': None,
        'heat_of_formation': None
    }
    
    try:
        with open(out_file_path, 'r') as f:
            content = f.read()
        
        # 最終座標の抽出
        coord_section = re.search(r'CARTESIAN COORDINATES\s*\n\s*NO\.\s+ATOM\s+X\s+Y\s+Z\s*\n(.*?)\n\s*\n', content, re.DOTALL)
        if coord_section:
            coord_lines = coord_section.group(1).strip().split('\n')
            for line in coord_lines:
                parts = line.strip().split()
                if len(parts) >= 5:
                    atom = parts[1]
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    data['atoms'].append(atom)
                    data['coordinates'].append([x, y, z])
        
        # 軌道エネルギーの抽出
        energy_section = re.search(r'EIGENVALUES\s*\n(.*?)(?=\n\s*[A-Z]|\n\s*$)', content, re.DOTALL)
        if energy_section:
            energy_text = energy_section.group(1)
            energies = re.findall(r'[-+]?\d+\.\d+', energy_text)
            data['orbital_energies'] = [float(e) for e in energies]
        
        # 全エネルギーの抽出
        energy_match = re.search(r'TOTAL ENERGY\s*=\s*([-+]?\d+\.\d+)', content)
        if energy_match:
            data['total_energy'] = float(energy_match.group(1))
        
        # 生成熱の抽出
        hof_match = re.search(r'HEAT OF FORMATION\s*=\s*([-+]?\d+\.\d+)', content)
        if hof_match:
            data['heat_of_formation'] = float(hof_match.group(1))
        
        # 双極子モーメントの抽出
        dipole_match = re.search(r'DIPOLE\s*=\s*([-+]?\d+\.\d+)', content)
        if dipole_match:
            data['dipole_moment'] = float(dipole_match.group(1))
        
        # HOMO軌道インデックスの推定
        if data['orbital_energies']:
            occupied_orbitals = [i for i, e in enumerate(data['orbital_energies']) if e < 0]
            if occupied_orbitals:
                data['homo_index'] = max(occupied_orbitals)
                
    except Exception as e:
        st.error(f"outファイルの解析中にエラーが発生しました: {e}")
        
    return data


def create_molecular_orbital_plot(orbital_energies, homo_index=None):
    """
    Create molecular orbital energy diagram
    
    Args:
        orbital_energies: List of orbital energies
        homo_index: HOMO orbital index
        
    Returns:
        matplotlib figure
    """
    fig, ax = plt.subplots(figsize=(8, 10))
    
    # Sort orbitals by energy (ascending order)
    sorted_indices = np.argsort(orbital_energies)
    sorted_energies = np.array(orbital_energies)[sorted_indices]
    
    # Create vertical energy diagram
    # X position for all orbitals (single column)
    x_center = 0.5
    line_width = 0.3
    
    # Plot energy levels as horizontal lines
    for i, energy in enumerate(sorted_energies):
        # Color coding: red for occupied (negative), blue for unoccupied (positive)
        color = 'red' if energy < 0 else 'blue'
        alpha = 0.8 if energy < 0 else 0.6
        
        # Draw horizontal line for each orbital
        ax.plot([x_center - line_width/2, x_center + line_width/2], 
                [energy, energy], 
                color=color, linewidth=3, alpha=alpha)
    
    # HOMO/LUMO markers
    if homo_index is not None:
        homo_energy = orbital_energies[homo_index]
        ax.axhline(y=homo_energy, color='darkred', linestyle='--', alpha=0.9, 
                  label='HOMO', linewidth=2)
        
        # Add HOMO label with orbital number
        ax.text(x_center + line_width/2 + 0.05, homo_energy, f'HOMO ({homo_index+1})', 
                fontsize=10, verticalalignment='center', fontweight='bold', color='darkred')
        
        # Find LUMO
        lumo_candidates = [e for e in orbital_energies if e > homo_energy]
        if lumo_candidates:
            lumo_energy = min(lumo_candidates)
            lumo_index = orbital_energies.index(lumo_energy)
            ax.axhline(y=lumo_energy, color='darkblue', linestyle='--', alpha=0.9, 
                      label='LUMO', linewidth=2)
            
            # Add LUMO label with orbital number
            ax.text(x_center + line_width/2 + 0.05, lumo_energy, f'LUMO ({lumo_index+1})', 
                    fontsize=10, verticalalignment='center', fontweight='bold', color='darkblue')
            
            # Calculate and display HOMO-LUMO gap
            gap = lumo_energy - homo_energy
            ax.text(0.7, (homo_energy + lumo_energy)/2, 
                   f'Gap: {gap:.2f} eV', 
                   fontsize=10, bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7))
    
    # Formatting
    ax.set_xlim(0, 1)
    ax.set_ylabel('Energy (eV)', fontsize=12)
    ax.set_title('Molecular Orbital Energy Diagram', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    ax.legend(loc='upper right')
    
    # Remove x-axis ticks and labels
    ax.set_xticks([])
    ax.set_xlabel('')
    
    # Add energy scale reference
    ax.text(0.05, ax.get_ylim()[1] * 0.95, 'Unstable', fontsize=10, fontweight='bold', 
           verticalalignment='top', color='blue')
    ax.text(0.05, ax.get_ylim()[0] * 0.95, 'Stable', fontsize=10, fontweight='bold', 
           verticalalignment='bottom', color='red')
    
    plt.tight_layout()
    return fig


def visualize_molecule_3d(atoms, coordinates, title="分子構造"):
    """
    分子の3D構造を可視化
    
    Args:
        atoms: 原子記号のリスト
        coordinates: 座標のリスト
        title: タイトル
        
    Returns:
        py3Dmolビューア
    """
    # XYZ形式の文字列を作成
    xyz_lines = [f"{len(atoms)}"]
    xyz_lines.append(title)
    
    for atom, coord in zip(atoms, coordinates):
        xyz_lines.append(f"{atom} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}")
    
    xyz_data = "\n".join(xyz_lines)
    
    # py3Dmolで可視化
    view = py3Dmol.view(width=600, height=400)
    view.addModel(xyz_data, 'xyz')
    view.setStyle({'stick': {'radius': 0.1}, 'sphere': {'scale': 0.3}})
    view.zoomTo()
    
    return view


def get_file_description(file_type):
    """
    ファイル種別の説明を取得
    
    Args:
        file_type: ファイルの拡張子（out, arc等）
        
    Returns:
        str: ファイルの説明
    """
    descriptions = {
        'out': 'メイン出力ファイル（計算結果詳細）',
        'arc': 'アーカイブファイル（構造最適化軌跡）',
        'mgf': 'Gaussian形式分子軌道ファイル',
        'mol': 'MDL Molファイル',
        'den': '電子密度ファイル',
        'esp': '静電ポテンシャルファイル',
        'gpt': 'グラフポイントファイル',
        'po': '分子軌道プロットファイル',
        'aux': '補助ファイル',
        'log': 'ログファイル'
    }
    return descriptions.get(file_type, '計算関連ファイル')


def process_mo_calculation(data, mo_theory, mo_charge, mo_multiplicity, 
                         include_graph, include_html, include_allvec):
    """
    分子軌道計算を実行し、結果を処理
    
    Args:
        data: 分子データ
        mo_theory: 計算理論
        mo_charge: 分子電荷
        mo_multiplicity: スピン多重度
        include_graph: GRAPH MOオプション
        include_html: HTML可視化オプション
        include_allvec: ALLVEC オプション
        
    Returns:
        dict: 計算結果
    """
    try:
        # XYZ形式の文字列を作成
        xyz_lines = [f"{len(data['atoms'])}"]
        xyz_lines.append("Generated from MOPAC output")
        
        for atom, coord in zip(data['atoms'], data['coordinates']):
            xyz_lines.append(f"{atom} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}")
        
        xyz_data = "\n".join(xyz_lines)
        
        # MoleculeHandlerとMopacCalculatorを作成
        try:
            calculator = get_mopac_calculator(xyz_data)
            
            # メソッドの存在確認
            if not hasattr(calculator, 'molecular_orbital_calculation'):
                return {'success': False, 'error': 'MopacCalculatorクラスにmolecular_orbital_calculationメソッドが見つかりません。'}
            
        except Exception as calc_error:
            # 代替方法でインスタンスを作成
            molecule_handler = MoleculeHandler(xyz_data, input_type="xyz")
            if molecule_handler.mol is None:
                return {'success': False, 'error': '分子構造の処理に失敗しました。'}
            calculator = MopacCalculator(molecule_handler)
        
        # 分子軌道計算を実行
        mo_result = calculator.molecular_orbital_calculation(
            theory=mo_theory,
            charge=mo_charge,
            multiplicity=mo_multiplicity,
            include_graph=include_graph,
            include_vectors=True,
            precise=True,
            title="MO_calculation_from_viz",
            # 可視化特化オプション
            include_html=include_html,
            include_allvec=include_allvec
        )
        
        return mo_result
        
    except Exception as e:
        return {'success': False, 'error': str(e)}


def create_orbital_dataframe(orbital_energies, homo_index=None):
    """
    軌道エネルギーのDataFrameを作成
    
    Args:
        orbital_energies: 軌道エネルギーのリスト
        homo_index: HOMO軌道のインデックス
        
    Returns:
        pandas.DataFrame: 軌道エネルギー表
    """
    orbital_df = pd.DataFrame({
        '軌道番号': range(1, len(orbital_energies) + 1),
        'エネルギー (eV)': orbital_energies,
        '占有状態': ['占有' if e < 0 else '非占有' for e in orbital_energies]
    })
    
    # HOMO/LUMOをハイライト
    if homo_index is not None:
        orbital_df.loc[homo_index, '軌道番号'] = f"{homo_index+1} (HOMO)"
        
        lumo_candidates = [e for e in orbital_energies if e > orbital_energies[homo_index]]
        if lumo_candidates:
            lumo_energy = min(lumo_candidates)
            lumo_index = orbital_energies.index(lumo_energy)
            orbital_df.loc[lumo_index, '軌道番号'] = f"{lumo_index+1} (LUMO)"
    
    return orbital_df
