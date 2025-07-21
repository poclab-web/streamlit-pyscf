import streamlit as st
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.interchange import Interchange
import openmm
from openmm import unit, app
from openmm.unit import kelvin, picoseconds
import py3Dmol
import stmol
from utils.module import load_css
from logic.molecule_handler import MoleculeHandler
from utils.openmm_ui import display_openmm_status, require_openmm, display_forcefield_selector

# カスタムCSSを適用
load_css("config/styles.css")

# Function to display 3D structure using py3Dmol
def show_3d_structure(mol_block):
    try:        
        viewer = py3Dmol.view(width=400, height=400)
        viewer.addModel(mol_block, "mol")  
        viewer.setStyle({"stick": {}})
        viewer.zoomTo()
        stmol.showmol(viewer, height=400)
    except Exception as e:
        st.warning(f"3D構造の生成に失敗しました: {e}")

def smiles_to_single_point_calculation(smiles: str, force_field_name="openff_unconstrained-2.0.0.offxml", temperature=300, debug=False):
    """OpenMMを使用してSMILESから1点エネルギー計算を実行"""
    if debug:
        st.info(f"デバッグモード: フォースフィールド={force_field_name}, 温度={temperature}K")
    
    # RDKitで初期構造を生成
    rdkit_mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(rdkit_mol)
    AllChem.UFFOptimizeMolecule(rdkit_mol)

    # シンプルな1点計算アプローチ
    try:
        if debug:
            st.info("1点エネルギー計算を開始...")
        
        # OpenFF Moleculeを作成
        off_mol = Molecule.from_rdkit(rdkit_mol, allow_undefined_stereo=True)
        
        if debug:
            st.info(f"OpenFF分子作成成功: 原子数={off_mol.n_atoms}")
        
        # トポロジーとフォースフィールドの設定
        topology = off_mol.to_topology()
        forcefield = ForceField(force_field_name)
        
        if debug:
            st.info(f"フォースフィールド読み込み成功: {force_field_name}")
        
        # OpenMMシステムを作成
        omm_system = forcefield.create_openmm_system(topology)
        
        if debug:
            st.info("OpenMMシステム作成成功")
        
        # 座標を取得（OpenMM用に変換）
        try:
            # RDKitのコンフォーマーから座標を取得（オングストローム単位）
            conf = rdkit_mol.GetConformer()
            positions_angstrom = []
            for i in range(rdkit_mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                positions_angstrom.append([pos.x, pos.y, pos.z])
            
            # オングストロームからナノメートルに変換
            positions_nm = np.array(positions_angstrom) / 10.0
            
            if debug:
                st.info(f"座標取得成功: 形状={positions_nm.shape}")
            
            # OpenMM Vec3リストに変換
            from openmm import Vec3
            positions = [Vec3(float(pos[0]), float(pos[1]), float(pos[2])) * unit.nanometer for pos in positions_nm]
            
            if debug:
                st.info(f"OpenMM座標変換成功: 要素数={len(positions)}")
            
        except Exception as e:
            if debug:
                st.error(f"座標変換失敗: {e}")
            raise
        
        # シンプルなコンテキスト作成（最小限）
        integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
        simulation = app.Simulation(topology.to_openmm(), omm_system, integrator)
        
        if debug:
            st.info("シミュレーション作成成功")
        
        # 座標を設定
        simulation.context.setPositions(positions)
        
        if debug:
            st.info("座標設定成功")
        
        # エネルギーのみを計算（座標は変更しない）
        state = simulation.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
        
        if debug:
            st.success(f"1点計算成功: エネルギー={energy:.2f} kcal/mol")
        
        return rdkit_mol, energy
        
    except Exception as e:
        if debug:
            st.error(f"1点計算失敗: {e}")
        raise Exception(f"1点エネルギー計算に失敗しました: {e}")

def compare_structures(original_mol, optimized_mol):
    """最適化前後の構造を比較表示"""
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("最適化前 (UFF)")
        original_block = Chem.MolToMolBlock(original_mol)
        show_3d_structure(original_block)
    
    with col2:
        st.subheader("最適化後 (OpenMM)")
        optimized_block = Chem.MolToMolBlock(optimized_mol)
        show_3d_structure(optimized_block)

def mol_to_3d_viewer(mol):
    """互換性のための関数（従来のコードとの互換性を保つ）"""
    mol_block = Chem.MolToMolBlock(mol)
    viewer = py3Dmol.view(width=400, height=400)
    viewer.addModel(mol_block, 'mol')
    viewer.setStyle({'stick': {}})
    viewer.zoomTo()
    return viewer

# --- Streamlit UI ---
st.title("💠 OpenMM 1点エネルギー計算アプリ")
st.markdown("SMILES を入力して、OpenMM で1点エネルギー計算を行います。")

# システム状態の確認と表示（新しいユーティリティを使用）
system_status = require_openmm()

# 計算設定
st.subheader("計算設定")
col1, col2 = st.columns(2)

with col1:
    force_field_option = display_forcefield_selector()
    temperature = st.slider("温度 (K)", 200, 400, 300, 10)

with col2:
    show_comparison = st.checkbox("初期構造との比較表示", value=True)
    debug_mode = st.checkbox("デバッグモード", value=True)  # デフォルトでオン

smiles = st.text_input("SMILESを入力してください", "CC(=O)Nc1ccc(O)cc1")  # アセトアミノフェン
run = st.button("1点計算実行")

if run and smiles:
    with st.spinner("1点エネルギー計算中..."):
        try:
            # MoleculeHandlerを使用して初期構造を生成
            handler = MoleculeHandler(smiles, input_type="smiles")
            original_mol = Chem.Mol(handler.mol)  # コピーを作成
            
            # OpenMMで1点計算
            mol, energy = smiles_to_single_point_calculation(smiles, force_field_option, temperature, debug_mode)
            st.success(f"1点計算エネルギー: {energy:.2f} kcal/mol")

            if show_comparison:
                col1, col2 = st.columns(2)
                
                with col1:
                    st.subheader("初期構造 (UFF最適化)")
                    original_block = Chem.MolToMolBlock(original_mol)
                    show_3d_structure(original_block)
                
                with col2:
                    st.subheader("同じ構造 (OpenMMエネルギー)")
                    mol_block = Chem.MolToMolBlock(mol)
                    show_3d_structure(mol_block)
            else:
                # 構造表示
                st.subheader("分子構造")
                mol_block = Chem.MolToMolBlock(mol)
                show_3d_structure(mol_block)
            
            # XYZ座標の表示
            st.subheader("XYZ座標")
            handler_result = MoleculeHandler(mol, input_type="rdkit")
            xyz_coordinates = handler_result.get_xyz_coordinates()
            xyz_text = "\n".join([f"{atom:<2} {x:>10.6f} {y:>10.6f} {z:>10.6f}" 
                                for atom, x, y, z in xyz_coordinates])
            st.text(xyz_text)
            
            # ダウンロードオプション
            st.subheader("ダウンロード")
            col1, col2, col3 = st.columns(3)
            
            with col1:
                xyz_data = handler_result.to_pyscf_input()
                st.download_button(
                    label="XYZ形式でダウンロード",
                    data=xyz_data,
                    file_name=f"{smiles}_single_point.xyz",
                    mime="text/plain"
                )
            
            with col2:
                mol_block = Chem.MolToMolBlock(mol)
                st.download_button(
                    label="SDF形式でダウンロード",
                    data=mol_block,
                    file_name=f"{smiles}_single_point.sdf",
                    mime="chemical/x-mdl-sdfile"
                )
            
            with col3:
                mopac_input = handler_result.to_mopac_input(
                    title=f"OpenMM single point: {smiles}",
                    keywords="PM7 PRECISE"
                )
                st.download_button(
                    label="MOPAC形式でダウンロード",
                    data=mopac_input,
                    file_name=f"{smiles}_single_point.mop",
                    mime="text/plain"
                )

        except Exception as e:
            st.error(f"エラーが発生しました: {e}")
            st.info("SMILESが正しいか、またはライブラリが正しくインストールされているか確認してください。")

elif run and not smiles:
    st.warning("SMILESを入力してください。")
