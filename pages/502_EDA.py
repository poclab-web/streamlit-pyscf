"""
個々のフラグメント（分子 A, 分子 B）のエネルギーを計算。
フラグメントが結合した複合体（A + B）のエネルギーを計算。
相互作用エネルギーを計算：
"""

import streamlit as st
from utils.module import load_css
import numpy as np

import py3Dmol  # 3D可視化用ライブラリ
import stmol

from logic.molecule_handler import MoleculeHandler
from controllers.energydecompositionanalysis import get_hf_energy_decomposition

from logic.calculation import theory_options, basis_set_options, run_quantum_calculation



# カスタムCSSを適用
load_css("config/styles.css")

st.title("エネルギー分解解析 (Energy Decomposition Analysis)")
st.markdown("分子間相互作用エネルギーを詳細に分解して解析します。")


# 分子入力セクション
st.header("分子入力")
col1, col2 = st.columns(2)

with col1:
    st.subheader("分子A")
    smiles_input1 = st.text_input("分子AのSMILES", "N", key="smiles1")
    
    st.write("🔧 配座生成設定")
    force_field1 = st.selectbox("Force Field", ["MMFF", "UFF"], key="force_field1")
    num_conformers1 = st.number_input("Conformers数", value=100, min_value=1, max_value=10000, key="num_conformers1")

with col2:
    st.subheader("分子B") 
    smiles_input2 = st.text_input("分子BのSMILES", "F", key="smiles2")
    
    st.write("🔧 配座生成設定")
    force_field2 = st.selectbox("Force Field", ["MMFF", "UFF"], key="force_field2")
    num_conformers2 = st.number_input("Conformers数", value=100, min_value=1, max_value=10000, key="num_conformers2")

# 分子配置制御オプション
st.header("分子配置制御")

# 配置方法の選択
placement_method = st.selectbox(
    "配置方法",
    ["最近接原子間距離", "C-H π相互作用", "特定原子間距離"],
    key="placement_method"
)

if placement_method == "最近接原子間距離":
    target_distance = st.slider("最近接原子間距離 (Å)", min_value=1.0, max_value=10.0, value=3.0, step=0.1)
    
elif placement_method == "C-H π相互作用":
    st.info("💡 C-H π相互作用の計算に適した配置を行います（分子Aがπ系、分子BがC-H結合）")
    
    # プリセット設定
    preset = st.selectbox(
        "プリセット設定",
        ["カスタム", "典型的なC-H π", "T字型配置", "平行配置"],
        key="ch_pi_preset"
    )
    
    # デフォルト値を設定
    default_distance = 3.5
    default_approach = 0.0
    default_rotation = 0.0
    
    if preset == "典型的なC-H π":
        default_distance = 3.5
        default_approach = 0.0
        default_rotation = 0.0
    elif preset == "T字型配置":
        default_distance = 3.0
        default_approach = 90.0
        default_rotation = 0.0
    elif preset == "平行配置":
        default_distance = 4.0
        default_approach = 90.0
        default_rotation = 90.0
    
    col1, col2 = st.columns(2)
    with col1:
        target_distance = st.slider("重心間距離 (Å)", min_value=2.0, max_value=8.0, value=default_distance, step=0.1)
        approach_angle = st.slider("接近角度 (°)", min_value=0.0, max_value=90.0, value=default_approach, step=5.0,
                                 help="0°: π面に垂直（典型的なC-H π相互作用）、90°: π面に平行")
    
    with col2:
        rotation_angle = st.slider("回転角度 (°)", min_value=0.0, max_value=360.0, value=default_rotation, step=15.0,
                                 help="π系分子の周りでの回転")
        
        # 配置の視覚的説明
        st.markdown("""
        **配置の説明:**
        - **重心間距離**: π系分子の重心からC-H分子までの距離
        - **接近角度**: 0°でπ面に垂直（典型的なC-H π）、90°でπ面に平行
        - **回転角度**: π系分子の周りでの回転位置
        """)
        
    # 推奨設定の表示
    st.markdown("""
    **推奨設定:**
    - **メタン-ベンゼン**: 距離 3.5 Å, 角度 0°（垂直）
    - **エタン-ベンゼン**: 距離 3.2 Å, 角度 0°（垂直）
    - **T字型配置**: 距離 3.0 Å, 角度 90°（平行）
    """)
        
elif placement_method == "特定原子間距離":
    st.info("💡 特定の原子間距離を制御して分子を配置します")
    
    col1, col2 = st.columns(2)
    with col1:
        atom_idx1 = st.number_input("分子Aの原子インデックス", min_value=0, value=0, step=1)
        atom_idx2 = st.number_input("分子Bの原子インデックス", min_value=0, value=0, step=1)
    
    with col2:
        target_distance = st.slider("原子間距離 (Å)", min_value=1.0, max_value=10.0, value=3.0, step=0.1)
        
        st.markdown("""
        **使用方法:**
        1. 分子構造生成後に原子座標を確認
        2. 制御したい原子のインデックスを入力
        3. 目標距離を設定
        """)

# 分子構造の生成と表示
st.header("分子構造の生成")

# 分子構造の生成と表示
st.header("分子構造の生成")

# 分子構造生成の状態管理
if 'molecules_generated' not in st.session_state:
    st.session_state.molecules_generated = False

# 分子構造生成ボタン
if st.button("分子構造を生成", type="secondary") or st.session_state.molecules_generated:
    try:
        with st.spinner("分子構造を生成中..."):
            # 各分子ごとに配座探索・最適化
            st.write("🔍 分子Aの処理を開始...")
            handler1 = MoleculeHandler(smiles_input1, input_type="smiles")
            st.write(f"分子A: {handler1.mol.GetNumAtoms()} 原子")
            
            if handler1.mol.GetNumAtoms() > 1:  # 単原子でない場合のみ配座生成
                try:
                    st.write(f"分子A: {num_conformers1} 配座を{force_field1}で生成中...")
                    conformers = handler1.generate_conformers(num_conformers=num_conformers1, forcefield=force_field1)
                    if conformers:  # 配座生成が成功した場合
                        try:
                            handler1.keep_lowest_energy_conformer()
                            st.success(f"分子A: {len(conformers)} 配座を生成し、最低エネルギー配座を選択しました")
                        except RuntimeError as e:
                            if "No conformer energies found" in str(e):
                                st.warning("分子A: エネルギー情報が見つかりません。最初の配座を使用します")
                                # 最初の配座を使用
                                if handler1.mol.GetNumConformers() > 0:
                                    try:
                                        first_conf = handler1.mol.GetConformer(0)
                                        handler1.mol.RemoveAllConformers()
                                        handler1.mol.AddConformer(first_conf, assignId=True)
                                    except:
                                        # 配座の取得に失敗した場合は3D構造を再生成
                                        st.warning("分子A: 配座情報の修復中...")
                                        from rdkit.Chem import AllChem
                                        AllChem.EmbedMolecule(handler1.mol)
                                        AllChem.UFFOptimizeMolecule(handler1.mol)
                            else:
                                raise e
                    else:
                        st.warning("分子A: 配座生成に失敗しました。初期構造を使用します")
                except Exception as e:
                    st.warning(f"分子A: 配座生成でエラーが発生しました ({e})。初期構造を使用します")
                    # 配座が全くない場合は3D構造を生成
                    if handler1.mol.GetNumConformers() == 0:
                        st.warning("分子A: 3D構造を再生成中...")
                        from rdkit.Chem import AllChem
                        AllChem.EmbedMolecule(handler1.mol)
                        AllChem.UFFOptimizeMolecule(handler1.mol)
            else:
                st.info("分子A: 単原子分子のため配座生成をスキップします")
            mol1 = handler1.mol

            st.write("🔍 分子Bの処理を開始...")
            handler2 = MoleculeHandler(smiles_input2, input_type="smiles")
            st.write(f"分子B: {handler2.mol.GetNumAtoms()} 原子")
            
            if handler2.mol.GetNumAtoms() > 1:  # 単原子でない場合のみ配座生成
                try:
                    st.write(f"分子B: {num_conformers2} 配座を{force_field2}で生成中...")
                    conformers = handler2.generate_conformers(num_conformers=num_conformers2, forcefield=force_field2)
                    if conformers:  # 配座生成が成功した場合
                        try:
                            handler2.keep_lowest_energy_conformer()
                            st.success(f"分子B: {len(conformers)} 配座を生成し、最低エネルギー配座を選択しました")
                        except RuntimeError as e:
                            if "No conformer energies found" in str(e):
                                st.warning("分子B: エネルギー情報が見つかりません。最初の配座を使用します")
                                # 最初の配座を使用
                                if handler2.mol.GetNumConformers() > 0:
                                    try:
                                        first_conf = handler2.mol.GetConformer(0)
                                        handler2.mol.RemoveAllConformers()
                                        handler2.mol.AddConformer(first_conf, assignId=True)
                                    except:
                                        # 配座の取得に失敗した場合は3D構造を再生成
                                        st.warning("分子B: 配座情報の修復中...")
                                        from rdkit.Chem import AllChem
                                        AllChem.EmbedMolecule(handler2.mol)
                                        AllChem.UFFOptimizeMolecule(handler2.mol)
                            else:
                                raise e
                    else:
                        st.warning("分子B: 配座生成に失敗しました。初期構造を使用します")
                except Exception as e:
                    st.warning(f"分子B: 配座生成でエラーが発生しました ({e})。初期構造を使用します")
                    # 配座が全くない場合は3D構造を生成
                    if handler2.mol.GetNumConformers() == 0:
                        st.warning("分子B: 3D構造を再生成中...")
                        from rdkit.Chem import AllChem
                        AllChem.EmbedMolecule(handler2.mol)
                        AllChem.UFFOptimizeMolecule(handler2.mol)
            else:
                st.info("分子B: 単原子分子のため配座生成をスキップします")
            mol2 = handler2.mol

            # 最近接原子間距離でmol2を配置
            st.write("🔍 分子を配置中...")
            
            # 分子の配座状態を確認
            st.write(f"分子A配座数: {mol1.GetNumConformers()}, 分子B配座数: {mol2.GetNumConformers()}")
            
            # 配座が存在しない場合は3D構造を生成
            if mol1.GetNumConformers() == 0:
                st.warning("分子A: 配座が見つかりません。3D構造を生成中...")
                from rdkit.Chem import AllChem
                AllChem.EmbedMolecule(mol1)
                AllChem.UFFOptimizeMolecule(mol1)
                
            if mol2.GetNumConformers() == 0:
                st.warning("分子B: 配座が見つかりません。3D構造を生成中...")
                from rdkit.Chem import AllChem
                AllChem.EmbedMolecule(mol2)
                AllChem.UFFOptimizeMolecule(mol2)
            
            # 配置方法に応じて分子を配置
            try:
                if placement_method == "最近接原子間距離":
                    mol2_placed = MoleculeHandler.place_mol_by_closest_distance(mol1, mol2, target_distance=target_distance)
                    st.success(f"分子Bを分子Aから{target_distance} Å離れた位置に配置しました")
                    
                elif placement_method == "C-H π相互作用":
                    # 変数が定義されていない場合のデフォルト値
                    if 'approach_angle' not in locals():
                        approach_angle = 0.0
                    if 'rotation_angle' not in locals():
                        rotation_angle = 0.0
                        
                    mol2_placed = MoleculeHandler.place_mol_for_ch_pi_interaction(
                        mol1, mol2, 
                        target_distance=target_distance,
                        approach_angle=approach_angle,
                        rotation_angle=rotation_angle
                    )
                    st.success(f"C-H π相互作用用の配置を完了しました")
                    st.info(f"重心間距離: {target_distance} Å, 接近角度: {approach_angle}°, 回転角度: {rotation_angle}°")
                    
                elif placement_method == "特定原子間距離":
                    # 変数が定義されていない場合のデフォルト値
                    if 'atom_idx1' not in locals():
                        atom_idx1 = 0
                    if 'atom_idx2' not in locals():
                        atom_idx2 = 0
                        
                    # 原子インデックスの範囲チェック
                    if atom_idx1 >= mol1.GetNumAtoms():
                        st.error(f"分子Aの原子インデックス {atom_idx1} が範囲外です (0-{mol1.GetNumAtoms()-1})")
                        mol2_placed = mol2
                    elif atom_idx2 >= mol2.GetNumAtoms():
                        st.error(f"分子Bの原子インデックス {atom_idx2} が範囲外です (0-{mol2.GetNumAtoms()-1})")
                        mol2_placed = mol2
                    else:
                        mol2_placed = MoleculeHandler.place_mol_by_specific_atoms(
                            mol1, mol2, 
                            atom_idx1=atom_idx1, 
                            atom_idx2=atom_idx2,
                            target_distance=target_distance
                        )
                        st.success(f"原子 {atom_idx1} と原子 {atom_idx2} 間の距離を{target_distance} Åに設定しました")
                
            except Exception as e:
                st.warning(f"分子の配置でエラーが発生しました ({e})。元の位置を使用します")
                mol2_placed = mol2
            
            from rdkit import Chem
            combo = Chem.CombineMols(mol1, mol2_placed)

            # MoleculeHandlerの作成前に配座の状態を確認
            st.write("🔍 PySCF入力形式に変換中...")
            
            try:
                handler = MoleculeHandler(combo, input_type="rdkit")
                handler_1 = MoleculeHandler(mol1, input_type="rdkit")
                handler_2 = MoleculeHandler(mol2_placed, input_type="rdkit")
                
                # PySCF形式の原子座標を取得（改行区切り → セミコロン区切りに変換）
                atom_coords_A = handler_1.to_pyscf_input().replace('\n', '; ')
                atom_coords_B = handler_2.to_pyscf_input().replace('\n', '; ')
                atom_coords_AB = handler.to_pyscf_input().replace('\n', '; ')
                
                st.success("PySCF入力形式への変換が完了しました")
                
            except Exception as e:
                st.error(f"PySCF入力形式への変換でエラーが発生しました: {e}")
                raise e
            
            # セッション状態に保存
            st.session_state.molecules_generated = True
            st.session_state.atom_coords_A = atom_coords_A
            st.session_state.atom_coords_B = atom_coords_B  
            st.session_state.atom_coords_AB = atom_coords_AB
            st.session_state.handler = handler
            st.session_state.current_placement_method = placement_method
            st.session_state.mol1 = mol1
            st.session_state.mol2_placed = mol2_placed
            
        st.success("分子構造の生成が完了しました")
        
        # 配座探索設定の表示
        st.subheader("配座探索設定")
        st.info("💡 配座生成時に選択した分子力場による構造最適化が自動的に実行されます")
        
        conf_col1, conf_col2 = st.columns(2)
        
        with conf_col1:
            st.info(f"**分子A**: {force_field1} force field, {num_conformers1} conformers")
            
        with conf_col2:
            st.info(f"**分子B**: {force_field2} force field, {num_conformers2} conformers")
        
        # 分子構造を3つのカラムで表示
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.subheader("分子A")
            with st.expander("原子座標 (PySCF形式)"):
                st.code(st.session_state.atom_coords_A.replace('; ', '\n'))
            
            # 原子インデックス情報を表示
            current_placement = st.session_state.get('current_placement_method', placement_method)
            if current_placement == "特定原子間距離":
                with st.expander("原子インデックス情報"):
                    atom_info_A = []
                    for i in range(mol1.GetNumAtoms()):
                        atom = mol1.GetAtomWithIdx(i)
                        atom_info_A.append(f"{i}: {atom.GetSymbol()}")
                    st.text("\n".join(atom_info_A))
        
        with col2:
            st.subheader("分子B")
            with st.expander("原子座標 (PySCF形式)"):
                st.code(st.session_state.atom_coords_B.replace('; ', '\n'))
            
            # 原子インデックス情報を表示
            current_placement = st.session_state.get('current_placement_method', placement_method)
            if current_placement == "特定原子間距離":
                with st.expander("原子インデックス情報"):
                    atom_info_B = []
                    for i in range(mol2_placed.GetNumAtoms()):
                        atom = mol2_placed.GetAtomWithIdx(i)
                        atom_info_B.append(f"{i}: {atom.GetSymbol()}")
                    st.text("\n".join(atom_info_B))
        
        with col3:
            st.subheader("複合体AB")
            with st.expander("原子座標 (PySCF形式)"):
                st.code(st.session_state.atom_coords_AB.replace('; ', '\n'))
            
            # 配置情報を表示
            with st.expander("配置情報"):
                current_placement = st.session_state.get('current_placement_method', placement_method)
                if current_placement == "最近接原子間距離":
                    st.text(f"最近接原子間距離: {target_distance} Å")
                elif current_placement == "C-H π相互作用":
                    st.text(f"重心間距離: {target_distance} Å")
                    st.text(f"接近角度: {approach_angle}°")
                    st.text(f"回転角度: {rotation_angle}°")
                elif current_placement == "特定原子間距離":
                    st.text(f"原子 {atom_idx1} - 原子 {atom_idx2}")
                    st.text(f"距離: {target_distance} Å")

    except Exception as e:
        st.error(f"分子の初期化に失敗しました: {e}")
        st.session_state.molecules_generated = False
        st.stop()
else:
    st.info("分子構造を生成するには、上のボタンをクリックしてください。")

# 3D構造表示
st.header("分子の3D構造")

if st.session_state.get('molecules_generated', False):
    try:
        mol_block = st.session_state.handler.generate_3d_molblock()
        viewer = py3Dmol.view(width=600, height=400)
        viewer.addModel(mol_block, "mol")
        viewer.setStyle({"stick": {}})
        viewer.zoomTo()
        stmol.showmol(viewer, height=400)
    except Exception as e:
        st.warning(f"3D構造の表示ができませんでした: {e}")
else:
    st.info("3D構造を表示するには、まず分子構造を生成してください。")

# エネルギー分解解析の実行
st.header("エネルギー分解解析")

# サイドバーで計算設定
st.header("計算設定")

# 理論レベルと基底関数の選択
theory = st.selectbox("Theory", theory_options, key="theory_selectbox")
basis_set = st.selectbox("Basis Set", basis_set_options, key="basis_set_selectbox")

charge = st.number_input("電荷", value=0, step=1)
spin = st.number_input("スピン多重度 (2S)", value=0, step=1)

# 並列計算の有無を選択
st.header("Parallel Calculation Options")
import multiprocessing

logical_cores = multiprocessing.cpu_count()
try:
    physical_cores = multiprocessing.cpu_count(logical=False)
except TypeError:
    physical_cores = logical_cores // 2  # Fallback

st.write(f"使用しているパソコンの論理コア数: {logical_cores}")
st.write(f"使用しているパソコンの物理コア数: {physical_cores}")

# 物理コアが3以上なら並列計算を推奨（分子A、分子B、複合体ABの3つの計算のため）
recommend_parallel = physical_cores >= 3

parallel_option = st.checkbox(
    "並列計算を有効にする（推奨）" if recommend_parallel else "並列計算を有効にする",
    value=recommend_parallel,
    key="parallel_option_checkbox"
)

if recommend_parallel:
    st.info("物理コア数が3以上のため、並列計算が推奨されます。（分子A、分子B、複合体ABの3つの計算を並列実行）")
elif physical_cores >= 2:
    st.info("物理コア数が2以上です。並列計算により計算時間を短縮できる可能性があります。")
else:
    st.warning("物理コア数が少ないため、並列計算は非推奨です。")

# 計算方法と参考文献の表示
with st.expander("計算方法と参考文献を表示", expanded=False):
    st.markdown("### 🧪 Method for Energy Decomposition Analysis")
    st.markdown(
        "**Computational Details**  \n"
        "Molecular structures were processed using RDKit [1] for initial 3D coordinate generation and conformational search.  \n"
        f"Conformational search was performed using the {force_field1} force field for molecule A and {force_field2} force field for molecule B.  \n"
        f"Single-point energy calculations were performed at the **{theory}/{basis_set}** level using PySCF [2].  \n"
        "The Energy Decomposition Analysis (EDA) provides detailed insights into intermolecular interactions by decomposing the total electronic energy into:  \n"
        "- **Nuclear repulsion energy (E_nuc)**: Electrostatic repulsion between nuclei  \n"
        "- **Core Hamiltonian energy (E_core)**: Electron-nuclear attraction and kinetic energy  \n"
        "- **Coulomb interaction energy (E_J)**: Classical electron-electron repulsion  \n"
        "- **Exchange energy (E_K)**: Quantum mechanical exchange interaction  \n"
        "- **Total electronic energy (E_elec)**: Sum of all electronic contributions  \n"
        "The interaction energy is calculated as: **ΔE_int = E_AB - (E_A + E_B)**  \n"
        "where E_AB is the energy of the complex and E_A, E_B are the energies of isolated fragments.  \n"
        "Each energy component is decomposed to understand the physical origins of intermolecular interactions.  \n"
        "Energy values are provided in Hartree (Ha) and converted to kcal/mol for comparative analysis (1 Ha = 627.509 kcal/mol).  \n"
        "This decomposition analysis is essential for understanding non-covalent interactions, hydrogen bonding, and van der Waals forces.  \n"
        "For C-H π interactions, the tool provides specialized molecular placement options:  \n"
        "- **Centroid distance control**: Distance between the π-system centroid and the C-H molecule  \n"
        "- **Approach angle**: 0° for perpendicular approach (typical C-H π geometry), 90° for parallel approach  \n"
        "- **Rotation angle**: Rotation around the π-system to explore different interaction orientations  \n"
        "These geometric parameters are crucial for accurately modeling C-H π interactions in molecular complexes."
    )
    st.markdown("---")
    st.markdown(
        "**References**  \n"
        "[1] Landrum, G. RDKit: Open-source cheminformatics. [https://www.rdkit.org](https://www.rdkit.org)  \n"
        "[2] Sun, Q. *et al.* PySCF: The Python-based Simulations of Chemistry Framework. **WIREs Comput Mol Sci** *2018*, **8**, e1340. DOI: [10.1002/wcms.1340](https://doi.org/10.1002/wcms.1340)  \n"
        "[3] PocLab streamlit-pyscf: Quantum chemistry web interface. [https://github.com/poclab-web/streamlit-pyscf](https://github.com/poclab-web/streamlit-pyscf)  \n"
        "[4] Szabo, A.; Ostlund, N. S. *Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory*; Dover Publications: New York, 1996.  \n"
        "[5] Helgaker, T.; Jørgensen, P.; Olsen, J. *Molecular Electronic-Structure Theory*; Wiley: Chichester, 2000.  \n"
        "[6] Kitaura, K.; Morokuma, K. A new energy decomposition scheme for molecular interactions within the Hartree-Fock approximation. **Int. J. Quantum Chem.** *1976*, **10**, 325-340. DOI: [10.1002/qua.560100211](https://doi.org/10.1002/qua.560100211)  \n"
        "[7] Nishio, M. *et al.* CH/π hydrogen bonds in crystals. **CrystEngComm** *2004*, **6**, 130-158. DOI: [10.1039/B313104A](https://doi.org/10.1039/B313104A)  \n"
        "[8] Takahashi, O. *et al.* Relevance of weak hydrogen bonds in the conformation of biological molecules and in the stabilization of supramolecular structures. **Chem. Rev.** *2010*, **110**, 6049-6076. DOI: [10.1021/cr100072x](https://doi.org/10.1021/cr100072x)"
    )

if st.button("計算実行", type="primary"):
    if not st.session_state.get('molecules_generated', False):
        st.warning("まず分子構造を生成してください。")
    elif smiles_input1 and smiles_input2:
        try:
            # セッション状態から座標を取得
            atom_coords_A = st.session_state.atom_coords_A
            atom_coords_B = st.session_state.atom_coords_B
            atom_coords_AB = st.session_state.atom_coords_AB
            
            # 計算パラメータの準備
            compound_names = [f"Mol_A_{smiles_input1}", f"Mol_B_{smiles_input2}", f"Complex_{smiles_input1}_{smiles_input2}"]
            smiles_list = [smiles_input1, smiles_input2, f"{smiles_input1}.{smiles_input2}"]
            atom_inputs = [atom_coords_A, atom_coords_B, atom_coords_AB]
            
            # プログレスバーと状態表示
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            status_text.text("エネルギー分解解析を開始しています...")
            progress_bar.progress(10)
            
            # エネルギー分解解析の実行
            with st.spinner("量子化学計算を実行中..."):
                eda_results = get_hf_energy_decomposition(
                    compound_names=compound_names,
                    smiles_list=smiles_list,
                    atom_inputs=atom_inputs,
                    basis=basis_set,
                    theory=theory,
                    charge=charge,
                    spin=spin,
                    parallel_option=parallel_option
                )
            
            progress_bar.progress(90)
            status_text.text("結果を処理中...")
            
            # 結果の表示
            st.success("計算が完了しました！")
            progress_bar.progress(100)
            status_text.text("完了")
            
            # 計算情報の表示
            st.subheader("計算情報")
            calc_info = eda_results.get("calculation_info", {})
            info_col1, info_col2, info_col3 = st.columns(3)
            
            with info_col1:
                st.metric("理論レベル", calc_info.get('theory', 'Unknown'))
                st.metric("基底関数", calc_info.get('basis', 'Unknown'))
                
            with info_col2:
                st.metric("電荷", calc_info.get('charge', 'Unknown'))
                st.metric("スピン", calc_info.get('spin', 'Unknown'))
                
            with info_col3:
                st.metric("並列処理", "有効" if parallel_option else "無効")
                st.metric("計算ID", calc_info.get('timestamp', 'Unknown')[:10] + "...")
            
            # エラーチェック
            energy_a = eda_results["molecule_A"]
            energy_b = eda_results["molecule_B"]
            energy_ab = eda_results["complex_AB"]
            
            if "error" in energy_a or "error" in energy_b or "error" in energy_ab:
                st.error("計算中にエラーが発生しました:")
                if "error" in energy_a:
                    st.error(f"分子A: {energy_a['error']}")
                if "error" in energy_b:
                    st.error(f"分子B: {energy_b['error']}")
                if "error" in energy_ab:
                    st.error(f"複合体AB: {energy_ab['error']}")
            else:
                # 収束チェック
                converged_a = energy_a.get("converged", False)
                converged_b = energy_b.get("converged", False)
                converged_ab = energy_ab.get("converged", False)
                
                if not (converged_a and converged_b and converged_ab):
                    st.warning("一部の計算が収束していません")
                    conv_col1, conv_col2, conv_col3 = st.columns(3)
                    with conv_col1:
                        st.metric("分子A収束", "✓" if converged_a else "✗")
                    with conv_col2:
                        st.metric("分子B収束", "✓" if converged_b else "✗")
                    with conv_col3:
                        st.metric("複合体AB収束", "✓" if converged_ab else "✗")
                
                # エネルギー分解結果の表示
                if "energy_decomposition" in eda_results:
                    decomp = eda_results["energy_decomposition"]
                    
                    st.subheader("個別エネルギー")
                    energy_keys = ["energy", "E_nuc", "E_core", "E_J", "E_K", "E_elec"]
                    
                    # テーブル形式で表示
                    import pandas as pd
                    
                    individual_data = []
                    for key in energy_keys:
                        if f"{key}_A" in decomp and f"{key}_B" in decomp and f"{key}_AB" in decomp:
                            display_key = "E_total" if key == "energy" else key
                            individual_data.append({
                                "エネルギー成分": display_key,
                                "分子A (Ha)": f"{decomp[f'{key}_A']:+.6f}",
                                "分子B (Ha)": f"{decomp[f'{key}_B']:+.6f}",
                                "複合体AB (Ha)": f"{decomp[f'{key}_AB']:+.6f}"
                            })
                    
                    if individual_data:
                        df_individual = pd.DataFrame(individual_data)
                        st.dataframe(df_individual, use_container_width=True)
                    
                    st.subheader("相互作用エネルギー")
                    interaction_data = []
                    for key in energy_keys:
                        delta_key = f"Δ{key}"
                        kcal_key = f"Δ{key}_kcal_mol"
                        if delta_key in decomp and kcal_key in decomp:
                            display_key = "ΔE_total" if key == "energy" else f"Δ{key}"
                            if isinstance(decomp[delta_key], (int, float, np.number)):
                                interaction_data.append({
                                    "エネルギー成分": display_key,
                                    "相互作用エネルギー (Ha)": f"{decomp[delta_key]:+.6f}",
                                    "相互作用エネルギー (kcal/mol)": f"{decomp[kcal_key]:+.2f}"
                                })
                            else:
                                interaction_data.append({
                                    "エネルギー成分": display_key,
                                    "相互作用エネルギー (Ha)": str(decomp[delta_key]),
                                    "相互作用エネルギー (kcal/mol)": "N/A"
                                })
                    
                    if interaction_data:
                        df_interaction = pd.DataFrame(interaction_data)
                        st.dataframe(df_interaction, use_container_width=True)
                        
                        # 重要な結果をハイライト
                        if len(interaction_data) > 0:
                            total_interaction = interaction_data[0]  # E_totalの相互作用エネルギー
                            st.info(f"**全相互作用エネルギー**: {total_interaction['相互作用エネルギー (Ha)']} Ha = {total_interaction['相互作用エネルギー (kcal/mol)']} kcal/mol")
                    
                    # 結果のダウンロード
                    st.subheader("結果のダウンロード")
                    
                    # JSONファイルとして結果をダウンロード
                    import json
                    result_json = json.dumps(eda_results, indent=2, default=str)
                    st.download_button(
                        label="結果をJSONでダウンロード",
                        data=result_json,
                        file_name=f"eda_results_{calc_info.get('timestamp', 'unknown')}.json",
                        mime="application/json"
                    )
                    
                    # CSVファイルとして相互作用エネルギーをダウンロード
                    if interaction_data:
                        csv_data = pd.DataFrame(interaction_data).to_csv(index=False)
                        st.download_button(
                            label="相互作用エネルギーをCSVでダウンロード",
                            data=csv_data,
                            file_name=f"interaction_energies_{calc_info.get('timestamp', 'unknown')}.csv",
                            mime="text/csv"
                        )
                        
        except Exception as e:
            st.error(f"計算中にエラーが発生しました: {str(e)}")
            st.exception(e)
    else:
        st.warning("分子A、分子B両方のSMILESを入力してください")