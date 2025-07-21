"""
Streamlit用OpenMM設定UI

OpenMMのインストール状況確認とフォースフィールド設定のためのStreamlitコンポーネント
"""

import streamlit as st
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


def check_openmm_status():
    """OpenMMのバージョンと使用可能性をチェック"""
    status = {
        "openmm_available": False,
        "openmm_version": None,
        "openff_available": False,
        "openff_version": None,
        "error_messages": []
    }
    
    # OpenMMのチェック
    try:
        import openmm
        status["openmm_available"] = True
        status["openmm_version"] = openmm.version.version
    except ImportError as e:
        status["error_messages"].append(f"OpenMM import error: {e}")
    
    # OpenFF Toolkitのチェック
    try:
        from openff.toolkit import __version__ as openff_version
        status["openff_available"] = True
        status["openff_version"] = openff_version
    except ImportError as e:
        status["error_messages"].append(f"OpenFF Toolkit import error: {e}")
    
    return status


def test_openmm_functionality():
    """OpenMMの基本機能をテスト"""
    try:
        # 必要なライブラリをインポート
        import openmm
        from openmm import unit, app
        from openff.toolkit.topology import Molecule
        from openff.toolkit.typing.engines.smirnoff import ForceField
        
        # 簡単なテスト分子でOpenMMの動作確認
        test_smiles = "C"  # メタン
        rdkit_mol = Chem.AddHs(Chem.MolFromSmiles(test_smiles))
        AllChem.EmbedMolecule(rdkit_mol)
        AllChem.UFFOptimizeMolecule(rdkit_mol)

        off_mol = Molecule.from_rdkit(rdkit_mol, allow_undefined_stereo=True)
        topology = off_mol.to_topology()
        forcefield = ForceField("openff_unconstrained-2.0.0.offxml")
        omm_system = forcefield.create_openmm_system(topology)
        
        # インテグレータとシミュレーションの作成テスト
        integrator = openmm.LangevinIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds)
        simulation = app.Simulation(topology.to_openmm(), omm_system, integrator)
        
        return True, "OpenMM動作テスト成功"
    except Exception as e:
        return False, f"OpenMM動作テスト失敗: {str(e)}"


def test_single_point_calculation():
    """1点計算のテスト"""
    try:
        # 必要なライブラリをインポート
        import openmm
        from openmm import unit, app, Vec3
        from openff.toolkit.topology import Molecule
        from openff.toolkit.typing.engines.smirnoff import ForceField
        
        # 簡単なテスト分子（水分子）
        test_smiles = "O"
        rdkit_mol = Chem.AddHs(Chem.MolFromSmiles(test_smiles))
        AllChem.EmbedMolecule(rdkit_mol)
        AllChem.UFFOptimizeMolecule(rdkit_mol)

        # OpenFF Moleculeを作成
        off_mol = Molecule.from_rdkit(rdkit_mol, allow_undefined_stereo=True)
        topology = off_mol.to_topology()
        forcefield = ForceField("openff_unconstrained-2.0.0.offxml")
        omm_system = forcefield.create_openmm_system(topology)
        
        # 座標を取得
        conf = rdkit_mol.GetConformer()
        positions_angstrom = []
        for i in range(rdkit_mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            positions_angstrom.append([pos.x, pos.y, pos.z])
        
        # オングストロームからナノメートルに変換
        positions_nm = np.array(positions_angstrom) / 10.0
        positions = [Vec3(float(pos[0]), float(pos[1]), float(pos[2])) * unit.nanometer for pos in positions_nm]
        
        # シミュレーション作成
        integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
        simulation = app.Simulation(topology.to_openmm(), omm_system, integrator)
        simulation.context.setPositions(positions)
        
        # エネルギー計算
        state = simulation.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
        
        return True, f"1点計算テスト成功: エネルギー = {energy:.2f} kcal/mol"
    except Exception as e:
        return False, f"1点計算テスト失敗: {str(e)}"


def display_openmm_status(show_config_section=True):
    """
    OpenMMのインストール状況を表示し、必要に応じて設定UIを表示する
    
    Args:
        show_config_section (bool): 設定セクションを表示するかどうか
        
    Returns:
        dict: OpenMMの状況を示す辞書
    """
    st.subheader("OpenMM Installation Status")
    status = check_openmm_status()

    # 基本ライブラリの状況表示
    col1, col2 = st.columns(2)
    
    with col1:
        # OpenMM
        if status["openmm_available"]:
            st.success(f"✅ OpenMM v{status['openmm_version']}")
        else:
            st.error("❌ OpenMM 利用不可")
    
    with col2:
        # OpenFF
        if status["openff_available"]:
            st.success(f"✅ OpenFF v{status['openff_version']}")
        else:
            st.error("❌ OpenFF 利用不可")

    # エラーメッセージがある場合
    if status["error_messages"]:
        st.subheader("エラー詳細")
        for msg in status["error_messages"]:
            st.error(msg)

    # 全てのライブラリが利用可能な場合のテスト
    all_available = all([
        status["openmm_available"], 
        status["openff_available"]
    ])

    if all_available and show_config_section:
        # 動作テストセクション
        st.subheader("🧪 OpenMM 動作テスト")
        
        col1, col2 = st.columns(2)
        
        with col1:
            if st.button("基本動作テスト", help="OpenMMの基本的なシステム作成をテスト"):
                with st.spinner("基本動作テスト中..."):
                    test_success, test_message = test_openmm_functionality()
                    if test_success:
                        st.success(test_message)
                    else:
                        st.error(test_message)
        
        with col2:
            if st.button("1点計算テスト", help="実際の1点エネルギー計算をテスト"):
                with st.spinner("1点計算テスト中..."):
                    test_success, test_message = test_single_point_calculation()
                    if test_success:
                        st.success(test_message)
                    else:
                        st.error(test_message)

    elif not all_available and show_config_section:
        # インストール手順の表示
        st.subheader("🔧 インストール手順")
        st.markdown("""
        **必要なライブラリをインストールするには:**
        """)
        
        if not status["openmm_available"]:
            st.code("conda install -c conda-forge openmm", language="bash")
            st.markdown("- OpenMM: 分子動力学シミュレーション")
            
        if not status["openff_available"]:
            st.code("conda install -c conda-forge openff-toolkit", language="bash")
            st.markdown("- OpenFF Toolkit: オープンフォースフィールド")

        st.markdown("""
        **推奨インストール手順:**
        1. 上記のコマンドを順番に実行
        2. ページを再読み込みして状況を確認
        3. 動作テストを実行して正常性を確認
        """)

    return status


def require_openmm():
    """
    OpenMMが必要なページで使用する関数
    OpenMMが利用できない場合はページの実行を停止する
    
    Returns:
        dict: OpenMMが利用可能な場合のステータス（実際には利用不可の場合は停止する）
    """
    status = display_openmm_status()
    
    if not all([status["openmm_available"], 
                status["openff_available"]]):
        st.stop()
    
    return status


def get_available_forcefields():
    """利用可能なフォースフィールドのリストを取得"""
    forcefields = [
        "openff_unconstrained-2.0.0.offxml",
        "openff-2.0.0.offxml",
        "openff-1.3.0.offxml",
        "openff_unconstrained-1.3.0.offxml"
    ]
    return forcefields


def display_forcefield_selector(default_ff="openff_unconstrained-2.0.0.offxml"):
    """フォースフィールド選択ウィジェットを表示"""
    available_ffs = get_available_forcefields()
    
    selected_ff = st.selectbox(
        "フォースフィールド",
        available_ffs,
        index=available_ffs.index(default_ff) if default_ff in available_ffs else 0,
        help="計算に使用するフォースフィールドを選択してください"
    )
    
    return selected_ff
