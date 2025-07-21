"""
Streamlit用PySCF設定UI

PySCFのインストール状況確認と基本機能テストのためのStreamlitコンポーネント
"""

import streamlit as st
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


def check_pyscf_status():
    """PySCFのバージョンと使用可能性をチェック"""
    status = {
        "pyscf_available": False,
        "pyscf_version": None,
        "pyscf_extensions": {},
        "available_modules": {},
        "error_messages": []
    }
    
    # PySCFコアのチェック
    try:
        import pyscf
        status["pyscf_available"] = True
        status["pyscf_version"] = pyscf.__version__
    except ImportError as e:
        status["error_messages"].append(f"PySCF import error: {e}")
    
    # PySCF拡張モジュールのチェック
    extensions_to_check = {
        "properties": "pyscf.properties",
        "qsdopt": "pyscf.qsdopt", 
        "dftd3": "dftd3",
        "ase": "ase"
    }
    
    for ext_name, module_name in extensions_to_check.items():
        try:
            module = __import__(module_name)
            version = getattr(module, '__version__', 'unknown')
            status["pyscf_extensions"][ext_name] = {
                "available": True,
                "version": version
            }
        except ImportError:
            status["pyscf_extensions"][ext_name] = {
                "available": False,
                "version": None
            }
    
    # 主要なPySCFモジュールのチェック
    if status["pyscf_available"]:
        modules_to_check = [
            "pyscf.scf",
            "pyscf.dft",
            "pyscf.mp",
            "pyscf.cc",
            "pyscf.mcscf",
            "pyscf.tddft",
            "pyscf.solvent",
            "pyscf.prop.polarizability",
            "pyscf.hessian",
            "pyscf.grad",
            "pyscf.geomopt",
            "pyscf.symm"
        ]
        
        for module_name in modules_to_check:
            try:
                __import__(module_name)
                status["available_modules"][module_name] = True
            except ImportError:
                status["available_modules"][module_name] = False
    
    return status


def test_pyscf_basic_functionality():
    """PySCFの基本機能をテスト"""
    try:
        import pyscf
        from pyscf import gto, scf
        
        # 簡単なテスト分子（水分子）
        mol = gto.Mole()
        mol.atom = '''
        O 0.0 0.0 0.0
        H 0.0 0.0 0.96
        H 0.927 0.0 -0.24
        '''
        mol.basis = 'sto-3g'
        mol.build()
        
        # HF計算のテスト
        mf = scf.RHF(mol)
        mf.verbose = 0  # ログを抑制
        energy = mf.kernel()
        
        return True, f"基本HF計算テスト成功: エネルギー = {energy:.6f} Hartree"
    except Exception as e:
        return False, f"基本機能テスト失敗: {str(e)}"


def test_pyscf_dft_functionality():
    """PySCFのDFT機能をテスト"""
    try:
        import pyscf
        from pyscf import gto, dft
        
        # メタン分子でのDFTテスト
        mol = gto.Mole()
        mol.atom = '''
        C 0.0 0.0 0.0
        H 0.0 0.0 1.1
        H 1.038 0.0 -0.367
        H -0.519 0.899 -0.367
        H -0.519 -0.899 -0.367
        '''
        mol.basis = '6-31g'
        mol.build()
        
        # B3LYP計算
        mf = dft.RKS(mol)
        mf.xc = 'b3lyp'
        mf.verbose = 0
        energy = mf.kernel()
        
        return True, f"DFT(B3LYP)計算テスト成功: エネルギー = {energy:.6f} Hartree"
    except Exception as e:
        return False, f"DFT機能テスト失敗: {str(e)}"


def test_pyscf_geometry_optimization():
    """PySCF構造最適化機能をテスト"""
    try:
        import pyscf
        from pyscf import gto, scf
        from pyscf.geomopt.geometric_solver import optimize
        
        # 水分子の構造最適化テスト
        mol = gto.Mole()
        mol.atom = '''
        O 0.0 0.0 0.0
        H 0.0 0.0 1.0
        H 0.9 0.0 -0.2
        '''
        mol.basis = 'sto-3g'
        mol.build()
        
        mf = scf.RHF(mol)
        mf.verbose = 0
        
        # 構造最適化
        mol_opt = optimize(mf, maxsteps=5)  # ステップ数を制限してテスト
        
        return True, "構造最適化テスト成功"
    except Exception as e:
        return False, f"構造最適化テスト失敗: {str(e)}"


def test_pyscf_solvent_calculation():
    """PySCF溶媒効果計算をテスト"""
    try:
        import pyscf
        from pyscf import gto, scf, solvent
        
        # 水分子の溶媒効果計算
        mol = gto.Mole()
        mol.atom = '''
        O 0.0 0.0 0.0
        H 0.0 0.0 0.96
        H 0.927 0.0 -0.24
        '''
        mol.basis = 'sto-3g'
        mol.build()
        
        # PCM溶媒効果を含むSCF計算
        mf = scf.RHF(mol).PCM()
        mf.with_solvent.eps = 78.39  # 水の誘電率
        mf.verbose = 0
        energy = mf.kernel()
        
        return True, f"溶媒効果計算テスト成功: エネルギー = {energy:.6f} Hartree"
    except Exception as e:
        return False, f"溶媒効果計算テスト失敗: {str(e)}"


def display_pyscf_status(show_config_section=True):
    """
    PySCFのインストール状況を表示し、必要に応じて設定UIを表示する
    
    Args:
        show_config_section (bool): 設定セクションを表示するかどうか
        
    Returns:
        dict: PySCFの状況を示す辞書
    """
    st.subheader("PySCF Installation Status")
    status = check_pyscf_status()

    # 基本ライブラリの状況表示
    col1, col2 = st.columns(2)
    
    with col1:
        # PySCFコア
        if status["pyscf_available"]:
            st.success(f"✅ PySCF v{status['pyscf_version']}")
        else:
            st.error("❌ PySCF 利用不可")
    
    with col2:
        # 拡張モジュールの概要
        available_extensions = sum(1 for ext in status["pyscf_extensions"].values() if ext["available"])
        total_extensions = len(status["pyscf_extensions"])
        st.info(f"🔧 拡張モジュール: {available_extensions}/{total_extensions}")

    # 拡張モジュールの詳細表示
    if status["pyscf_available"]:
        with st.expander("拡張モジュール詳細"):
            for ext_name, ext_info in status["pyscf_extensions"].items():
                if ext_info["available"]:
                    st.success(f"✅ {ext_name} v{ext_info['version']}")
                else:
                    st.warning(f"⚠️ {ext_name} 利用不可")
        
        # 主要モジュールの状況
        with st.expander("PySCFモジュール詳細"):
            col1, col2 = st.columns(2)
            module_items = list(status["available_modules"].items())
            mid_point = len(module_items) // 2
            
            with col1:
                for module_name, available in module_items[:mid_point]:
                    display_name = module_name.replace("pyscf.", "")
                    if available:
                        st.success(f"✅ {display_name}")
                    else:
                        st.error(f"❌ {display_name}")
            
            with col2:
                for module_name, available in module_items[mid_point:]:
                    display_name = module_name.replace("pyscf.", "")
                    if available:
                        st.success(f"✅ {display_name}")
                    else:
                        st.error(f"❌ {display_name}")

    # エラーメッセージがある場合
    if status["error_messages"]:
        st.subheader("エラー詳細")
        for msg in status["error_messages"]:
            st.error(msg)

    # PySCFが利用可能な場合のテスト
    if status["pyscf_available"] and show_config_section:
        # 動作テストセクション
        st.subheader("🧪 PySCF 動作テスト")
        
        col1, col2 = st.columns(2)
        
        with col1:
            if st.button("基本HF計算テスト", help="PySCFの基本的なHartree-Fock計算をテスト"):
                with st.spinner("基本HF計算テスト中..."):
                    test_success, test_message = test_pyscf_basic_functionality()
                    if test_success:
                        st.success(test_message)
                    else:
                        st.error(test_message)
        
            if st.button("構造最適化テスト", help="PySCFの構造最適化機能をテスト"):
                with st.spinner("構造最適化テスト中..."):
                    test_success, test_message = test_pyscf_geometry_optimization()
                    if test_success:
                        st.success(test_message)
                    else:
                        st.error(test_message)
        
        with col2:
            if st.button("DFT計算テスト", help="PySCFのDFT計算機能をテスト"):
                with st.spinner("DFT計算テスト中..."):
                    test_success, test_message = test_pyscf_dft_functionality()
                    if test_success:
                        st.success(test_message)
                    else:
                        st.error(test_message)
        
            if st.button("溶媒効果計算テスト", help="PySCFの溶媒効果計算をテスト"):
                with st.spinner("溶媒効果計算テスト中..."):
                    test_success, test_message = test_pyscf_solvent_calculation()
                    if test_success:
                        st.success(test_message)
                    else:
                        st.error(test_message)

    elif not status["pyscf_available"] and show_config_section:
        # インストール手順の表示
        st.subheader("🔧 インストール手順")
        st.markdown("""
        **PySCFをインストールするには:**
        """)
        
        st.code("pip install --prefer-binary pyscf", language="bash")
        st.markdown("- PySCF: Python量子化学計算フレームワーク")
        
        st.markdown("""
        **拡張モジュールのインストール:**
        """)
        st.code("""pip install dftd3
pip install ase
pip install git+https://github.com/pyscf/properties
pip install git+https://github.com/pyscf/qsdopt""", language="bash")

        st.markdown("""
        **推奨インストール手順:**
        1. 上記のコマンドを順番に実行
        2. ページを再読み込みして状況を確認
        3. 動作テストを実行して正常性を確認
        """)

    return status


def require_pyscf():
    """
    PySCFが必要なページで使用する関数
    PySCFが利用できない場合はページの実行を停止する
    
    Returns:
        dict: PySCFが利用可能な場合のステータス（実際には利用不可の場合は停止する）
    """
    status = display_pyscf_status()
    
    if not status["pyscf_available"]:
        st.stop()
    
    return status


def get_available_basis_sets():
    """利用可能な基底関数系のリストを取得"""
    basis_sets = [
        # 最小基底
        "sto-3g",
        # 分割価電子基底
        "3-21g",
        "6-31g",
        "6-31g(d)",
        "6-31g(d,p)",
        "6-31+g(d)",
        "6-31+g(d,p)",
        "6-311g(d,p)",
        "6-311+g(d,p)",
        # Dunning基底
        "cc-pvdz",
        "cc-pvtz",
        "cc-pvqz",
        "aug-cc-pvdz",
        "aug-cc-pvtz",
        # def2基底
        "def2-svp",
        "def2-svpd",
        "def2-tzvp",
        "def2-tzvpd",
        "def2-qzvp"
    ]
    return basis_sets


def get_available_functionals():
    """利用可能なDFT汎関数のリストを取得"""
    functionals = [
        # LDA
        "lda,vwn",
        # GGA
        "pbe",
        "blyp",
        "bp86",
        # ハイブリッド
        "b3lyp",
        "pbe0",
        "m06",
        "m06-2x",
        # 長距離補正
        "wb97x-d",
        "camb3lyp",
        # Meta-GGA
        "tpss",
        "m11"
    ]
    return functionals


def display_basis_selector(default_basis="6-31g(d,p)"):
    """基底関数系選択ウィジェットを表示"""
    available_basis = get_available_basis_sets()
    
    selected_basis = st.selectbox(
        "基底関数系",
        available_basis,
        index=available_basis.index(default_basis) if default_basis in available_basis else 0,
        help="計算に使用する基底関数系を選択してください"
    )
    
    return selected_basis


def display_functional_selector(default_functional="b3lyp"):
    """DFT汎関数選択ウィジェットを表示"""
    available_functionals = get_available_functionals()
    
    selected_functional = st.selectbox(
        "DFT汎関数",
        available_functionals,
        index=available_functionals.index(default_functional) if default_functional in available_functionals else 0,
        help="DFT計算に使用する交換相関汎関数を選択してください"
    )
    
    return selected_functional


def display_calculation_options():
    """計算オプションの設定UIを表示"""
    st.subheader("⚙️ 計算設定")
    
    col1, col2 = st.columns(2)
    
    with col1:
        charge = st.number_input(
            "電荷",
            value=0,
            step=1,
            help="分子の電荷を指定してください"
        )
        
        multiplicity = st.number_input(
            "多重度",
            value=1,
            min_value=1,
            step=1,
            help="スピン多重度を指定してください（1=一重項、2=二重項、3=三重項）"
        )
    
    with col2:
        symmetry = st.checkbox(
            "対称性を使用",
            value=True,
            help="分子対称性を利用して計算を高速化"
        )
        
        convergence = st.selectbox(
            "収束判定基準",
            ["normal", "tight", "loose"],
            index=0,
            help="SCF収束の厳密さを選択"
        )
    
    return {
        "charge": charge,
        "multiplicity": multiplicity,
        "symmetry": symmetry,
        "convergence": convergence
    }
