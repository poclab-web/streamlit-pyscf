import concurrent.futures
from pyscf import gto, scf
import numpy as np

from logic.calculation import run_quantum_calculation


def run_ced(parallel_option, compound_name1, compound_name2, atom_input1, atom_input2, smiles1, smiles2, basis_set, theory, charge, spin, opt_theory=None, opt_basis_set=None, solvent_model=None, eps=None, symmetry=False):
    import time
    import random
    
    # 一意のタイムスタンプを生成してファイル名の重複を回避
    timestamp = str(int(time.time() * 1000))
    random_suffix = str(random.randint(1000, 9999))
    
    # 各分子に固有の識別子を追加
    compound_name1_unique = f"{compound_name1}_mol1_{timestamp}_{random_suffix}"
    compound_name2_unique = f"{compound_name2}_mol2_{timestamp}_{random_suffix}"
    
    if parallel_option:
        # 並列実行
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
            future_A = executor.submit(
                run_quantum_calculation, compound_name1_unique, smiles1, atom_input1,
                basis_set, theory, charge, spin, opt_theory, opt_basis_set, solvent_model, eps, symmetry
            )
            future_B = executor.submit(
                run_quantum_calculation, compound_name2_unique, smiles2, atom_input2,
                basis_set, theory, charge, spin, opt_theory, opt_basis_set, solvent_model, eps, symmetry
            )
            eda_A = future_A.result()
            eda_B = future_B.result()
    else:
        # 通常実行
        eda_A = run_quantum_calculation(compound_name1_unique, smiles1, atom_input1, basis_set, theory, charge, spin, opt_theory, opt_basis_set, solvent_model, eps, symmetry)
        eda_B = run_quantum_calculation(compound_name2_unique, smiles2, atom_input2, basis_set, theory, charge, spin, opt_theory, opt_basis_set, solvent_model, eps, symmetry)

    # 戻り値が辞書であることを確認
    if not isinstance(eda_A, dict) or not isinstance(eda_B, dict):
        raise TypeError(f"Expected dict, got eda_A: {type(eda_A)}, eda_B: {type(eda_B)}")
    
    # エネルギー分解の対象となる数値キーのみを抽出
    energy_keys = ['energy', 'E_nuc', 'E_core', 'E_J', 'E_K', 'E_elec']
    
    delta = {}
    
    # 個別の分子のエネルギーも結果に含める
    if 'energy' in eda_A:
        delta['energy_A'] = eda_A['energy']
    if 'energy' in eda_B:
        delta['energy_B'] = eda_B['energy']
    
    for key in energy_keys:
        if key in eda_A and key in eda_B:
            if isinstance(eda_A[key], (int, float, np.number)) and isinstance(eda_B[key], (int, float, np.number)):
                delta[key] = eda_A[key] - eda_B[key]
                # kcal/mol単位でも表示
                delta[f"{key}_kcal_mol"] = (eda_A[key] - eda_B[key]) * 627.509
            else:
                delta[key] = f"Cannot subtract: {type(eda_A[key])} - {type(eda_B[key])}"
        else:
            delta[key] = f"Key '{key}' not found in one or both results"
    
    # 収束情報も含める
    if 'converged' in eda_A and 'converged' in eda_B:
        delta['converged_A'] = eda_A['converged']
        delta['converged_B'] = eda_B['converged']

    return delta


def get_hf_energy_decomposition(compound_names, smiles_list, atom_inputs, basis="sto-3g", theory="RHF", charge=0, spin=0, solvent_model=None, eps=None, symmetry=False, parallel_option=True):
    """
    3つの分子（A、B、AB複合体）のエネルギー分解を実行する関数
    
    Parameters:
    compound_names: 化合物名のリスト [compound_A, compound_B, compound_AB]
    smiles_list: SMILESのリスト [smiles_A, smiles_B, smiles_AB] 
    atom_inputs: 原子座標のリスト [atom_A, atom_B, atom_AB]
    basis: 基底関数セット
    theory: 理論レベル（RHF, UHF, RKS, UKS等）
    charge: 電荷
    spin: スピン多重度（2S）
    solvent_model: 溶媒モデル（PCM等）
    eps: 誘電率
    symmetry: 対称性の使用可否
    parallel_option: 並列実行するかどうか
    """
    
    if len(compound_names) != 3 or len(smiles_list) != 3 or len(atom_inputs) != 3:
        raise ValueError("compound_names, smiles_list, atom_inputsはそれぞれ3つの要素が必要です")
    
    import time
    import random
    
    # 一意のタイムスタンプを生成してファイル名の重複を回避
    timestamp = str(int(time.time() * 1000))
    random_suffix = str(random.randint(1000, 9999))
    
    # 各分子に固有の識別子を追加
    compound_names_unique = [
        f"{compound_names[0]}_eda_{timestamp}_{random_suffix}_0",
        f"{compound_names[1]}_eda_{timestamp}_{random_suffix}_1", 
        f"{compound_names[2]}_eda_{timestamp}_{random_suffix}_2"
    ]
    
    results = []
    
    if parallel_option:
        # 並列実行
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
            futures = []
            for i, (compound_name, smiles, atom_input) in enumerate(zip(compound_names_unique, smiles_list, atom_inputs)):
                future = executor.submit(
                    run_quantum_calculation,
                    compound_name, smiles, atom_input,
                    basis, theory, charge, spin,
                    opt_theory=None, opt_basis_set=None,
                    solvent_model=solvent_model, eps=eps, symmetry=symmetry
                )
                futures.append(future)
            
            # 結果を順番に取得
            for i, future in enumerate(futures):
                try:
                    result = future.result()                    
                    if isinstance(result, dict) and 'energy' in result:
                        results.append(result)
                    else:
                        results.append({"error": f"Invalid result format for {compound_names[i]}", "converged": False})
                        
                except Exception as e:
                    print(f"エラー: 分子 {i+1} ({compound_names[i]}) の並列計算中にエラーが発生しました: {str(e)}")
                    results.append({
                        "error": str(e),
                        "converged": False,
                        "compound_name": compound_names[i]
                    })
    else:
        # 通常実行（順次処理）
        for i, (compound_name, smiles, atom_input) in enumerate(zip(compound_names_unique, smiles_list, atom_inputs)):
            try:
                result = run_quantum_calculation(
                    compound_name, smiles, atom_input,
                    basis, theory, charge, spin, 
                    opt_theory=None, opt_basis_set=None, 
                    solvent_model=solvent_model, eps=eps, symmetry=symmetry
                )
                
                if isinstance(result, dict) and 'energy' in result:
                    results.append(result)
                else:
                    results.append({"error": f"Invalid result format for {compound_names[i]}", "converged": False})
                    
            except Exception as e:
                results.append({
                    "error": str(e),
                    "converged": False,
                    "compound_name": compound_names[i]
                })
    
    # 結果が3つ揃っているかチェック
    if len(results) != 3:
        raise RuntimeError("3つの分子すべての計算が完了しませんでした")
    
    # エラーチェック
    for i, result in enumerate(results):
        if "error" in result:
            print(f"分子 {i+1} でエラーが発生: {result['error']}")
    
    # エネルギー分解の計算を追加
    eda_A = results[0]  # 分子A
    eda_B = results[1]  # 分子B
    eda_AB = results[2] # 複合体AB
    
    # エネルギー分解の対象となる数値キーのみを抽出
    energy_keys = ['energy', 'E_nuc', 'E_core', 'E_J', 'E_K', 'E_elec']
    
    delta = {}
    
    # 個別の分子のエネルギーも結果に含める
    if 'energy' in eda_A:
        delta['energy_A'] = eda_A['energy']
    if 'energy' in eda_B:
        delta['energy_B'] = eda_B['energy']
    if 'energy' in eda_AB:
        delta['energy_AB'] = eda_AB['energy']
    
    # 相互作用エネルギーの計算: E_AB - (E_A + E_B)
    for key in energy_keys:
        if key in eda_A and key in eda_B and key in eda_AB:
            if isinstance(eda_A[key], (int, float, np.number)) and isinstance(eda_B[key], (int, float, np.number)) and isinstance(eda_AB[key], (int, float, np.number)):
                # 相互作用エネルギー = 複合体 - (分子A + 分子B)
                interaction_energy = eda_AB[key] - (eda_A[key] + eda_B[key])
                delta[f"Δ{key}"] = interaction_energy
                # kcal/mol単位でも表示
                delta[f"Δ{key}_kcal_mol"] = interaction_energy * 627.509
                
                # 個別のエネルギー成分も保存
                delta[f"{key}_A"] = eda_A[key]
                delta[f"{key}_B"] = eda_B[key]
                delta[f"{key}_AB"] = eda_AB[key]
            else:
                delta[f"Δ{key}"] = f"Cannot calculate: incompatible types"
        else:
            delta[f"Δ{key}"] = f"Key '{key}' not found in one or more results"
    
    # 収束情報も含める
    if 'converged' in eda_A and 'converged' in eda_B and 'converged' in eda_AB:
        delta['converged_A'] = eda_A['converged']
        delta['converged_B'] = eda_B['converged']
        delta['converged_AB'] = eda_AB['converged']
    
    return {
        "molecule_A": results[0],
        "molecule_B": results[1], 
        "complex_AB": results[2],
        "energy_decomposition": delta,  # エネルギー分解結果を追加
        "calculation_info": {
            "theory": theory,
            "basis": basis,
            "charge": charge,
            "spin": spin,
            "solvent_model": solvent_model,
            "eps": eps,
            "symmetry": symmetry,
            "timestamp": timestamp
        }
    }


if __name__ == "__main__":
    # テスト用の分子データ
    compound_names = ["HF", "NH3", "HF_NH3_complex"]
    smiles_list = ["F", "N", "F.N"]  # 簡単なSMILES（実際の複合体では正確でない場合があります）

    # 原子座標
    hf_coords = "H 0 0 0; F 0 0 0.9"
    nh3_coords = "N 0 0 0; H 0 0.9 0; H 0.9 0 0; H 0 0 -0.9"
    complex_coords = """
    N 0.0000 0.0000 0.0000
    H 0.9000 0.0000 0.0000
    H 0.0000 0.9000 0.0000
    H 0.0000 0.0000 -0.9000
    F 0.0000 0.0000 2.5000
    H 0.0000 0.0000 3.4000
    """

    atom_inputs = [hf_coords, nh3_coords, complex_coords]

    # 計算の実行（並列処理を有効化）
    print("エネルギー分解解析を開始します...")
    eda_results = get_hf_energy_decomposition(
        compound_names=compound_names,
        smiles_list=smiles_list, 
        atom_inputs=atom_inputs,
        basis="sto-3g",
        theory="HF",
        charge=0,
        spin=0,
        parallel_option=True  # 並列処理を有効化
    )

    # 差分計算関数も新しいデータ構造に対応
    def print_energy_diff(eda_results):
        print("===== 相互作用エネルギーの分解 =====")
        
        # 計算情報の表示
        calc_info = eda_results.get("calculation_info", {})
        print(f"理論レベル: {calc_info.get('theory', 'Unknown')}")
        print(f"基底関数: {calc_info.get('basis', 'Unknown')}")
        print(f"計算ID: {calc_info.get('timestamp', 'Unknown')}")
        print()
        
        energy_a = eda_results["molecule_A"]
        energy_b = eda_results["molecule_B"] 
        energy_ab = eda_results["complex_AB"]
        
        # エラーチェック
        if "error" in energy_a or "error" in energy_b or "error" in energy_ab:
            print("計算中にエラーが発生しました:")
            if "error" in energy_a:
                print(f"分子A: {energy_a['error']}")
            if "error" in energy_b:
                print(f"分子B: {energy_b['error']}")
            if "error" in energy_ab:
                print(f"複合体AB: {energy_ab['error']}")
            return
        
        # 収束チェック
        converged_a = energy_a.get("converged", False)
        converged_b = energy_b.get("converged", False) 
        converged_ab = energy_ab.get("converged", False)
        
        if not (converged_a and converged_b and converged_ab):
            print("警告: 一部の計算が収束していません")
            print(f"分子A収束: {converged_a}")
            print(f"分子B収束: {converged_b}")
            print(f"複合体AB収束: {converged_ab}")
        
        # エネルギー分解結果の表示
        if "energy_decomposition" in eda_results:
            decomp = eda_results["energy_decomposition"]
            print("\n=== 個別エネルギー ===")
            energy_keys = ["energy", "E_nuc", "E_core", "E_J", "E_K", "E_elec"]
            
            for key in energy_keys:
                if f"{key}_A" in decomp and f"{key}_B" in decomp and f"{key}_AB" in decomp:
                    display_key = "E_total" if key == "energy" else key
                    print(f"{display_key:10s} : A={decomp[f'{key}_A']:+.6f} Ha, B={decomp[f'{key}_B']:+.6f} Ha, AB={decomp[f'{key}_AB']:+.6f} Ha")
            
            print("\n=== 相互作用エネルギー ===")
            for key in energy_keys:
                delta_key = f"Δ{key}"
                kcal_key = f"Δ{key}_kcal_mol"
                if delta_key in decomp and kcal_key in decomp:
                    display_key = "ΔE_total" if key == "energy" else f"Δ{key}"
                    if isinstance(decomp[delta_key], (int, float, np.number)):
                        print(f"{display_key:10s} : {decomp[delta_key]:+.6f} Ha = {decomp[kcal_key]:+.2f} kcal/mol")
                    else:
                        print(f"{display_key:10s} : {decomp[delta_key]}")
        else:
            # フォールバック：従来の方法で計算
            print("\n=== 相互作用エネルギー（従来計算） ===")
            energy_keys = ["energy", "E_nuc", "E_core", "E_J", "E_K", "E_elec"]
            
            for key in energy_keys:
                if key in energy_a and key in energy_b and key in energy_ab:
                    delta = energy_ab[key] - (energy_a[key] + energy_b[key])
                    kcal = delta * 627.509
                    display_key = "E_total" if key == "energy" else key
                    print(f"{display_key:10s} : ΔE = {delta:+.6f} Ha = {kcal:+.2f} kcal/mol")
                else:
                    print(f"{key:10s} : データが不足しています")

    # 結果の表示
    print_energy_diff(eda_results)
