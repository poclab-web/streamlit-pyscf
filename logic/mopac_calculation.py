# mopacがインストールされていることを確認
import subprocess
import os
import tempfile
import shutil
from pathlib import Path
import re
from rdkit import Chem

# MOPACで選べる計算理論
theory_options = [
    "PM7",     # 最新の半経験的法。高精度で汎用的。
    "PM6",     # 以前の主力モデル。PM7より少し粗い。
    "AM1",     # 古典的なモデル（古いが軽量）。
    "MNDO",    # 最も基本的な手法（教材向き）。
]

# 精度や構造最適化の挙動を調整するキーワード
precision_options = [
    "PRECISE",    # 高精度モード（厳密な収束、やや重い）
    "GEO-OK",     # 座標異常があっても実行を続ける
    "RECALC=10",  # 10ステップごとにヘッセ行列再計算（TS探索に有効）
    "1SCF",       # 構造最適化なし、1回のSCFのみ（エネルギーだけ欲しいとき）
    "CHARGE=1",   # 分子の電荷（例: カチオン）
    "CHARGE=-1",  # 分子の電荷（例: アニオン）
    "SINGLET",    # 一重項（既定値だが明示したい場合）
    "TRIPLET",    # 三重項（励起状態やラジカルに必要）
]

# 結果の出力内容に関わるキーワード
output_options = [
    "THERMO",     # 熱化学量（Gibbsエネルギー・エンタルピーなど）
    "FORCE",      # 振動解析を実行（IRスペクトルにも必要）
    "VECTORS",    # 分子軌道ベクトル（波動関数）を出力
    "GRAPH MO",   # 分子軌道を3Dグラフィック形式で出力
    "GRADIENTS",  # 原子ごとのエネルギー勾配出力
    "MOZYME",     # 局在軌道法（大規模系に対応、タンパク質など）
    "NOINTER",    # 中間出力を抑制（軽量化、ログを短く）
    "DEBUG",      # 詳細な内部情報を出力（トラブル時用）
]

def check_mopac_installation():
    """
    MOPACのインストール状況をチェックする
    
    Returns:
        dict: インストール状況の情報
            - installed (bool): インストールされているかどうか
            - path (str): 実行ファイルのパス（見つからない場合はNone）
            - error (str): エラーメッセージ（エラーがない場合は空文字列）
    """
    result = {
        'installed': False,
        'path': None,
        'error': ''
    }
    
    # 可能なMOPAC実行ファイル名
    possible_names = ['mopac', 'MOPAC2016.exe', 'MOPAC2023.exe', 'MOPAC']
    
    # 実行ファイルを探す
    mopac_path = None
    for name in possible_names:
        path = shutil.which(name)
        if path:
            mopac_path = path
            break
    
    if not mopac_path:
        result['error'] = 'MOPAC executable not found in PATH'
        return result
    
    result['path'] = mopac_path
    result['installed'] = True
    
    return result


class MopacCalculator:
    """
    MoleculeHandlerから受け取った分子データを使ってMOPAC計算を実行するクラス
    """
    
    def __init__(self, molecule_handler, work_dir=None):
        """
        初期化
        
        Args:
            molecule_handler: MoleculeHandlerのインスタンス
            work_dir: 作業ディレクトリ（指定しない場合はInChIKeyベースのディレクトリを作成）
        """
        self.molecule_handler = molecule_handler
        
        # work_dirが指定されていない場合は、InChIKeyからディレクトリ名を作成
        if work_dir is None:
            inchikey = self._get_inchikey_from_molecule()
            work_dir = Path("data") / inchikey / "mopac_work"
        
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)
        
        # MOPACの実行可能ファイルのパスを確認
        self.mopac_executable = self._find_mopac_executable()
        if not self.mopac_executable:
            raise RuntimeError("MOPAC executable not found. Please install MOPAC.")
    
    def _get_inchikey_from_molecule(self):
        """
        分子からInChIKeyを取得する
        
        Returns:
            str: InChIKey（取得できない場合はフォールバック値）
        """
        try:
            from rdkit import Chem
            if self.molecule_handler.mol:
                inchikey = Chem.MolToInchiKey(self.molecule_handler.mol)
                # ファイルシステムで安全な文字のみ使用
                return inchikey.replace('/', '_').replace('\\', '_')
            else:
                return "unknown_molecule"
        except Exception:
            return "unknown_molecule"
    
    def _get_inchikey_from_title(self, title):
        """
        タイトル行からInChIKeyを抽出する（InChIKeyの形式に一致する場合）
        
        Args:
            title (str): タイトル文字列
            
        Returns:
            str: InChIKey（見つからない場合は分子から取得）
        """
        # InChIKeyの形式をチェック（XXX-XXX-X の形式）
        import re
        inchikey_pattern = r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$'
        
        if re.match(inchikey_pattern, title.strip()):
            return title.strip().replace('/', '_').replace('\\', '_')
        else:
            # タイトルがInChIKeyでない場合は分子から取得
            return self._get_inchikey_from_molecule()
    
    def set_work_dir_from_title(self, title):
        """
        タイトルからInChIKeyを抽出してディレクトリを設定する
        
        Args:
            title (str): MOPAC入力ファイルのタイトル行
        """
        inchikey = self._get_inchikey_from_title(title)
        self.work_dir = Path("data") / inchikey / "mopac_work"
        self.work_dir.mkdir(parents=True, exist_ok=True)
    
    def _find_mopac_executable(self):
        """MOPACの実行可能ファイルを探す"""
        mopac_status = check_mopac_installation()
        if mopac_status['installed']:
            return os.path.basename(mopac_status['path'])
        return None
    
    def create_input_file(self, 
                         theory="PM7", 
                         keywords=None, 
                         charge=0, 
                         multiplicity=1,
                         title=None,
                         filename="temp"):
        """
        MOPAC入力ファイルを作成
        
        Args:
            theory: 計算理論 (PM7, PM6, AM1, MNDO)
            keywords: 追加キーワードのリスト
            charge: 分子の電荷
            multiplicity: スピン多重度 (1=一重項, 2=二重項, 3=三重項)
            title: 計算のタイトル（Noneの場合はInChIKeyを使用）
            filename: ファイル名（拡張子なし）
            
        Returns:
            Path: 作成されたMOPファイルのパス
        """
        if not self.molecule_handler.mol:
            raise ValueError("Molecule is not initialized in MoleculeHandler.")
        
        # キーワード行を構築
        keyword_line_parts = [theory.upper()]
        
        # 電荷の設定
        if charge != 0:
            keyword_line_parts.append(f"CHARGE={charge}")
        
        # スピン多重度の設定
        if multiplicity == 1:
            keyword_line_parts.append("SINGLET")
        elif multiplicity == 2:
            keyword_line_parts.append("DOUBLET")
        elif multiplicity == 3:
            keyword_line_parts.append("TRIPLET")
        elif multiplicity > 3:
            keyword_line_parts.append(f"MS={(multiplicity-1)/2}")
        
        # 追加キーワード
        if keywords:
            keyword_line_parts.extend(keywords)
        
        keyword_line = " ".join(keyword_line_parts)
        
        # MOPAC入力ファイルを作成
        mopac_input = self.molecule_handler.to_mopac_input(
            title=title, 
            keywords=keyword_line
        )
        
        # タイトルからディレクトリを設定（必要に応じて）
        if title:
            self.set_work_dir_from_title(title)
        
        mop_file = self.work_dir / f"{filename}.mop"
        
        with open(mop_file, "w") as f:
            f.write(mopac_input)
        
        return mop_file
    
    def run_calculation(self, input_file):
        """
        MOPAC計算を実行
        
        Args:
            input_file: 入力ファイルのパス
            
        Returns:
            dict: 計算結果の情報
        """
        input_path = Path(input_file)
        output_file = input_path.with_suffix('.out')
        arc_file = input_path.with_suffix('.arc')
        
        try:
            # MOPAC実行 - ファイル名から拡張子を除いて実行
            input_name = input_path.stem  # 拡張子なしのファイル名
            cmd = [self.mopac_executable, input_name]
            
            result = subprocess.run(
                cmd, 
                cwd=self.work_dir,
                capture_output=True, 
                text=True, 
                timeout=300  # 5分でタイムアウト
            )
            
            # 結果を解析
            # MOPACの成功判定: return code 0 かつ 出力ファイルが存在する
            mopac_success = (result.returncode == 0 and output_file.exists())
            
            calculation_result = {
                'success': mopac_success,
                'input_file': str(input_path),
                'output_file': str(output_file) if output_file.exists() else None,
                'arc_file': str(arc_file) if arc_file.exists() else None,
                'work_directory': str(self.work_dir),
                'output_files': self._get_output_files_info(input_path),
                'stdout': result.stdout,
                'stderr': result.stderr,
                'return_code': result.returncode,
                'command_executed': ' '.join(cmd)
            }
            
            # 出力ファイルが存在する場合、エネルギーなどの情報を抽出
            if output_file.exists():
                calculation_result.update(self._parse_output_file(output_file))
            
            return calculation_result
            
        except subprocess.TimeoutExpired:
            return {
                'success': False,
                'error': 'Calculation timed out',
                'input_file': str(input_path)
            }
        except Exception as e:
            return {
                'success': False,
                'error': str(e),
                'input_file': str(input_path)
            }
    
    def _get_output_files_info(self, input_path):
        """
        計算で生成された全ての出力ファイルの情報を取得
        
        Args:
            input_path: 入力ファイルのパス
            
        Returns:
            dict: 出力ファイルの詳細情報
        """
        input_path = Path(input_path)
        base_name = input_path.stem
        
        # MOPACが生成する可能性のあるファイル拡張子
        output_extensions = [
            '.out',    # メイン出力ファイル
            '.arc',    # アーカイブファイル（最適化の軌跡）
            '.log',    # ログファイル
            '.aux',    # 補助ファイル
            '.end',    # 終了ファイル
            '.mgf',    # Gaussian形式のファイル
            '.mol',    # MDL Molファイル
            '.pdb',    # PDBファイル
            '.xyz',    # XYZ座標ファイル
            '.den',    # 密度ファイル
            '.esp',    # 静電ポテンシャルファイル
            '.pol',    # 分極ファイル
        ]
        
        output_files = {}
        
        for ext in output_extensions:
            file_path = input_path.with_suffix(ext)
            if file_path.exists():
                file_stat = file_path.stat()
                output_files[ext[1:]] = {  # 拡張子から先頭の'.'を除く
                    'path': str(file_path),
                    'absolute_path': str(file_path.absolute()),
                    'size': file_stat.st_size,
                    'size_human': self._format_file_size(file_stat.st_size),
                    'modified_time': file_stat.st_mtime,
                    'exists': True
                }
            else:
                output_files[ext[1:]] = {
                    'path': str(file_path),
                    'absolute_path': str(file_path.absolute()),
                    'exists': False
                }
        
        return output_files
    
    def _format_file_size(self, size_bytes):
        """
        ファイルサイズを人間が読みやすい形式に変換
        
        Args:
            size_bytes: バイト単位のサイズ
            
        Returns:
            str: 読みやすい形式のサイズ
        """
        if size_bytes == 0:
            return "0 B"
        
        size_names = ["B", "KB", "MB", "GB"]
        i = 0
        size = float(size_bytes)
        
        while size >= 1024.0 and i < len(size_names) - 1:
            size /= 1024.0
            i += 1
        
        return f"{size:.1f} {size_names[i]}"
    
    def get_calculation_summary(self, result):
        """
        計算結果のサマリーを取得
        
        Args:
            result: run_calculationの結果
            
        Returns:
            dict: 計算結果のサマリー
        """
        summary = {
            'calculation_status': 'success' if result.get('success', False) else 'failed',
            'work_directory': result.get('work_directory', ''),
            'input_file': result.get('input_file', ''),
            'main_output_file': result.get('output_file', ''),
            'archive_file': result.get('arc_file', ''),
            'generated_files': [],
            'key_results': {}
        }
        
        # 生成されたファイルのリスト
        output_files = result.get('output_files', {})
        for file_type, file_info in output_files.items():
            if file_info.get('exists', False):
                summary['generated_files'].append({
                    'type': file_type,
                    'path': file_info['absolute_path'],
                    'size': file_info['size_human']
                })
        
        # 主要な計算結果
        if 'heat_of_formation' in result:
            summary['key_results']['heat_of_formation'] = f"{result['heat_of_formation']:.3f} kcal/mol"
        
        if 'electronic_energy' in result:
            summary['key_results']['electronic_energy'] = f"{result['electronic_energy']:.3f} eV"
        
        if 'dipole_moment' in result:
            summary['key_results']['dipole_moment'] = f"{result['dipole_moment']:.3f} Debye"
        
        if 'homo_energy' in result and 'lumo_energy' in result:
            summary['key_results']['homo_energy'] = f"{result['homo_energy']:.3f} eV"
            summary['key_results']['lumo_energy'] = f"{result['lumo_energy']:.3f} eV"
            summary['key_results']['homo_lumo_gap'] = f"{result['homo_lumo_gap']:.3f} eV"
        
        # エラー情報
        if not result.get('success', False):
            summary['error_info'] = result.get('error', result.get('stderr', ''))
        
        return summary
    
    def _parse_output_file(self, output_file):
        """
        MOPAC出力ファイルから重要な情報を抽出
        
        Args:
            output_file: 出力ファイルのパス
            
        Returns:
            dict: 抽出された計算結果
        """
        result = {}
        
        try:
            with open(output_file, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            
            # デバッグ用: ファイルサイズと最初の数行を記録
            result['file_size'] = len(content)
            result['first_lines'] = content.split('\n')[:10]
            
            # より柔軟なパターンで最終エネルギーを抽出
            energy_patterns = [
                r'FINAL HEAT OF FORMATION\s*=\s*([-+]?\d*\.?\d+)\s*KCAL/?MOL',
                r'HEAT OF FORMATION\s*=\s*([-+]?\d*\.?\d+)\s*KCAL/?MOL',
                r'FORMATION\s*=\s*([-+]?\d*\.?\d+)\s*KCAL/?MOL'
            ]
            
            for pattern in energy_patterns:
                energy_match = re.search(pattern, content, re.IGNORECASE)
                if energy_match:
                    result['heat_of_formation'] = float(energy_match.group(1))
                    break
            
            # 電子エネルギーを抽出
            electronic_patterns = [
                r'ELECTRONIC ENERGY\s*=\s*([-+]?\d*\.?\d+)\s*EV',
                r'TOTAL ENERGY\s*=\s*([-+]?\d*\.?\d+)\s*EV',
                r'ENERGY\s*=\s*([-+]?\d*\.?\d+)\s*EV'
            ]
            
            for pattern in electronic_patterns:
                electronic_match = re.search(pattern, content, re.IGNORECASE)
                if electronic_match:
                    result['electronic_energy'] = float(electronic_match.group(1))
                    break
            
            # 双極子モーメントを抽出
            dipole_patterns = [
                r'DIPOLE\s*=\s*([-+]?\d*\.?\d+)\s*DEBYE',
                r'DIPOLE MOMENT\s*=\s*([-+]?\d*\.?\d+)\s*DEBYE'
            ]
            
            for pattern in dipole_patterns:
                dipole_match = re.search(pattern, content, re.IGNORECASE)
                if dipole_match:
                    result['dipole_moment'] = float(dipole_match.group(1))
                    break
            
            # HOMO/LUMOエネルギーを抽出
            homo_patterns = [
                r'HOMO LUMO ENERGIES \(EV\)\s*=\s*([-+]?\d*\.?\d+)\s*([-+]?\d*\.?\d+)',
                r'HOMO\s*=\s*([-+]?\d*\.?\d+)\s*EV.*?LUMO\s*=\s*([-+]?\d*\.?\d+)\s*EV',
                r'LUMO\s*=\s*([-+]?\d*\.?\d+)\s*EV.*?HOMO\s*=\s*([-+]?\d*\.?\d+)\s*EV'
            ]
            
            for pattern in homo_patterns:
                homo_match = re.search(pattern, content, re.IGNORECASE | re.DOTALL)
                if homo_match:
                    if 'LUMO.*HOMO' in pattern:
                        # LUMOが先の場合
                        result['lumo_energy'] = float(homo_match.group(1))
                        result['homo_energy'] = float(homo_match.group(2))
                    else:
                        # HOMOが先の場合
                        result['homo_energy'] = float(homo_match.group(1))
                        result['lumo_energy'] = float(homo_match.group(2))
                    result['homo_lumo_gap'] = result['lumo_energy'] - result['homo_energy']
                    break
            
            # 計算が正常終了したかチェック
            termination_patterns = [
                "CALCULATION DONE",
                "NORMAL TERMINATION", 
                "JOB COMPLETED",
                "MOPAC DONE"
            ]
            
            result['normal_termination'] = any(pattern in content.upper() for pattern in termination_patterns)
            
            # エラーメッセージがあるかチェック
            if "ERROR" in content.upper():
                error_lines = [line.strip() for line in content.split('\n') 
                             if 'ERROR' in line.upper()]
                result['errors'] = error_lines
            
        except Exception as e:
            result['parsing_error'] = str(e)
        
        return result
    
    def optimize_geometry(self, 
                         theory="PM7", 
                         charge=0, 
                         multiplicity=1,
                         precise=True,
                         keywords=None,
                         title=None):
        """
        構造最適化を実行
        
        Args:
            theory: 計算理論
            charge: 分子の電荷
            multiplicity: スピン多重度
            precise: 高精度モードを使用するか
            keywords: 追加のキーワードリスト
            title: 計算のタイトル
            
        Returns:
            dict: 計算結果
        """
        base_keywords = []
        if precise:
            base_keywords.append("PRECISE")
        
        # 追加キーワードをマージ
        if keywords:
            base_keywords.extend(keywords)
        
        input_file = self.create_input_file(
            theory=theory,
            keywords=base_keywords,
            charge=charge,
            multiplicity=multiplicity,
            title=title,
            filename="geometry_opt"
        )
        
        result = self.run_calculation(input_file)
        result['summary'] = self.get_calculation_summary(result)
        return result
    
    def single_point_energy(self, 
                           theory="PM7", 
                           charge=0, 
                           multiplicity=1,
                           keywords=None,
                           title=None):
        """
        シングルポイントエネルギー計算を実行
        
        Args:
            theory: 計算理論
            charge: 分子の電荷
            multiplicity: スピン多重度
            keywords: 追加のキーワードリスト
            title: 計算のタイトル
            
        Returns:
            dict: 計算結果
        """
        base_keywords = ["1SCF"]  # 構造最適化なし
        
        # 追加キーワードをマージ
        if keywords:
            base_keywords.extend(keywords)
        
        input_file = self.create_input_file(
            theory=theory,
            keywords=base_keywords,
            charge=charge,
            multiplicity=multiplicity,
            title=title,
            filename="single_point"
        )
        
        result = self.run_calculation(input_file)
        result['summary'] = self.get_calculation_summary(result)
        return result
    
    def frequency_calculation(self, 
                            theory="PM7", 
                            charge=0, 
                            multiplicity=1,
                            title=None):
        """
        振動解析を実行
        
        Args:
            theory: 計算理論
            charge: 分子の電荷
            multiplicity: スピン多重度
            title: 計算のタイトル
            
        Returns:
            dict: 計算結果
        """
        keywords = ["FORCE", "THERMO"]  # 振動解析と熱化学量計算
        
        input_file = self.create_input_file(
            theory=theory,
            keywords=keywords,
            charge=charge,
            multiplicity=multiplicity,
            title=title,
            filename="frequency"
        )
        
        result = self.run_calculation(input_file)
        result['summary'] = self.get_calculation_summary(result)
        return result
    
    def molecular_orbital_calculation(self, 
                                    theory="PM7", 
                                    charge=0, 
                                    multiplicity=1,
                                    include_graph=True,
                                    include_vectors=True,
                                    precise=True,
                                    title=None):
        """
        分子軌道計算を実行（VECTORS、GRAPH MOキーワード付き）
        
        Args:
            theory: 計算理論
            charge: 分子の電荷
            multiplicity: スピン多重度
            include_graph: GRAPH MOキーワードを含めるか
            include_vectors: VECTORSキーワードを含めるか
            precise: PRECISEキーワードを含めるか
            title: 計算のタイトル
            
        Returns:
            dict: 計算結果
        """
        keywords = ["1SCF"]  # 構造最適化なし、SCFのみ
        
        if precise:
            keywords.append("PRECISE")
        
        if include_vectors:
            keywords.append("VECTORS")
        
        if include_graph:
            keywords.append("GRAPH MO")
        
        input_file = self.create_input_file(
            theory=theory,
            keywords=keywords,
            charge=charge,
            multiplicity=multiplicity,
            title=title,
            filename="molecular_orbital"
        )
        
        result = self.run_calculation(input_file)
        result['summary'] = self.get_calculation_summary(result)
        
        # 分子軌道ファイルの情報を追加
        if result.get('success'):
            mo_files = self._find_molecular_orbital_files()
            result['molecular_orbital_files'] = mo_files
        
        return result
    
    def optimize_and_mo_calculation(self, 
                                   theory="PM7", 
                                   charge=0, 
                                   multiplicity=1,
                                   include_graph=True,
                                   include_vectors=True,
                                   precise=True,
                                   title=None):
        """
        構造最適化後に分子軌道計算を実行
        
        Args:
            theory: 計算理論
            charge: 分子の電荷
            multiplicity: スピン多重度
            include_graph: GRAPH MOキーワードを含めるか
            include_vectors: VECTORSキーワードを含めるか
            precise: PRECISEキーワードを含めるか
            title: 計算のタイトル
            
        Returns:
            dict: 両計算の結果を含む辞書
        """
        # まず構造最適化を実行
        opt_result = self.optimize_geometry(
            theory=theory,
            charge=charge,
            multiplicity=multiplicity,
            precise=precise,
            title=title
        )
        
        mo_result = {}
        
        # 最適化が成功した場合、分子軌道計算を実行
        if opt_result.get('success'):
            # 最適化された構造から新しいMoleculeHandlerを作成
            arc_file = opt_result.get('arc_file')
            if arc_file and os.path.exists(arc_file):
                optimized_mol = self.read_optimized_structure(arc_file)
                
                if optimized_mol:
                    # 最適化された構造で分子軌道計算を実行
                    from logic.molecule_handler import MoleculeHandler
                    opt_handler = MoleculeHandler(optimized_mol, input_type="rdkit")
                    opt_calculator = MopacCalculator(opt_handler, work_dir=str(self.work_dir))
                    
                    mo_result = opt_calculator.molecular_orbital_calculation(
                        theory=theory,
                        charge=charge,
                        multiplicity=multiplicity,
                        include_graph=include_graph,
                        include_vectors=include_vectors,
                        precise=precise,
                        title=f"{title}_MO" if title else "MO_calculation"
                    )
                else:
                    mo_result = {
                        'success': False,
                        'error': 'Failed to read optimized structure for MO calculation'
                    }
            else:
                mo_result = {
                    'success': False,
                    'error': 'Optimization did not produce ARC file for MO calculation'
                }
        else:
            mo_result = {
                'success': False,
                'error': 'Optimization failed, MO calculation skipped'
            }
        
        return {
            'optimization': opt_result,
            'molecular_orbital': mo_result,
            'combined_success': opt_result.get('success', False) and mo_result.get('success', False)
        }
    
    def _find_molecular_orbital_files(self):
        """
        分子軌道計算で生成されたファイルを探す
        
        Returns:
            dict: 分子軌道関連ファイルの情報
        """
        mo_files = {}
        
        # MOPACが生成する可能性のある分子軌道関連ファイル
        mo_extensions = [
            '.mgf',     # Gaussian形式の分子軌道ファイル
            '.mol',     # MDL Molファイル
            '.den',     # 密度ファイル
            '.esp',     # 静電ポテンシャルファイル
            '.gpt',     # グラフポイントファイル
            '.po',      # 分子軌道プロットファイル
        ]
        
        for ext in mo_extensions:
            pattern = self.work_dir / f"*{ext}"
            import glob
            files = glob.glob(str(pattern))
            
            if files:
                mo_files[ext[1:]] = []
                for file_path in files:
                    file_stat = os.stat(file_path)
                    mo_files[ext[1:]].append({
                        'path': file_path,
                        'size': file_stat.st_size,
                        'size_human': self._format_file_size(file_stat.st_size),
                        'modified_time': file_stat.st_mtime
                    })
        
        return mo_files
    
    def read_optimized_structure(self, arc_file_path):
        """
        ARCファイルから最適化された構造を読み込んでRDKit Molオブジェクトを作成
        
        Args:
            arc_file_path: ARCファイルのパス
            
        Returns:
            Mol: 最適化された構造のRDKit Molオブジェクト（失敗時はNone）
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import rdDetermineBonds
            
            with open(arc_file_path, 'r') as f:
                content = f.read()
            
            print(f"Debug: ARC file content length: {len(content)}")
            
            # ARCファイルを行に分割
            lines = content.strip().split('\n')
            print(f"Debug: Number of lines in ARC file: {len(lines)}")
            
            # ARCファイルの最後の構造を探す
            # MOPACのARCファイルは複数のステップの座標を含むので、最後のものを取得
            coord_sets = []
            current_coords = []
            
            for i, line in enumerate(lines):
                line = line.strip()
                if not line:
                    continue
                    
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        # ARCファイルの形式を判定
                        # MOPACのARCファイル形式: 元素記号 X座標 フラグ Y座標 フラグ Z座標 フラグ
                        # 例: C    -2.52311728 +1   1.71470856 +1  -0.02973795 +1
                        
                        if len(parts) >= 7 and parts[0].isalpha():
                            # MOPACのARC形式: 元素記号 X +1 Y +1 Z +1
                            atom_symbol = parts[0].strip()
                            x = float(parts[1])
                            y = float(parts[3])  # フラグをスキップしてY座標
                            z = float(parts[5])  # フラグをスキップしてZ座標
                        elif len(parts) >= 5 and parts[0].isdigit():
                            # 番号付きの場合
                            atom_symbol = parts[1].strip()
                            x = float(parts[2])
                            y = float(parts[3])
                            z = float(parts[4])
                        elif len(parts) >= 4:
                            # 番号なしの基本形式
                            atom_symbol = parts[0].strip()
                            x = float(parts[1])
                            y = float(parts[2])
                            z = float(parts[3])
                        else:
                            # 解析不可能な行
                            continue
                        
                        print(f"Debug ARC parsing: Line {i}: {line}")
                        print(f"Debug ARC parsing: Parts: {parts}")
                        print(f"Debug ARC parsing: Parsed coordinates: {atom_symbol} {x} {y} {z}")
                        
                        # 有効な原子記号かチェック
                        if atom_symbol.isalpha() and len(atom_symbol) <= 2:
                            current_coords.append((atom_symbol, x, y, z))
                            print(f"Debug ARC: Added {atom_symbol} ({x:.6f}, {y:.6f}, {z:.6f}) to current_coords")
                        else:
                            # 座標セットの終了と判断
                            if current_coords:
                                coord_sets.append(current_coords)
                                current_coords = []
                    except (ValueError, IndexError) as e:
                        print(f"Debug ARC parsing: Error parsing line {i}: {line} - {e}")
                        # 座標として解釈できない行
                        if current_coords:
                            coord_sets.append(current_coords)
                            current_coords = []
                        continue
            
            # 最後の座標セットを追加
            if current_coords:
                coord_sets.append(current_coords)
            
            if not coord_sets:
                print("Debug: No coordinate sets found in ARC file")
                return None
            
            # 最後の座標セット（最適化された構造）を使用
            final_coords = coord_sets[-1]
            print(f"Debug: Found {len(coord_sets)} coordinate sets, using the last one with {len(final_coords)} atoms")
            print(f"Debug ARC: Final coordinates to be used:")
            for i, (symbol, x, y, z) in enumerate(final_coords):
                print(f"  {i}: {symbol} ({x:.6f}, {y:.6f}, {z:.6f})")
            
            if not final_coords:
                return None
            
            # RDKit Molオブジェクトを作成
            mol = Chem.RWMol()
            
            # 有効な原子のみを追加
            valid_coords = []
            for atom_symbol, x, y, z in final_coords:
                atom_symbol = atom_symbol.strip().upper()
                if atom_symbol in ['H', 'C', 'N', 'O', 'P', 'S', 'F', 'CL', 'BR', 'I', 'B', 'SI', 'AL']:
                    atom = Chem.Atom(atom_symbol)
                    mol.AddAtom(atom)
                    valid_coords.append((x, y, z))
                    print(f"Debug: Added atom {atom_symbol} at ({x:.3f}, {y:.3f}, {z:.3f})")
                else:
                    print(f"Debug: Skipped invalid atom symbol: {atom_symbol}")
            
            if mol.GetNumAtoms() == 0:
                print("Debug: No valid atoms added to molecule")
                return None
            
            # コンフォーマーを追加
            conf = Chem.Conformer(mol.GetNumAtoms())
            for i, (x, y, z) in enumerate(valid_coords):
                print(f"Debug SetAtomPosition: Atom {i}: Setting position to ({x:.6f}, {y:.6f}, {z:.6f})")
                conf.SetAtomPosition(i, (x, y, z))
                # 設定後の確認
                set_pos = conf.GetAtomPosition(i)
                print(f"Debug SetAtomPosition: Atom {i}: Actual position after setting: ({set_pos.x:.6f}, {set_pos.y:.6f}, {set_pos.z:.6f})")
            
            mol.AddConformer(conf)
            
            # 結合を推定
            mol = mol.GetMol()
            rdDetermineBonds.DetermineBonds(mol)
            
            print(f"Debug: Successfully created molecule with {mol.GetNumAtoms()} atoms")
            return mol
            
        except Exception as e:
            print(f"Error reading optimized structure from {arc_file_path}: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def read_structure_from_output(self, output_file_path):
        """
        出力ファイルから最適化された構造を読み込んでRDKit Molオブジェクトを作成
        
        Args:
            output_file_path: 出力ファイルのパス
            
        Returns:
            Mol: 最適化された構造のRDKit Molオブジェクト（失敗時はNone）
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import rdDetermineBonds
            
            with open(output_file_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            
            print(f"Debug: Output file content length: {len(content)}")
            
            # "CARTESIAN COORDINATES"セクションを探す
            lines = content.split('\n')
            coord_lines = []
            capturing = False
            
            # 最後の座標セクションを探すため、逆順で処理
            for i in range(len(lines) - 1, -1, -1):
                line = lines[i].strip()
                
                # 座標セクションの終了を検出（逆順なので）
                if 'CARTESIAN COORDINATES' in line.upper():
                    # この行から前向きに座標を読む
                    capturing = True
                    coord_lines = []
                    
                    # このセクションから座標を抽出
                    for j in range(i + 1, len(lines)):
                        coord_line = lines[j].strip()
                        if not coord_line:
                            continue
                            
                        parts = coord_line.split()
                        if len(parts) >= 4:
                            try:
                                # MOPACの出力形式を判定
                                # 「番号 元素記号 X Y Z」の形式か「元素記号 X Y Z」の形式かをチェック
                                if len(parts) >= 5 and parts[0].isdigit():
                                    # 番号付きの場合: 1 C -2.523117279 1.714708561 -0.029737947
                                    atom_symbol = parts[1].strip()
                                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                                else:
                                    # 番号なしの場合: C -2.523117279 1.714708561 -0.029737947
                                    atom_symbol = parts[0].strip()
                                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                                
                                print(f"Debug OUTPUT parsing: Line: {coord_line}")
                                print(f"Debug OUTPUT parsing: Parts: {parts}")
                                print(f"Debug OUTPUT parsing: Parsed coordinates: {atom_symbol} {x} {y} {z}")
                                
                                # 有効な元素記号かチェック
                                if atom_symbol.isalpha() and len(atom_symbol) <= 2:
                                    coord_lines.append((atom_symbol, x, y, z))
                                    print(f"Debug OUTPUT: Added {atom_symbol} ({x:.6f}, {y:.6f}, {z:.6f}) to coord_lines")
                                else:
                                    # 座標セクションの終了
                                    break
                            except (ValueError, IndexError) as e:
                                print(f"Debug OUTPUT parsing: Error parsing line: {coord_line} - {e}")
                                # 座標セクションの終了
                                break
                        else:
                            # 空行や短い行は無視
                            if len(parts) == 0:
                                continue
                            else:
                                # 座標セクションの終了
                                break
                    
                    # 最後の（最新の）座標セクションが見つかったので終了
                    break
            
            if not coord_lines:
                print("Debug: No coordinate section found in output file")
                return None
            
            print(f"Debug: Found {len(coord_lines)} atoms in output file")
            print(f"Debug OUTPUT: Final coordinates to be used:")
            for i, (symbol, x, y, z) in enumerate(coord_lines):
                print(f"  {i}: {symbol} ({x:.6f}, {y:.6f}, {z:.6f})")
            
            # RDKit Molオブジェクトを作成
            mol = Chem.RWMol()
            
            # 有効な原子のみを追加
            valid_coords = []
            for atom_symbol, x, y, z in coord_lines:
                atom_symbol = atom_symbol.strip().upper()
                if atom_symbol in ['H', 'C', 'N', 'O', 'P', 'S', 'F', 'CL', 'BR', 'I', 'B', 'SI', 'AL']:
                    atom = Chem.Atom(atom_symbol)
                    mol.AddAtom(atom)
                    valid_coords.append((x, y, z))
                    print(f"Debug: Added atom {atom_symbol} at ({x:.3f}, {y:.3f}, {z:.3f})")
                else:
                    print(f"Debug: Skipped invalid atom symbol: {atom_symbol}")
            
            if mol.GetNumAtoms() == 0:
                print("Debug: No valid atoms added to molecule from output file")
                return None
            
            # コンフォーマーを追加
            conf = Chem.Conformer(mol.GetNumAtoms())
            for i, (x, y, z) in enumerate(valid_coords):
                conf.SetAtomPosition(i, (x, y, z))
            
            mol.AddConformer(conf)
            
            # 結合を推定
            mol = mol.GetMol()
            rdDetermineBonds.DetermineBonds(mol)
            
            print(f"Debug: Successfully created molecule from output with {mol.GetNumAtoms()} atoms")
            return mol
            
        except Exception as e:
            print(f"Error reading structure from output file {output_file_path}: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def cleanup(self):
        """作業ディレクトリを削除"""
        if self.work_dir.exists():
            shutil.rmtree(self.work_dir)

    def optimize_multiple_conformers(self,
                                    theory="PM7",
                                    charge=0,
                                    multiplicity=1,
                                    keywords=None,
                                    title="Multi-conformer optimization",
                                    max_conformers=None,
                                    num_conformers=10,
                                    max_iterations=200,
                                    prune_rms_threshold=0.5,
                                    forcefield="MMFF",
                                    include_mo_calculation=False):
        """
        複数の配座を生成してMOPAC構造最適化を実行
        
        Args:
            theory: 計算理論
            charge: 分子の電荷
            multiplicity: スピン多重度
            keywords: 追加キーワード
            title: 計算のタイトル
            max_conformers: 計算する最大配座数（Noneで全配座）
            num_conformers: 生成する配座数
            max_iterations: 配座生成の最大反復回数
            prune_rms_threshold: 配座の類似性判定閾値
            forcefield: 分子力場（UFF or MMFF）
            include_mo_calculation: 最適化後に分子軌道計算を実行するか
            
        Returns:
            dict: 各配座の計算結果を含む辞書
        """
        # MoleculeHandlerを使って配座を生成
        print(f"Generating {num_conformers} conformers using {forcefield} force field...")
        try:
            conformers = self.molecule_handler.generate_conformers(
                num_conformers=num_conformers,
                max_iterations=max_iterations,
                prune_rms_threshold=prune_rms_threshold,
                forcefield=forcefield
            )
            
            if not conformers:
                raise ValueError("No conformers were generated")
            
            print(f"Successfully generated {len(conformers)} conformers")
            
        except Exception as e:
            raise ValueError(f"Failed to generate conformers: {e}")
        
        mol = self.molecule_handler.mol
        if not mol or mol.GetNumConformers() == 0:
            raise ValueError("No conformers available after generation")
        
        num_available_conformers = mol.GetNumConformers()
        if max_conformers is not None:
            num_available_conformers = min(num_available_conformers, max_conformers)
        
        print(f"Proceeding with MOPAC optimization for {num_available_conformers} conformers...")
        
        results = {}
        conformer_energies = []
        
        for conf_id in range(num_available_conformers):
            try:
                print(f"Processing conformer {conf_id + 1}/{num_available_conformers}...")
                
                # 特定の配座のみを使用するためのMoleculeHandlerを作成
                conf_mol = Chem.Mol(mol)
                conf_mol.RemoveAllConformers()
                
                # 対象配座をコピー
                original_conf = mol.GetConformer(conf_id)
                new_conf = Chem.Conformer(conf_mol.GetNumAtoms())
                for i in range(conf_mol.GetNumAtoms()):
                    pos = original_conf.GetAtomPosition(i)
                    new_conf.SetAtomPosition(i, pos)
                conf_mol.AddConformer(new_conf, assignId=True)
                
                # 一時的なMoleculeHandlerを作成
                from logic.molecule_handler import MoleculeHandler
                temp_handler = MoleculeHandler(conf_mol, input_type="rdkit")
                
                # 配座専用の作業ディレクトリを作成
                conf_work_dir = self.work_dir / f"conformer_{conf_id}"
                conf_work_dir.mkdir(parents=True, exist_ok=True)
                
                # 配座専用のMopacCalculatorを作成
                conf_calculator = MopacCalculator(temp_handler, work_dir=str(conf_work_dir))
                
                # 最適化を実行
                print(f"Running MOPAC optimization for conformer {conf_id}...")
                result = conf_calculator.optimize_geometry(
                    theory=theory,
                    charge=charge,
                    multiplicity=multiplicity,
                    keywords=keywords,
                    title=f"{title}_conformer_{conf_id}"
                )
                
                # 結果を保存
                result['conformer_id'] = conf_id
                result['initial_energy'] = None
                result['initial_forcefield'] = forcefield
                
                # 初期エネルギー（分子力場）が利用可能な場合は取得
                if mol.HasProp(f"Energy_{conf_id}"):
                    try:
                        result['initial_energy'] = float(mol.GetProp(f"Energy_{conf_id}"))
                        print(f"Conformer {conf_id} - Initial energy ({forcefield}): {result['initial_energy']:.6f} kcal/mol")
                    except Exception as e:
                        print(f"Warning: Could not parse initial energy for conformer {conf_id}: {e}")
                
                # MOPAC最適化結果
                if result.get('success') and 'heat_of_formation' in result:
                    print(f"Conformer {conf_id} - MOPAC optimization successful")
                    print(f"Conformer {conf_id} - Heat of formation: {result['heat_of_formation']:.6f} kcal/mol")
                    
                    # エネルギー変化を計算
                    if result['initial_energy'] is not None:
                        energy_change = result['heat_of_formation'] - result['initial_energy']
                        result['energy_change'] = energy_change
                        print(f"Conformer {conf_id} - Energy change: {energy_change:.6f} kcal/mol")
                    
                    # 分子軌道計算の実行（オプション）
                    if include_mo_calculation:
                        print(f"Running molecular orbital calculation for conformer {conf_id}...")
                        try:
                            # 最適化された構造を読み込み
                            arc_file = result.get('arc_file')
                            if arc_file and os.path.exists(arc_file):
                                optimized_mol = conf_calculator.read_optimized_structure(arc_file)
                                
                                if optimized_mol:
                                    # 最適化された構造で分子軌道計算
                                    opt_handler = MoleculeHandler(optimized_mol, input_type="rdkit")
                                    mo_calculator = MopacCalculator(opt_handler, work_dir=str(conf_work_dir))
                                    
                                    mo_result = mo_calculator.molecular_orbital_calculation(
                                        theory=theory,
                                        charge=charge,
                                        multiplicity=multiplicity,
                                        include_graph=True,
                                        include_vectors=True,
                                        title=f"{title}_conformer_{conf_id}_MO"
                                    )
                                    
                                    result['molecular_orbital'] = mo_result
                                    if mo_result.get('success'):
                                        print(f"Conformer {conf_id} - MO calculation successful")
                                    else:
                                        print(f"Conformer {conf_id} - MO calculation failed")
                                else:
                                    print(f"Conformer {conf_id} - Failed to read optimized structure for MO calculation")
                            else:
                                print(f"Conformer {conf_id} - No ARC file found for MO calculation")
                        
                        except Exception as mo_error:
                            print(f"Conformer {conf_id} - MO calculation error: {mo_error}")
                            result['molecular_orbital'] = {
                                'success': False,
                                'error': str(mo_error)
                            }
                    
                else:
                    print(f"Conformer {conf_id} - MOPAC optimization failed")
                
                results[conf_id] = result
                
                # エネルギー情報を収集（成功した計算のみ）
                if result.get('success') and 'heat_of_formation' in result:
                    conformer_energies.append({
                        'conformer_id': conf_id,
                        'heat_of_formation': result['heat_of_formation'],
                        'initial_energy': result['initial_energy'],
                        'energy_change': result.get('energy_change'),
                        'forcefield': forcefield
                    })
                
            except Exception as e:
                # 個別の配座で失敗した場合もログに記録
                print(f"Error processing conformer {conf_id}: {e}")
                results[conf_id] = {
                    'conformer_id': conf_id,
                    'success': False,
                    'error': str(e),
                    'initial_energy': None,
                    'initial_forcefield': forcefield
                }
                
                # 初期エネルギーは取得を試行
                if mol.HasProp(f"Energy_{conf_id}"):
                    try:
                        results[conf_id]['initial_energy'] = float(mol.GetProp(f"Energy_{conf_id}"))
                    except:
                        pass
        
        # 全体的な結果サマリーを作成
        summary_result = {
            'total_conformers': num_available_conformers,
            'successful_calculations': sum(1 for r in results.values() if r.get('success', False)),
            'failed_calculations': sum(1 for r in results.values() if not r.get('success', False)),
            'conformer_results': results,
            'energy_ranking': sorted(conformer_energies, key=lambda x: x['heat_of_formation']),
            'best_conformer': None,
            'energy_range': None,
            'forcefield_used': forcefield,
            'generation_parameters': {
                'num_conformers_requested': num_conformers,
                'max_iterations': max_iterations,
                'prune_rms_threshold': prune_rms_threshold,
                'forcefield': forcefield
            }
        }
        
        # 最安定配座を特定
        if conformer_energies:
            best_conf = min(conformer_energies, key=lambda x: x['heat_of_formation'])
            summary_result['best_conformer'] = best_conf
            
            energies = [conf['heat_of_formation'] for conf in conformer_energies]
            summary_result['energy_range'] = {
                'min': min(energies),
                'max': max(energies),
                'span': max(energies) - min(energies)
            }
            
            print(f"Best conformer: {best_conf['conformer_id']} with heat of formation {best_conf['heat_of_formation']:.6f} kcal/mol")
            print(f"Energy range: {summary_result['energy_range']['span']:.6f} kcal/mol")
        
        return summary_result


def run_mopac_calculation(molecule_handler, 
                         calculation_type="optimize",
                         theory="PM7",
                         charge=0,
                         multiplicity=1,
                         work_dir=None,
                         include_mo=False):
    """
    MoleculeHandlerを使ってMOPAC計算を実行する便利関数
    
    Args:
        molecule_handler: MoleculeHandlerのインスタンス
        calculation_type: 計算の種類 ("optimize", "single_point", "frequency", "molecular_orbital", "optimize_and_mo")
        theory: 計算理論
        charge: 分子の電荷
        multiplicity: スピン多重度
        work_dir: 作業ディレクトリ
        include_mo: 分子軌道計算を含めるか（optimize時のみ有効）
        
    Returns:
        dict: 計算結果
    """
    calculator = MopacCalculator(molecule_handler, work_dir=work_dir)
    
    try:
        if calculation_type == "optimize":
            if include_mo:
                result = calculator.optimize_and_mo_calculation(
                    theory=theory, 
                    charge=charge, 
                    multiplicity=multiplicity
                )
            else:
                result = calculator.optimize_geometry(
                    theory=theory, 
                    charge=charge, 
                    multiplicity=multiplicity
                )
        elif calculation_type == "single_point":
            result = calculator.single_point_energy(
                theory=theory, 
                charge=charge, 
                multiplicity=multiplicity
            )
        elif calculation_type == "frequency":
            result = calculator.frequency_calculation(
                theory=theory, 
                charge=charge, 
                multiplicity=multiplicity
            )
        elif calculation_type == "molecular_orbital":
            result = calculator.molecular_orbital_calculation(
                theory=theory, 
                charge=charge, 
                multiplicity=multiplicity
            )
        elif calculation_type == "optimize_and_mo":
            result = calculator.optimize_and_mo_calculation(
                theory=theory, 
                charge=charge, 
                multiplicity=multiplicity
            )
        else:
            raise ValueError(f"Unsupported calculation type: {calculation_type}")
        
        return result
        
    finally:
        # 作業ディレクトリが一時的に作成された場合は削除
        if work_dir is None:
            calculator.cleanup()

