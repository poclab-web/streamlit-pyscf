import os
import ast
import streamlit as st

def load_readme(file_path="README.md"):
    with open(file_path, "r", encoding="utf-8") as file:
        return file.read()

def get_module_docstrings(directory):
    """
    指定されたディレクトリ内のPythonモジュールファイルを探索し、
    そのファイル名とdocstringを取得します。

    Parameters:
        directory (str): モジュールが含まれるディレクトリのパス。

    Returns:
        list of tuple: ファイル名とそのdocstringのペアのリスト。
    """
    docstrings = []
    for filename in os.listdir(directory):
        if filename.endswith(".py"):  # Pythonファイルのみを対象
            filepath = os.path.join(directory, filename)
            with open(filepath, "r", encoding="utf-8") as file:
                try:
                    module = ast.parse(file.read())
                    docstring = ast.get_docstring(module)
                    docstrings.append((filename, docstring or "No docstring found"))
                except Exception as e:
                    docstrings.append((filename, f"Error parsing file: {e}"))
    return docstrings


# モジュール名から不要な接頭辞（例: '1_'や'2_'）を削除
def clean_module_name(filename):
    """
    モジュール名から接頭辞を削除し、拡張子を除いた名前を返します。

    Parameters:
        filename (str): ファイル名。

    Returns:
        str: 整形されたモジュール名。
    """
    module_name = filename.replace(".py", "")  # 拡張子を除去
    module_name = module_name.lstrip("0123456789_")  # 接頭辞の数字とアンダースコアを削除
    return module_name