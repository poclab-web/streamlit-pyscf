import os
import ast
import streamlit as st

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


# Streamlit アプリの設定
st.title("Conputational Chemistry Tool")

# `pages`フォルダの存在確認
pages_directory = "pages"  # `pages`フォルダのパス
if os.path.exists(pages_directory) and os.path.isdir(pages_directory):
    module_docstrings = get_module_docstrings(pages_directory)
    sorted_files = sorted(module_docstrings)

    if sorted_files:
        st.header("モジュール一覧")
        for i, (filename, docstring) in enumerate(sorted_files, start=1): 
            # モジュール名を整形
            clean_name = clean_module_name(filename)
            # 各モジュール名を作成
            st.markdown(f"### {i} - {clean_name}")

            # Docstring を折りたたみ可能に表示
            with st.expander("Docstring を表示"):
                st.code(docstring, language="python")

    else:
        st.info("`pages` フォルダには Python モジュールが見つかりませんでした。")
else:
    st.error(f"ディレクトリ '{pages_directory}' が存在しないか、有効なディレクトリではありません。")
