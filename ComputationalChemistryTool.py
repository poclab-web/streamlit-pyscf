import os
import ast
import streamlit as st

from utils.module import get_module_docstrings, clean_module_name, load_readme

# Streamlit App
st.title("Conputational Chemistry Tool")


st.markdown("## Overview")

# README.mdの内容を表示
readme_content = load_readme()
st.markdown(readme_content, unsafe_allow_html=True)

# `pages`フォルダの存在確認
pages_directory = "pages"  # `pages`フォルダのパス
if os.path.exists(pages_directory) and os.path.isdir(pages_directory):
    module_docstrings = get_module_docstrings(pages_directory)
        # __init__.py を除外
    # __init__.py を除外（リスト形式）
    module_docstrings = [
        (filename, docstring)
        for filename, docstring in module_docstrings
        if filename != '__init__.py'
    ]
    sorted_files = sorted(module_docstrings)

    if sorted_files:
        st.markdown("---")
        st.header("Streamlit Pages")
        st.markdown("以下のページが、あります。")
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


# フォルダパスを固定
logic_directory = "logic"  # `logic`フォルダのパス

# フォルダの存在確認
if os.path.exists(logic_directory) and os.path.isdir(logic_directory):
    # フォルダ内のPythonファイルをリストアップ
    files = [f for f in os.listdir(logic_directory) if f.endswith('.py') and f != '__init__.py']
    
    if files:
        st.markdown("---")
        st.header("Logic Pages")
        st.markdown("以下は、各種コードの中身です。")
        for file in files:
            # ファイルの中身を折りたたんで表示
            with st.expander(f"{file}"):
                file_path = os.path.join(logic_directory, file)
                with open(file_path, "r", encoding="utf-8") as f:
                    file_content = f.read()
                st.code(file_content, language="python")
    else:
        st.warning("logicフォルダ内にPythonファイルが見つかりません。")
else:
    st.error("logicフォルダが存在しません。")

