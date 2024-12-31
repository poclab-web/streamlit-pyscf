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
    sorted_files = sorted(module_docstrings)

    if sorted_files:
        st.markdown("---")
        st.header("Streamlit Pages")
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


# `logic`フォルダの存在確認
logic_directory = "logic"  # `logic`フォルダのパス
if os.path.exists(logic_directory) and os.path.isdir(logic_directory):
    module_docstrings = get_module_docstrings(logic_directory)
    sorted_files = sorted(module_docstrings)

    if sorted_files:
        st.markdown("---")
        st.header("logic Pages")
        for filename, docstring in sorted_files: 
            # モジュール名を整形
            clean_name = clean_module_name(filename)
            # 各モジュール名を作成
            st.markdown(f"### {clean_name}")

            # Docstring を折りたたみ可能に表示
            with st.expander("Docstring を表示"):
                st.code(docstring, language="python")

    else:
        st.info("`logic` フォルダには Python モジュールが見つかりませんでした。")
else:
    st.error(f"ディレクトリ '{logic_directory}' が存在しないか、有効なディレクトリではありません。")
