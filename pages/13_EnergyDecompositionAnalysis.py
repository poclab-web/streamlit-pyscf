"""
個々のフラグメント（分子 A, 分子 B）のエネルギーを計算。
フラグメントが結合した複合体（A + B）のエネルギーを計算。
相互作用エネルギーを計算：
"""

import streamlit as st
from utils.module import load_css

# カスタムCSSを適用
load_css("config/styles.css")

