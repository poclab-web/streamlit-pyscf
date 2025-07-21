import os
import sys
sys.path.insert(0, os.path.abspath(".."))  # streamlit-pyscf/ をパスに追加

project = "streamlit-pyscf"
copyright = "2025, poclab"
author = "poclab"
release = "0.1.0"
version = "0.1"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.githubpages",
    "myst_parser",  # Markdown対応
    "sphinx_copybutton",  # コピーボタン
    "sphinx.ext.autosummary",
]

# AutoDoc設定
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': True,
    'exclude-members': '__weakref__'
}

# Napoleon設定
napoleon_google_docstring = True
napoleon_numpy_docstring = True

# MyST設定
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "tasklist",
    "html_image",
]

templates_path = ["_templates"]
exclude_patterns = []

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

# テーマオプション
html_theme_options = {
    'analytics_id': '',
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

# Intersphinx設定
intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'rdkit': ('https://www.rdkit.org/docs/', None),
}
