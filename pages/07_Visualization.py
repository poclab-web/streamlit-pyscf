"""
軌道の計算と表示を行うもの
TODO: 改修中
計算後のcheck pointを呼び出して行う。
静電ポテンシャル（ESP）の等値面も表示可能。
"""

import os
import tempfile

import streamlit as st
import streamlit.components.v1 as components
import py3Dmol
import stmol

from pyscf import gto, scf, tools, lib
from pyscf.tools import cubegen
import numpy as np
import io
from rdkit import Chem
from rdkit.Chem import Draw

def list_chk_files(data_dir="data"):
    chk_paths = []
    for root, _, files in os.walk(data_dir):
        for f in files:
            if f.endswith(".chk"):
                chk_paths.append(os.path.join(root, f))
    return chk_paths

def load_mo_info(chk_path):
    try:
        mol = lib.chkfile.load_mol(chk_path)
        mf = scf.RHF(mol)
        mf.__dict__.update(lib.chkfile.load(chk_path, 'scf'))
        mo_occ = np.array(mf.mo_occ)
        homo_idx = np.where(mo_occ > 0)[0][-1]
        lumo_idx = homo_idx + 1
        return mol, mf.mo_coeff, mo_occ, homo_idx, lumo_idx, mf
    except Exception as e:
        st.error(f"❌ 読み込み失敗: {e}")
        return None, None, None, None, None

def generate_xyz_from_mol(mol):
    coords = mol.atom_coords() * 0.529177
    xyz_lines = []
    for i in range(mol.natm):
        symbol = mol.atom_symbol(i)
        x, y, z = coords[i]
        xyz_lines.append(f"{symbol} {x:.5f} {y:.5f} {z:.5f}")
    xyz_str = f"{mol.natm}\nGenerated from mol.atom_coords()\n" + "\n".join(xyz_lines)
    return xyz_str

def read_cube_values(cube_path):
    with open(cube_path) as f:
        lines = f.readlines()
    n_atoms = int(lines[2].split()[0])
    data_lines = lines[6 + n_atoms:]
    values = np.array([float(v) for line in data_lines for v in line.strip().split()])
    return values

def save_cube_from_data(filename, cube_data):
    try:
        with open(filename, "w") as f:
            f.write(cube_data)
        st.success(f"✅ {filename} を保存しました。")
    except Exception as e:
        st.error(f"❌ 保存に失敗しました: {e}")

def mol_to_rdkit_mol(mol):
    """pyscfのmolからxyz文字列を作り、RDKit Molに変換"""
    xyz = generate_xyz_from_mol(mol)
    rdkit_mol = Chem.MolFromXYZBlock(xyz)
    if rdkit_mol is None:
        # xyzからMolが作れない場合はSMILES経由なども検討
        return None
    return rdkit_mol

def draw_rdkit_mol_with_numbers(rdkit_mol):
    """原子番号付きで2D画像を生成"""
    if rdkit_mol is None:
        return None
    # 原子番号ラベルを付与
    atom_nums = {i: str(atom.GetAtomicNum()) for i, atom in enumerate(rdkit_mol.GetAtoms())}
    img = Draw.MolToImage(rdkit_mol, size=(400, 300), legend="", 
                          highlightAtomColors=None, 
                          highlightAtoms=None, 
                          kekulize=True, 
                          wedgeBonds=True, 
                          atomLabels=atom_nums)
    return img

st.title("🧮 任意の軌道を選んで可視化 & 保存")

chk_files = list_chk_files("data")
if chk_files:
    selected_chk = st.selectbox("📁 .chk ファイルを選んでください", chk_files)

    if st.button("📊 軌道情報を読み込む"):
        mol, mo_coeff, mo_occ, homo_idx, lumo_idx, mf = load_mo_info(selected_chk)
        if mol is not None:
            st.session_state["mol"] = mol
            st.session_state["mf"] = mf
            st.session_state["mo_coeff"] = mo_coeff
            st.session_state["selected_chk"] = selected_chk
            st.session_state["homo_idx"] = homo_idx
            st.session_state["lumo_idx"] = lumo_idx
            st.session_state["orbital_idx"] = homo_idx
            st.session_state.pop("cube_data", None)
            st.session_state.pop("cube_values", None)
            st.session_state.pop("cached_mode", None)

if "mol" in st.session_state and "mo_coeff" in st.session_state:
    mol = st.session_state["mol"]
    mf = st.session_state["mf"]
    mo_coeff = st.session_state["mo_coeff"]
    selected_chk = st.session_state["selected_chk"]
    homo_idx = st.session_state["homo_idx"]
    lumo_idx = st.session_state["lumo_idx"]
    orbital_idx = st.number_input("表示する軌道番号", 0, mo_coeff.shape[1] - 1,
                                  value=st.session_state.get("orbital_idx", homo_idx))
    st.session_state["orbital_idx"] = orbital_idx

    mode = st.radio("表示対象", ["軌道", "電子密度", "静電ポテンシャル"])
    mode_map = {"軌道": "orbital", "電子密度": "density", "静電ポテンシャル": "esp"}
    internal_mode = mode_map[mode]


    if ("cube_data" not in st.session_state or
        st.session_state.get("cached_mode") != internal_mode or
        (internal_mode == "orbital" and st.session_state.get("cached_orbital_idx") != orbital_idx)):

        with st.spinner(f"{mode} cube を生成中..."):
            with tempfile.NamedTemporaryFile(delete=False, suffix=".cube") as tmp:
                cube_path = tmp.name
                if internal_mode == "orbital":
                    tools.cubegen.orbital(mol, cube_path, mo_coeff[:, orbital_idx])
                elif internal_mode == "density":
                    tools.cubegen.density(mol, cube_path, mf.make_rdm1())
                elif internal_mode == "esp":
                    cubegen.mep(mol, cube_path, mf.make_rdm1())

                with open(cube_path, "r") as f:
                    st.session_state["cube_data"] = f.read()
                st.session_state["cube_values"] = read_cube_values(cube_path)
                st.session_state["cached_mode"] = internal_mode
                st.session_state["cached_orbital_idx"] = orbital_idx

    cube_data = st.session_state["cube_data"]
    values = st.session_state["cube_values"]

    if internal_mode == "density":
        # Mulliken電荷の計算と表示
        st.subheader("Mulliken 電荷")
        
        try:
            pop, chg = mf.mulliken_pop()
            atom_syms = [mol.atom_symbol(i) for i in range(mol.natm)]
            import pandas as pd
            df = pd.DataFrame({
                "atom": atom_syms,
                "Mulliken charge": chg
            })
            st.table(df)
        except Exception as e:
            st.error(f"Mulliken電荷の計算に失敗しました: {e}")

    if internal_mode == "orbital":
        suggested_isoval = float(np.percentile(np.abs(values), 95))

        # --- フロンティア軌道の原子ごと電子密度を計算・表示 ---
        st.subheader("選択軌道の原子ごとのフロンティア電子密度寄与")
        ao_slices = mol.aoslice_by_atom()
        homo_contrib = []
        for i_atom in range(mol.natm):
            p0, p1 = ao_slices[i_atom][2], ao_slices[i_atom][3]
            contrib = float((mo_coeff[p0:p1, orbital_idx] ** 2).sum())
            homo_contrib.append(contrib)
        # 正規化（全体で1になるように）
        total = sum(homo_contrib)
        contrib_norm = [c / total for c in homo_contrib]
        import pandas as pd
        df_contrib = pd.DataFrame({
            "atom": [mol.atom_symbol(i) for i in range(mol.natm)],
            "contribution": contrib_norm
        })
        st.table(df_contrib)
    elif internal_mode == "density":
        suggested_isoval = float(np.percentile(values, 99))
    elif internal_mode == "esp":
        min_val = np.percentile(values, 5)
        max_val = np.percentile(values, 95)
        suggested_isoval = float((abs(min_val) + abs(max_val)) / 2)
    else:
        suggested_isoval = 0.02

    st.markdown(f"**データ範囲**: min={values.min():.4f}, max={values.max():.4f}, 推奨 isoval ≈ {suggested_isoval:.4f}")

    opacity = st.slider("透明度 (opacity)", 0.0, 1.0, 0.6, step=0.05)
    color_pos = st.color_picker("正の等値面の色", "#FF0000")
    color_neg = st.color_picker("負の等値面の色", "#0000FF")

    xyz_data = generate_xyz_from_mol(mol)
    viewer = py3Dmol.view(width=500, height=500)
    viewer.setBackgroundColor("white")
    viewer.addModel(xyz_data, "xyz")
    viewer.setStyle({"stick": {}, "sphere": {"scale": 0.3}})

    if internal_mode == "esp":
        pos_isoval = st.slider("正の等値面 (isoval)", 0.001, 0.1, min(suggested_isoval, 0.05), step=0.001)
        neg_isoval = st.slider("負の等値面 (isoval)", 0.001, 0.1, min(suggested_isoval, 0.05), step=0.001)
        viewer.addVolumetricData(cube_data, "cube", {"isoval": pos_isoval, "color": color_pos, "opacity": opacity})
        viewer.addVolumetricData(cube_data, "cube", {"isoval": -neg_isoval, "color": color_neg, "opacity": opacity})
    else:
        isoval = st.slider("等値面 (isoval)", 0.001, 0.1, min(suggested_isoval, 0.05), step=0.001)
        viewer.addVolumetricData(cube_data, "cube", {"isoval": isoval, "color": color_pos, "opacity": opacity})
        if internal_mode == "orbital":
            viewer.addVolumetricData(cube_data, "cube", {"isoval": -isoval, "color": color_neg, "opacity": opacity})
            n_orb = st.session_state["mo_coeff"].shape[1]
            st.info(f"HOMO: {st.session_state['homo_idx']}, LUMO: {st.session_state['lumo_idx']}, 軌道数: {n_orb}")

            # 現在選択中の軌道情報を表示
            current_occ = st.session_state["mf"].mo_occ[orbital_idx]
            st.markdown(
                f"""
                **現在の軌道番号:** {orbital_idx}  
                **占有数 (occupation):** {current_occ:.2f}  
                **エネルギー (a.u.):** {st.session_state['mf'].mo_energy[orbital_idx]:.6f}  
                **エネルギー (eV):** {(st.session_state['mf'].mo_energy[orbital_idx] * 27.211386245988):.3f}
                """
            )

    # --- 原子番号ラベル表示の有無を選択 ---
    show_atom_numbers = st.checkbox("原子番号（インデックス）を3D表示する", value=False)
    # --- 原子番号ラベルを追加（オプション） ---
    coords = mol.atom_coords() * 0.529177  # xyzと同じ単位
    if show_atom_numbers:
        for i in range(mol.natm):
            x, y, z = coords[i]
            viewer.addLabel(str(i), {'position': {'x': float(x), 'y': float(y), 'z': float(z)},
                                    'backgroundColor': 'white', 'fontColor': 'black', 'fontSize': 16})

    viewer.zoomTo()
    stmol.showmol(viewer, height=500)

    if st.button("💾 表示中のcubeを保存"):
        filename = selected_chk.replace(".chk", f"_{internal_mode}.cube")
        save_cube_from_data(filename, cube_data)


else:
    st.warning("🔍 .chk ファイルが見つかりません。")