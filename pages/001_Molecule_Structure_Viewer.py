import streamlit as st
from streamlit_ketcher import st_ketcher
from rdkit import Chem
from rdkit.Chem import Draw
import py3Dmol
import stmol
from logic.molecule_handler import MoleculeHandler
from utils.module import load_css

# カスタムCSSを適用
load_css("config/styles.css")

# ページタイトル
st.title("🧪 分子構造描画・表示ツール")
st.markdown("構造式エディタで分子を描画し、2D/3D構造や各種分子情報を確認できます。")

# 分子構造描画エディタ
st.subheader("📝 分子構造描画")
smiles = st_ketcher()

# ユーザーがまだ構造を入力していない場合
if not smiles:
    st.warning("構造が描画されていません。構造を描画してApplyをクリックしてください。")
    st.stop()

# SMILESの表示
st.subheader("📊 分子情報")
st.write("**入力されたSMILES:**")
st.code(smiles)

try:
    # MoleculeHandlerを使用して分子を処理
    handler = MoleculeHandler(smiles, input_type="smiles")
    mol = handler.mol
    
    if mol is None:
        st.error("SMILESの解析に失敗しました。正しいSMILES記法で入力してください。")
        st.stop()
    
    # 分子の基本情報を表示
    col_info1, col_info2, col_info3 = st.columns(3)
    
    with col_info1:
        st.metric("原子数", mol.GetNumAtoms())
    
    with col_info2:
        st.metric("結合数", mol.GetNumBonds())
    
    with col_info3:
        try:
            from rdkit.Chem import rdMolDescriptors
            molecular_weight = rdMolDescriptors.CalcExactMolWt(mol)
            st.metric("分子量", f"{molecular_weight:.2f}")
        except:
            st.metric("分子量", "N/A")

    # InChI情報
    st.write("**InChI:**")
    st.code(Chem.MolToInchi(mol))

    st.write("**InChIKey:**")
    st.code(Chem.MolToInchiKey(mol))

    # 2D/3D構造表示
    st.subheader("🔍 分子構造表示")
    col1, col2 = st.columns(2)

    # Display 2D structure in the first column
    with col1:
        st.markdown("### 2D構造")
        try:
            img = Draw.MolToImage(mol, size=(400, 400))
            st.image(img, caption="2D分子構造", use_container_width=True)
        except Exception as e:
            st.error(f"2D構造の生成に失敗しました: {e}")

    # Display 3D structure in the second column
    with col2:
        st.markdown("### 3D構造")
        try:
            # 3D構造をMOLブロック形式で取得
            mol_block = handler.generate_3d_molblock()
            if mol_block:
                viewer = py3Dmol.view(width=400, height=400)
                viewer.addModel(mol_block, "mol")
                viewer.setStyle({"stick": {}})
                viewer.zoomTo()
                stmol.showmol(viewer, height=400)
            else:
                st.warning("3D構造を生成できませんでした")
        except Exception as e:
            st.error(f"3D構造の表示に失敗しました: {e}")
    
    # pyscf input
    st.subheader("PySCF Input")
    if handler:
        try:
            pyscf_input = handler.to_pyscf_input()
            st.code(pyscf_input, language="text")
        except Exception as e:
            st.error(f"PySCF入力の生成に失敗しました: {e}")

    # mopac input
    st.subheader("MOPAC Input")
    if handler:
        try:
            mopac_input = handler.to_mopac_input()
            st.code(mopac_input, language="text")
        except Exception as e:
            st.error(f"MOPAC入力の生成に失敗しました: {e}")

    # SDF情報
    st.subheader("📄 SDF(Structure-Data File)")
    try:
        mol_block = handler.generate_3d_molblock()
        st.code(mol_block, language="text")
        
        # ダウンロードボタン
        st.download_button(
            label="📥 SDFファイルをダウンロード",
            data=mol_block,
            file_name=f"{Chem.MolToInchiKey(mol)}.sdf",
            mime="chemical/x-mdl-sdfile"
        )
    except Exception as e:
        st.error(f"SDF生成に失敗しました: {e}")

except Exception as e:
    st.error(f"分子構造の処理中にエラーが発生しました: {e}")
    st.info("正しいSMILES記法で分子構造を入力してください。例: CCO (エタノール), c1ccccc1 (ベンゼン)")