"""
複数計算についてのまとめについて行う。
"""

import streamlit as st
import pandas as pd
import os
import json
import glob
import re

from utils.module import load_css

from logic.output_handler import (
    parse_folder_and_filename,
    parse_filename,
    extract_homo_lumo_scf_from_out,
    convert_energy_units
)

# dataフォルダ内のフォルダ名を取得して表示
data_path = "data"  # dataフォルダのパス

if os.path.exists(data_path):
    folder_names = [f for f in os.listdir(data_path) if os.path.isdir(os.path.join(data_path, f))]

    # 各フォルダ内のcalculation_log.jsonと画像を読み取る
    folder_smiles_list = []
    for folder in folder_names:
        json_path = os.path.join(data_path, folder, "calculation_log.json")
        image_path = os.path.join(data_path, folder, "molecule_2d.png")  # 画像ファイル名を指定

        if os.path.exists(json_path):
            with open(json_path, "r") as f:
                data = json.load(f)
                smiles = data[0]["smiles"]
                # 画像が存在する場合はパスを追加、存在しない場合は空文字列
                image = image_path if os.path.exists(image_path) else ""
                folder_smiles_list.append({
                    "Folder": folder,
                    "SMILES": smiles,
                    "Image": image,
                    "JsonPath": json_path
                })
        else:
            st.write(f"'{json_path}' が存在しません。")

    # フォルダ名、SMILES、画像をデータフレームとして表示
    if folder_smiles_list:
        # 画像を表示するかどうかを選択
        show_details = st.checkbox("Show Details", value=False)

        df = pd.DataFrame(folder_smiles_list)
        
        if show_details:
            # データフレームを表示
            st.write("Folder, SMILES, Structure:")
            for _, row in df.iterrows():
                with st.expander(f"Folder: {row['Folder']}"):
                    st.write(f"SMILES: {row['SMILES']}")
                    st.image(row["Image"], use_container_width=False, width=200)
                    
                    # JSONファイルの内容をデータフレームとして表示
                    with open(row["JsonPath"], "r") as f:
                        json_data = json.load(f)
                        # "time": "End_Time" のエントリのみフィルタリング
                        filtered_data = [entry for entry in json_data if entry.get("time") == "Normal_End_Time" or entry.get("time") == "Error_End_Time"]
                        if filtered_data:
                            # parameters, result, file_info を展開
                            expanded_data = []
                            for entry in filtered_data:
                                base_data = {key: entry[key] for key in entry if key not in ["parameters", "result", "file_info"]}
                                parameters = entry.get("parameters", {})
                                result = entry.get("result", {})
                                file_info = entry.get("file_info", {})
                                
                                # 各項目を展開して1つの辞書にまとめる
                                expanded_entry = {**base_data, **parameters, **result, **file_info}
                                expanded_data.append(expanded_entry)
                            
                            # 展開したデータをデータフレームに変換
                            json_df = pd.DataFrame(expanded_data)
                            st.dataframe(json_df)
                        else:
                            st.write("該当するデータがありません。")
        else:
            st.dataframe(df[["Folder", "SMILES"]], use_container_width=True)
    else:
        st.write("SMILES データが見つかりませんでした。")
else:
    st.write(f"'{data_path}' フォルダが存在しません。")

if st.checkbox("Show OUT File Information with HOMO, LUMO, and SCF Energies", value=False):
    out_files = glob.glob(os.path.join(data_path, "**", "*.out"), recursive=True)
    combined_data = []

    for out_file in out_files:
        # ファイルパスからInChIKeyとファイル名を抽出
        inchikey, filename = parse_folder_and_filename(out_file)

        # ファイル名を解析して情報を抽出（6要素に変更）
        theory, basis_set, opt_theory, opt_basis_set, index, extension = parse_filename(filename)

        # .outファイルからHOMO、LUMO、SCFエネルギーを取得
        out_homo_energy, out_lumo_energy, out_scf_energy = extract_homo_lumo_scf_from_out(out_file)

        # エネルギーをeVに変換（HOMOとLUMOのみ）
        out_homo_energy_ev = convert_energy_units(out_homo_energy, unit="eV") if out_homo_energy is not None else None
        out_lumo_energy_ev = convert_energy_units(out_lumo_energy, unit="eV") if out_lumo_energy is not None else None

        # SCFエネルギーをkcal/molに変換
        out_scf_energy_kcal = convert_energy_units(out_scf_energy, unit="kcal/mol") if out_scf_energy is not None else None

        # データを統合
        combined_data.append({
            "InChIKey": inchikey,
            "Theory": theory,
            "Basis Set": basis_set,
            "Opt Theory": opt_theory,
            "Opt Basis Set": opt_basis_set,
            "Index": index,
            "HOMO Energy (eV)": out_homo_energy_ev,
            "LUMO Energy (eV)": out_lumo_energy_ev,
            "SCF Energy (Hartree)": out_scf_energy,
            "SCF Energy (kcal/mol)": out_scf_energy_kcal
        })

    if combined_data:
        combined_df = pd.DataFrame(combined_data)

        # グループ化して差分を計算
        combined_df["Energy Difference (kcal/mol)"] = combined_df.groupby(
            ["InChIKey", "Theory", "Basis Set", "Opt Theory", "Opt Basis Set"]
        )["SCF Energy (kcal/mol)"].transform(lambda x: x - x.min())

        st.subheader("OUT File Information with HOMO, LUMO, and SCF Energies")
        st.dataframe(combined_df)
    else:
        st.write("No OUT files found.")


