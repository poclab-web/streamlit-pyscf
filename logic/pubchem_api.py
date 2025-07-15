"""
pubchemapi
"""

import requests

class PubChemSpectraClient:
    def __init__(self, inchikey):
        self.inchikey = inchikey
        self.cid = self._get_cid_from_inchikey()
        self.view_data = self._get_view_data() if self.cid else None

    def _get_cid_from_inchikey(self):
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{self.inchikey}/cids/JSON"
        res = requests.get(url)
        if res.status_code != 200:
            print("CID取得失敗")
            return None
        return res.json().get("IdentifierList", {}).get("CID", [None])[0]

    def _get_view_data(self):
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{self.cid}/JSON"
        res = requests.get(url)
        if res.status_code != 200:
            print("pug_view取得失敗")
            return None
        return res.json()

    def _find_section(self, heading):
        """指定されたTOCHeading（例: 'UV/Visible Spectrum'）を検索"""
        def search(sections):
            for section in sections:
                if section.get("TOCHeading", "").lower() == heading.lower():
                    return section
                if "Section" in section:
                    result = search(section["Section"])
                    if result:
                        return result
            return None

        if self.view_data:
            return search(self.view_data.get("Record", {}).get("Section", []))
        return None

    def get_uv_spectrum(self):
        return self._find_section("UV/Visible Spectrum")

    def get_ir_spectrum(self):
        return self._find_section("Infrared Spectra")

    def get_nmr_spectrum(self):
        return self._find_section("1D NMR Spectra")


if __name__ == "__main__":
    inchikey = "CZDYPVPMEAXLPK-UHFFFAOYSA-N"  # 例: アセトフェノンのInChIKey
    client = PubChemSpectraClient(inchikey)

    uv_spectrum = client.get_uv_spectrum()

    if uv_spectrum:
        print("UV/Visible Spectrum:")
        print(uv_spectrum)
    else:
        print("UV/Visible Spectrum not found.")

    ir_spectrum = client.get_ir_spectrum()
    if ir_spectrum:
        print("Infrared Spectrum:")
        print(ir_spectrum)
    else:
        print("Infrared Spectrum not found.")

    nmr_spectrum = client.get_nmr_spectrum()
    if nmr_spectrum:
        print("1D NMR Spectrum:")
        print(nmr_spectrum)
    else:
        print("1D NMR Spectrum not found.")