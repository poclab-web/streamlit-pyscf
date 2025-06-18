
# 量子化学計算の構造最適化のときの設定値を定義

conv_preset_values = {
    "Loose":  {"energy": 1.0e-4, "grms": 1.0e-3, "gmax": 3.0e-3, "drms": 4.0e-3, "dmax": 6.0e-3},
    "Normal": {"energy": 1.0e-5, "grms": 5.0e-4, "gmax": 1.5e-3, "drms": 2.0e-3, "dmax": 3.0e-3},
    "Tight":  {"energy": 1.0e-6, "grms": 3.0e-4, "gmax": 1.2e-3, "drms": 1.2e-3, "dmax": 1.8e-3},
}