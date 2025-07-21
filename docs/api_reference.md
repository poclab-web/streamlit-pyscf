# API Reference

## Logic Modules

### Calculation Module

```{eval-rst}
.. automodule:: logic.calculation
   :members:
   :undoc-members:
   :show-inheritance:
```

主要な量子化学計算を実行します：
- 一点エネルギー計算
- 構造最適化
- 振動解析
- 励起状態計算（TD-DFT）
- NMR計算

### Molecule Handler

```{eval-rst}
.. automodule:: logic.molecule_handler
   :members:
   :undoc-members:
   :show-inheritance:
```

分子構造の入力・変換・出力を管理します：
- SMILES文字列の解析
- XYZ座標の読み込み
- 配座探索
- 2D/3D可視化

### Database Module

```{eval-rst}
.. automodule:: logic.database
   :members:
   :undoc-members:
   :show-inheritance:
```

計算結果のデータベース管理機能：
- 分子データの保存・検索
- エネルギー・振動解析結果の管理
- 励起状態データの管理

### MOPAC Calculator

```{eval-rst}
.. automodule:: logic.mopac_calculation
   :members:
   :undoc-members:
   :show-inheritance:
```

MOPAC計算エンジンとの連携機能：
- 半経験的量子化学計算
- 分子軌道可視化
- 構造最適化

## Controllers

### Energy Decomposition Analysis

```{eval-rst}
.. automodule:: controllers.energydecompositionanalysis
   :members:
   :undoc-members:
   :show-inheritance:
```

### Excited State Calculation

```{eval-rst}
.. automodule:: controllers.excited_state_calculation
   :members:
   :undoc-members:
   :show-inheritance:
```

### Fragment Calculation

```{eval-rst}
.. automodule:: controllers.fragment_calculation
   :members:
   :undoc-members:
   :show-inheritance:
```

## Configuration

### User Preferences

```{eval-rst}
.. automodule:: config.user_preferences
   :members:
   :undoc-members:
   :show-inheritance:
```

### External Software Configuration

```{eval-rst}
.. automodule:: config.external_software_config
   :members:
   :undoc-members:
   :show-inheritance:
```
