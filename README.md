# streamlit-pyscf

![Python](https://img.shields.io/badge/python-3.7+-blue.svg)
![PySCF](https://img.shields.io/badge/PySCF-2.x-green.svg)
![Streamlit](https://img.shields.io/badge/Streamlit-1.x-red.svg)

**A web application that makes advanced quantum chemistry calculations accessible to everyone using Streamlit and PySCF.**

PySCF (Python-based Simulations of Chemistry Framework) is a library for performing advanced quantum chemistry calculations. This project provides an intuitive web interface to utilize PySCF's powerful capabilities.

## 🌟 Features

- **No GUI Required**: Intuitive interface running in web browsers
- **Comprehensive Calculation Functions**: Structure optimization, vibrational analysis, excited states, NMR prediction, and more
- **Visualization Support**: 3D display of molecular structures and calculation results
- **Database Management**: Automatic saving and search functionality for calculation results
- **No Programming Required**: Usable with chemical knowledge alone

## 🚀 Demo

Online version available without installation (local environment recommended):

**[https://conputational-chemistry-tool.streamlit.app/](https://conputational-chemistry-tool.streamlit.app/)**

## 📦 Feature List

### Basic Calculations
- **Single Point Calculation**: Energy, orbital energies, dipole moment
- **Structure Optimization**: Molecular geometry optimization
- **Vibrational Analysis**: Frequencies, thermodynamic quantities, IR/Raman spectra
- **Conformational Search**: Structure exploration using molecular force fields

### Spectroscopic Properties
- **NMR Prediction**: Theoretical calculation of chemical shifts
- **UV-Vis Prediction**: Time-dependent DFT calculation of excited states
- **IR/Raman Spectra**: Prediction of vibrational spectra

### Advanced Analysis
- **pKa Calculation**: Theoretical prediction of proton affinity
- **Bond Dissociation Energy**: Evaluation of chemical bond strength
- **Energy Decomposition Analysis**: Detailed analysis of intermolecular interactions
- **Solvent Effects**: Solvation calculations using PCM/ddCOSMO models

### Reaction Analysis
- **Transition State Search**: Activation barriers for chemical reactions
- **IRC Calculation**: Reaction path tracking
- **NEB Calculation**: Exploration of complex reaction pathways

## 📋 Requirements

- Python 3.7 or higher
- Recommended Memory: 8GB or more
- Recommended CPU: Multi-core (for faster calculations)

## 🛠️ Installation

### 1. Clone the Repository
```bash
git clone https://github.com/poclab-web/streamlit-pyscf.git
cd streamlit-pyscf
```

### 2. Create Virtual Environment (Recommended)
```bash
# Using Conda
conda create -n streamlit-pyscf python=3.11
conda activate streamlit-pyscf

# Using venv
python -m venv streamlit-pyscf
source streamlit-pyscf/bin/activate  # Linux/Mac
# streamlit-pyscf\Scripts\activate  # Windows
```

### 3. Install Dependencies
```bash
pip install -r requirements.txt
```

**Key Dependencies:**
- `pyscf`: Quantum chemistry calculation engine
- `streamlit`: Web application framework
- `rdkit`: Cheminformatics library
- `matplotlib`, `plotly`: Data visualization
- `py3Dmol`, `stmol`: 3D molecular display

## 🚀 Usage

### Basic Launch
```bash
streamlit run ComputationalChemistryTool.py
```

The browser will automatically open and you can access the application at `http://localhost:8501`.

### Basic Usage Steps
1. **Molecular Input**: Enter molecule using SMILES notation or XYZ coordinates
2. **Calculation Settings**: Select theory level (HF, DFT, etc.) and basis set
3. **Run Calculation**: Click the calculation button to execute
4. **View Results**: Check and download visualized results

### Supported Theory Levels
- **Hartree-Fock (HF)**
- **Density Functional Theory (DFT)**: B3LYP, PBE, M06-2X, etc.
- **Correlation Theory**: MP2
- **Dispersion Correction**: D3 correction

### Supported Basis Sets
- STO-3G, 3-21G, 6-31G(d,p)
- cc-pVDZ, cc-pVTZ, aug-cc-pVDZ
- def2-SVP, def2-TZVP

## 📁 Project Structure

```
streamlit-pyscf/
├── ComputationalChemistryTool.py   # Main application
├── pages/                          # Feature pages
│   ├── 01_GeneralCalculation.py    # General calculations
│   ├── 02_StructureCheck.py        # Structure validation
│   ├── 03_Singlepointcalculation.py # Single point calculations
│   ├── 04_Optimization.py          # Structure optimization
│   ├── 06_OPTandFreq.py           # Optimization + vibrational analysis
│   ├── 11_Pka.py                  # pKa calculations
│   ├── 12_UV_SpectrumPrediction.py # UV-Vis spectrum
│   └── ...
├── logic/                          # Calculation logic
│   ├── calculation.py              # Quantum chemistry calculations
│   ├── molecule_handler.py         # Molecular data processing
│   ├── visualization.py            # Result visualization
│   └── database.py                 # Database management
├── config/                         # Configuration files
├── data/                          # Calculation data storage
└── utils/                         # Utilities
```

## 💡 Usage Examples

### Example 1: Water Molecule Structure Optimization
1. Select "04_Optimization" page
2. Enter SMILES: `O`
3. Select Theory Level: `B3LYP`, Basis Set: `6-31G(d,p)`
4. Click "Run Optimization"

### Example 2: Benzene UV Spectrum Prediction
1. Select "12_UV_SpectrumPrediction" page
2. Enter SMILES: `c1ccccc1`
3. Set Number of Excited States: `10`
4. Click "Run Calculation"

## 🔬 Calculation Examples and Benchmarks

### Estimated Calculation Times
| Molecule Size | Theory Level | Basis Set | Single Point | Structure Optimization |
|---------------|--------------|-----------|--------------|----------------------|
| H2O (3 atoms) | B3LYP | 6-31G(d,p) | ~1 sec | ~30 sec |
| Benzene (12 atoms) | B3LYP | 6-31G(d,p) | ~10 sec | ~5 min |
| Naphthalene (18 atoms) | B3LYP | 6-31G(d,p) | ~30 sec | ~15 min |

*Calculation times vary significantly depending on the environment

### Accuracy Information
- **Structure**: Error < 2% compared to experimental values (for many organic molecules)
- **Frequencies**: Error < 10% compared to experimental values (harmonic approximation)
- **NMR Chemical Shifts**: Useful for qualitative predictions
- **Excitation Energies**: Depends on TD-DFT accuracy

## 🐛 Troubleshooting

### Common Issues

**Q: Calculation fails**
- Check if molecular structure is valid
- Try smaller basis sets
- Verify charge and spin multiplicity settings

**Q: Memory shortage error**
- Use smaller basis sets
- Reduce molecule size
- Increase available memory

**Q: Convergence issues**
- Change initial structure
- Relax convergence criteria
- Try different theory levels

## 🤝 Contributing

We welcome contributions to the project!

### Development Environment Setup
```bash
# Install development dependencies
pip install -r requirements-dev.txt

# Install pre-commit hooks
pre-commit install
```

### Contribution Guidelines
1. Fork this repository
2. Create a feature branch: `git checkout -b feature/new-feature`
3. Commit your changes: `git commit -m "Add new feature"`
4. Push the branch: `git push origin feature/new-feature`
5. Create a Pull Request

### Bug Reports and Feature Requests
Please report on the Issues page:
**[https://github.com/poclab-web/streamlit-pyscf/issues](https://github.com/poclab-web/streamlit-pyscf/issues)**

## 📖 Documentation

### Detailed Documentation
- [Architecture Overview](docs/architecture.md)
- [API Reference](docs/api.md)
- [Developer Guide](docs/development.md)

### Learning Resources
- [PySCF Official Documentation](https://pyscf.org/)
- [Introduction to Quantum Chemistry Calculations](https://example.com/quantum-chemistry)
- [Fundamentals of DFT Theory](https://example.com/dft-basics)

## 🏆 Citation

If you use this software in your research, please cite:

```bibtex
@software{streamlit_pyscf,
  title = {streamlit-pyscf: Web-based Quantum Chemistry Interface},
  author = {PocLab},
  url = {https://github.com/poclab-web/streamlit-pyscf},
  year = {2024}
}
```

### Citing Dependent Libraries
- **PySCF**: Sun, Q. et al. *WIREs Comput Mol Sci* **2018**, *8*, e1340.
- **RDKit**: Landrum, G. RDKit: Open-source cheminformatics.

## 🛡️ Security

If you discover security issues, please email the developers directly rather than posting public issues.

## 👥 Development Team

- **PocLab** - Project Lead
- See [Contributors](https://github.com/poclab-web/streamlit-pyscf/contributors) for the full list of contributors

## 🔗 Related Links

- [PySCF Official Site](https://pyscf.org/)
- [Streamlit Official Site](https://streamlit.io/)
- [RDKit Official Site](https://www.rdkit.org/)
- [Quantum Chemistry Calculation Best Practices](https://example.com/qc-best-practices)

## 📊 Statistics

![GitHub stars](https://img.shields.io/github/stars/poclab-web/streamlit-pyscf?style=social)
![GitHub forks](https://img.shields.io/github/forks/poclab-web/streamlit-pyscf?style=social)
![GitHub issues](https://img.shields.io/github/issues/poclab-web/streamlit-pyscf)
![GitHub last commit](https://img.shields.io/github/last-commit/poclab-web/streamlit-pyscf)

---

**🔬 Making advanced quantum chemistry calculations more accessible.**