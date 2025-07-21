# Computational Chemistry Tool (streamlit-pyscf)

![Python](https://img.shields.io/badge/python-3.7+-blue.svg)
![PySCF](https://img.shields.io/badge/PySCF-2.x-green.svg)
![Gaussian](https://img.shields.io/badge/Gaussian-16/09-blue.svg)
![OpenMM](https://img.shields.io/badge/OpenMM-8.x-orange.svg)
![xTB](https://img.shields.io/badge/xTB-6.x-purple.svg)
![MOPAC](https://img.shields.io/badge/MOPAC-2016-red.svg)
![Streamlit](https://img.shields.io/badge/Streamlit-1.x-red.svg)

**A comprehensive computational chemistry platform integrating diverse calculation engines**

This project integrates various computational chemistry programs including PySCF, Gaussian, OpenMM, MOPAC, and xTB, providing a wide range of calculation methods from molecular force fields to quantum chemistry calculations through an intuitive web interface.

## 🌟 Features

- **Integrated Platform**: Access multiple computational chemistry programs through a single interface
- **No GUI Required**: Intuitive interface running in web browsers
- **Multi-tier Calculation Methods**: Hierarchical approach from molecular force fields → semi-empirical → quantum chemistry
- **Comprehensive Calculation Functions**: Structure optimization, vibrational analysis, excited states, NMR prediction, molecular dynamics, and more
- **Visualization Support**: 3D display of molecular structures and calculation results
- **Database Management**: Automatic saving and search functionality for calculation results
- **No Programming Required**: Usable with chemical knowledge alone
- **Customizable**: Page display settings and user preferences

## 🚀 Demo

Online version available without installation (local environment recommended):

**[https://conputational-chemistry-tool.streamlit.app/](https://conputational-chemistry-tool.streamlit.app/)**

## 📦 Feature List

### 🔧 Preprocessing & Utilities
- **Structure Check**: Molecular structure validation and verification
- **Format Conversion**: Conversion between various molecular file formats

### 🔀 Integrated Workflows
- **General Calculation**: Automated workflows integrating multiple calculations

### 🧲 Molecular Force Field Calculations
- **Conformational Search**: Conformational analysis using RDKit
- **OpenMM Molecular Dynamics**: High-speed molecular dynamics simulations
  - Classical force field calculations (AMBER, CHARMM, OPLS, etc.)
  - Solvation effect simulations
  - Energy minimization

### 🧬 Semi-empirical Methods
- **MOPAC Calculations**: High-speed calculations using PM6/PM7
  - Single point energy calculations
  - Structure optimization
  - Vibrational analysis
  - Result visualization
- **xTB Calculations**: Extended tight-binding calculations
  - GFN1-xTB, GFN2-xTB
  - Solvation effects (GBSA, ALPB)
  - Conformational analysis

### 🧪 Quantum Chemistry Calculations
- **PySCF Calculations**: Open-source quantum chemistry package
  - Single point energy: Energy, orbital energies, dipole moment
  - Structure optimization: Molecular geometry optimization
  - Vibrational analysis: Frequencies, thermodynamic quantities, IR/Raman spectra
- **Gaussian Calculations**: Commercial quantum chemistry package (license required)
  - High-precision DFT calculations
  - Correlation energy methods (MP2, CCSD(T), etc.)
  - Large molecular system calculations

### 🔍 Visualization and Analysis
- **Molecular Orbital Visualization**: 3D display of molecular structures and orbitals
- **Energy Decomposition Analysis (EDA)**: Detailed analysis of intermolecular interactions
- **Fragment Decomposition**: Fragment analysis from xyz coordinates
- **Conformational Energy Decomposition**: Energy analysis of conformational isomers

### ⚡ Property Calculations
- **Ionization Potential**: Theoretical calculation of electron affinity
- **Solvation Effects**: Solvation calculations using PCM/ddCOSMO models
- **Bond Dissociation Energy**: Evaluation of chemical bond strength
- **pKa Calculation**: Theoretical prediction of proton affinity

### 📊 Spectroscopic Calculations
- **IR/Raman Spectra**: Prediction of vibrational spectra
- **NMR Prediction**: Theoretical calculation of chemical shifts
- **Polarizability Calculation**: Molecular polarizability
- **UV-Vis Prediction**: Time-dependent DFT calculation of excited states

### 🔄 Transition State Calculations
- **NEB Calculation**: Exploration of complex reaction pathways
- **Transition State Search**: Activation barriers for chemical reactions
- **IRC Calculation**: Reaction path tracking

### ⚙️ System & Settings
- **Database**: Management and search of calculation results
- **Result Summarization**: Statistics and visualization of calculation data
- **System Settings**: Application configuration
- **Page Settings**: Customization of displayed pages

## 📋 Requirements

### Essential Requirements
- Python 3.7 or higher
- Recommended Memory: 8GB or more
- Recommended CPU: Multi-core (for faster calculations)

### Calculation Engine Requirements
- **PySCF**: Automatically installed in Python environment
- **OpenMM**: Automatically installed in Python environment
- **xTB**: Executable file path configuration required
- **MOPAC**: Executable file path configuration required
- **Gaussian**: Commercial license and execution environment setup required (optional)

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
- `pyscf`: Quantum chemistry calculation engine (open source)
- `streamlit`: Web application framework
- `rdkit`: Cheminformatics library
- `matplotlib`, `plotly`: Data visualization
- `py3Dmol`, `stmol`: 3D molecular display
- `openmm`: Molecular dynamics simulation
- `geometric`, `pyberny`: Structure optimization
- `ase`: Atomic Simulation Environment

### 4. External Program Configuration (Optional)
Install and configure the following programs as needed:

```bash
# xTB (semi-empirical calculations)
# Download and install from https://github.com/grimme-lab/xtb

# MOPAC (semi-empirical calculations)  
# Download and install from http://openmopac.net/

# Gaussian (quantum chemistry calculations, license required)
# Install from official site after license purchase
```

Configure program paths in `config/external_software_config.py`.

##  Usage

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

### Supported Calculation Methods and Theory Levels

#### Molecular Force Fields
- **OpenMM Force Fields**: AMBER, CHARMM, OPLS-AA, GAFF, etc.

#### Semi-empirical Methods
- **MOPAC**: PM6, PM7, PM3, AM1, MNDO
- **xTB**: GFN1-xTB, GFN2-xTB, GFN-FF

#### Quantum Chemistry Calculations
- **Hartree-Fock (HF)**
- **Density Functional Theory (DFT)**: B3LYP, PBE, M06-2X, ωB97X-D, etc.
- **Correlation Theory**: MP2, CCSD, CCSD(T) (available in Gaussian)
- **Dispersion Correction**: D3, D4 correction

### Supported Basis Sets
- STO-3G, 3-21G, 6-31G(d,p)
- cc-pVDZ, cc-pVTZ, aug-cc-pVDZ
- def2-SVP, def2-TZVP
## 📁 Project Structure

```
streamlit-pyscf/
├── ComputationalChemistryTool.py   # Main application
├── pages/                          # Feature pages
│   ├── 001_StructureCheck.py       # Structure check
│   ├── 101_GeneralCalculation.py   # General calculations
│   ├── 201_RDKit_ConfSearch.py     # RDKit conformational search
│   ├── 211_openmm_singlepoint.py   # OpenMM calculations
│   ├── 301-304_mopac_*.py          # MOPAC calculation suite
│   ├── 311_xtb_calculation_suite.py# xTB calculations
│   ├── 401-403_*.py                # Quantum chemistry calculations
│   ├── 501-504_*.py                # Visualization & analysis
│   ├── 601-604_*.py                # Property calculations
│   ├── 701-704_*.py                # Spectroscopic calculations
│   ├── 801-803_*.py                # Transition state calculations
│   └── 901-904_*.py                # System settings
├── logic/                          # Calculation logic
│   ├── calculation.py              # Quantum chemistry calculations
│   ├── molecule_handler.py         # Molecular data processing
│   ├── visualization.py            # Result visualization
│   └── database.py                 # Database management
├── config/                         # Configuration files
│   ├── page_visibility.json        # Page display settings
│   ├── user_preferences.py         # User preferences
│   └── external_software_config.py # External software configuration
├── data/                          # Calculation data storage
├── controllers/                   # Calculation controllers
├── utils/                         # Utilities
└── tests/                         # Test code
```

## 💡 Usage Examples

### Example 1: Hierarchical Calculation Approach
1. **Initial structure optimization with molecular force fields**
   - Generate initial conformations with "Molecular Force Field > RDKit Conformational Search"
   - Perform rough optimization with "Molecular Force Field > OpenMM Calculation"
2. **Refinement with semi-empirical methods**
   - Optimize with medium accuracy using "Semi-empirical > xTB Calculation"
3. **High-precision calculation with quantum chemistry**
   - Final high-precision calculation with "Quantum Chemistry > PySCF Structure Optimization"

### Example 2: Water Molecule Structure Optimization (PySCF)
1. Select "Quantum Chemistry > Structure Optimization" page
2. SMILES input: `O`
3. Select Theory Level: `B3LYP`, Basis Set: `6-31G(d,p)`
4. Click "Run Optimization"

### Example 3: Benzene Molecular Dynamics Simulation (OpenMM)
1. Select "Molecular Force Field > OpenMM Calculation" page
2. SMILES input: `c1ccccc1`
3. Set Force Field: `GAFF`, Temperature: `300K`
4. Click "Run MD"

### Example 4: Organic Molecule Semi-empirical Calculation (xTB)
1. Select "Semi-empirical > xTB Calculation" page
2. SMILES input: `CCO` (ethanol)
3. Select Method: `GFN2-xTB`, Solvent: `Water`
4. Click "Run Calculation"

### Example 5: Intermolecular Interaction EDA Analysis
1. Select "Visualization & Analysis > Energy Decomposition Analysis" page
2. Input complex molecule
3. Configure fragment division
4. Run EDA analysis

## 🔬 Calculation Examples and Benchmarks

### Estimated Calculation Time (protein amino acid residue-sized molecules)

| Calculation Method | Theory Level | Basis Set/Settings | Single Point | Structure Optimization |
|-------------------|--------------|-------------------|--------------|----------------------|
| **Molecular Force Field** | OpenMM | GAFF | ~0.1 sec | ~1 sec |
| **Semi-empirical** | xTB | GFN2-xTB | ~0.5 sec | ~5 sec |
| **Semi-empirical** | MOPAC | PM7 | ~1 sec | ~10 sec |
| **Quantum Chemistry** | PySCF | B3LYP/6-31G(d,p) | ~10 sec | ~5 min |
| **Quantum Chemistry** | Gaussian | B3LYP/6-31G(d,p) | ~5 sec | ~3 min |

### Larger Molecular Systems (Benzene = 12 atoms)

| Calculation Method | Theory Level | Settings | Single Point | Structure Optimization |
|-------------------|--------------|----------|--------------|----------------------|
| **Molecular Force Field** | OpenMM | GAFF | ~0.1 sec | ~2 sec |
| **Semi-empirical** | xTB | GFN2-xTB | ~1 sec | ~10 sec |
| **Semi-empirical** | MOPAC | PM7 | ~2 sec | ~30 sec |
| **Quantum Chemistry** | PySCF | B3LYP/6-31G(d,p) | ~30 sec | ~15 min |

*Calculation times vary significantly depending on the environment

### Accuracy vs. Cost Trade-off
| Calculation Method | Computational Cost | Structural Accuracy | Spectral Accuracy | Application Range |
|-------------------|-------------------|-------------------|------------------|------------------|
| **Molecular Force Field** | Very Low | Medium | × | Large systems, MD |
| **Semi-empirical** | Low | Medium-High | Medium | Medium systems, screening |
| **DFT** | Medium | High | High | General chemical systems |
| **ab initio** | High | Very High | Very High | High-precision systems |

## 🐛 Troubleshooting

### Common Issues

**Q: Calculation fails**
- Check if molecular structure is valid
- Try faster calculation methods first (molecular force field → semi-empirical → quantum chemistry)
- Try smaller basis sets (for quantum chemistry calculations)
- Verify charge and spin multiplicity settings

**Q: Memory shortage error**
- Use faster calculation methods (OpenMM, xTB, MOPAC)
- Use smaller basis sets (for quantum chemistry calculations)
- Reduce molecule size
- Increase available memory

**Q: External program not found**
- Check path configuration in `config/external_software_config.py`
- Verify that the program is correctly installed
- Check environment variable settings

**Q: Convergence issues**
- Optimize initial structure with coarser calculation methods
- Change initial structure
- Relax convergence criteria
- Try different theory levels

**Q: Gaussian license error**
- Verify that you have a valid license
- Check Gaussian environment configuration
- Consider alternative calculations with PySCF

## 🎨 Page Customization

The application allows customization of displayed pages:

1. Access "System & Settings > Page Settings"
2. Configure page visibility by category
3. Different settings can be configured for each user

Available Categories:
- 🔧 Preprocessing & Utilities
- 🔀 Integrated Workflows  
- 🧲 Molecular Force Fields
- 🧬 Semi-empirical Methods
- 🧪 Quantum Chemistry Calculations
- 🔍 Visualization & Analysis
- ⚡ Property Calculations
- 📊 Spectroscopic Calculations
- 🔄 Transition State Calculations
- ⚙️ System & Settings

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
@software{computational_chemistry_tool,
  title = {Computational Chemistry Tool: Integrated Multi-Engine Platform},
  author = {PocLab},
  url = {https://github.com/poclab-web/streamlit-pyscf},
  year = {2024}
}
```

### Citing Dependent Libraries and Programs
- **PySCF**: Sun, Q. et al. *WIREs Comput Mol Sci* **2018**, *8*, e1340.
- **OpenMM**: Eastman, P. et al. *PLoS Comput Biol* **2017**, *13*, e1005659.
- **xTB**: Bannwarth, C. et al. *J Chem Theory Comput* **2019**, *15*, 1652-1671.
- **MOPAC**: Stewart, J. J. P. *J Mol Model* **2013**, *19*, 1-32.
- **Gaussian**: Frisch, M. J. et al. Gaussian 16, Revision C.01, Gaussian, Inc., Wallingford CT, 2016.
- **RDKit**: Landrum, G. RDKit: Open-source cheminformatics.

## 🛡️ Security

If you discover security issues, please email the developers directly rather than posting public issues.

## 👥 Development Team

- **PocLab** - Project Lead
- See [Contributors](https://github.com/poclab-web/streamlit-pyscf/contributors) for the full list of contributors

## 🔗 Related Links

### Integrated Programs
- [PySCF Official Site](https://pyscf.org/)
- [OpenMM Official Site](http://openmm.org/)
- [xTB Official Site](https://xtb-docs.readthedocs.io/)
- [MOPAC Official Site](http://openmopac.net/)
- [Gaussian Official Site](https://gaussian.com/)

### Frameworks & Libraries
- [Streamlit Official Site](https://streamlit.io/)
- [RDKit Official Site](https://www.rdkit.org/)
- [ASE Official Site](https://wiki.fysik.dtu.dk/ase/)

### Learning Resources
- [Quantum Chemistry Calculation Best Practices](https://example.com/qc-best-practices)
- [Introduction to Molecular Dynamics](https://example.com/md-basics)
- [Computational Chemistry Method Comparison](https://example.com/computational-methods)

## 📊 Statistics

![GitHub stars](https://img.shields.io/github/stars/poclab-web/streamlit-pyscf?style=social)
![GitHub forks](https://img.shields.io/github/forks/poclab-web/streamlit-pyscf?style=social)
![GitHub issues](https://img.shields.io/github/issues/poclab-web/streamlit-pyscf)
![GitHub last commit](https://img.shields.io/github/last-commit/poclab-web/streamlit-pyscf)

---

**🔬 A comprehensive computational chemistry platform from molecular force fields to quantum chemistry.**