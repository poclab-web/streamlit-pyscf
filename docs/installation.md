# Installation

## Requirements

- Python 3.8+
- PySCF
- RDKit
- Streamlit

## Install from source

```bash
git clone https://github.com/poclab-web/streamlit-pyscf.git
cd streamlit-pyscf
pip install -r requirements.txt
```

## Running the application

```bash
streamlit run ComputationalChemistryTool.py
```

## Docker installation (optional)

```bash
docker build -t streamlit-pyscf .
docker run -p 8501:8501 streamlit-pyscf
```
