# QSAR Molecule Scoring

This repository contains a set of Python scripts for filtering, analyzing, and scoring bioactive molecules using SMILES data and IC50 values. The workflow is based on cheminformatics tools like RDKit and pandas.

## Features

- Filters molecular data based on activity values and units
- Calculates molecular descriptors (TPSA, MW, LogP, etc.)
- Simulates basic pharmacokinetic properties (MDCK, BBB, hERG)
- Scores molecules based on custom criteria and QED
- Visualizes top candidate structures

## File Structure

```

QSAR-Molecule-Scoring/
├── data/
│   └── archivo\_filtrado.csv
├── scripts/
│   ├── 01\_filtrado\_datos.py
│   ├── 02\_visualizar\_top3\_moleculas.py
│   ├── 03\_top3\_con\_propiedades\_simuladas.py
│   ├── 04\_mejor\_molecula\_con\_propiedades\_simuladas.py
│   ├── 05\_evaluacion\_completa\_con\_scores\_simulados.py
│   ├── 06\_evaluacion\_basica\_con\_score.py
│   └── 07\_evaluacion\_con\_qed.py

````

## Requirements

Install dependencies with:

```bash
pip install -r requirements.txt
````

## Usage

1. Place your molecular dataset in the `data/` folder.
2. Run the scripts in sequence or independently to filter and evaluate molecules.
3. Output includes visualizations and scoring information.

## License

MIT License

```

---
