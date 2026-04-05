# Mol-Morpher

Molecule morphing using RDKit.

The main project workflow is the notebook:

- scripts/multi_beam.ipynb

This notebook is the core of the repository. You insert two SMILES strings in the notebook input section and it generates intermediate molecules between them.

The other scripts are test and utility helpers:

- scripts/smiles_to_json.py
- scripts/weight_sensitivity.py

## Project layout

```text
Mol-Morpher/
  data/
    input/
      Smiles.csv
    output/
      Data.jsonl
      explored_molecules.csv
      weight_search_full/
        run_config.json
        weight_leaderboard.csv
        best_weight.json
        <one folder per weight setting>/
          pair_metrics.csv
          summary.json
  scripts/
    multi_beam.ipynb
    smiles_to_json.py
    sascorer.py
    weight_sensitivity.py
    fpscores.pkl.gz
  docs/
    FlowChart.png
    ReadME.pdf
  requirements.txt
```

## Main workflow

### 1) Environment setup

Use Python 3.12+.

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### 2) Run the main notebook

Open scripts/multi_beam.ipynb in VS Code or Jupyter.

Insert your two SMILES in the notebook input section, then run the notebook cells.

Main output from this core flow is written under data/output, including explored_molecules.csv.

## Testing and utility scripts

These are not the main project deliverable. They are used for testing, diagnostics, and batch experiments.

### 1) scripts/smiles_to_json.py (testing utility)

- Reads data/input/Smiles.csv
- Canonicalizes/filters SMILES
- Builds source-target pairs
- Writes data/output/Data.jsonl

Run:

```powershell
python scripts/smiles_to_json.py
```

### 2) scripts/weight_sensitivity.py (testing utility)

Used to test how alpha/beta/gamma weight changes affect outputs and whether results differ across weights.

- Loads algorithm functions from scripts/multi_beam.ipynb
- Runs a weight grid
- Writes diagnostics to data/output/weight_search_full

Run:

```powershell
python scripts/weight_sensitivity.py
```

Optional runtime cap for faster test runs:

```powershell
$env:MOL_MORPHER_MAX_PAIRS = "200"
python scripts/weight_sensitivity.py
Remove-Item Env:MOL_MORPHER_MAX_PAIRS
```
