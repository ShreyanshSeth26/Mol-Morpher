import pandas as pd
import json
from rdkit import Chem
from itertools import product

def load_smiles(csv_path, column_name='smiles'):
    df = pd.read_csv(csv_path)
    smiles_list = []
    for smi in df[column_name].dropna().astype(str):
        mol = Chem.MolFromSmiles(smi)
        if mol:
            smiles_list.append(Chem.MolToSmiles(mol))
    return smiles_list

def count_atoms(smi):
    mol = Chem.MolFromSmiles(smi)
    return mol.GetNumHeavyAtoms() if mol else 0

def generate_smiles_pairs(
    smiles_list,
    size_metric='length',
    min_diff=10
):
    # Precompute sizes
    if size_metric == 'atoms':
        sizes = {s: count_atoms(s) for s in smiles_list}
    else:
        sizes = {s: len(s) for s in smiles_list}

    for s1, s2 in product(smiles_list, smiles_list):
        if s1 == s2: continue
        if sizes[s2] - sizes[s1] >= min_diff:
            yield s1, s2

if __name__ == '__main__':
    CSV_PATH = r'X:\Projects\IP\old scripts\Smiles.csv'
    OUTPUT_PATH = "Data.jsonl"

    smiles = load_smiles(CSV_PATH, column_name='smiles')

    with open(OUTPUT_PATH, 'w', encoding='utf-8') as fout:
        for s1, s2 in generate_smiles_pairs(smiles, size_metric='length', min_diff=80): #change min diff. in sizew between smiles A and B
            record = {"src": s1, "tgt": s2}
            fout.write(json.dumps(record) + "\n")

    print(f"Written all pairs to {OUTPUT_PATH}")