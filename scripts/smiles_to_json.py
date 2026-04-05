import pandas as pd 
import json 
from rdkit import Chem ,RDLogger 
from itertools import product 
from pathlib import Path 

def _resolve_column (df ,column_name ='smiles'):
    if column_name in df .columns :
        return column_name 

    target =column_name .strip ().lower ()
    for col in df .columns :
        if str (col ).strip ().lower ()==target :
            return col 

    raise KeyError (f"Column '{column_name }' not found in CSV")


def load_smiles (csv_path ,column_name ='smiles'):
    df =pd .read_csv (csv_path ,dtype =str )
    col =_resolve_column (df ,column_name =column_name )

    smiles_list =[]
    seen =set ()
    invalid_count =0 
    skipped_text_count =0 


    RDLogger .DisableLog ('rdApp.error')
    try :
        for raw in df [col ].dropna ().astype (str ):
            smi =raw .strip ()
            if not smi :
                continue 

            low =smi .lower ()
            if low in {'smiles','smiles not found'}or any (ch .isspace ()for ch in smi ):
                skipped_text_count +=1 
                continue 

            mol =Chem .MolFromSmiles (smi )
            if not mol :
                invalid_count +=1 
                continue 

            canon =Chem .MolToSmiles (mol ,canonical =True ,isomericSmiles =True )
            if canon not in seen :
                seen .add (canon )
                smiles_list .append (canon )
    finally :
        RDLogger .EnableLog ('rdApp.error')

    return smiles_list ,invalid_count ,skipped_text_count 

def count_atoms (smi ):
    mol =Chem .MolFromSmiles (smi )
    return mol .GetNumHeavyAtoms ()if mol else 0 

def generate_smiles_pairs (
smiles_list ,
size_metric ='length',
min_diff =10 
):

    if size_metric =='atoms':
        sizes ={s :count_atoms (s )for s in smiles_list }
    else :
        sizes ={s :len (s )for s in smiles_list }

    for s1 ,s2 in product (smiles_list ,smiles_list ):
        if s1 ==s2 :continue 
        if sizes [s2 ]-sizes [s1 ]>=min_diff :
            yield s1 ,s2 

if __name__ =='__main__':
    root_dir =Path (__file__ ).resolve ().parent .parent 
    CSV_PATH =root_dir /"data"/"input"/"Smiles.csv"
    OUTPUT_PATH =root_dir /"data"/"output"/"Data.jsonl"
    OUTPUT_PATH .parent .mkdir (parents =True ,exist_ok =True )
    smiles ,invalid_count ,skipped_text_count =load_smiles (CSV_PATH ,column_name ='smiles')
    print (
    f"Loaded {len (smiles )} valid unique SMILES "
    f"(skipped text rows: {skipped_text_count }, invalid molecules: {invalid_count })"
)

    with OUTPUT_PATH .open ('w',encoding ='utf-8')as fout :
        for s1 ,s2 in generate_smiles_pairs (smiles ,size_metric ='length',min_diff =80 ):
            record ={"src":s1 ,"tgt":s2 }
            fout .write (json .dumps (record )+"\n")

    print (f"Written all pairs to {OUTPUT_PATH }")