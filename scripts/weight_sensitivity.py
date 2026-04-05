from __future__ import annotations 
import csv 
import json 
import os 
import random 
import statistics 
import sys 
import time 
from dataclasses import dataclass 
from itertools import product 
from pathlib import Path 
from rdkit import Chem ,RDLogger 

SCRIPT_DIR =Path (__file__ ).resolve ().parent 
PROJECT_ROOT =SCRIPT_DIR .parent 
NOTEBOOK_PATH =SCRIPT_DIR /"multi_beam.ipynb"
SMILES_CSV_PATH =PROJECT_ROOT /"data"/"input"/"Smiles.csv"
SMILES_COLUMN ="smiles"
OUTPUT_DIR =PROJECT_ROOT /"data"/"output"/"weight_search_full"

PAIR_MODE ="all_ordered"
SIZE_METRIC ="length"
MIN_SIZE_DIFF =80 
MAX_PAIRS :int |None =5 
MAX_PAIRS_ENV_VAR ="MOL_MORPHER_MAX_PAIRS"

if os .environ .get (MAX_PAIRS_ENV_VAR ):
    try :
        MAX_PAIRS =int (os .environ [MAX_PAIRS_ENV_VAR ])
    except ValueError as exc :
        raise ValueError (f"{MAX_PAIRS_ENV_VAR } must be an integer")from exc 

if MAX_PAIRS is not None and MAX_PAIRS <=0 :
    raise ValueError ("MAX_PAIRS must be a positive integer or None")

WEIGHT_GRID_MODE ="focused"
WEIGHT_CANDIDATES ={
"focused":{
"alpha":[0.8 ,1.0 ,1.2 ],
"beta":[0.05 ,0.1 ,0.2 ],
"gamma":[0.3 ,0.5 ,0.7 ],
},
"refined":{
"alpha":[0.85 ,0.95 ,1.05 ,1.15 ],
"beta":[0.06 ,0.08 ,0.1 ,0.12 ,0.15 ],
"gamma":[0.35 ,0.45 ,0.55 ,0.65 ],
},
"expanded":{
"alpha":[0.7 ,0.85 ,1.0 ,1.15 ,1.3 ],
"beta":[0.03 ,0.05 ,0.1 ,0.15 ,0.25 ],
"gamma":[0.2 ,0.35 ,0.5 ,0.65 ,0.8 ],
},
}
if WEIGHT_GRID_MODE not in WEIGHT_CANDIDATES :
    raise ValueError (f"Unknown WEIGHT_GRID_MODE: {WEIGHT_GRID_MODE }")

ALPHA_VALUES =WEIGHT_CANDIDATES [WEIGHT_GRID_MODE ]["alpha"]
BETA_VALUES =WEIGHT_CANDIDATES [WEIGHT_GRID_MODE ]["beta"]
GAMMA_VALUES =WEIGHT_CANDIDATES [WEIGHT_GRID_MODE ]["gamma"]

RUN_BIDIRECTIONAL =True 
RESUME =True 
WRITE_INTERMEDIATES =False 
PROGRESS_EVERY =20 
PAIR_RANDOM_SEED =42 
QUIET_PIPELINE_PRINTS =True 

SEARCH_OVERRIDES :dict [str ,float |int ]={}

RANKING_METRICS =[
("coverage",0.20 ,True ),
("hit_rate",0.35 ,True ),
("mean_best_sim_target",0.20 ,True ),
("mean_sim_target",0.10 ,True ),
("mean_sa",0.10 ,False ),
("error_rate",0.05 ,False ),
]

@dataclass 
class PairStats :
    rows :int =0 
    unique_smiles :set [str ]|None =None 
    sim_target :list [float ]|None =None 
    sim_source :list [float ]|None =None 
    sa_score :list [float ]|None =None 
    comp_score :list [float ]|None =None 

    def __post_init__ (self )->None :
        if self .unique_smiles is None :
            self .unique_smiles =set ()
        if self .sim_target is None :
            self .sim_target =[]
        if self .sim_source is None :
            self .sim_source =[]
        if self .sa_score is None :
            self .sa_score =[]
        if self .comp_score is None :
            self .comp_score =[]

    def add_row (self ,row :dict )->None :
        self .rows +=1 
        smi =row .get ("smiles","")
        if smi :
            self .unique_smiles .add (smi )
        self .sim_target .append (_to_float (row .get ("sim_target",0.0 )))
        self .sim_source .append (_to_float (row .get ("sim_source",0.0 )))
        self .sa_score .append (_to_float (row .get ("sa_score",0.0 )))
        self .comp_score .append (_to_float (row .get ("comp_score",0.0 )))

    def as_metrics (self ,threshold :float )->dict [str ,float |int |bool ]:
        if self .rows ==0 :
            return {
            "rows":0 ,
            "unique_smiles":0 ,
            "max_sim_target":0.0 ,
            "mean_sim_target":0.0 ,
            "mean_sim_source":0.0 ,
            "mean_sa":0.0 ,
            "mean_comp":0.0 ,
            "hit_threshold":False ,
            }

        max_sim_t =max (self .sim_target )
        return {
        "rows":self .rows ,
        "unique_smiles":len (self .unique_smiles ),
        "max_sim_target":max_sim_t ,
        "mean_sim_target":statistics .mean (self .sim_target ),
        "mean_sim_source":statistics .mean (self .sim_source ),
        "mean_sa":statistics .mean (self .sa_score ),
        "mean_comp":statistics .mean (self .comp_score ),
        "hit_threshold":max_sim_t >=threshold ,
        }


class MetadataWriter :
    def __init__ (
    self ,
    pair_stats :PairStats ,
    sink_writer :csv .DictWriter |None ,
    weight_label :str ,
    pair_id :int ,
    source :str ,
    target :str ,
    direction :str ,
    )->None :
        self .pair_stats =pair_stats 
        self .sink_writer =sink_writer 
        self .weight_label =weight_label 
        self .pair_id =pair_id 
        self .source =source 
        self .target =target 
        self .direction =direction 

    def writerow (self ,row :dict )->None :
        self .pair_stats .add_row (row )
        if self .sink_writer is None :
            return 

        out ={
        "weight_label":self .weight_label ,
        "pair_id":self .pair_id ,
        "source":self .source ,
        "target":self .target ,
        "direction":self .direction ,
        "smiles":row .get ("smiles",""),
        "sa_score":row .get ("sa_score",""),
        "sim_source":row .get ("sim_source",""),
        "sim_target":row .get ("sim_target",""),
        "comp_score":row .get ("comp_score",""),
        "path_id":row .get ("path_id",""),
        }
        self .sink_writer .writerow (out )

def _to_float (value :float |str )->float :
    try :
        return float (value )
    except Exception :
        return 0.0 

def _parse_bool (value :str |bool )->bool :
    return str (value ).strip ().lower ()in {"1","true","yes","y"}

def _resolve_column_name (fieldnames :list [str ],column_name :str )->str :
    if column_name in fieldnames :
        return column_name 

    target =column_name .strip ().lower ()
    for col in fieldnames :
        if str (col ).strip ().lower ()==target :
            return col 

    raise KeyError (f"Column '{column_name }' not found in CSV")

def weight_label (alpha :float ,beta :float ,gamma :float )->str :
    def fmt (x :float )->str :
        return f"{x :.3f}".rstrip ("0").rstrip (".").replace (".","p")

    return f"a{fmt (alpha )}_b{fmt (beta )}_g{fmt (gamma )}"

def load_pipeline_namespace ()->dict :
    if str (SCRIPT_DIR )not in sys .path :
        sys .path .insert (0 ,str (SCRIPT_DIR ))

    with NOTEBOOK_PATH .open ("r",encoding ="utf-8")as f :
        nb =json .load (f )

    source_code =None 
    for cell in nb .get ("cells",[]):
        if cell .get ("cell_type")!="code":
            continue 
        src ="\n".join (cell .get ("source",[]))
        if "def morph_direction"in src and "def beam_search"in src :
            source_code =src 
            break 

    if not source_code :
        raise RuntimeError ("Could not find main algorithm cell in multi_beam.ipynb")

    ns ={"__name__":"weight_search_runner"}
    exec (source_code ,ns )

    ns ["display"]=lambda *args ,**kwargs :None 
    if QUIET_PIPELINE_PRINTS :
        ns ["print"]=lambda *args ,**kwargs :None 

    for key ,value in SEARCH_OVERRIDES .items ():
        if key not in ns :
            raise KeyError (f"Invalid SEARCH_OVERRIDES key: {key }")
        ns [key ]=value 

    return ns 


def load_smiles_list (csv_path :Path ,column_name :str )->list [str ]:
    smiles :list [str ]=[]
    seen :set [str ]=set ()
    invalid =0 
    skipped_text =0 

    if not csv_path .exists ():
        raise FileNotFoundError (f"SMILES CSV not found: {csv_path }")

    with csv_path .open ("r",newline ="",encoding ="utf-8")as f :
        reader =csv .DictReader (f )
        fieldnames =list (reader .fieldnames or [])
        col =_resolve_column_name (fieldnames ,column_name )


        RDLogger .DisableLog ("rdApp.error")
        try :
            for row in reader :
                raw =(row .get (col )or "").strip ()
                if not raw :
                    continue 

                low =raw .lower ()
                if low in {"smiles","smiles not found"}or any (ch .isspace ()for ch in raw ):
                    skipped_text +=1 
                    continue 

                mol =Chem .MolFromSmiles (raw )
                if mol is None :
                    invalid +=1 
                    continue 
                can =Chem .MolToSmiles (mol ,canonical =True ,isomericSmiles =True )
                if can not in seen :
                    seen .add (can )
                    smiles .append (can )
        finally :
            RDLogger .EnableLog ("rdApp.error")

    print (
    f"Loaded {len (smiles )} unique valid SMILES "
    f"(invalid skipped: {invalid }, text rows skipped: {skipped_text })"
    )
    return smiles 

def precompute_sizes (smiles :list [str ])->list [int ]:
    if SIZE_METRIC =="length":
        return [len (s )for s in smiles ]

    if SIZE_METRIC =="atoms":
        out =[]
        for s in smiles :
            mol =Chem .MolFromSmiles (s )
            out .append (mol .GetNumHeavyAtoms ()if mol is not None else 0 )
        return out 

    raise ValueError ("SIZE_METRIC must be 'length' or 'atoms'")

def iter_pairs (smiles :list [str ],sizes :list [int ]):
    emitted =0 
    pair_id =1 
    n =len (smiles )

    for i in range (n ):
        src =smiles [i ]
        src_size =sizes [i ]
        for j in range (n ):
            if i ==j :
                continue 

            tgt =smiles [j ]
            if PAIR_MODE =="size_filtered":
                if sizes [j ]-src_size <MIN_SIZE_DIFF :
                    continue 
            elif PAIR_MODE !="all_ordered":
                raise ValueError ("PAIR_MODE must be 'all_ordered' or 'size_filtered'")

            yield pair_id ,src ,tgt 
            emitted +=1 
            pair_id +=1 
            if MAX_PAIRS is not None and emitted >=MAX_PAIRS :
                return 


def count_pairs (smiles :list [str ],sizes :list [int ])->int :
    return sum (1 for _ in iter_pairs (smiles ,sizes ))

def run_one_pair (
ns :dict ,
alpha :float ,
beta :float ,
gamma :float ,
pair_id :int ,
source :str ,
target :str ,
threshold :float ,
label :str ,
intermediate_writer :csv .DictWriter |None ,
)->dict :
    seed =PAIR_RANDOM_SEED +pair_id 
    random .seed (seed )
    ns ["random"].seed (seed )

    ns ["ALPHA"]=alpha 
    ns ["BETA"]=beta 
    ns ["GAMMA"]=gamma 

    stats =PairStats ()
    start =time .perf_counter ()

    stop_forward =""
    stop_reverse =""
    status ="ok"
    error_message =""

    try :
        seen :set [str ]=set ()
        fw =MetadataWriter (
        pair_stats =stats ,
        sink_writer =intermediate_writer ,
        weight_label =label ,
        pair_id =pair_id ,
        source =source ,
        target =target ,
        direction ="forward",
        )
        stop_forward ,next_pid =ns ["morph_direction"](
        source ,
        target ,
        fw ,
        seen ,
        path_id_start =1 ,
        )

        if RUN_BIDIRECTIONAL :
            rw =MetadataWriter (
            pair_stats =stats ,
            sink_writer =intermediate_writer ,
            weight_label =label ,
            pair_id =pair_id ,
            source =source ,
            target =target ,
            direction ="reverse",
            )
            stop_reverse ,_ =ns ["morph_direction"](
            target ,
            source ,
            rw ,
            seen ,
            path_id_start =next_pid ,
            )

    except Exception as exc :
        status ="error"
        error_message =f"{type (exc ).__name__ }: {exc }"

    elapsed =time .perf_counter ()-start 
    m =stats .as_metrics (threshold )

    if status !="ok":
        m ["hit_threshold"]=False 

    return {
    "weight_label":label ,
    "pair_id":pair_id ,
    "source":source ,
    "target":target ,
    "alpha":alpha ,
    "beta":beta ,
    "gamma":gamma ,
    "status":status ,
    "error_message":error_message ,
    "stop_forward":stop_forward ,
    "stop_reverse":stop_reverse ,
    "rows":m ["rows"],
    "unique_smiles":m ["unique_smiles"],
    "max_sim_target":f"{m ['max_sim_target']:.6f}",
    "mean_sim_target":f"{m ['mean_sim_target']:.6f}",
    "mean_sim_source":f"{m ['mean_sim_source']:.6f}",
    "mean_sa":f"{m ['mean_sa']:.6f}",
    "mean_comp":f"{m ['mean_comp']:.6f}",
    "hit_threshold":m ["hit_threshold"],
    "seconds":f"{elapsed :.6f}",
    }

def read_processed_pair_ids (pair_metrics_path :Path )->set [int ]:
    done :set [int ]=set ()
    with pair_metrics_path .open ("r",newline ="",encoding ="utf-8")as f :
        reader =csv .DictReader (f )
        for row in reader :
            try :
                done .add (int (row ["pair_id"]))
            except Exception :
                continue 
    return done 

def summarize_weight (pair_metrics_path :Path ,expected_pairs :int )->dict [str ,float |int ]:
    completed =0 
    ok_pairs =0 
    error_pairs =0 
    hit_pairs =0 

    best_sim_list :list [float ]=[]
    mean_sim_t_list :list [float ]=[]
    mean_sim_s_list :list [float ]=[]
    mean_sa_list :list [float ]=[]
    mean_comp_list :list [float ]=[]
    rows_per_pair :list [int ]=[]
    seconds_per_pair :list [float ]=[]

    with pair_metrics_path .open ("r",newline ="",encoding ="utf-8")as f :
        reader =csv .DictReader (f )
        for row in reader :
            completed +=1 
            status =(row .get ("status")or "").strip ().lower ()
            if status =="ok":
                ok_pairs +=1 
                if _parse_bool (row .get ("hit_threshold","false")):
                    hit_pairs +=1 
                best_sim_list .append (_to_float (row .get ("max_sim_target",0.0 )))
                mean_sim_t_list .append (_to_float (row .get ("mean_sim_target",0.0 )))
                mean_sim_s_list .append (_to_float (row .get ("mean_sim_source",0.0 )))
                mean_sa_list .append (_to_float (row .get ("mean_sa",0.0 )))
                mean_comp_list .append (_to_float (row .get ("mean_comp",0.0 )))
                rows_per_pair .append (int (_to_float (row .get ("rows",0 ))))
                seconds_per_pair .append (_to_float (row .get ("seconds",0.0 )))
            else :
                error_pairs +=1 
                seconds_per_pair .append (_to_float (row .get ("seconds",0.0 )))

    denom_total =expected_pairs if expected_pairs >0 else 1 

    return {
    "expected_pairs":expected_pairs ,
    "completed_pairs":completed ,
    "coverage":completed /denom_total ,
    "ok_pairs":ok_pairs ,
    "error_pairs":error_pairs ,
    "error_rate":error_pairs /denom_total ,
    "hit_pairs":hit_pairs ,
    "hit_rate":hit_pairs /denom_total ,
    "mean_best_sim_target":statistics .mean (best_sim_list )if best_sim_list else 0.0 ,
    "mean_sim_target":statistics .mean (mean_sim_t_list )if mean_sim_t_list else 0.0 ,
    "mean_sim_source":statistics .mean (mean_sim_s_list )if mean_sim_s_list else 0.0 ,
    "mean_sa":statistics .mean (mean_sa_list )if mean_sa_list else 0.0 ,
    "mean_comp":statistics .mean (mean_comp_list )if mean_comp_list else 0.0 ,
    "mean_rows_per_pair":statistics .mean (rows_per_pair )if rows_per_pair else 0.0 ,
    "total_rows":sum (rows_per_pair ),
    "total_seconds":sum (seconds_per_pair ),
    "mean_seconds_per_pair":statistics .mean (seconds_per_pair )if seconds_per_pair else 0.0 ,
    }

def rank_weight_summaries (summaries :list [dict ])->list [dict ]:
    if not summaries :
        return []

    extrema :dict [str ,tuple [float ,float ]]={}
    for metric ,_coef ,_maximize in RANKING_METRICS :
        values =[float (s .get (metric ,0.0 ))for s in summaries ]
        extrema [metric ]=(min (values ),max (values ))

    ranked =[]
    for row in summaries :
        score =0.0 
        for metric ,coef ,maximize in RANKING_METRICS :
            value =float (row .get (metric ,0.0 ))
            lo ,hi =extrema [metric ]
            if hi ==lo :
                norm =1.0 
            else :
                norm =(value -lo )/(hi -lo )
            if not maximize :
                norm =1.0 -norm 
            score +=coef *norm 

        out =dict (row )
        out ["rank_score"]=score 
        ranked .append (out )

    ranked .sort (key =lambda x :x ["rank_score"],reverse =True )
    for i ,row in enumerate (ranked ,start =1 ):
        row ["rank"]=i 
    return ranked 


def save_config (smiles_count :int ,pair_count :int ,weights :list [tuple [float ,float ,float ]])->None :
    cfg ={
    "smiles_csv_path":str (SMILES_CSV_PATH ),
    "smiles_column":SMILES_COLUMN ,
    "notebook_path":str (NOTEBOOK_PATH ),
    "output_dir":str (OUTPUT_DIR ),
    "weight_grid_mode":WEIGHT_GRID_MODE ,
    "pair_mode":PAIR_MODE ,
    "size_metric":SIZE_METRIC ,
    "min_size_diff":MIN_SIZE_DIFF ,
    "max_pairs":MAX_PAIRS ,
    "run_bidirectional":RUN_BIDIRECTIONAL ,
    "resume":RESUME ,
    "write_intermediates":WRITE_INTERMEDIATES ,
    "pair_random_seed":PAIR_RANDOM_SEED ,
    "quiet_pipeline_prints":QUIET_PIPELINE_PRINTS ,
    "search_overrides":SEARCH_OVERRIDES ,
    "weight_grid":[{"alpha":a ,"beta":b ,"gamma":g }for a ,b ,g in weights ],
    "smiles_count":smiles_count ,
    "pair_count":pair_count ,
    "ranking_metrics":[
    {"metric":m ,"coefficient":c ,"maximize":mx }
    for m ,c ,mx in RANKING_METRICS 
    ],
    }
    with (OUTPUT_DIR /"run_config.json").open ("w",encoding ="utf-8")as f :
        json .dump (cfg ,f ,indent =2 )


def main ()->None :
    start_all =time .perf_counter ()

    OUTPUT_DIR .mkdir (parents =True ,exist_ok =True )

    smiles =load_smiles_list (SMILES_CSV_PATH ,SMILES_COLUMN )
    if len (smiles )<2 :
        raise RuntimeError ("Need at least two valid SMILES for pair runs")

    sizes =precompute_sizes (smiles )
    pair_count =count_pairs (smiles ,sizes )

    weights =list (product (ALPHA_VALUES ,BETA_VALUES ,GAMMA_VALUES ))
    if not weights :
        raise RuntimeError ("Weight grid is empty")

    print (f"Pair mode: {PAIR_MODE }")
    print (f"Total candidate pairs: {pair_count }")
    if MAX_PAIRS is None :
        print ("Pair cap: all pairs")
    else :
        print (f"Pair cap: first {MAX_PAIRS } pairs")
    print (f"Total weight settings: {len (weights )}")
    total_pair_weight_runs =pair_count *len (weights )
    print (f"Total pair-weight evaluations: {total_pair_weight_runs }")
    if pair_count >100000 :
        print ("WARNING: Pair count is very large and can take a long time.")
    if total_pair_weight_runs >1000000 :
        print ("WARNING: Pair-weight evaluations are very large. Consider MAX_PAIRS or focused grid.")

    save_config (len (smiles ),pair_count ,weights )

    ns =load_pipeline_namespace ()
    threshold =float (ns .get ("SIM_THRESHOLD",95.0 ))

    all_summaries :list [dict ]=[]

    pair_metrics_fields =[
    "weight_label",
    "pair_id",
    "source",
    "target",
    "alpha",
    "beta",
    "gamma",
    "status",
    "error_message",
    "stop_forward",
    "stop_reverse",
    "rows",
    "unique_smiles",
    "max_sim_target",
    "mean_sim_target",
    "mean_sim_source",
    "mean_sa",
    "mean_comp",
    "hit_threshold",
    "seconds",
    ]

    interm_fields =[
    "weight_label",
    "pair_id",
    "source",
    "target",
    "direction",
    "smiles",
    "sa_score",
    "sim_source",
    "sim_target",
    "comp_score",
    "path_id",
    ]

    for idx ,(alpha ,beta ,gamma )in enumerate (weights ,start =1 ):
        label =weight_label (alpha ,beta ,gamma )
        print ("-"*80 )
        print (
        f"[{idx }/{len (weights )}] Running weight {label } "
        f"(alpha={alpha }, beta={beta }, gamma={gamma })"
        )

        weight_dir =OUTPUT_DIR /label 
        weight_dir .mkdir (exist_ok =True )

        pair_metrics_path =weight_dir /"pair_metrics.csv"
        interm_path =weight_dir /"intermediates.csv"

        processed_ids :set [int ]=set ()
        if RESUME and pair_metrics_path .exists ():
            processed_ids =read_processed_pair_ids (pair_metrics_path )
            print (f"Resuming: {len (processed_ids )} already processed pairs")

        pair_metrics_mode ="a"if pair_metrics_path .exists ()and RESUME else "w"
        interm_mode ="a"if interm_path .exists ()and RESUME else "w"

        processed_now =0 
        started =time .perf_counter ()

        with pair_metrics_path .open (pair_metrics_mode ,newline ="",encoding ="utf-8")as pair_f :
            pair_writer =csv .DictWriter (pair_f ,fieldnames =pair_metrics_fields )
            if pair_metrics_mode =="w":
                pair_writer .writeheader ()

            if WRITE_INTERMEDIATES :
                interm_f =interm_path .open (interm_mode ,newline ="",encoding ="utf-8")
                interm_writer =csv .DictWriter (interm_f ,fieldnames =interm_fields )
                if interm_mode =="w":
                    interm_writer .writeheader ()
            else :
                interm_f =None 
                interm_writer =None 

            try :
                for pair_id ,src ,tgt in iter_pairs (smiles ,sizes ):
                    if pair_id in processed_ids :
                        continue 

                    metrics =run_one_pair (
                    ns =ns ,
                    alpha =alpha ,
                    beta =beta ,
                    gamma =gamma ,
                    pair_id =pair_id ,
                    source =src ,
                    target =tgt ,
                    threshold =threshold ,
                    label =label ,
                    intermediate_writer =interm_writer ,
                    )
                    pair_writer .writerow (metrics )
                    processed_now +=1 

                    if processed_now %PROGRESS_EVERY ==0 :
                        elapsed =time .perf_counter ()-started 
                        done_total =len (processed_ids )+processed_now 
                        print (
                        f"{label }: processed {done_total }/{pair_count } pairs "
                        f"({elapsed /60.0 :.1f} min in this session)"
                        )
            finally :
                if interm_f is not None :
                    interm_f .close ()

        summary =summarize_weight (pair_metrics_path ,pair_count )
        summary .update (
        {
        "weight_label":label ,
        "alpha":alpha ,
        "beta":beta ,
        "gamma":gamma ,
        }
        )

        with (weight_dir /"summary.json").open ("w",encoding ="utf-8")as f :
            json .dump (summary ,f ,indent =2 )

        all_summaries .append (summary )
        print (
        f"Done {label }: coverage={summary ['coverage']:.3f}, "
        f"hit_rate={summary ['hit_rate']:.3f}, "
        f"mean_best_sim_t={summary ['mean_best_sim_target']:.3f}, "
        f"mean_sa={summary ['mean_sa']:.3f}, "
        f"errors={summary ['error_pairs']}"
        )

    ranked =rank_weight_summaries (all_summaries )

    leaderboard_path =OUTPUT_DIR /"weight_leaderboard.csv"
    if ranked :
        leaderboard_fields =list (ranked [0 ].keys ())
        with leaderboard_path .open ("w",newline ="",encoding ="utf-8")as f :
            writer =csv .DictWriter (f ,fieldnames =leaderboard_fields )
            writer .writeheader ()
            for row in ranked :
                writer .writerow (row )

        best =ranked [0 ]
        with (OUTPUT_DIR /"best_weight.json").open ("w",encoding ="utf-8")as f :
            json .dump (best ,f ,indent =2 )

        print ("="*80 )
        print ("Best default candidate based on configured ranking:")
        print (
        f"weight={best ['weight_label']} "
        f"alpha={best ['alpha']} beta={best ['beta']} gamma={best ['gamma']}"
        )
        print (
        f"rank_score={best ['rank_score']:.4f} "
        f"hit_rate={best ['hit_rate']:.4f} "
        f"mean_best_sim_target={best ['mean_best_sim_target']:.4f} "
        f"mean_sa={best ['mean_sa']:.4f}"
        )
    else :
        print ("No summaries found. Nothing to rank.")

    elapsed_all =time .perf_counter ()-start_all 
    print ("="*80 )
    print (f"Finished full weight search in {elapsed_all /60.0 :.2f} minutes")
    print (f"Output directory: {OUTPUT_DIR }")

if __name__ =="__main__":
    main ()