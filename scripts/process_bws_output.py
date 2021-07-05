import pandas as pd 
import sys
import numpy as np
from glob import glob
import os
from os.path import basename
import ast
from itertools import product
import re
from string import Formatter
import shutil
from pathlib import Path

source_format = "/home/clustor/ma/g/gunsim/bws_output/bws_{}*"
target_format = "/home/clustor/ma/g/gunsim/bws_cache/{}"
version_param  = "V1-9"
clear_data = True

def m_columns(cols): return [c for c in cols if c[0:1]=="M" and c[0:2]!= "M0"]
def ensure_dir(d): 
    if not os.path.exists(d): os.makedirs(d)
def make_key(t):  
    d = dict(t)
    s= "H"+str(d["H"])+"D"+str(int(d["D"]))+"B"+str(int(d["BCs"]))
    return s.replace(".","")

def read_df_custom(f):
    lines = []
    #print("opening file ", f)
    with open(f) as fl:
        for i,l in enumerate(fl):
            if l.startswith('#'):
                continue
            l = l.replace("\n", "").replace('\n', "")
            lines.append(l.split('\t'))
    df = pd.DataFrame([l for l in lines[1:-1]], columns = lines[0] ).astype(float)
    return df

def _parse_params(l):
    """
    Assumes parameter line is a valid dictionary/json like format and parses it 
    otherwise just returns the full line
    """
    d = {}
    if l.replace(" ", "").startswith("#Parameters"):
        l = l.replace("#Parameters=", "")
        l = l.replace("# Parameters=", "")
        l=l.replace(" ","").replace("(","").replace(")","")            
        try:
            d = ast.literal_eval(l)
            #if we needed to map names we would do it here
            return d
        
        except:
            try:
                #try look for terms
                terms = l.replace("{","").replace("}","").split(",")
                for t in terms:
                    a,b = t.split(":")
                    a = a.lstrip().rstrip().replace(" ","")
                    #print("parsing terms",a)
                    #we use the non-spaced tokens
                    cl_no_space = [s.replace(" ","") for s in ["Hopping rate", "Realisations", "Chunk size"]]
                    if a in cl_no_space: 
                        d[a] = b
            except:
                return {"unparsed", p}
            
    return d
        
    
def get_file_list(vkey, format_string):
    search_path = format_string.format(vkey)
    print("searching for files", search_path)
    for f in glob(search_path):
        yield f
        
        
def read_info(file, param_line="#Parameters",sep='\t'):
    with open(file) as _f:
        d = {"version": "default"}
        #TODO: extract version and other info from file header
        
        tokens = file.split("-") #should use a smarter regex maybe    
        if len(tokens) > 1: d["version"]= tokens[0].split("_")[-1]+"-"+tokens[1]
        for i,line in enumerate(_f):
            if line.startswith(param_line):
                d.update(_parse_params(line))                    
            if i == 10:#some reasonable comment scanning
                break
        d["file"] = file
        return d  #{"file":f, "L" : None, "D": None, "BCs" : None, "seed" : None, "N":None, "h":None, "T"} #info


def process(version_param, source_format,target_format):
    for f in get_file_list(version_param, source_format): 
        info = read_info(f)
        #print(info)
        H = info["Hoppingrate"]
        df = read_df_custom(f)
        for k,g in df.groupby(["L","D","BCs"]): 
            key = tuple(zip(["L","D","BCs", "H"], list(k)+[H]))
            print("saving data for", key, "to the cache location")
            moment_cols = m_columns(g.columns)
            for m in moment_cols: g[m] /= g["M0"]
            mn = g.groupby("t").mean()[moment_cols]
            st = g.groupby("t").std()[moment_cols] / np.sqrt(g.groupby("t").count())[moment_cols]
            g = mn.join(st, rsuffix="error")
            group_name=make_key(key)
            directory = (target_format+"/{}/").format(version_param,group_name)
            ensure_dir(directory)
            g.to_csv(directory+"L"+str(int(key[0][1]))+".csv")

    print("data saved to cache, see", target_format.format(version_param))
    print("combing data in all folders")
    for directory in glob(target_format.format(version_param)+"/*"):
        print("merging", directory)
        data = {}
        for f in glob(directory+"/L*.*"):
            L = int(basename(f).split(".")[0][1:])
            df= pd.read_csv(f).set_index("t")
            L = "{:0>5}".format(L)
            sorted_cols = list(df.columns.values)
            sorted_cols.sort()
            df = df[sorted_cols]
            #print(df.head(5))
            df.columns = [L+c for c in df.columns]
            data[int(L)] = df
            P = Path(directory)
        ordered_l = list(data.keys())
        ordered_l.sort()
        if len(ordered_l) > 0:
            print("saving concatenation of", len(ordered_l), "files")
            adf = pd.concat([data[i] for i in ordered_l],axis=1)
            filename = directory+"/{}_ALL{}.csv".format(version_param,P.parts[-1])
            adf.to_csv(filename, header=False)
            with open (filename+".headers", "w") as hf:
                for i,c in enumerate(adf.columns):  hf.write("{:0>2}:{}\t".format(i+1, c))
    
    print("done")


if __name__ == "__main__":
    l= len(sys.argv)
    if l > 2:
        source_format = "C:/Users/sirsh/Documents/Code/C/BWS_GITLAB/BranchingWienerProcess/data/bws_output/bws_{}*"
        target_format = "C:/Users/sirsh/Documents/Code/C/BWS_GITLAB/BranchingWienerProcess/data/bws_cache/{}"
    if l >1:  process(sys.argv[1],source_format=source_format,target_format=target_format)
    else:  print("Please supply a version parameter with format (example) v1-10")
    
        