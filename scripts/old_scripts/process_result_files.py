source_format = "/home/clustor/ma/g/gunsim/bws_output/bws_{}*"

settings ={
    "cache_location":"/home/clustor/ma/g/gunsim/bws_output/cache", 
    "param_space_format" : "{0:i5}H{1:f}D{2:i}B{3:i}",
    "grouping_keys" : ["L", "H", "D","B"],
    "cols_to_drop" : ["M0", "chunk"],
    "combiner" :{
        "reduce_on" : ["L"],
        "filters" : { "L" : [15,31] , "D" : [2,3,4,5] },
        "enabled" : True
    },
    "header_coerce_list" : ["Hopping rate", "Realisations", "Chunk size"]  ,    
}

#todo - we should validate the inout files against the expectations above early on and give useful warnings!
test_lin = 'C:/Users/sirsh/Documents/Code/C/BWS_GITLAB/BranchingWienerProcess/BWS_VS_PROJECT/BWS/out.data'
#source_format = "C:/Users/sirsh/Documents/Code/C/BWS_GITLAB/BranchingWienerProcess/data/bws_output/bws_{}*"
H_NAME= "Hoppingrate"

import sys
import pandas as pd 
import numpy as np
from glob import glob
import os
from os.path import basename
import ast
from itertools import product
import re
from string import Formatter
  
class MyFormatter(Formatter):
    def format_field(self, value, format_spec):
        if value == None:return "__"
        #print(format_spec)
        if format_spec == 'i':  # Truncate and render as int
            return str(int(value))
        if  format_spec[0] == 'i' and format_spec[1:].isdigit():
            fstring = "{:0"+format_spec[1:]+"d}"
            #print(fstring)
            return fstring.format(int(value))
        if format_spec == 'f':  # Truncate and render as int
            value = int(str(round(float(value),2)).replace("0.","").replace(".0",""))
            fstring = "{:<02}"
            return fstring.format(int(value))
        #default
        return super(MyFormatter, self).format_field(value, format_spec)

    
def process(key, skip_to_combine=False):
    if skip_to_combine == False:
        print("generating atoms for",key)
        l = list(get_file_list(key))

        if len(l) ==0:
            print("there are no files to process when matching", key)
            return

        for f in l: reduce_file(f)
   
    print("")
    print("combining ", key)
    
    combiner(key)
    
    print("done!")
    
def combiner(version_key):
    #using settings, find a regex for file name and its parameters
    #determine various filters 
    #find out what things we want to keep in file e.g. L
    #for each such file load it and concatenate it horz
    #save file according to format using a format that 
    
    _dir = os.path.join(settings["cache_location"], "atoms")
    _dir = os.path.join(_dir, version_key)
    glob_dir = os.path.join(_dir, "*.*")
    print("combining files in ",glob_dir)
    
    file_from_key_desc = lambda f: os.path.join(_dir,f)
    combiner_settings = settings["combiner"]
    reduction = []
    if "reduce_on" in combiner_settings: reduction = combiner_settings["reduce_on"]   
    ordered_keys = settings["grouping_keys"]
    template = {}
    for k in ordered_keys: template[k] = None
    grouping_keys = [g for g in ordered_keys if g not in reduction]
    
    print("splittings groups in", grouping_keys)
    #get groups which are made of matching things that are not in the reduction
    
    #key the files, group by grouping keys
    #iterate groups and get lists - concat them
    #save as a file with a wildcard
    
    files = [f for f in glob(glob_dir)]
    
    print("processing ", len(files), "files")
    keys = []
    for f in files:
        #use just the file name without path or extension
        d = parameters_from_key(basename(f).split(".")[0])
        d["file"] =f
        keys.append(d)
    
    #create a map of the file parameters
    df = pd.DataFrame(keys)

    print("grouping on keys",grouping_keys)
    #split them into groups
    grp = df.groupby(grouping_keys)
    for i,g in grp:
        #create a template for the entire parameter space
        _t = template.copy()
        #for a template of parameter space, take the ordered value from the grouping
        for ordinal, k in enumerate(grouping_keys): _t[k] = i[ordinal]
            
        print("combining on group", _t)
        #these null fields are the ones that are not in the group e.g. {L:None, H=0.95, B=0...}
        null_keys = []
        for k in _t.keys():
            if _t[k] == None: null_keys.append(k)
                
        files = g.file.values
        #print(files)
        dfs = []
        for _f in files:
            dkk = parameters_from_key(basename(_f).split(".")[0])
            #for tihs file name, figure out what is not in the template key
            print("null keys are ", null_keys, "and file keys are", dkk)
            
            qual = "" #make a key for columns with the qualifier
            for k in null_keys: 
                hack = k+dkk[k]
                #this is a terrible hack just for L formatting - need to rethink all formatting
                #if hack.startswith("L"):  hack = str(str(hack.replace("L", "")))
                qual+= hack
                
            _df = pd.read_csv(_f).set_index("t")
            _df.columns = [qual+"_"+c for i, c in enumerate(_df.columns)]
            dfs.append(_df)
        #qualify the columns based on what group they came from
        result = pd.concat(dfs,axis=1)
        #result = result.set_index(result.t)#.drop("t", 1)
        print("sorting columns on moment ordering")
        cols = result.columns.values
        cols = sorted(cols, key = lambda c : str(c))
        result = result[cols]   
        #append an integer col number AFTER choosing an ordering - RENAME
        result.columns =  [str(i) +"_"+c for i,c in enumerate(result.columns)]
        
        file_key = [_t[k] for k in ordered_keys]
        file_name = MyFormatter().format(settings["param_space_format"], *file_key)
        
        print("generating tail for large t")
        
        #this is atemp hack : think about how this large-t generalises. I guess one want have to choose the primary axis only
        df = result.tail(1).T
        df.columns = ["inf"]
        df["M"] = "sigma"
        df.loc[df.index.to_series().str.endswith("mu"), "M"] = "mu"
        #SYS is an optional BUt it might be necessary to do more complex disentanglement from mu/sigs
        #if more dims
        #THESE MAGIC NUMBERS 1 and 2 are the key parts for sys and moment. This shit is SOOOO danderous. 
        #Need to formalise the KEY system and encapsulate in a class ASAP!!
        df["sys"] = df.index.to_series().str.split("_").str.get(1)
        df["Moment"] = df.index.to_series().str.split("_").str.get(2)
        df = df.pivot("sys","Moment", "inf")
        

        print("saving files to ...", os.path.join(settings["cache_location"]))
        df.to_csv(os.path.join(os.path.join(settings["cache_location"]), file_name + "_tail.csv"), sep='\t')
        #reset index to put t inplace and then used the regular index as row counter called 'index'
        result = result.reset_index()
        # SAVE WITHOUT IDEX COLUMN - T has been pop'd
        result.to_csv(os.path.join(os.path.join(settings["cache_location"]), file_name + ".csv"), sep='\t', index=False)

    return pd.DataFrame(keys)
    
    

def parameters_from_key(s):
    """
    the key MUST be a combination of words and numbers: AA12B12C987D93
    the config knows which things are floats based on the format string
    and this can be used to recover the floating point number
    """
    return dict(re.findall(r"([a-z]+)([0-9]+)", s, re.I))


def get_file_list(vkey, format_string =source_format):
    search_path = format_string.format(vkey)
    print("searching for files", search_path)
    for f in glob(search_path):
        yield f
        
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
                    cl_no_space = [s.replace(" ","") for s in settings["header_coerce_list"]]
                    if a in cl_no_space: 
                        d[a] = b
            except:
                return {"unparsed", p}
            
    return d
        
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
    
def _format_keys(keys):
    ordered_keys = settings["grouping_keys"]
    vals = [keys[k] for k in ordered_keys if k in keys]
    print("formating keys ", settings["param_space_format"],vals)
    return MyFormatter().format(settings["param_space_format"], *vals)
    
def iter_coordinate_frame(df):
    print("checking parameter space of dataframe for space ", str(settings["grouping_keys"]))
    space = {}
    
    for c in settings["grouping_keys"]:
        if c in df.columns:
            space [c] = list(df[c].unique())
            
    print("found space ",space, "\nproducing cartesian product parameter space data...")
    
    keys = list(space.keys())
    lists = [space[k] for k in keys]
    P = list(product(*lists))
    print("We have coordinates ("+str(keys)+"): "+str(P))   
    print("slicing dataframes for each of these coords...")
    for i,c in enumerate(P):
        map_keys = dict(zip(keys,list(c)))
        print( str(i)+":", str(map_keys))   
        print("generated a key: ",_format_keys(map_keys))
        yield _format_keys(map_keys), df
    
def reduce_file(f):
    """
    reduce accross chunks on other categorical columns. capture statistics for Moments M0+
    because we might have multiple keys per dataframe, for simplicity with split them here and reduce to one coord
    """
    _df = None
    
    metadata = read_info(f)  
    print("reducing file",f)
    print("HEADER:")
    print(metadata)
    print("*************************")
    #df = pd.read_csv(f, comment="#",sep='\t',dtype=float)   
    df = read_df_custom(f)
    #add or modifiy columns on the data datafram
    #TODO: broken abstraction: we want to use different names for our params - in generall 1/2 chars upper case preferred
    df["H"] = metadata[H_NAME]
    #rename from what we called them in the output
    df = df.rename(columns={"BCs":"B"})
    
    for key, _df in iter_coordinate_frame(df):
        _df = reduce_df(_df)
        _dir = os.path.join(settings["cache_location"],"atoms")
        if metadata["version"] != None: _dir = os.path.join(_dir, metadata["version"])
        _try_make_dir(_dir)
        _file = os.path.join(_dir, key+".csv")
        print("saving file", _file)
        #method to write comments in header
        #http://stackoverflow.com/questions/29233496/write-comments-in-csv-file-with-pandas
        _df.to_csv(_file)
    #returns both the dataframe and the version for accumulation
    return _df, metadata["version"]

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
    
def reduce_df(df):   
    moments = ["M"+str(i) for i in range(1,9)]
    print("normalising moments ",str(moments))
    for m in moments: df[m] /= df["M0"]
    print("checking grouping keys {} are in dataframe".format(str(settings["grouping_keys"])))
    gkeys = [g for g in settings["grouping_keys"] if g in df.columns]
    #add t to the list of things to group by 
    grp = df.groupby(gkeys + ["t"])
    df_mean = grp.mean().reset_index().set_index("t")
    #sample a chunk value
    N = grp.count()["M0"].values[0]
    #standard error
    df_std = grp.std() / np.sqrt(N)
    #setup the join index
    df_std =df_std.reset_index().set_index("t")
    print("dropping columns; moments ",str(settings["cols_to_drop"]))
    for c in settings["cols_to_drop"]:
        if c in df_mean.columns: df_mean = df_mean.drop(c,1)
        if c in df_std.columns: df_std = df_std.drop(c,1)
            
    print("dropping columns in the grouping, post grouping ",gkeys)
    for c in gkeys:
        if c in df_mean.columns: df_mean = df_mean.drop(c,1)
        if c in df_std.columns: df_std = df_std.drop(c,1)

    result = df_mean.join(df_std, lsuffix="mu", rsuffix="sigma")
    
    print("sorting columns on moment ordering")
    cols = result.columns.values
    cols.sort()

    result = result[cols]
    return result

def collect_hist(key, cache=False):
    #go through the files and save /histograms or return them for plotting with their key
    pass


def plot_helper(key, data=None):
    """
    if data is none, try to load from cache
    plotting type can be added to key/dictionary
    """
    pass

def _try_make_dir(dir_name):
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
        

if __name__ == "__main__":
    l= len(sys.argv) 
    if l >1:
        process(sys.argv[1],False)
    else:
        print("Please supply a version parameter with format (example) v1-10")
    
        
    