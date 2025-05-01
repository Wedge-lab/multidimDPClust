#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import sys

# python 02_GenerateMutWT.py  DW1-4 0 4


case = sys.argv[1]
run_info_path=sys.argv[2]
outpath=sys.argv[3]


Input_path = outpath + "/" + case + "/DPinput/"

union_SNV = pd.read_csv(Input_path + case + "_loci_combined.txt",sep="\t",header=None)
union_SNV.columns = ["chr", "pos","ref","alt"]
union_SNV.index = union_SNV.apply(lambda x : x["chr"] + "_" + str(x["pos"]), axis=1)
union_SNV.drop_duplicates(subset=["chr","pos"],inplace=True)
samplenames = pd.read_csv(run_info_path + "/run_info.txt",sep="\t",index_col=0).loc[case].Tumour.values

mutCount = pd.DataFrame(index = union_SNV.index)
WTcount = pd.DataFrame(index = union_SNV.index)
for name  in samplenames:
    print(name)
    AF = pd.read_csv(Input_path + name + "_aC_out",sep="\t")
    AF.columns = ["chr","pos","A","C","G","T","depth"]
    AF.index = AF.apply(lambda x : x["chr"] + "_" + str(x["pos"]), axis=1)
    get_from_ref = union_SNV[union_SNV["ref"].apply(lambda x: len(x)==1)]
    info = AF.loc[get_from_ref.index].drop_duplicates()
    refcount = pd.concat([info,get_from_ref[["ref","alt"]]],axis=1)
    refcount["ref_count"] = refcount.apply(lambda x: x[x["ref"]],axis=1)
    refcount["alt_count"] = refcount.apply(lambda x: int(x["depth"])- int(x["ref_count"]),axis=1)

    # get_from_alt = union_SNV[union_SNV["ref"].apply(lambda x: len(x)!=1)]
    # info1 = AF.loc[get_from_alt.index].drop_duplicates()
    # altcount = pd.concat([info1,get_from_alt[["ref","alt"]]],axis=1)
    # altcount["alt_count"] = altcount.apply(lambda x: x[x["alt"]],axis=1)
    # altcount["ref_count"] = altcount.apply(lambda x: int(x["depth"])- int(x["alt_count"]),axis=1)
    # count = pd.concat([refcount,altcount])
    count = refcount[['chr', 'pos', 'ref', 'alt', 'ref_count','alt_count']]
    
    WT = count[['ref_count']]
    mut = count[['alt_count']]
    WTcount[name] = WT
    mutCount[name] = mut
    
    
WTcount.reset_index(inplace = True)
WTcount = WTcount.rename(columns={"index":"ID"})
mutCount.reset_index(inplace = True)
mutCount = mutCount.rename(columns={"index":"ID"})
mutCount.fillna(0).to_csv(Input_path + case + "_mutCount.txt",sep="\t", index= None)
WTcount.fillna(0).to_csv(Input_path + case + "_WTCount.txt",sep="\t", index= None)