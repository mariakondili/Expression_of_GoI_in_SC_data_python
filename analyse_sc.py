#!/usr/bin/env python

## Subject: Run the script for a selection of Genes,given in a .tsv with Name and Description
##          Calculate Mean_expr and SEM per cell-type along replicates
##          Then: Run then Rscript "create_histo.R", which produces Histograms per cell-type for the genes-of-Interest,along Timepoints.

import os
import numpy as np
import pandas as pd
import argparse
from scipy import stats #> for .sem
import glob
import re


##> Read genes-list from user:
parser = argparse.ArgumentParser(description = 'give a list of genes of interest')
parser.add_argument("-g", "--genelist", required=True,type=argparse.FileType('r'),help="a list of genes (in 1 column) to be searched in the dataset for expression")
parser.add_argument("-o", "--outdir", required=True,help="output directory to save mean-expr. tables.")
parser.add_argument("-d", "--datafile", required=True, type=str,help="the file containing all expression data")
parser.add_argument("-m", "--metadata", required=True, type=str, help="the file of metadata explaining columns of expr.data")

args = parser.parse_args()
GoI_file = args.genelist
GoI = pd.read_table(GoI_file,header=0,sep="\t")
## header= "mm_gene_name" -> in line 0
GoI = GoI.loc[0] #> if >1 cols, use only 1st -> on mouse IDs.

print("Number of genes given: ", len(GoI))
print("The Genes are : \n")
print(GoI)


md = pd.read_table( args.metadata, header=0,sep="\t",index_col=0 )
## md.columns =  sampleID	nUMI	nGene	percent_mito	injury	cell_annotation

## Read the columns of FAPs in the normalized.data.table
data = pd.read_table(args.datafile, header=0,sep="\t",index_col=0)

import pickle
with open('sc_data.pkl', 'wb') as outdataf:
      pickle.dump(data, outdataf, pickle.HIGHEST_PROTOCOL)

#with open('sc_data.pkl', 'rb') as indataf:
#    data = pickle.load(indataf)

data.shape[1]
data.shape[0]

gene_names = data.index
found =  [g in gene_names for g in GoI]
genesfound = dict(zip(GoI,found))

line_genes = [l for l,g in enumerate(gene_names) if g in genesfound and genesfound[g] == True]

outdir = args.outdir
if not(os.path.exists(outdir)):
    os.makedirs(outdir) ## Attention,previous tables will be written-over ,if different genes !


merged_expr = pd.DataFrame()

## Find ALL CELL-types :
cell_types = md["cell_annotation"].unique()

for ct in cell_types :

    print(f"{ct}  will be processed...")
    cell_ids = [md.index[i] for i in range(md.shape[0])  if md.iloc[i]["cell_annotation"] == ct ]
    print(f"{len(cell_ids)} cells are included in {ct}.")
    ## Choose only columns of the cell-type
    tab_celltype    = data[cell_ids]
    gene_cnt_group  = tab_celltype.iloc[line_genes]
    day0_cols = gene_cnt_group.filter(regex="D0_*").columns
    day2_cols = gene_cnt_group.filter(regex="D2_*").columns
    day5_cols = gene_cnt_group.filter(regex="D5_*").columns
    day7_cols = gene_cnt_group.filter(regex="D7_*").columns
    ## Save smaller tables with one cell-type and all timepoints + mean of timepoint over replicates
    ct_day0_Mean_expr = gene_cnt_group[day0_cols].mean(axis=1)
    ct_day2_Mean_expr = gene_cnt_group[day2_cols].mean(axis=1)
    ct_day5_Mean_expr = gene_cnt_group[day5_cols].mean(axis=1)
    ct_day7_Mean_expr = gene_cnt_group[day7_cols].mean(axis=1)
    ct_day0_sem = gene_cnt_group[day0_cols].sem(axis=1)
    ct_day2_sem = gene_cnt_group[day2_cols].sem(axis=1)
    ct_day5_sem = gene_cnt_group[day5_cols].sem(axis=1)
    ct_day7_sem = gene_cnt_group[day7_cols].sem(axis=1)

    ##--- Save a table of Means + SEM
    ct = ct.replace("/","_").replace(" ","_")
    mean_expr_ct  = pd.concat([ct_day0_Mean_expr,ct_day2_Mean_expr, ct_day5_Mean_expr,ct_day7_Mean_expr],axis=1, join="inner", keys= ["D0","D2","D5","D7"])
    mean_expr_ct["GeneName"] = mean_expr_ct.index
    mean_expr_ct["CellType"] = ct
    reshaped_expr_ct = pd.melt(mean_expr_ct, id_vars=["GeneName","CellType"],var_name="TimePoint",value_name="MeanExpr")
    #mean_expr_ct.to_csv(outdir+"mean_expr_per_day_"+ct+".csv",header=True, sep="\t", na_rep = "NA",index=True,encoding="utf-8")
    #print("Mean-expr.table written!")
    sem_expr_ct  = pd.concat([ct_day0_sem, ct_day2_sem, ct_day5_sem,ct_day7_sem],axis=1, join="inner", keys= ["D0","D2","D5","D7"])
    sem_expr_ct["gene_name"] = sem_expr_ct.index
    sem_expr_ct["cell_type"] = ct
    reshaped_sem_ct = pd.melt(sem_expr_ct, id_vars=["GeneName","CellType"],var_name="TimePoint",value_name="SEM")
    #sem_expr_ct.to_csv(outdir+"sem_expr_per_day_"+ct+".csv",header=True, sep="\t", na_rep = "NA",index=True,encoding="utf-8")
    #print("SEM of mean expr.table written!")

    ###--- APPEND to existing TABLES( after adding one column for SEM):
    reshaped_expr_ct["SEM"] = reshaped_sem_ct["SEM"]
    merged_expr = merged_expr.append(reshaped_expr_ct,ignore_index =True)
# end.for


###
### WRITE FINAL MERGED TABLE TO A FILE :
###
os.makedirs(outdir+"Merged/")
outfilepath = os.path.basename(GoI_file)
merged_expr.to_csv(outdir+"Merged/MeanExpr_All_CellTypes_"+ outfilepath, header=True, sep="\t", na_rep = "NA",index=True,encoding="utf-8")


###____SAVE ALL VARIABLES IN ENV ____##
# https://stackoverflow.com/questions/2960864/how-to-save-all-the-variables-in-the-current-python-session
