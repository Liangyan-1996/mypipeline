#! /public/home/liunangroup/liangyan/software/miniconda3/bin/python

# example:
# python stringtieOut2TPM.py 2-cell-rep1.txt 2-cell-rep2.txt
# output file: 2-cell_tpm.txt

import sys
import pandas as pd

def GetTPM(exp_file):
    df = pd.read_table(exp_file,sep="\t",header=0)
    clean = df.loc[:,["Gene ID","TPM"]].groupby("Gene ID").apply(lambda x:x.sum(),include_groups=False)
    return(clean)

def get_id():
    gtf = pd.read_table("/public/home/liunangroup/liangyan/Genome/ensembl/hg38/Homo_sapiens.GRCh38.110.chr.gtf",sep="\t",header=None,comment="#")
    # gtf = pd.read_table("/public/home/liunangroup/liangyan/Genome/ensembl/mm39/Mus_musculus.GRCm39.110.chr.gtf",sep="\t",header=None,comment="#")
    gtf = gtf.loc[gtf.iloc[:,2]=="gene",[0,3,4,6,8]]
    gtf["Gene ID"] = gtf.loc[:,8].str.extract(r'gene_id "(.*?)"')
    gtf["Gene Name"] = gtf.loc[:,8].str.extract(r'gene_name "(.*?)"')
    gtf = gtf.iloc[:,[-2,-1,0,1,2,3]]
    gtf.columns = ["Gene ID","Gene Name","Reference","Start","End","Strand"]
    return(gtf.set_index("Gene ID"))

filelist = sys.argv[1:]

id = get_id()
df = pd.concat(map(GetTPM,filelist),axis="columns")
df.columns = filelist
df = pd.concat([id,df],axis="columns")
df["TPM_mean"] = df.iloc[:,(-1 * len(filelist)):].mean(axis="columns")

df.to_csv("1.expr_TPM.txt",index=True,sep="\t",header=True,float_format = '%.2f')
