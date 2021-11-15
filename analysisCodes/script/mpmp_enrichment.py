import pandas as pd

import pickle
import os
import pandas as pd
from scipy.stats import hypergeom ## for hypergeomtric distributions
import numpy as np
import statsmodels.stats.multitest as mtest ## for multiple hypothesis testing


def apply_enrich(pathway_to_genes,targetList,M,outfile):
    ## groupby along Phenotyepes
    #First argumnet: temperory dataframe
    #Second argument: list of male-specific genes
    #Third argumnet: number of genes in background or population
    #Fourth argumnet: file where we store enrichment analysis

    pvals=[] # list for storing pvalues
    n_list=[] # store number of genes in cluster or pathway
    x_list=[] ## store number of target genes in cluster or pathway
    genes=[] ## we store genes target geens
    pathways=[] ## store pathway

    ## group clusters
    for k,v in pathway_to_genes.items():
        ## print(k, len(v))
        pathways.append(k) ## store cluster number
        clust_genes=set(v) ## genes in the cluster

        sex_gene_in_clust= set(targetList)& clust_genes ## how many sex-specific genes are found in cluster
        ## how many are in phenoype
        M=M ## total genes
        N=len(targetList) ## sex-specific genes
        n=len(clust_genes) ## total number of cluster genes

        rv = hypergeom(M, n, N) ## hypergeometric  distributions
        x=len(sex_gene_in_clust) ## number of sex-specific genes in cluster

        ## for p-value we need to get summation of all probalities greater than > x
        pval= 1-rv.cdf(x)
        pvals.append(pval)
        n_list.append(n)
        x_list.append(x)
        genes.append(', '.join(sex_gene_in_clust))
    fdr=mtest.multipletests(pvals, alpha=0.05, method='fdr_bh')
    enrich={'Pathways':pathways, 'Pvals':pvals,'Pathway size':n_list,'Number of genes found in pathway':x_list,'FDR':fdr[1],'Genes':genes}
    ## create enrichment table
    enrich_df=pd.DataFrame.from_dict(enrich)
    enrich_soreted=enrich_df.sort_values(by=['Pvals'],ascending=True)
    enrich_soreted.to_csv(outfile,sep='\t',index=None)
    return enrich_soreted



## load mpmp data
mpmp_data=pickle.load(open('/Users/vpandey/projects/enrichmentData/pbe/mpmp_pbe.pickle','rb'))
pathway_to_genes=mpmp_data[0]
print(len(pathway_to_genes))
## remove some of the categories
mpmp_remove=pd.read_excel('/Users/vpandey/projects/enrichmentData/pbe/MPMP.xlsx',sheet_name='rCS',header=None)

for k in mpmp_remove[0].to_list():
    if k in pathway_to_genes:
       del pathway_to_genes[k]

print(len(pathway_to_genes))


##
## read input file
input_df=pd.read_csv("/Users/vpandey/projects/githubs/Fertility_screen/preFinals/male_only_genes.txt",sep= "\t",header=None)
genelist=input_df.loc[:,0].to_list()


## number of background genes it could be total genes or screen genens
flat_list = [item for k,sublist in pathway_to_genes.items() for item in sublist]
M=len(set(flat_list)) ## number of genes are screened
outfile="/Users/vpandey/projects/githubs/Fertility_screen/preFinals/male_mpmp_enrich.txt"
outdf=apply_enrich(pathway_to_genes,genelist,M,outfile)


## for female

input_df=pd.read_csv("/Users/vpandey/projects/githubs/Fertility_screen/preFinals/female_only_genes.txt",sep= "\t",header=None)
genelist=input_df.loc[:,0].to_list()
outfile="/Users/vpandey/projects/githubs/Fertility_screen/preFinals/female_mpmp_enrich.txt"
outdf=apply_enrich(pathway_to_genes,genelist,M,outfile)


##
input_df=pd.read_csv("//Users/vpandey/projects/githubs/Fertility_screen/preFinals/female_male_genes.txt",sep= "\t",header=None)
genelist=input_df.loc[:,0].to_list()
outfile="/Users/vpandey/projects/githubs/Fertility_screen/preFinals/female_male_mpmp_enrich.txt"
outdf=apply_enrich(pathway_to_genes,genelist,M,outfile)
