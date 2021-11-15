import pandas as pd
import subprocess
import numpy as np
import pickle
import time
import pickle
from goatools.go_enrichment import GOEnrichmentStudy
from goatools import obo_parser

# load file form final table
fertility_df=pd.read_excel('/Users/vpandey/projects/githubs/Fertility_screen/analysisCodes/Phenotype_call_final_080920.xlsx',sheet_name='data')
go_df=pd.read_csv('/Users/vpandey/projects/gitlabs/singleGeneDel/featuresData/GenesByTaxon_Summary.txt',sep='\t')
pathway_df=pd.read_csv('/Users/vpandey/projects/enrichmentData/pbe/PathwaysByGeneIds_Summary.txt',sep='\t')

gene_pathway={}
for i,genes in enumerate(pathway_df['Genes'].to_list()):
    pathway=pathway_df['Pathway'][i]
    gene_list=genes.split('|')

    for gene in gene_list:
        if gene not in gene_pathway.keys():
            gene_pathway[gene]=[]
        gene_pathway[gene].append(pathway)

gene_pathway2={}
for gene,pathways in gene_pathway.items():
    gene_pathway2[gene]='|'.join(pathways)



### Prepare for MPMP
## get orthologus genes pfa nad pbe
df=pd.read_csv('/Users/vpandey/projects/githubs/Fertility_screen/analysisCodes/GenesPFAtoPbe.txt','\t')
pfaTopbe={}
aa=df['Input Ortholog(s)'].str.split(',').to_list()
for i,item in enumerate(aa):
    for pfid in item:
        if pfid not in pfaTopbe.keys():
             pfaTopbe[pfid]=[]
        pfaTopbe[pfid].append(df.loc[i,'Gene ID'])

## read MPMP
mpmp=pd.read_csv('/Users/vpandey/projects/githubs/Fertility_screen/analysisCodes/mpmp.tsv','\t')

gene_mpmp={}
mpmp_gene={}
for idx in mpmp.index:
    pfid=mpmp.loc[idx,'PFID New Name']
    map_name=mpmp.loc[idx,'Map Name']
    # test id in orthology
    if pfid in pfaTopbe.keys():
        for pbeid in pfaTopbe[pfid]:
            if pbeid not in gene_mpmp.keys():
                gene_mpmp[pbeid]=[]
            if map_name not in mpmp_gene.keys():
                mpmp_gene[map_name]=[]
            mpmp_gene[map_name].append(pbeid)
            gene_mpmp[pbeid].append(map_name)


gene_mpmp2={}
for gene,pathways in gene_mpmp.items():
    gene_mpmp2[gene]='|'.join(pathways)



#pickle.dump([mpmp_gene,gene_mpmp],open('/Users/vpandey/projects/enrichmentData/pbe/mpmp_pbe.pickle','wb'))




##

female_cut=-2
male_cut=-2

#female_reduced=fertility_df[(fertility_df['GCKO2_diff_max']<=female_cut) and (fertility_df['GCKO2_diff_min']<female_cut) and (not(fertility_df['GCKO2_pheno']=='No Data'))]
fertility_df['GCKO2_pheno_new']='No power'
fertility_df['g145480_pheno_new']='No power'


fertility_df['MPMP pathway']=np.nan
fertility_df['Pathway']=np.nan

columns=['Computed GO Component IDs', 'Computed GO Components',
       'Computed GO Function IDs', 'Computed GO Functions',
       'Computed GO Process IDs', 'Computed GO Processes',
       'Curated GO Component IDs', 'Curated GO Components',
       'Curated GO Function IDs', 'Curated GO Functions',
       'Curated GO Process IDs', 'Curated GO Processes']
for col in columns:
    fertility_df[col]=np.nan




for ind in fertility_df.index:
    # get information from go term
    pbid=fertility_df.loc[ind,'pbanka_id']
    # new pheno types
    # if (fertility_df.loc[ind,'GCKO2_diff_max']<female_cut) & (fertility_df.loc[ind,'GCKO2_pheno']!='No Data'):
    #     fertility_df.loc[ind,'GCKO2_pheno_new']='Reduced'
    # elif (fertility_df.loc[ind,'GCKO2_diff_min']>=female_cut) & (fertility_df.loc[ind,'GCKO2_pheno']!='No Data'):
    #     fertility_df.loc[ind,'GCKO2_pheno_new']='Not reduced'
    # elif (fertility_df.loc[ind,'GCKO2_pheno']=='No Data'):
    #     fertility_df.loc[ind,'GCKO2_pheno_new']='No data'


    if (fertility_df.loc[ind,'GCKO2_diff_max']<female_cut) & (fertility_df.loc[ind,'GCKO2_pheno']!='No Data'):
        fertility_df.loc[ind,'GCKO2_pheno_new']='Reduced'

    elif (fertility_df.loc[ind,'GCKO2_pheno']=='No Data'):
        fertility_df.loc[ind,'GCKO2_pheno_new']='No data'
    else:
        fertility_df.loc[ind,'GCKO2_pheno_new']='Not reduced'




    # if (fertility_df.loc[ind,'GCKO2_RGR']<female_cut) & (fertility_df.loc[ind,'GCKO2_pheno']!='No Data'):
    #     fertility_df.loc[ind,'GCKO2_pheno_new']='Reduced'
    # elif (fertility_df.loc[ind,'GCKO2_RGR']>=female_cut) & (fertility_df.loc[ind,'GCKO2_pheno']!='No Data'):
    #     fertility_df.loc[ind,'GCKO2_pheno_new']='Not reduced'
    # elif (fertility_df.loc[ind,'GCKO2_pheno']=='No Data'):
    #     fertility_df.loc[ind,'GCKO2_pheno_new']='No data'

    # new pheno types
    # if (fertility_df.loc[ind,'g145480_diff_max']<=male_cut) & (fertility_df.loc[ind,'g145480_pheno']!='No Data'):
    #     fertility_df.loc[ind,'g145480_pheno_new']='Reduced'
    # elif (fertility_df.loc[ind,'g145480_diff_min']>=male_cut) & (fertility_df.loc[ind,'g145480_pheno']!='No Data'):
    #     fertility_df.loc[ind,'g145480_pheno_new']='Not reduced'
    # elif (fertility_df.loc[ind,'g145480_pheno']=='No Data'):
    #     fertility_df.loc[ind,'g145480_pheno_new']='No data'


    if (fertility_df.loc[ind,'g145480_diff_max']<male_cut) & (fertility_df.loc[ind,'g145480_pheno']!='No Data'):
        fertility_df.loc[ind,'g145480_pheno_new']='Reduced'
    elif (fertility_df.loc[ind,'g145480_pheno']=='No Data'):
        fertility_df.loc[ind,'g145480_pheno_new']='No data'
    else:
        fertility_df.loc[ind,'g145480_pheno_new']='Not reduced'


    # ## new pheno types
    # if (fertility_df.loc[ind,'g145480_RGR']<male_cut) & (fertility_df.loc[ind,'g145480_pheno']!='No Data'):
    #     fertility_df.loc[ind,'g145480_pheno_new']='Reduced'
    # elif (fertility_df.loc[ind,'g145480_RGR']>=male_cut) & (fertility_df.loc[ind,'g145480_pheno']!='No Data'):
    #     fertility_df.loc[ind,'g145480_pheno_new']='Not reduced'
    # elif (fertility_df.loc[ind,'g145480_pheno']=='No Data'):
    #     fertility_df.loc[ind,'g145480_pheno_new']='No data'

    ## for MPMP pathway
    if pbid in gene_mpmp2.keys():
        fertility_df.loc[ind,'MPMP pathway']=gene_mpmp2[pbid]


    ## for metabolic pathway
    if pbid in gene_pathway2.keys():
        fertility_df.loc[ind,'Pathway']=gene_pathway2[pbid]




    ## GO terms
    tmp=go_df[go_df['Gene ID']==pbid]
    if not tmp.empty:
        for col in columns:
            fertility_df.loc[ind,col]=tmp[col].to_list()[0]


#fertility_df.to_excel('/Users/vpandey/projects/githubs/Fertility_screen/analysisCodes/Phenotype_call_final_100621.xlsx',sheet_name='data',index=None)



##
female_reduced=fertility_df[fertility_df['GCKO2_pheno_new']=='Reduced']
female_reduced_genes=female_reduced['pbanka_id'].unique().tolist()
female_notreduced=fertility_df[fertility_df['GCKO2_pheno_new']=='Not reduced']
female_notreduced_genes=female_notreduced['pbanka_id'].unique().tolist()

male_reduced=fertility_df[fertility_df['g145480_pheno_new']=='Reduced']
male_reduced_genes=male_reduced['pbanka_id'].unique().tolist()
male_notreduced=fertility_df[fertility_df['g145480_pheno_new']=='Not reduced']
male_notreduced_genes=male_notreduced['pbanka_id'].unique().tolist()

#
female_only=set(female_reduced_genes)&set(male_notreduced_genes)
male_only=set(male_reduced_genes) & set(female_notreduced_genes)
female_male=set(male_reduced_genes) & set(female_reduced_genes)


##
# fm=pd.read_excel('/Users/vpandey/Downloads/Phenotype_call_final_100621_sexspecificlists.xlsx',sheet_name='fm')
# fo=pd.read_excel('/Users/vpandey/Downloads/Phenotype_call_final_100621_sexspecificlists.xlsx',sheet_name='fo')
# mo=pd.read_excel('/Users/vpandey/Downloads/Phenotype_call_final_100621_sexspecificlists.xlsx',sheet_name='mo')
# ##


# female_reduced=fertility_df[(fertility_df['GCKO2_diff_max']<=female_cut) & (fertility_df['GCKO2_pheno']!='No Data')]
# male_reduced=fertility_df[(fertility_df['g145480_diff_max']<=male_cut) & (fertility_df['g145480_pheno']!='No Data')]
#
# female_notreduced=fertility_df[(fertility_df['GCKO2_diff_min']>=female_cut) & (fertility_df['GCKO2_pheno']!='No Data')]
# male_notreduced=fertility_df[(fertility_df['g145480_diff_min']>=male_cut) & (fertility_df['g145480_pheno']!='No Data')]




### add some columns to fertility df

###

xx=fertility_df[~(fertility_df['g145480_pheno_new']=='No data')]
geneIds=xx['pbanka_id'].to_list()
textfile = open("/Users/vpandey/projects/enrichmentData/pbe/fertility_screen_male.txt", "w")
for element in geneIds:
    textfile.write(element + "\n")
textfile.close()


xx=fertility_df[~(fertility_df['GCKO2_pheno_new']=='No data')]
geneIds=xx['pbanka_id'].to_list()
textfile = open("/Users/vpandey/projects/enrichmentData/pbe/fertility_screen_female.txt", "w")
for element in geneIds:
    textfile.write(element + "\n")
textfile.close()




geneIds=male_only
textfile = open("/Users/vpandey/projects/enrichmentData/pbe/male_only_genes.txt", "w")
for element in geneIds:
    textfile.write(element + "\n")
textfile.close()


geneIds=female_only
textfile = open("/Users/vpandey/projects/enrichmentData/pbe/female_only_genes.txt", "w")
for element in geneIds:
    textfile.write(element + "\n")
textfile.close()


geneIds=female_male
textfile = open("/Users/vpandey/projects/enrichmentData/pbe/female_male_genes.txt", "w")
for element in geneIds:
    textfile.write(element + "\n")
textfile.close()



## downlaod go basic
# go_asso=pickle.load(open('/Users/vpandey/projects/enrichmentData/pbe/gene_go_association_all.pickle','rb'))
# pop = go_asso[0].keys()
#
# assoc = {}
#
go = obo_parser.GODag('/Users/vpandey/projects/githubs/goatools/data/go-basic.obo')
#
# for k,v in go_asso[0].items():
#     if k not in assoc:
#         assoc[k] = set()
#     for vv in v:
#         assoc[k].add(vv)

# study=male_only
# methods = ["bonferroni", "fdr"]
# methods = ["bonferroni", "sidak", "holm", "fdr"]
# g = GOEnrichmentStudy(pop, assoc, go,
#                          propagate_counts=True,
#                          alpha=0.05,
#                          methods=methods)
# g_res = g.run_study(study)
# result = []
# for x in g_res:
#     result.append([x.goterm.id, x.goterm.name,x.p_fdr,x.p_bonferroni])



time.sleep(0.5)
out = subprocess.run('python scripts/find_enrichment.py ../../enrichmentData/pbe/male_only_genes.txt ../../enrichmentData/pbe/fertility_screen_male.txt ../../enrichmentData/pbe/gene_go_association_combined_BO.txt --pval=0.5  --ns=BP --method=fdr_bh --pval_field=fdr_bh --obo ../../enrichmentData/pbe/go-basic.obo --outfile=../../enrichmentData/pbe/results_GO_male_BO.xlsx', shell=True)


import pdb; pdb.set_trace()


time.sleep(0.5)
out = subprocess.run('python scripts/find_enrichment.py ../../enrichmentData/pbe/male_only_genes.txt ../../enrichmentData/pbe/pbe_all_genes.txt ../../enrichmentData/pbe/gene_go_association_combined.txt --pval=0.5 --method=fdr_bh --pval_field=fdr_bh --obo ../../enrichmentData/pbe/go-basic.obo --outfile=../../enrichmentData/pbe/results_GO_male.xlsx', shell=True)
time.sleep(0.5)
out = subprocess.run('python scripts/find_enrichment.py ../../enrichmentData/pbe/female_only_genes.txt ../../enrichmentData/pbe/pbe_all_genes.txt ../../enrichmentData/pbe/gene_go_association_combined.txt --pval=0.5 --method=fdr_bh --pval_field=fdr_bh --obo ../../enrichmentData/pbe/go-basic.obo --outfile=../../enrichmentData/pbe/results_GO_female.xlsx', shell=True)
time.sleep(0.5)
out = subprocess.run('python scripts/find_enrichment.py ../../enrichmentData/pbe/female_male_genes.txt ../../enrichmentData/pbe/pbe_all_genes.txt ../../enrichmentData/pbe/gene_go_association_combined.txt --pval=0.5 --method=fdr_bh --pval_field=fdr_bh --obo ../../enrichmentData/pbe/go-basic.obo --outfile=../../enrichmentData/pbe/results_GO_male_female.xlsx', shell=True)
time.sleep(0.5)
out = subprocess.run('python scripts/find_enrichment.py ../../enrichmentData/pbe/male_only_genes.txt ../../enrichmentData/pbe/fertility_screen_male.txt ../../enrichmentData/pbe/gene_go_association_combined.txt --pval=0.5 --method=fdr_bh --pval_field=fdr_bh --obo ../../enrichmentData/pbe/go-basic.obo --outfile=../../enrichmentData/pbe/results_GO_male_screen.xlsx', shell=True)
time.sleep(0.5)
out = subprocess.run('python scripts/find_enrichment.py ../../enrichmentData/pbe/female_only_genes.txt ../../enrichmentData/pbe/fertility_screen_female.txt ../../enrichmentData/pbe/gene_go_association_combined.txt --pval=0.5 --method=fdr_bh --pval_field=fdr_bh --obo ../../enrichmentData/pbe/go-basic.obo --outfile=../../enrichmentData/pbe/results_GO_female_screen.xlsx', shell=True)
time.sleep(0.5)
out = subprocess.run('python scripts/find_enrichment.py ../../enrichmentData/pbe/female_male_genes.txt ../../enrichmentData/pbe/fertility_screen_male.txt ../../enrichmentData/pbe/gene_go_association_combined.txt --pval=0.8 --method=fdr_bh --pval_field=fdr_bh --obo ../../enrichmentData/pbe/go-basic.obo --outfile=../../enrichmentData/pbe/results_GO_male_female_screen.xlsx', shell=True)
