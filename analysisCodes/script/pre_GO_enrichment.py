##
import pandas as pd
## read GO terms form PlasmoDB 1) curated and computed GO terms
## GO analysis works for 3 types 1) biological process
## 2) Molecular function
## 3) Components

# read gene ontology file for Plasmodb


def get_gene_GO_asso(df,columns,file=None):
    '''identify gene to GO association'''
    print('start')
    gene_go_asso={}

    for col in columns:
        xdf=df[~df[col].isna()]
        asso_df=xdf[['Gene ID',col]]
        for ind,gene in enumerate(asso_df['Gene ID'].to_list()):
            go_list=asso_df.loc[asso_df.index[ind],col].split(';')
            if gene not in gene_go_asso.keys():
                gene_go_asso[gene]=[]
            for go in go_list:
                if go not in gene_go_asso[gene]:
                    gene_go_asso[gene].append(go)
    if file:
        with open(file, 'w') as f:
            for key, value in gene_go_asso.items():
                vv=';'.join(value)
                f.write('%s\t%s\n' % (key, vv))


    return gene_go_asso


def find_goid_name(df,combined_columns,combined_columns_name):
    go_id_name={};

    for i,item in enumerate(combined_columns):

        xdf=df[~df[item].isna()]
        for ind in xdf.index:
            id=xdf.loc[ind,item].split(';')
            name=xdf.loc[ind,combined_columns_name[i]].split(';')
            for go in zip(id,name):
                # if go[0] not in go_id_name.keys():
                #     go_id_name[go[0]]=[]
                go_id_name[go[0]]=go[1]
    return go_id_name




df=pd.read_csv('/Users/vpandey/projects/gitlabs/singleGeneDel/featuresData/GenesByTaxon_Summary.txt',sep='\t')
# make 7 files for pbe enrichment analysis


geneIds=df['Gene ID'].unique().tolist()
textfile = open("/Users/vpandey/projects/enrichmentData/pbe/pbe_all_genes.txt", "w")
for element in geneIds:
    textfile.write(element + "\n")
textfile.close()

## first make 3 files
#1)computed
#2)curated
#3)combined

anotate_df=df[~(df['Computed GO Component IDs'].isna() & df['Computed GO Function IDs'].isna() & df['Computed GO Process IDs'].isna()
& df['Curated GO Component IDs'].isna() & df['Curated GO Function IDs'].isna() & df['Curated GO Process IDs'].isna())]

geneIds=anotate_df['Gene ID'].unique().tolist()
textfile = open("/Users/vpandey/projects/enrichmentData/pbe/pbe_all_genes.txt", "w")
for element in geneIds:
    textfile.write(element + "\n")
textfile.close()

computed_columns=['Computed GO Component IDs','Computed GO Function IDs','Computed GO Process IDs']
curated_columns=['Curated GO Component IDs','Curated GO Function IDs','Curated GO Process IDs']
combined_columns=computed_columns+curated_columns

columns1=['Computed GO Components','Computed GO Functions','Computed GO Processes']
columns2=['Curated GO Components','Curated GO Functions','Curated GO Processes']
combined_columns_name=columns1+columns2



types=['CC','MF','BP']

# sperate GO terms
assso=[]
go_id_names=[]
for i,type in enumerate(types):

    xx=gene_go_association_combined=get_gene_GO_asso(anotate_df,[computed_columns[i],curated_columns[i]],file='/Users/vpandey/projects/enrichmentData/pbe/gene_go_association_combined.txt')

    assso.append(xx)

    yy=find_goid_name(df,[computed_columns[i],curated_columns[i]],[columns1[i],columns2[i]])
    go_id_names.append(yy)


import pickle
pickle.dump([assso,go_id_names,types],open('/Users/vpandey/projects/enrichmentData/pbe/gene_go_association_all.pickle','wb'))



BO=['Computed GO Process IDs','Curated GO Process IDs']
gene_go_association_combined=get_gene_GO_asso(anotate_df,BO,file='/Users/vpandey/projects/enrichmentData/pbe/gene_go_association_combined_BO.txt')

geneIds=gene_go_association_combined.keys()
textfile = open("/Users/vpandey/projects/enrichmentData/pbe/pbe_all_genes_BO.txt", "w")
for element in geneIds:
    textfile.write(element + "\n")
textfile.close()



import pdb; pdb.set_trace()
gene_go_association_computed=get_gene_GO_asso(anotate_df,computed_columns,file='/Users/vpandey/projects/enrichmentData/pbe/gene_go_association_computed.txt')
gene_go_association_combined=get_gene_GO_asso(anotate_df,combined_columns,file='/Users/vpandey/projects/enrichmentData/pbe/gene_go_association_combined.txt')
gene_go_association_curated=get_gene_GO_asso(anotate_df,curated_columns,file='/Users/vpandey/projects/enrichmentData/pbe/gene_go_association_curated.txt')



columns1=['Computed GO Components','Computed GO Functions','Computed GO Processes']
columns2=['Curated GO Components','Curated GO Functions','Curated GO Processes']
combined_columns_name=columns1+columns2

## gete GO process to GO ids
go_id_name={};

for i,item in enumerate(combined_columns):

    xdf=df[~df[item].isna()]


    for ind in xdf.index:
        id=xdf.loc[ind,item].split(';')
        name=xdf.loc[ind,combined_columns_name[i]].split(';')
        for go in zip(id,name):
            # if go[0] not in go_id_name.keys():
            #     go_id_name[go[0]]=[]
            go_id_name[go[0]]=go[1]

#pickle.dump([gene_go_association_combined,go_id_name,gene_go_association_curated,gene_go_association_computed,['combined','go_id_to_name_dict','curated','computed']],open('/Users/vpandey/projects/enrichmentData/pbe/gene_go_association_all.pickle','wb'))




##
geneIds=gene_go_association_combined.keys()
textfile = open("/Users/vpandey/projects/enrichmentData/pbe/pbe_all_genes.txt", "w")
for element in geneIds:
    textfile.write(element + "\n")
textfile.close()

### gene to go association

with open('/Users/vpandey/projects/enrichmentData/pbe/gene_go_association_combined.txt', 'w') as file:
    for k,v in gene_go_association_combined.items():
     file.write("%s\t%s\n"%(k,';'.join(v)))






# asso_df=xdf[['Gene ID','Computed GO Component IDs']]
# asso_df.to_csv('/Users/vpandey/projects/enrichmentData/pbe/testasso.txt',sep='\t',index=None)
