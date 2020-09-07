# merge all screening data

## read fertility screen data
import pandas as pd
import numpy as np
from sqlalchemy import create_engine
###
fertility_screen_file='/Users/vikash/git-hub/Fertility_screen/Figures/Phenotype_call.txt'
fertility_df=pd.read_csv(fertility_screen_file,sep='\t')
## Gametocyte screen data
gam_screen_file='/Users/vikash/git-hub/Fertility_screen/Figures/Gametocyte_screen.txt'
gam_df=pd.read_csv(gam_screen_file,sep='\t')




## read combined df from the database
POSTGRES = {
    'user': 'vpandey',
    'pw': 'om16042020',
    'db': 'pbe_db',
    'host': 'localhost',
    'port': '5432',
}

DATABASE_URI='postgresql+psycopg2://%(user)s:%(pw)s@%(host)s:%(port)s/%(db)s' % POSTGRES
engine = create_engine(DATABASE_URI)
geneanot = pd.read_sql_query('select * from "geneanot"',con=engine)
pheno_data = pd.read_sql_query('select * from "phenodata"',con=engine)
fertility_columns=fertility_df.columns[5:] ### for fertility screening
gam_columns=['Male gametcytes', 'difference.x', 'differencesd.x', 'p.x',
       'diffmax.x', 'diffmin.x', 'phenotype.x', 'Condition.y',
       'Female gametocytes', 'differencesd.y', 'p.y', 'diffmax.y', 'diffmin.y',
       'power.y', 'minmax']
pheno_columns= pheno_data.columns[2:]

all_columns=[['pbanka_id','old_pbanka_id','gene_description','gene_shortname'],fertility_columns,gam_columns,pheno_columns]
flat_columns = [item for sublist in all_columns for item in sublist]


found_fertility=set(geneanot['pbankaNewID'])&set(fertility_df['pbanka_id'])
not_found_fertility=set(geneanot['pbankaNewID'])-set(fertility_df['pbanka_id'])
extra=set(fertility_df['pbanka_id'])-set(geneanot['pbankaNewID'])

### now we want to make table
new_to_old=dict(zip(geneanot['pbankaNewID'],geneanot['pbankaOldID']))
### get columns Right

##




def map_based_on_given_genes(final_df,list_of_genes,new_to_old,fertility_columns,gam_columns,pheno_columns):
    ''' map gene to all screening data'''
    ind=len(final_df.index)### this is the index
    for new_gene_id in list_of_genes:

        t_gene_df=geneanot[geneanot['pbankaNewID']==new_gene_id]
        if t_gene_df.shape[0]==1:
            final_df.loc[ind,'pbanka_id']=t_gene_df['pbankaNewID'].to_list()[0]
            final_df.loc[ind,'old_pbanka_id']=t_gene_df['pbankaOldID'].to_list()[0]
            final_df.loc[ind,'gene_description']=t_gene_df['description'].to_list()[0]
            final_df.loc[ind,'gene_shortname']=t_gene_df['shortName'].to_list()[0]
        elif t_gene_df.shape[0]>1:
            print('same gene anotated at two place')
        else:
            ## we can not find this gene in gloabl databse because of naming
            #p230p-tag
            final_df.loc[ind,'pbanka_id']=new_gene_id
            final_df.loc[ind,'old_pbanka_id']=new_gene_id
            final_df.loc[ind,'gene_description']='NA'
            final_df.loc[ind,'gene_shortname']='NA'

        ## now put data from the fertility screen
        t_fertility=fertility_df[fertility_df['pbanka_id']==new_gene_id]

        if len(t_fertility.index)>0:
            try:
                final_df.loc[ind,fertility_columns]=t_fertility.loc[t_fertility.index,fertility_columns].squeeze()
            except:
                print('may be two indexes are comming')
        else:

            final_df.loc[ind,fertility_columns]=np.nan


        ## gametocytes screen is made with old ids
        if new_gene_id in new_to_old.keys():
            old_id=new_to_old[new_gene_id]
        else:
            old_id=new_gene_id

        t_gam=gam_df[gam_df['gene']==old_id]

        if len(t_gam.index)>0:
            try:
                final_df.loc[ind,gam_columns]=t_gam.loc[t_gam.index,gam_columns].squeeze()
            except:
                print('may be two indexes are comming')
        else:

            final_df.loc[ind,gam_columns]=np.nan


        ## liver and blood screen
        t_pheno=pheno_data[pheno_data['name']==new_gene_id]

        if len(t_pheno.index)>0:
            try:
                final_df.loc[ind,pheno_columns]=t_pheno.loc[t_pheno.index,pheno_columns].squeeze()
            except:
                print('may be two indexes are comming')
        else:

            final_df.loc[ind,pheno_columns]=np.nan

        ind=ind+1
    return final_df


### create a function for mapping to all data
final_df=pd.DataFrame( columns=flat_columns)
final_df1=map_based_on_given_genes(final_df,fertility_df['pbanka_id'].to_list(),new_to_old,fertility_columns,gam_columns,pheno_columns)

final_df2=map_based_on_given_genes(final_df1,list(not_found_fertility),new_to_old,fertility_columns,gam_columns,pheno_columns)


all_columns=[fertility_columns,gam_columns,pheno_columns]
NA_columns = [item for sublist in all_columns for item in sublist]

xx=final_df2[NA_columns]
NAidx=xx.index[xx.isnull().all(1)]
notNAidx=xx.index[~xx.isnull().all(1)]
result_df=final_df2.loc[np.concatenate([notNAidx,NAidx]),:].copy()
result_df.to_csv('/Users/vikash/git-hub/Fertility_screen/Figures/all_screening_data.txt',sep='\t')
res_df=final_df2.loc[np.concatenate([notNAidx]),:].copy()
res_df.to_csv('/Users/vikash/git-hub/pbeDB/data/all_screening_data_without_nan.txt',sep='\t',index=None)
