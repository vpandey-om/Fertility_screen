import pandas as pd
import re
import pickle
old_to_new=pickle.load(open('/Users/vikash/git-hub/Fertility_screen/data/prevTonew_PBANKA.pickle','rb'))

new_to_old = dict([(value, key) for key, value in old_to_new.items()])

phenocall_df=pd.read_csv('/Users/vikash/git-hub/Fertility_screen/Figures/Phenotype_call_final.txt',sep='\t')
phenocall_df=phenocall_df.drop(columns=['Unnamed: 0'])
phenocall_df['PbGEM-ID']='NA'

## read mapping vector
barcode=pd.read_csv('/Users/vikash/git-hub/barseq-obilab/barcode_gene_file.csv',sep=',')
## misssing plasmogem id from claire

#plsamogem=pd.read_csv('/Users/vikash/git-hub/Fertility_screen/data/map_plasmogem_missing_claire.txt',sep='\t')
plsamogem=pd.read_csv('/Users/vikash/git-hub/Fertility_screen/data/all_input.txt',sep='\t')

## old to new
plsamogem['new_gene_id']='NA'

for idx,pid in enumerate(plsamogem['gene_id'].to_list()):
    if pid in old_to_new.keys():
        plsamogem.loc[idx,'new_gene_id']=old_to_new[pid]






for idx,pid in enumerate(phenocall_df['pbanka_id'].to_list()):
    ## test is plasmogem
    tmp=plsamogem[plsamogem['new_gene_id']==pid]
    if not tmp.empty:
        phenocall_df.loc[idx,'PbGEM-ID']=tmp['PbGEM-ID'].to_list()[0]
    else:
        ## try to see in barcode
        npid=pid
        if pid in new_to_old.keys():
            npid=new_to_old[pid]



        df1 = barcode[barcode['gene'].str.contains(npid)]
         ## test size
        if df1.shape[0]>1:
            tmp_list=[]
            for item in df1['gene'].to_list():
                probable_gems=re.findall('(PbGEM-\d+)', item)
                for pr in probable_gems:
                    tmp_list.append(pr)
            phenocall_df.loc[idx,'PbGEM-ID']=','.join(tmp_list)

        elif df1.shape[0]==1:
            probable_gems=re.findall('(PbGEM-\d+)', df1['gene'].to_list()[0])

            phenocall_df.loc[idx,'PbGEM-ID']=','.join(probable_gems)
        else:
            continue

phenocall_df.to_csv('/Users/vikash/git-hub/Fertility_screen/Figures/Phenotype_call_plasmogem_id.txt',sep='\t')
