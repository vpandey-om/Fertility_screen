import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import itertools
import numpy as np
import pickle

prev_to_new=pickle.load(open('/Users/vpandey/projects/githubs/Fertility_screen/data/prevTonew_PBANKA.pickle','rb'))
new_to_prev = dict((v,k) for k,v in prev_to_new.items())

def get_subset_data(df,col_name,gene_list,types=['Female','Male','Female and male']):
    ## find common genes
    res_df=df.copy()
    list_tmp=[]
    for i,type in enumerate(types):
        for gene in gene_list[i]:
            tmp=df[df[col_name]==gene]

            if not tmp.empty:
                tmp['Sex']='NA'
                tmp.loc[tmp.index,'Sex']=type
                list_tmp.append(tmp)
    tmp_df=pd.concat(list_tmp)

    return tmp_df

def get_subset_data_2IDs(df,col_name,col_name2,gene_list,types=['Female','Male','Female and male']):
    ## find common genes
    res_df=df.copy()
    list_tmp=[]
    for i,type in enumerate(types):
        for gene in gene_list[i]:
            tmp=df[df[col_name]==gene]

            if not tmp.empty:
                tmp['Sex']='NA'
                tmp.loc[tmp.index,'Sex']=type
                list_tmp.append(tmp)
            else:

                if gene in new_to_prev.keys():
                    old_gene=new_to_prev[gene]
                    tmp=df[df[col_name2]==old_gene]
                    if not tmp.empty:
                        tmp['Sex']='NA'
                        tmp.loc[tmp.index,'Sex']=type
                        list_tmp.append(tmp)


    tmp_df=pd.concat(list_tmp)

    return tmp_df


def compute_ttest(df,col='FG RGR', types=['Female','Male','Female and male']):
    ''' We are going to compute t test  '''
    ##
    # get list
    RGR_values=[]
    for t in types:
        type_0=df[df['Sex']==t].copy()
        type_0=type_0[~type_0[col].isnull()].copy()
        RGR_values.append(type_0[col].values)
    ## now do pair wise t test
    combi=[item for item in itertools.combinations(np.arange(len(types)), 2)]
    df=pd.DataFrame()
    comp1=[]
    comp2=[]
    pvals=[]
    tvals=[]
    for item in combi:
        print(" for %s\n" %(col))
        print("%s vs %s" %(types[item[0]],types[item[1]]))
        tval,pval=stats.ttest_ind(RGR_values[item[0]], RGR_values[item[1]],equal_var = False)
        print("t-statistic=%g, p-value=%g" %(tval,pval))
        comp1.append(types[item[0]])
        comp2.append(types[item[1]])
        pvals.append(pval)
        tvals.append(tval)
    df['Condition 1']=comp1
    df['Condition 2']=comp2
    df[col.replace(' RGR','')+'_p_values']=pvals
    df[col.replace(' RGR','')+'_t_statistic']=tvals
    return df

folder='/Users/vpandey/projects/githubs/Fertility_screen/preFinals/'

df=pd.read_csv('/Users/vpandey/projects/githubs/Fertility_screen/preFinals/female_only_genes.txt',sep='\t',header=None)
female_only=df[0].tolist()
female_only_old=[new_to_prev[item] for item in female_only ]

df=pd.read_csv('/Users/vpandey/projects/githubs/Fertility_screen/preFinals/male_only_genes.txt',sep='\t',header=None)
male_only=df[0].tolist()
male_only_old=[new_to_prev[item]for item in male_only if item in new_to_prev.keys()]

df=pd.read_csv('/Users/vpandey/projects/githubs/Fertility_screen/preFinals/female_male_genes.txt',sep='\t',header=None)
male_female=df[0].tolist()

male_female_old=[new_to_prev[item] for item in male_female ]

all_screen=pd.read_excel('/Users/vpandey/projects/githubs/Fertility_screen/analysisCodes/all_stage_phenodata.xlsx',sheet_name='pheno_data')

fertility_screen=pd.read_excel('/Users/vpandey/projects/githubs/Fertility_screen/preFinals/Phenotype_call_final_100621.xlsx',sheet_name='data')
gam_screen=pd.read_excel('/Users/vpandey/projects/githubs/Fertility_screen/analysisCodes/Gametocyte_screen_results_preprint_final.xlsx',sheet_name='Combined phenoypes')

gam_df=get_subset_data_2IDs(gam_screen,'current_version_ID', 'gene',[female_only,male_only,male_female],['Female','Male','Female and male'])
gam_df['pbanka_id']=gam_df['current_version_ID']

hue_order=['Not reduced','Reduced']
joint_df=pd.merge(gam_df,fertility_screen, on='pbanka_id',how='left')
joint_df['gam']='NA'
joint_df.loc[(joint_df['GFP_result']=='Reduced'),'gam']='Male'
joint_df.loc[(joint_df['RFP_result']=='Reduced'),'gam']='Female'
joint_df.loc[(joint_df['RFP_result']=='Reduced')& (joint_df['GFP_result']=='Reduced'),'gam']='Female and male'


joint_df_male=joint_df[(joint_df['GFP_result']=='Reduced')]


rel=sns.scatterplot(data=joint_df_male, x="GFP_log2change", y="g145480_RGR",hue="g145480_pheno_new",hue_order=hue_order).set(title='Redcued at male gametocyte')
plt.savefig('male.pdf')
plt.close()
joint_df_female=joint_df[joint_df['RFP_result']=='Reduced']
rel=sns.scatterplot(data=joint_df_female, x="RFP_log2change", y="GCKO2_RGR",hue="GCKO2_pheno_new",hue_order=hue_order).set(title='Redcued at female gametocyte')
plt.savefig('female.pdf')
plt.close()





joint_df_required=joint_df[(joint_df['gam']=='Male')| (joint_df['gam']=='Female')| (joint_df['gam']=='Female and male')]

df_MG=compute_ttest(joint_df_required,col="g145480_RGR")
colors=['#F7F7F7', '#CCCCCC', '#969696']
sns.stripplot(x="gam", y="g145480_RGR", data=joint_df_required,color="black", edgecolor="gray")

sns_plot = sns.violinplot(x="gam", y="g145480_RGR", data=joint_df_required,inner='box',
                    order=['Female','Male','Female and male'],palette=colors) ## 'quartile' 'box'
figure = sns_plot.get_figure()
figure.savefig(folder+'maleFertility_vs_gam.pdf', dpi=300)
plt.close()



df_MG=compute_ttest(joint_df_required,col="GCKO2_RGR")
colors=['#F7F7F7', '#CCCCCC', '#969696']
sns.stripplot(x="gam", y="GCKO2_RGR", data=joint_df_required,color="black", edgecolor="gray")

sns_plot = sns.violinplot(x="gam", y="GCKO2_RGR", data=joint_df_required,inner='box',
                    order=['Female','Male','Female and male'],palette=colors) ## 'quartile' 'box'
figure = sns_plot.get_figure()
figure.savefig(folder+'femaleFertility_vs_gam.pdf', dpi=300)
plt.close()


joint_df_required.to_excel('male_female_gam.xlsx')



import pdb; pdb.set_trace()

## for gametocyte analysis

res_df=get_subset_data(all_screen,'name', [female_only,male_only,male_female],['Female','Male','Female and male'])


### Do analysis with phenotyeps
tmp=res_df[res_df['O P']=='reduced']

plt.clf()

df_FG=compute_ttest(tmp,col='O RGR')
colors=['#F7F7F7', '#CCCCCC', '#969696']
sns.stripplot(x="Sex", y="O RGR", data=tmp,color="black", edgecolor="gray")

sns_plot = sns.violinplot(x="Sex", y="O RGR", data=tmp,inner='box',
                    order=['Female','Male','Female and male'],palette=colors) ## 'quartile' 'box'
figure = sns_plot.get_figure()
figure.savefig(folder+'occyst_pheno.pdf', dpi=300)

import pdb; pdb.set_trace()









###




## For female gametocyte

plt.clf()

df_FG=compute_ttest(res_df,col='FG RGR')
colors=['#F7F7F7', '#CCCCCC', '#969696']
sns.stripplot(x="Sex", y="FG RGR", data=res_df,color="black", edgecolor="gray")

sns_plot = sns.violinplot(x="Sex", y="FG RGR", data=res_df,inner='box',
                    order=['Female','Male','Female and male'],palette=colors) ## 'quartile' 'box'
figure = sns_plot.get_figure()
figure.savefig(folder+'FemaleG.pdf', dpi=300)


## For female gametocyte

plt.clf()

df_MG=compute_ttest(res_df,col='MG RGR')
colors=['#F7F7F7', '#CCCCCC', '#969696']
sns.stripplot(x="Sex", y="MG RGR", data=res_df,color="black", edgecolor="gray")

sns_plot = sns.violinplot(x="Sex", y="MG RGR", data=res_df,inner='box',
                    order=['Female','Male','Female and male'],palette=colors) ## 'quartile' 'box'
figure = sns_plot.get_figure()
figure.savefig(folder+'MaleG.pdf', dpi=300)


## for occyst
plt.clf()

df_O=compute_ttest(res_df,col='O RGR')
colors=['#F7F7F7', '#CCCCCC', '#969696']
sns.stripplot(x="Sex", y="O RGR", data=res_df,color="black", edgecolor="gray")

sns_plot = sns.violinplot(x="Sex", y="O RGR", data=res_df,inner='box',
                    order=['Female','Male','Female and male'],palette=colors) ## 'quartile' 'box'
figure = sns_plot.get_figure()
figure.savefig(folder+'occyst.pdf', dpi=300)

two_df=df_MG.merge(df_FG, how='left', on=['Condition 1','Condition 2'])
final_df=two_df.merge(df_O, how='left', on=['Condition 1','Condition 2'])
final_df.to_csv(folder+'t_stats.txt',sep='\t')
