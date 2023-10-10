import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import itertools
import numpy as np
import pickle
# import pylab
# pylab.rcParams['xtick.major.pad']='30'
# pylab.rcParams['ytick.major.pad']='4'
prev_to_new=pickle.load(open('/Users/vpandey/projects/githubs/Fertility_screen/data/prevTonew_PBANKA.pickle','rb'))
new_to_prev = dict((v,k) for k,v in prev_to_new.items())

colors=['#AF58BA', '#009ADE', '#FFC61E']
colorv=['#FFFFFF', '#FFFFFF', '#FFFFFF']
# color_gam=['#636363', '#636363', '#636363']
# color_gam=['#66C2A5','#FC8D62']
color_gam_dict={'Reduced':'#FC8D62','Not reduced':'#66C2A5'}
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

folder='/Users/vpandey/projects/githubs/Fertility_screen_2/preFinals/'

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

df_MG=compute_ttest(joint_df_required,col='GFP_log2change')

### plot for male
male_tmp=joint_df_required.copy()
rename_sex={'Female':'Female','Male':'Male','Female and male':'Both sexes'}
male_tmp['gam']=male_tmp['gam'].replace(rename_sex)
male_tmp['Gametocyte']=male_tmp['gam'].copy()
male_tmp['Relative male fertility']=male_tmp['g145480_RGR'].copy()
male_tmp['Phenotypes']=male_tmp['g145480_pheno_new'].copy()



# Relative female fertility (y), Female/Male/Both sexes (x)
# sns.set(font="Arial",font_scale =2,rc = {'figure.figsize':(6,8)})
sns.set(font="Arial",font_scale =2)
sns.set_style("white")


# sns.set_context("talk", font_scale=1.2)
plt.figure()

# sns.stripplot(x="Gametocyte", y="Relative male fertility", data=male_tmp,palette=color_gam_dict,  hue='Phenotypes',edgecolor="gray",size=8)

sns.violinplot(x="Gametocyte", y="Relative male fertility", data=male_tmp,inner='box',palette=colorv,order=['Female','Male','Both sexes']) ## 'quartile' 'box'

# sns.swarmplot(y="Relative male fertility",
#                 x="Gametocyte",
#                 data=male_tmp,
#                    hue='Phenotypes',edgecolor="gray",size=8,palette=color_gam_dict,order=['Female','Male','Both sexes'])

sns.stripplot(y="Relative male fertility",
                x="Gametocyte",
                data=male_tmp,
                   hue='Phenotypes',edgecolor="gray",size=8,palette=color_gam_dict,order=['Female','Male','Both sexes'])



# figure = sns_plot.get_figure()

plt.legend().remove()
plt.tick_params(axis='x', rotation=0)
# plt.xlabel('Gametocyte', fontsize=20,labelpad=10)
plt.xlabel('')
plt.ylabel('')
plt.ylim((-18,8))

plt.savefig(folder+'maleFertility_vs_gam.pdf',
            format='pdf',dpi=300)


# # plt.ylabel('Relative male fertility', fontsize=20)
# figure.savefig(folder+'maleFertility_vs_gam.pdf', dpi=300)

plt.close()


### plot for fmale
female_tmp=joint_df_required.copy()
rename_sex={'Female':'Female','Male':'Male','Female and male':'Both sexes'}
female_tmp['gam']=female_tmp['gam'].replace(rename_sex)
female_tmp['Gametocyte']=female_tmp['gam'].copy()
female_tmp['Relative female fertility']=female_tmp['GCKO2_RGR'].copy()
female_tmp['Phenotypes']=female_tmp['GCKO2_pheno_new'].copy()



df_FG=compute_ttest(joint_df_required,col='RFP_log2change')
# colors=['#F7F7F7', '#CCCCCC', '#969696']

plt.figure()

sns.violinplot(x="Gametocyte", y="Relative female fertility", data=female_tmp,inner='box',
                    order=['Female','Male','Both sexes'],palette=colorv) ## 'quartile' 'box'
# sns.swarmplot(y="Relative female fertility",
#                 x="Gametocyte",
#                 data=female_tmp,
#                    hue='Phenotypes',edgecolor="gray",size=8,palette=color_gam_dict,order=['Female','Male','Both sexes'])
sns.stripplot(y="Relative female fertility",
                x="Gametocyte",
                data=female_tmp,
                   hue='Phenotypes',edgecolor="gray",size=8,palette=color_gam_dict,order=['Female','Male','Both sexes'])


#sns.stripplot(x="Gametocyte", y="Relative female fertility", data=female_tmp,palette=color_gam_dict, hue='Phenotypes',edgecolor="gray",size=8)


# figure = sns_plot.get_figure()
plt.legend().remove()
plt.tick_params(axis='x', rotation=0)
# plt.xlabel('Gametocyte', fontsize=20,labelpad=10)
plt.xlabel('')
plt.ylabel('')
plt.ylim((-18,8))
plt.savefig(folder+'femaleFertility_vs_gam.pdf',
            format='pdf',dpi=300)
# figure.savefig(folder+'femaleFertility_vs_gam.pdf', dpi=300)
plt.close()


joint_df_required.to_excel('male_female_gam.xlsx')





## for gametocyte analysis

res_df=get_subset_data(all_screen,'name', [female_only,male_only,male_female],['Female','Male','Female and male'])


### Do analysis with phenotyeps
tmp=res_df[res_df['O P']=='reduced']

plt.clf()


# colors=['#F7F7F7', '#CCCCCC', '#969696']

### do analysis for oocyst

oocys_df=tmp.copy()
rename_sex={'Female':'Female','Male':'Male','Female and male':'Both sexes'}
oocys_df['Sex']=oocys_df['Sex'].replace(rename_sex)
oocys_df['oo']=oocys_df['Sex'].copy()
oocys_df['Relative oocyst growth']=oocys_df["O RGR"].copy()
oocys_df['Phenotypes']=oocys_df["O P"].copy()

plt.figure()

sns.violinplot(x="oo", y="Relative oocyst growth", data=oocys_df,inner='box',
                    order=['Female','Male','Both sexes'],color='#FFFFFF') ## 'quartile' 'box'
# sns.swarmplot(y="Relative oocyst growth",
#                 x="oo",data=oocys_df,edgecolor="gray",size=4,color='#FC8D62',order=['Female','Male','Both sexes'])
sns.stripplot(y="Relative oocyst growth",
                x="oo",data=oocys_df,edgecolor="gray",size=8,color='#FC8D62',order=['Female','Male','Both sexes'])

#sns.stripplot(x="Gametocyte", y="Relative female fertility", data=female_tmp,palette=color_gam_dict, hue='Phenotypes',edgecolor="gray",size=8)


# figure = sns_plot.get_figure()
plt.legend().remove()
plt.tick_params(axis='x', rotation=0)
# plt.xlabel('Gametocyte', fontsize=20,labelpad=10)
plt.xlabel('')
plt.ylabel('')
plt.ylim((-18,8))
plt.savefig(folder+'occyst_pheno.pdf',
            format='pdf',dpi=300)

plt.close()

df_O=compute_ttest(res_df,col='O RGR')

two_df=df_MG.merge(df_FG, how='left', on=['Condition 1','Condition 2'])
final_df=two_df.merge(df_O, how='left', on=['Condition 1','Condition 2'])
final_df.to_csv(folder+'t_stats.txt',sep='\t')

import pdb; pdb.set_trace()

# sns.stripplot(x="Sex", y="O RGR", data=tmp,palette=colors, edgecolor="gray")
#
# sns_plot = sns.violinplot(x="Sex", y="O RGR", data=tmp,inner='box',
#                     order=['Female','Male','Female and male'],palette=colorv) ## 'quartile' 'box'
# figure = sns_plot.get_figure()
# figure.savefig(folder+'occyst_pheno.pdf', dpi=300)










###




# ## For female gametocyte
#
# plt.clf()
#
# df_FG=compute_ttest(res_df,col='FG RGR')
# colors=['#F7F7F7', '#CCCCCC', '#969696']
# sns.stripplot(x="Sex", y="FG RGR", data=res_df,color="black", edgecolor="gray")
#
# sns_plot = sns.violinplot(x="Sex", y="FG RGR", data=res_df,inner='box',
#                     order=['Female','Male','Female and male'],palette=colors) ## 'quartile' 'box'
# figure = sns_plot.get_figure()
# figure.savefig(folder+'FemaleG.pdf', dpi=300)
#
#
# ## For female gametocyte
#
# plt.clf()
#
# df_MG=compute_ttest(res_df,col='MG RGR')
# # colors=['#F7F7F7', '#CCCCCC', '#969696']
# sns.stripplot(x="Sex", y="MG RGR", data=res_df,color="black", edgecolor="gray")
#
# sns_plot = sns.violinplot(x="Sex", y="MG RGR", data=res_df,inner='box',
#                     order=['Female','Male','Female and male'],palette=colors) ## 'quartile' 'box'
# figure = sns_plot.get_figure()
# figure.savefig(folder+'MaleG.pdf', dpi=300)
#
#
# ## for occyst
# plt.clf()
#
# df_O=compute_ttest(res_df,col='O RGR')
# # colors=['#F7F7F7', '#CCCCCC', '#969696']
# sns.stripplot(x="Sex", y="O RGR", data=res_df,color="black", edgecolor="gray")
#
# sns_plot = sns.violinplot(x="Sex", y="O RGR", data=res_df,inner='box',
#                     order=['Female','Male','Female and male'],palette=colors) ## 'quartile' 'box'
# figure = sns_plot.get_figure()
# figure.savefig(folder+'occyst.pdf', dpi=300)
#
# two_df=df_MG.merge(df_FG, how='left', on=['Condition 1','Condition 2'])
# final_df=two_df.merge(df_O, how='left', on=['Condition 1','Condition 2'])
# final_df.to_csv(folder+'t_stats.txt',sep='\t')
