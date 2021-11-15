import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import itertools
import numpy as np
import pickle

from scipy.stats import hypergeom ## for hypergeomtric distributions
import statsmodels.stats.multitest as mtest ## for multiple hypothesis testing

### import our final xlsx file


def test_corr_cluster(df,col='g145480_RGR',outfile='test.pdf'):
    ## get clusters by groupby
    time=['6s', '12s', '18s', '60s']

    clusters=np.arange(len(df['cluster'].unique()))
    dictIdx=df.groupby(['time','cluster']).indices
    # tdf=pd.DataFrame(index=time,columns=clusters)
    corrs=np.zeros((len(time),len(clusters)))
    pvals=np.zeros((len(time),len(clusters)))
    geneNum=np.zeros((len(time),len(clusters)))
    logpvals=np.zeros((len(time),len(clusters)))
    for i,t in enumerate(time):
        for j,cidx in enumerate(clusters):
            if (t,cidx+1) in dictIdx.keys():
                tmp=df.iloc[dictIdx[(t,cidx+1)],:].copy()
                ## get correlations
                (corr,pval)=stats.spearmanr(tmp[col].values,tmp['phospho_data'].values )
                corrs[i,j]=corr
                pvals[i,j]=pval
                geneNum[i,j]=tmp.shape[0]
                if pval<=0.05:
                    logpvals[i,j]=1
                else:
                    logpvals[i,j]=0
                # tdf.iloc[i,j]=corr
    plt.clf()

    #sns_plot=sns.heatmap(corrs.T,xticklabels=time,yticklabels=clusters,center=0,annot=logpvals.T)
    sns_plot=sns.heatmap(corrs.T,xticklabels=time,yticklabels=clusters,center=0,annot=geneNum.T)
    plt.xlabel('time')
    plt.ylabel('Cluster number')
    figure = sns_plot.get_figure()
    figure.savefig(outfile, dpi=300)



def investigate_phospho_data(sexdf,col,sex_only,M,time,fold=2,outfile='test.xlsx',outfile1='test.pdf'):
    print('compute up and downregulated as well as phenotypes')
    ## groupby
    reqDict=sexdf.groupby([col+'_pheno_new','time']).indices
    reqDictTime=sexdf.groupby(['time']).indices
    up_genes=[]
    down_genes=[]
    pvals_up=[]
    pvals_down=[]
    up_genes_num=[]
    down_genes_num=[]
    up_genes_phosphosite=[]
    down_genes_phosphosite=[]
    up_genes_with_phospho_data=[]
    down_genes_with_phospho_data=[]
    up_corr=[]
    down_corr=[]
    for t in time:
        if ('Reduced',t) in reqDict.keys():
            ## get data
            tmp=sexdf.iloc[reqDict[('Reduced',t)],:].copy()
            # get upregulated genes

            tmp_up=tmp[tmp['phospho_data']>=np.log2(fold)]
            # get sex only
            tmp_up=tmp_up[tmp_up['pbanka_id'].isin(sex_only)]
            # get downregulated genes
            tmp_down=tmp[tmp['phospho_data']<=np.log2(1/fold)]
            tmp_down=tmp_down[tmp_down['pbanka_id'].isin(sex_only)]

            ## compute pvalues

            ## how many are in phenoype
            # M=1220 ## total genes
            N=len(sex_only) ## sex-specific genes
            tmp2=sexdf.iloc[reqDictTime[t],:].copy()
            n=len(tmp2['pbanka_id'].unique())#len(clust_genes) ## total number of cluster genes

            rv = hypergeom(M, n, N) ## hypergeometric  distributions
            x=len(tmp_up['pbanka_id'].unique()) ## number of sex-specific genes in cluster
            ## for p-value we need to get summation of all probalities greater than > x
            up_genes_num.append(x)
            pval_up= 1-rv.cdf(x)
            pvals_up.append(pval_up)
            x=len(tmp_down['pbanka_id'].unique())

            pval_down= 1-rv.cdf(x)
            pvals_down.append(pval_down)
            ## get p-value
            down_genes_num.append(x)
            up_genes_phosphosite.append(tmp_up.shape[0])
            down_genes_phosphosite.append(tmp_down.shape[0])




            ## correlation analysis for phospho_data

            (corr,pval)=stats.spearmanr(tmp_up[col+'_RGR'].values,tmp_up['phospho_data'].values )
            up_corr.append(corr)
            (corr,pval)=stats.spearmanr(tmp_down[col+'_RGR'].values,tmp_down['phospho_data'].values )
            down_corr.append(corr)
            xx=tmp_up['pbanka_id'].value_counts()
            yy=xx.to_string()
            up_genes_with_phospho_data.append(yy.replace('\n','|'))

            xx=tmp_down['pbanka_id'].value_counts()
            yy=xx.to_string()
            down_genes_with_phospho_data.append(yy.replace('\n','|'))

            up_genes.append((tmp_up.shape[0],len(tmp_up['pbanka_id'].unique()),tmp_up['pbanka_id'].value_counts(),tmp_up))
            down_genes.append((tmp_down.shape[0],len(tmp_down['pbanka_id'].unique()),tmp_down['pbanka_id'].value_counts(),tmp_down))

    ## result
    result={'Number of unique upregulated genes':up_genes_num,
    'Number of upregulated phosphosites':up_genes_phosphosite,
    'p-values for upregualted genes':pvals_up,
    'Number of unique downregulated genes':down_genes_num,
    'Number of downregualted phosphosites':down_genes_phosphosite,
    'p-values for downregualted genes':pvals_down,
    'Number of phosphosites per upregualted gene':up_genes_with_phospho_data,
    'Number of phosphosites per downregualted gene':down_genes_with_phospho_data}
    ## create enrichment table
    result_df=pd.DataFrame.from_dict(result)
    result_df.index=time
    result_df.to_excel(outfile,sheet_name='data')
    # sns heatmap
    A = np.zeros((4, 2))
    A[:,0]=-np.log10(pvals_up)
    A[:,1]=-np.log10(pvals_down)

    A2 = np.zeros((4, 2))
    A2[:,0]=pvals_up
    A2[:,1]=pvals_down


    sns_plot=sns.heatmap(A,xticklabels=['up','down'],yticklabels=time,cbar_kws={'label': '-log10(pval)'},annot=A2)
    plt.ylabel('time')
    plt.xlabel('Phospho regulation')
    figure = sns_plot.get_figure()
    figure.savefig(outfile1, dpi=300)

    return A2








def apply_cluster_membership(pathway_to_genes,targetList,M,cluster_size,outfile,tmpdf=None):
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
    n1_list=[] # store actual number of genes in cluster
    data=[]

    ## group clusters
    for k,v in pathway_to_genes.items():
        ## print(k, len(v))
        pathways.append(k) ## store cluster number
        clust_genes=set(v) ## genes in the cluster

        sex_gene_in_clust= set(targetList)& clust_genes ## how many sex-specific genes are found in cluster
        xx=tmpdf[tmpdf['pbanka_id'].isin(sex_gene_in_clust)]

        ## how many are in phenoype
        M=M ## total genes
        N=len(targetList) ## sex-specific genes
        n=len(clust_genes) ## total number of cluster genes
        n1_list.append(len(cluster_size[k]))

        rv = hypergeom(M, n, N) ## hypergeometric  distributions
        x=len(sex_gene_in_clust) ## number of sex-specific genes in cluster

        ## for p-value we need to get summation of all probalities greater than > x
        pval= 1-rv.cdf(x)
        pvals.append(pval)
        n_list.append(n)
        x_list.append(x)
        genes.append(', '.join(sex_gene_in_clust))
        data.append(xx.to_string())
    fdr=mtest.multipletests(pvals, alpha=0.05, method='fdr_bh')
    enrich={'Cluster':pathways, 'Pvals':pvals,'Cluster size (number of unique genes in phospho data)':n1_list,'Cluster size ( overlapping genes with fertility data)':n_list,'Number of phenotype genes in cluster':x_list,'FDR':fdr[1],'Genes':genes,'Data':data}
    ## create enrichment table
    enrich_df=pd.DataFrame.from_dict(enrich)
    enrich_soreted=enrich_df.sort_values(by=['Pvals'],ascending=True)
    # enrich_soreted.to_csv(outfile,sep='\t',index=None)
    enrich_soreted.to_excel(outfile,sheet_name='enrichment',index=None)
    return enrich_soreted

fertility_df=pd.read_excel('/Users/vpandey/projects/githubs/Fertility_screen/analysisCodes/Phenotype_call_final_100621.xlsx',sheet_name='data')

## import proteomics file data
phospho_data=pd.read_excel('/Users/vpandey/projects/githubs/Fertility_screen/analysisCodes/Invergo_et_al.xlsx',sheet_name='data')
## see comparisons

## id conversion
prev_to_new=pickle.load(open( '/Users/vpandey/projects/githubs/Fertility_screen/data/prevTonew_PBANKA.pickle','rb'))

new_ids=[]
for ind in phospho_data.index:
    protein_id=phospho_data.loc[ind,'protein']
    if protein_id in prev_to_new.keys():
        protein_id=prev_to_new[protein_id]
    new_ids.append(protein_id)
tmp_df=phospho_data.copy()
tmp_df['pbanka_id']=new_ids

combined_df=fertility_df.merge(tmp_df, how='inner', on=['pbanka_id'])
#combined_df.to_excel('/Users/vpandey/projects/githubs/Fertility_screen/analysisCodes/Merged_phospho.xlsx',sheet_name='data')
# get cluster to genes

selcted_df=combined_df[['pbanka_id','gene_product','GCKO2_RGR','g145480_RGR','GCKO2_pheno_new','g145480_pheno_new','0s', '6s', '12s', '18s', '60s','cluster']].copy()


############### cluster membership analysis
MaleData=fertility_df[~(fertility_df['g145480_pheno_new']=='No data')]
cluster_genes={}
cluster_size={}
for k,v in tmp_df.groupby(['cluster']).indices.items():
    s1=set(tmp_df.loc[v,'pbanka_id'].to_list())
    s2=set(MaleData['pbanka_id'].to_list())
    pbanka_ids=list(s1 & s2)
    cluster_size[k]=list(s1)
    cluster_genes[k]=pbanka_ids

M=MaleData.shape[0]
# apply enrichment for male candidate
tmp=fertility_df[fertility_df['g145480_pheno_new']=='Reduced']
# genelist=tmp['pbanka_id'].to_list()
input_df=pd.read_csv("/Users/vpandey/projects/githubs/Fertility_screen/preFinals/male_only_genes.txt",sep= "\t",header=None)
genelist1=input_df.loc[:,0].to_list()

outfile="/Users/vpandey/projects/githubs/Fertility_screen/preFinals/male_phospho_cluster_membership.xlsx"
outdf=apply_cluster_membership(cluster_genes,genelist1,M,cluster_size,outfile,tmp_df)

# apply enrichment for female

FeMaleData=fertility_df[~(fertility_df['GCKO2_pheno_new']=='No data')]
cluster_genes={}
cluster_size={}
for k,v in tmp_df.groupby(['cluster']).indices.items():
    s1=set(tmp_df.loc[v,'pbanka_id'].to_list())
    s2=set(FeMaleData['pbanka_id'].to_list())
    pbanka_ids=list(s1 & s2)
    cluster_size[k]=list(s1)
    cluster_genes[k]=pbanka_ids

M=FeMaleData.shape[0]


tmp=fertility_df[fertility_df['GCKO2_pheno_new']=='Reduced']
#genelist=tmp['pbanka_id'].to_list()
input_df=pd.read_csv("/Users/vpandey/projects/githubs/Fertility_screen/preFinals/female_only_genes.txt",sep= "\t",header=None)
genelist2=input_df.loc[:,0].to_list()
outfile="/Users/vpandey/projects/githubs/Fertility_screen/preFinals/female_phospho_cluster_membership.xlsx"
outdf=apply_cluster_membership(cluster_genes,genelist2,M,cluster_size,outfile,tmp_df)



#######

## tranform
tr_df=selcted_df.copy()
cols=['6s', '12s', '18s', '60s']

phos_time=[]
phos_data=[]
for col in cols:
    phos_time.append([col]*tr_df.shape[0])
    phos_data.append(selcted_df[col].values)
phos_time=np.array(phos_time)
phos_data=np.array(phos_data)
phos_time=phos_time.flatten()
phos_data=phos_data.flatten()

tr_df['phospho_data']=np.nan
tr_df['time']=np.nan

tr_df2=pd.concat([tr_df]*len(cols), ignore_index=True)
tr_df2['phospho_data']=phos_data
tr_df2['time']=phos_time

## remove nan
male_df=tr_df2[(~tr_df2['g145480_RGR'].isna())& (~tr_df2['phospho_data'].isna()) & (~(tr_df2['g145480_pheno_new']=='No data'))].copy()
female_df=tr_df2[(~tr_df2['GCKO2_RGR'].isna())& (~tr_df2['phospho_data'].isna()) & (~(tr_df2['GCKO2_pheno_new']=='No data'))].copy()


# ## test correlation
# male_df_redcued=male_df[male_df['g145480_pheno_new']=='Reduced']
# outpdf='/Users/vpandey/projects/githubs/Fertility_screen/analysisCodes/cluster_corr_male_reduced.pdf'
# test_corr_cluster(male_df,col='g145480_RGR',outfile=outpdf)
# import pdb; pdb.set_trace()
# plt.clf()
# sns_plot=sns.relplot(data=male_df, x="g145480_RGR", y="phospho_data", hue="time", row="cluster")
# sns_plot.savefig("output.png")



## get up regulated and downregulated phosphosites genes with time
outfile="/Users/vpandey/projects/githubs/Fertility_screen/preFinals/male_fertility_role_analysis.xlsx"
outfile1="/Users/vpandey/projects/githubs/Fertility_screen/preFinals/male_fertility_role_analysis.pdf"
time=['6s', '12s', '18s', '60s']
fold=2
plt.clf()
male_pvals=investigate_phospho_data(male_df,'g145480',genelist1,MaleData.shape[0],time,fold,outfile,outfile1)
outfile="/Users/vpandey/projects/githubs/Fertility_screen/preFinals/female_fertility_role_analysis.xlsx"
outfile1="/Users/vpandey/projects/githubs/Fertility_screen/preFinals/female_fertility_role_analysis.pdf"
plt.clf()
female_pvals=investigate_phospho_data(female_df,'GCKO2',genelist2,FeMaleData.shape[0],time,fold,outfile,outfile1)

A=np.concatenate((male_pvals,female_pvals),axis=1)
A2=-np.log10(A)
## combined male and female
plt.clf()
sns_plot=sns.heatmap(A2,xticklabels=['male_up','male_down','female_up','female_down'],yticklabels=time,cbar_kws={'label': '-log10(pvalue)'},annot=A)
plt.ylabel('time')
plt.xlabel('Sex-specific phosphosite regulation')
outfile1="/Users/vpandey/projects/githubs/Fertility_screen/preFinals/male_female_fertility_role_analysis.pdf"
figure = sns_plot.get_figure()
figure.savefig(outfile1, dpi=300)
import pdb; pdb.set_trace()




plt.clf()
sns_plot=sns.scatterplot(data=male_df, x='g145480_RGR', y="phospho_data", hue="time")
figure = sns_plot.get_figure()
figure.savefig('/Users/vpandey/projects/githubs/Fertility_screen/analysisCodes/male_phospho.pdf', dpi=300)

plt.clf()
sns_plot=sns.scatterplot(data=female_df, x='GCKO2_RGR', y="phospho_data", hue="time")
figure = sns_plot.get_figure()
figure.savefig('/Users/vpandey/projects/githubs/Fertility_screen/analysisCodes/female_phospho.pdf', dpi=300)

##

import pdb; pdb.set_trace()
