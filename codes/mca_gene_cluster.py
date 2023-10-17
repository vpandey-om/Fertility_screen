import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import hypergeom

pheno_names=['Not reduced','Reduced',]
pheno_colors=['#66C2A5','#FC8D62']

kNN_clusters=[4,5,19,1,15,12,8,13,16,7,20,11,3,14,18,17,9,10,2,6]
kNN_cluster_colors=["#3CB371", "#FFA500", "#EEC900", "#32CD32", "#98FB98",
        "#00CDCD", "#6495ED", "#8B008B", "#5F9EA0","#1E90FF", "#FA8072", "#00688B",
         "#A020F0", "#CDB5CD", "#20B2AA", "#0000FF", "#FFC0CB", "#FF0000","#A9A9A9", "#79CDCD"]

def get_data():
    # get x, y co-ordinates for the and mca gene cluster
    #mca_df=pd.read_excel('/Users/vpandey/projects/githubs/Fertility_screen/data/ginny_gene_mca.xlsx', sheet_name='mca')
    #fertility_df=pd.read_excel('/Users/vpandey/projects/githubs/Fertility_screen/preFinals/Phenotype_call_final_100621.xlsx',sheet_name='data')
    mca_df=pd.read_excel(snakemake.input[0], sheet_name='mca')
    fertility_df=pd.read_excel(snakemake.input[1],sheet_name='data')


    female_df=fertility_df[~(fertility_df['GCKO2_pheno_new']=='No data')]
    male_df=fertility_df[~(fertility_df['g145480_pheno_new']=='No data')]
    mca_df['Cluster']=mca_df['Cluster'].astype('category')
    male_mca_df=mca_df.copy()
    female_mca_df=mca_df.copy()
    male_mca_df=mca_df.copy()
    female_mca_df=mca_df.copy()

    ##
    male_mca_df['pheno']='NA'
    female_mca_df['pheno']='NA'
    male_mca_df['rgr']=np.nan
    female_mca_df['rgr']=np.nan


    for idx in male_mca_df.index:
        tmp=male_df[male_df['pbanka_id']==male_mca_df.loc[idx,'feature_symbol']]
        if not tmp.empty:
            male_mca_df.loc[idx,'pheno']=tmp['g145480_pheno_new'].to_list()[0]
            male_mca_df.loc[idx,'rgr']=tmp['g145480_RGR'].to_list()[0]
            male_mca_df.loc[idx,'sd']=tmp['g145480_sd'].to_list()[0]
        tmp=female_df[female_df['pbanka_id']==female_mca_df.loc[idx,'feature_symbol']]
        if not tmp.empty:
            female_mca_df.loc[idx,'pheno']=tmp['GCKO2_pheno_new'].to_list()[0]
            female_mca_df.loc[idx,'rgr']=tmp['GCKO2_RGR'].to_list()[0]
            female_mca_df.loc[idx,'sd']=tmp['GCKO2_sd'].to_list()[0]

    return mca_df,male_df,female_df,male_mca_df,female_mca_df


def plot_mca_cluster(df,cluster,colors):
    ''' plot all cluster with given colors '''
    df['cluster_color']=df['Cluster'].copy()
    color_dict = dict(zip(cluster, colors))
    df['cluster_color']=df['cluster_color'].replace(color_dict)
    # df.plot(kind='scatter',x='name',y='num_children',ax=ax)
    f, ax = plt.subplots(figsize=[6.4*1.2,6.4*1.2])
    handels=[]
    for i,c in enumerate(cluster):
        tmp=df[df['Cluster']==c]
        # x_m=tmp['knn_graph_X'].mean()
        # y_m=tmp['knn_graph_Y'].mean()
        points = ax.scatter(x=tmp['knn_graph_X'], y=tmp['knn_graph_Y'], c=tmp['cluster_color'], s=6, label=str(c))

        # plt.text(x_m,y_m, str(c), size=12)
        handels.append(points)
    # f.colorbar(points)
    plt.legend(handles=handels,title='Clusters', bbox_to_anchor=(1.01, 1))
    # plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
    plt.grid(False)
    plt.savefig("/Users/vpandey/projects/githubs/Fertility_screen_2/mca_figures/Cluster_mca.pdf")
    plt.close()

def plot_phenotype_on_mca(df,color_dict,filename):
    print('plot phenotype on mca cluster')
    df['cluster_color']=df['pheno'].replace(color_dict)
    plot_order=['NA','Not reduced','Reduced']
    # df.plot(kind='scatter',x='name',y='num_children',ax=ax)
    f, ax = plt.subplots(figsize=[6.4*1.2,6.4*1.2])
    handels=[]
    for p in plot_order:
        tmp=df[df['pheno']==p]
        # x_m=tmp['knn_graph_X'].mean()
        # y_m=tmp['knn_graph_Y'].mean()
        if p=='Reduced':
            points = ax.scatter(x=tmp['knn_graph_X'], y=tmp['knn_graph_Y'], c=tmp['cluster_color'], s=6, label=p)
        else:
            points = ax.scatter(x=tmp['knn_graph_X'], y=tmp['knn_graph_Y'], c=tmp['cluster_color'], s=6)


        # plt.text(x_m,y_m, str(c), size=12)
        handels.append(points)
    # f.colorbar(points)
    plt.legend(handles=handels,title='Phenotypes', bbox_to_anchor=(1.01, 1))
    # plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
    plt.grid(False)
    plt.savefig(filename)
    plt.close()


def plot_phenotype_on_mca_male_female(df,color_dict,filename,plot_order=['NA','Male', 'Female','Male and Female' ]):
    print('plot phenotype on mca cluster')
    df['cluster_color']=df['category'].replace(color_dict)
    # df.plot(kind='scatter',x='name',y='num_children',ax=ax)
    f, ax = plt.subplots(figsize=[6.4*1.2,6.4*1.2])
    handels=[]
    for p in plot_order:
        tmp=df[df['category']==p]
        # x_m=tmp['knn_graph_X'].mean()
        # y_m=tmp['knn_graph_Y'].mean()
        if not p=='NA':
            points = ax.scatter(x=tmp['knn_graph_X'], y=tmp['knn_graph_Y'], c=tmp['cluster_color'], s=6, label=p)
        else:
            points = ax.scatter(x=tmp['knn_graph_X'], y=tmp['knn_graph_Y'], c=tmp['cluster_color'], s=6)


        # plt.text(x_m,y_m, str(c), size=12)
        handels.append(points)
    # f.colorbar(points)
    plt.legend(handles=handels,title='Groups', bbox_to_anchor=(1.01, 1))
    # plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
    plt.grid(False)
    plt.savefig(filename)
    plt.close()


def plot_phenotype_one_sex(df,color_dict,filename,plot_order=['NA','Male', 'Female','Male and Female' ]):
    print('plot phenotype on mca cluster')
    df['cluster_color']=df['category'].replace(color_dict)
    # df.plot(kind='scatter',x='name',y='num_children',ax=ax)
    f, ax = plt.subplots(figsize=[6.4*1.2,6.4*1.2])
    handels=[]
    for p in plot_order:
        tmp=df[df['category']==p]
        # x_m=tmp['knn_graph_X'].mean()
        # y_m=tmp['knn_graph_Y'].mean()
        if not p=='NA':
            points = ax.scatter(x=tmp['knn_graph_X'], y=tmp['knn_graph_Y'], c=tmp['cluster_color'], s=6, label=p)
        else:
            points = ax.scatter(x=tmp['knn_graph_X'], y=tmp['knn_graph_Y'], c=tmp['cluster_color'], s=6)


        # plt.text(x_m,y_m, str(c), size=12)
        handels.append(points)
    # f.colorbar(points)
    plt.legend(handles=handels,title='Groups', bbox_to_anchor=(1.01, 1))
    # plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
    plt.grid(False)
    plt.savefig(filename)
    plt.close()


def enrich_mca_cluster(mca_df,screen_genens,listgenes,types):
    ''' Pefrorming cluster enrichment bassed on screen    '''

    clust=[i+1 for i in range(20)]

    for k in range(len(types)):
        pvals=[]
        N_list=[]
        x_list=[]
        genes=[]
        clusters=[]

        for c in clust:

            tmp=mca_df[mca_df['Cluster']==c]
            clusters.append(str(c))
            clust_genes=set(tmp['feature_symbol'].to_list())
            clust_genes_in_screen= set(screen_genens)& clust_genes
            clust_genes_in_pheno=set(listgenes[k])& clust_genes

            M=len(set(screen_genens))
            N=len(clust_genes_in_screen)
            n=len(set(listgenes[k]))
            rv = hypergeom(M, n, N)
            x=len(clust_genes_in_pheno)
            pval= 1-rv.cdf(x)
            pvals.append(pval)
            N_list.append(N)
            x_list.append(x)
            genes.append(', '.join(clust_genes_in_pheno))



        enrich={'Cluster':clusters, 'pvals':pvals,'cluster_size_in_screen':N_list,'# of genes found in phenotype':x_list,'genes':genes}
        enrich_df=pd.DataFrame.from_dict(enrich)
        enrich_soreted=enrich_df.sort_values(by=['pvals'],ascending=True)
        enrich_soreted.to_csv('/Users/vpandey/projects/githubs/Fertility_screen_2/mca_figures/Cluster_enrichment_with_phenotype_'+ types[k] +'.txt',sep='\t',index=None)



def make_male_female(mdf,fdf):
    final_df=mdf.copy()
    final_df['category']='NA'
    both_df=mdf[(mdf['pheno']=='Reduced') & (fdf['pheno']=='Reduced')]
    final_df.loc[both_df.index,'category']='Male and Female'
    male_df=mdf[(mdf['pheno']=='Reduced') & (~ (fdf['pheno']=='Reduced'))]
    final_df.loc[male_df.index,'category']='Male'
    female_df=mdf[(~(mdf['pheno']=='Reduced')) & (fdf['pheno']=='Reduced')]
    final_df.loc[female_df.index,'category']='Female'
    return final_df


def plot_male_female(mdf,fdf):
    ''' we can combine male and female both '''
    print('combining both')
    final_df=make_male_female(mdf,fdf)
    filename="/Users/vpandey/projects/githubs/Fertility_screen_2/mca_figures/Male_female_mca.pdf"
    color_dict={'NA':'#f5f5f5','Male and Female':'#FFC61E', 'Male':'#009ADE', 'Female':'#AF58BA'}
    plot_phenotype_on_mca_male_female(final_df,color_dict,filename)


def plot_male_female_both(df,male_only,female_only,male_female):
    ''' we can combine male and female both '''
    print('combining both')
    final_df=df.copy()
    final_df['category']='NA'
    tmp=final_df[final_df['feature_symbol'].isin(male_only)]
    final_df.loc[tmp.index,'category']='Male'
    tmp=final_df[final_df['feature_symbol'].isin(female_only)]
    final_df.loc[tmp.index,'category']='Female'
    tmp=final_df[final_df['feature_symbol'].isin(male_female)]
    final_df.loc[tmp.index,'category']='Male and Female'
    xx=final_df[final_df['feature_symbol'].isin(male_only)]
    print(set(male_only)-set(xx['feature_symbol'].to_list()))
    filename="/Users/vpandey/projects/githubs/Fertility_screen_2/mca_figures/Male_female_mca.pdf"
    # color_dict={'NA':'#f5f5f5','Male and Female':'#FFC61E', 'Male':'#009ADE', 'Female':'#AF58BA'}
    color_dict={'NA':'#f5f5f5','Male and Female':'#FFC61E', 'Male':'#009ADE', 'Female':'#AF58BA'}
    plot_phenotype_on_mca_male_female(final_df,color_dict,filename)






mca_df,male_df,female_df,male_mca_df,female_mca_df=get_data()


## enrichment analysis
# male_female_df=make_male_female(male_mca_df,female_mca_df)
# male_female_genes=male_female_df[male_female_df['category']=='Male and Female']['feature_symbol'].to_list()
# male_genes=male_female_df[male_female_df['category']=='Male']['feature_symbol'].to_list()
# female_genes=male_female_df[male_female_df['category']=='Female']['feature_symbol'].to_list()
# folder='/Users/vpandey/projects/githubs/Fertility_screen_2/preFinals/'

# df=pd.read_csv('/Users/vpandey/projects/githubs/Fertility_screen/preFinals/female_only_genes.txt',sep='\t',header=None)
# female_only=df[0].tolist()
#
#
# df=pd.read_csv('/Users/vpandey/projects/githubs/Fertility_screen/preFinals/male_only_genes.txt',sep='\t',header=None)
# male_only=df[0].tolist()
#
# df=pd.read_csv('/Users/vpandey/projects/githubs/Fertility_screen/preFinals/female_male_genes.txt',sep='\t',header=None)
# male_female=df[0].tolist()

df=pd.read_csv(snakemake.input[3],sep='\t',header=None)
female_only=df[0].tolist()


df=pd.read_csv(snakemake.input[2],sep='\t',header=None)
male_only=df[0].tolist()

# df=pd.read_csv('/Users/vpandey/projects/githubs/Fertility_screen/preFinals/female_male_genes.txt',sep='\t',header=None)
# male_female=df[0].tolist()
#
# screen_genens=list(set(male_df['pbanka_id'].to_list()) | set(female_df['pbanka_id'].to_list()))
#
# listgenes=[male_only,female_only,male_female]
# types=['Male','Female','Male_and_Female']
#
# # plot_male_female_both(mca_df,male_only,female_only,male_female)
#
#
# enrich_mca_cluster(mca_df,screen_genens,listgenes,types)
# 'Male':'#009ADE', 'Female':'#AF58BA'
### male specific
final_df=mca_df.copy()
final_df['category']='NA'
tmp=final_df[final_df['feature_symbol'].isin(male_only)]
final_df.loc[tmp.index,'category']='Male'
# plot_mca_cluster(mca_df,kNN_clusters,kNN_cluster_colors)
#filename="/Users/vpandey/projects/githubs/Fertility_screen_2/mca_figures/Cluster_mca_male_reduced.pdf"
filename=snakemake.output[0]
color_dict={'NA':'#f5f5f5', 'Male':'#009ADE'}
# color_dict={'NA':'#f5f5f5', 'Not reduced':'#f5f5f5', 'Reduced':'#FC8D62'}
plot_phenotype_one_sex(final_df,color_dict,filename,plot_order=['NA','Male' ])



### male specific
final_df=mca_df.copy()
final_df['category']='NA'
tmp=final_df[final_df['feature_symbol'].isin(female_only)]
final_df.loc[tmp.index,'category']='Female'
# plot_mca_cluster(mca_df,kNN_clusters,kNN_cluster_colors)
#filename="/Users/vpandey/projects/githubs/Fertility_screen_2/mca_figures/Cluster_mca_female_reduced.pdf"
filename=snakemake.output[1]
# color_dict={'NA':'#f5f5f5', 'Not reduced':'#f5f5f5', 'Reduced':'#FC8D62'}
color_dict={'NA':'#f5f5f5', 'Female':'#AF58BA'}
plot_phenotype_one_sex(final_df,color_dict,filename,plot_order=['NA','Female' ])





# #color_dict={'NA':'#f5f5f5', 'Not reduced':'#f5f5f5', 'Reduced':'#AF58BA'}
#
# # plot_phenotype_on_mca(female_mca_df,color_dict,filename)
# filename="/Users/vpandey/projects/githubs/Fertility_screen_2/mca_figures/Cluster_mca_male_reduced.pdf"
# color_dict={'NA':'#f5f5f5', 'Not reduced':'#f5f5f5', 'Reduced':'#FC8D62'}
# #color_dict={'NA':'#f5f5f5', 'Not reduced':'#f5f5f5', 'Reduced':'#009ADE'}
# plot_phenotype_on_mca(male_mca_df,color_dict,filename)
#
#
#
