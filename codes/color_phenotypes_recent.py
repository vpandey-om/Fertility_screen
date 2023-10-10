import pickle
import pandas as pd

## read excel file
import os
import sys
import plotly.express as px
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np
import statsmodels.stats.multitest as mtest
#import scanpy as sc
import math
code=os.getcwd()
upLevel=code.replace('codes','') ####### we are going to upper level of code directory
sys.path.insert(0,upLevel+'/data')
sys.path.insert(1, upLevel+'/Figures')

import matplotlib.pyplot as plt

import matplotlib
import matplotlib.cm as cm
from copy import copy
from matplotlib.colors import rgb2hex
from scipy.stats import hypergeom


# ['#687847', '#979b8d', '#706d59', '#6a734d', '#8b5e66', '#5e5051']




# print(sys.path)


# input files which will be needed for analysis of cliare screen data

data_folder=sys.path[0]

## output folder where we will write figures and output files
out_folder=sys.path[1]
prev_to_new=pickle.load(open(data_folder+'/prevTonew_PBANKA.pickle','rb'))
cross_df=pd.read_excel(open(data_folder+'/cross_pheno.xlsx','rb'), sheet_name='cross')

#import pdb; pdb.set_trace()
# ID conversion: we used plasmoDB to convert ID for P. Berghai
# mca_df=pd.read_excel(open(data_folder+'/ginny_gene_mca.xlsx','rb'), sheet_name='mca')
# obs=mca_df[['feature_symbol','product_description','Cluster']]
# obs.set_index('feature_symbol',inplace=True)
#
# fertility_df=pd.read_excel(open(data_folder+'/Phenotype_call_plasmogem_id_Vikash_271020.xlsx','rb'),sheet_name='fertility')
# female_df=fertility_df[~(fertility_df['GCKO2_pheno']=='No Data')]
# male_df=fertility_df[~(fertility_df['g145480_pheno']=='No Data')]
#
# var=pd.DataFrame(index=['knn_graph_X','knn_graph_Y'],columns=['co-ordinates'])
# var.loc['knn_graph_X','co-ordinates']='x'
# var.loc['knn_graph_Y','co-ordinates']='y'
# # adata = sc.AnnData(X=mca_df[['knn_graph_X','knn_graph_Y']].to_numpy(),obs=obs,var=var)
# mca_df['Cluster']=mca_df['Cluster'].astype('category')
#
# male_mca_df=mca_df.copy()
# female_mca_df=mca_df.copy()
#
# ##
# male_mca_df['pheno']='NA'
# female_mca_df['pheno']='NA'
# male_mca_df['rgr']=np.nan
# female_mca_df['rgr']=np.nan
#
#
# for idx in male_mca_df.index:
#     tmp=male_df[male_df['pbanka_id']==male_mca_df.loc[idx,'feature_symbol']]
#     if not tmp.empty:
#         male_mca_df.loc[idx,'pheno']=tmp['g145480_pheno'].to_list()[0]
#         male_mca_df.loc[idx,'rgr']=tmp['g145480_RGR'].to_list()[0]
#     tmp=female_df[female_df['pbanka_id']==female_mca_df.loc[idx,'feature_symbol']]
#     if not tmp.empty:
#         female_mca_df.loc[idx,'pheno']=tmp['GCKO2_pheno'].to_list()[0]
#         female_mca_df.loc[idx,'rgr']=tmp['GCKO2_RGR'].to_list()[0]
#
#
# #import pdb; pdb.set_trace()
# # color gradient

from cobra.io import load_matlab_model


def read_metabolic_genes_subsystem(genelist=None):

### load FBA model
    cobra_model= load_matlab_model(data_folder+'/ipbeblood.mat')
    genes=[g.id.replace('putative_','') for g in cobra_model.genes]
    if genelist ==None:
        genelist=genes
    ### read genes

    subsys_genes={}
    for rxn in cobra_model.reactions:
        if rxn.gene_reaction_rule!="":
            if rxn.subsystem not in subsys_genes.keys():
                subsys_genes[rxn.subsystem]=[]
            rxn_genes=[item.id.replace('putative_','') for item in rxn.genes]
            for item in rxn_genes:
                subsys_genes[rxn.subsystem].append(item)


    subsys_df=pd.DataFrame(index=subsys_genes.keys())

    for subsys,gs in subsys_genes.items():
        comm=set(genelist)&set(gs)
        if len(comm)>0:
            subsys_df.loc[subsys,'num']=len(comm)
            pheno_genes=','.join(comm)
            subsys_df.loc[subsys,'genes']=pheno_genes
    return genes,subsys_df



#
def get_category_colors(df):
    ''' get colors from dataframe'''
    clusters=[i+1 for i in range(20)]
    vals=[(i+1)/20 for i in range(20)]
    # cmap = matplotlib.cm.get_cmap('plasma')
    #
    #
    # colors=[]
    #
    # for v in vals:
    #     rgba=cmap(v)
    #     hex=rgb2hex(rgba[0:3])
    #     colors.append(hex)
    # sns.palplot(sns.diverging_palette(220, 20, n=7))
    cols=sns.color_palette("husl",len(clusters))
    colors=[rgb2hex(item) for item in cols]
    colors=[
    '#6EAB4B','#A6A7A7','#664992','#61AC75','#EBA43F','#8CC6CA','#5383BD','#718FC6','#EDC0CA','#D33B2F',
    '#2A678A','#64B8C0','#792B7F','#C8B6C8','#ACCB8F','#6D9D9E','#2F4C95','#53AEA9','#E6CA44','#E28278'
    ]


    f, ax = plt.subplots(figsize=[6.4*1.2,4.8*1.2])
    handels=[]
    for i,c in enumerate(clusters):
        tmp=df[df['Cluster']==c]
        x_m=tmp['knn_graph_X'].mean()
        y_m=tmp['knn_graph_Y'].mean()

        points = ax.scatter(x=tmp['knn_graph_X'], y=tmp['knn_graph_Y'], c=colors[i], s=6, label=str(c))

        # plt.text(x_m,y_m, str(c), size=12)
        handels.append(points)
    # f.colorbar(points)
    plt.legend(handles=handels,title='Clusters', bbox_to_anchor=(1.01, 1))
    # plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
    plt.grid(False)
    plt.savefig(out_folder+"/Cluster_mca.pdf")
    plt.close()



    # colors=cmap[df['Cluster']]




# sns_plot=sns.scatterplot(data=mca_df, x="knn_graph_X", y="knn_graph_Y", hue="Cluster",s=6)
# plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
# plt.savefig(out_folder+"/mca_gene_cluster.pdf")
# plt.close()
#
# # Create an array with the colors you want to use
# colorsIdx = {'No Power': '#cccccc', 'Not Reduced': '#b2e2e2','Reduced':'#d8b365'}
# # colors = ["#f7f7f7", "#b2e2e2",'#d8b365','#cccccc']
# colors = ["#f7f7f7", "#99d594",'#fc8d59','#ffffbf']
# # Set your custom color palette
# customPalette = sns.set_palette(sns.color_palette(colors))
#
#
# sns_plot=sns.scatterplot(data=male_mca_df, x="knn_graph_X", y="knn_graph_Y", hue="pheno",s=6)
# plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
# plt.savefig(out_folder+"/mca_gene_cluster_male.pdf")
# plt.close()
#
# sns_plot=sns.scatterplot(data=female_mca_df, x="knn_graph_X", y="knn_graph_Y", hue="pheno",s=6)
# plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
# plt.savefig(out_folder+"/mca_gene_cluster_female.pdf")
# plt.close()



def mask_dataframe_table_plot(df,fname,pheno=-5,no_pheno=-1):
    ## get nan index
    from matplotlib.font_manager import FontProperties

    fontP = FontProperties()
    fontP.set_size('xx-small')
    df_nan=df[df['rgr'].isna()]
    df_extreme=df[df['rgr']<pheno]
    df_under=df[df['rgr']>no_pheno]
    df_remain=df[~((df['rgr']>no_pheno) |(df['rgr']<pheno)  | (df['rgr'].isna())) ]
    ### now we will plot all points
    f, ax = plt.subplots(figsize=[6.4*1.2,4.8*1.2])
    p1=ax.scatter(df_nan["knn_graph_X"], df_nan["knn_graph_Y"], s=6,color='lightgrey', label='No data')
    # p2=ax.scatter(df_extreme["knn_graph_X"], df_extreme["knn_graph_Y"], s=6,color='r', label='RGR < '+str(pheno))
    # p3=ax.scatter(df_under["knn_graph_X"], df_under["knn_graph_Y"], s=6,color='g', label='RGR > '+str(no_pheno))
    p3=ax.scatter(df_under["knn_graph_X"], df_under["knn_graph_Y"], s=6,color='#EDF155', label='RGR > '+str(no_pheno))
    points = ax.scatter(x=df_remain['knn_graph_X'], y=df_remain['knn_graph_Y'], c=df_remain['rgr'], s=6, cmap='plasma')
    # f.colorbar(points)
    p2=ax.scatter(df_extreme["knn_graph_X"], df_extreme["knn_graph_Y"], s=6,color='#0B1283', label='RGR < '+str(pheno))

    # plt.legend(handles=[p1, p2,p3], title='Data labels', bbox_to_anchor=(1.1, 1), loc='upper left', prop=fontP)
    # plt.legend(handles=[p1, p2,p3], title='Data labels', bbox_to_anchor=(1.1, 1), loc='upper left')
    plt.legend(handles=[p1, p2,p3], title='RGR data', loc='best',frameon=True,prop=fontP)
    # points = ax.scatter(x= male_mca_df['knn_graph_X'], y=male_mca_df['knn_graph_Y'], c=colors, s=6, cmap='Spectral')
    #points = ax.scatter(x=df_remain['knn_graph_X'], y=df_remain['knn_graph_Y'], c=df_remain['rgr'], s=6, cmap='RdYlGn')

    # plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
    plt.grid(False)
    plt.savefig(out_folder+"/"+fname)

    plt.close()


    # plt.plot(x2*0.4, y2, 'o-', label='Points removed')
    # plt.plot(x*0.7, y3, 'o-', label='Masked values')
    # plt.plot(x*1.0, y4, 'o-', label='NaN values')
    #











# def get_color(lst):
#     minima = lst.min()
#     maxima = lst.max()
#     norm = matplotlib.colors.Normalize(vmin=minima, vmax=maxima, clip=True)
#     palette = copy(cm.Spectral)
#     # palette.set_over('r', 1.0)
#     # palette.set_under('g', 1.0)
#     palette.set_bad('lightgrey', 0.0)
#     mapper = cm.ScalarMappable(norm=norm, cmap=palette)
#     colors=[]
#     for v in lst:
#         colors.append(mapper.to_rgba(v))
#     return colors

# cmap = copy.copy(mpl.cm.get_cmap("Spectral"))

# palette = copy(plt.cm.Spectral)
# # palette.set_over('r', 1.0)
# # palette.set_under('g', 1.0)
# palette.set_bad('b', 0.0)

# cmap = plt.get_cmap('Spectral')
# cmap.set_bad(color = 'grey', alpha = 0)


def define_much_greater(mdf,fdf,cut_gr=-2,cut_less=-5):

    print('hi')
    # get male fertility genes
    cond_m=((fdf['rgr']>cut_gr) & (mdf['rgr']<cut_less))
    male_fertility_df=mdf[cond_m]
    #male_fertility_genes=male_fertility_df['feature_symbol']
    # female fertility genes
    cond_f=((fdf['rgr']<cut_less) & (mdf['rgr']>cut_gr))
    female_fertility_df=mdf[cond_f]
    #female_fertility_genes=female_fertility_df['feature_symbol']
    ## both
    cond_m_f=((fdf['rgr']<cut_less) & (mdf['rgr']<cut_less))
    both_df=mdf[cond_m_f]
    #both_genes=both_df['feature_symbol']
    cond_neutral=((fdf['rgr']>=cut_gr) & (mdf['rgr']>=cut_gr))
    nopheno_df=mdf[cond_neutral]

    ## get nan index
    from matplotlib.font_manager import FontProperties

    fontP = FontProperties()
    fontP.set_size('xx-small')
    df_nan=mdf[mdf['rgr'].isna()| (fdf['rgr'].isna())]
    # df_extreme=df[df['rgr']<pheno]
    # df_under=df[df['rgr']>no_pheno]
    df_remain=mdf[~(cond_m |cond_f | cond_m_f | cond_neutral |(mdf['rgr'].isna()) | (fdf['rgr'].isna())) ]

    ### now we will plot all points
    f, ax = plt.subplots(figsize=[6.4*1.2,4.8*1.2])
    p1=ax.scatter(df_nan["knn_graph_X"], df_nan["knn_graph_Y"], s=6,color='#f5f5f5', label='No data')
    # p1=ax.scatter(df_nan["knn_graph_X"], df_nan["knn_graph_Y"], s=6,color='#A0B1BA', label='No data')
    p5=ax.scatter(nopheno_df["knn_graph_X"], nopheno_df["knn_graph_Y"], s=6,color='#00CD6C', label='MF and FF ( '+'f > '+str(cut_gr)+',m > '+str(cut_gr)+' )' )
    # p6=ax.scatter(df_remain["knn_graph_X"], df_remain["knn_graph_Y"], s=6,color='#A0B1BA', label='others (slight reduce in a sex)' )
    p2=ax.scatter(female_fertility_df["knn_graph_X"], female_fertility_df["knn_graph_Y"], s=6,color='#AF58BA', label='FR( '+'f < '+str(cut_less)+',m > '+str(cut_gr)+' )' )
    p3=ax.scatter(male_fertility_df["knn_graph_X"], male_fertility_df["knn_graph_Y"], s=6,color='#009ADE', label='MR( '+'f > '+str(cut_gr)+',m < '+str(cut_less)+' )' )
    p4=ax.scatter(both_df["knn_graph_X"], both_df["knn_graph_Y"], s=6,color='#F28522', label='MR and FR( '+'f < '+str(cut_less)+',m < '+str(cut_less)+' )' )

    # p2=ax.scatter(df_extreme["knn_graph_X"], df_extreme["knn_graph_Y"], s=6,color='r', label='RGR < '+str(pheno))
    # p3=ax.scatter(df_under["knn_graph_X"], df_under["knn_graph_Y"], s=6,color='g', label='RGR > '+str(no_pheno))
    # p3=ax.scatter(df_under["knn_graph_X"], df_under["knn_graph_Y"], s=6,color='#EDF155', label='RGR > '+str(no_pheno))
    # points = ax.scatter(x=df_remain['knn_graph_X'], y=df_remain['knn_graph_Y'], c=df_remain['rgr'], s=6, cmap='plasma')
    # # f.colorbar(points)
    # p2=ax.scatter(df_extreme["knn_graph_X"], df_extreme["knn_graph_Y"], s=6,color='#0B1283', label='RGR < '+str(pheno))

    # plt.legend(handles=[p1, p2,p3], title='Data labels', bbox_to_anchor=(1.1, 1), loc='upper left', prop=fontP)
    # plt.legend(handles=[p1, p2,p3], title='Data labels', bbox_to_anchor=(1.1, 1), loc='upper left')
    plt.legend(handles=[p1, p5,p2,p3,p4], title='RGR data', loc='best',frameon=True,prop=fontP)
    # points = ax.scatter(x= male_mca_df['knn_graph_X'], y=male_mca_df['knn_graph_Y'], c=colors, s=6, cmap='Spectral')
    #points = ax.scatter(x=df_remain['knn_graph_X'], y=df_remain['knn_graph_Y'], c=df_remain['rgr'], s=6, cmap='RdYlGn')

    # plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
    plt.grid(False)
    plt.savefig(out_folder+"/"+'mca_cluster_gene'+str(abs(cut_gr))+str(abs(cut_less))+'.pdf')

    plt.close()


# cut_grs=[-2, -1.5, -1 ,0.5]
# cut_lesss=[-2,-2.5, -3, -3.5, -4, -4.5, -5, -5.5, -6, -6.5]


# mask_dataframe_table_plot(male_mca_df,'mca_gene_cluster_male_color_gradient_3.pdf',pheno=-3,no_pheno=-0.25)
# mask_dataframe_table_plot(female_mca_df,'mca_gene_cluster_female_color_gradient_3.pdf',pheno=-3,no_pheno=-0.25)

### get colors


def male_female_only(mdf,fdf,sex,cut_gr=-2,cut_less=-2):

    print('hi')
    # get male fertility genes
    cond_m=((fdf['rgr']>cut_gr) & (mdf['rgr']<cut_less))
    male_fertility_df=mdf[cond_m]
    #male_fertility_genes=male_fertility_df['feature_symbol']
    # female fertility genes
    cond_f=((fdf['rgr']<cut_less) & (mdf['rgr']>cut_gr))
    female_fertility_df=mdf[cond_f]
    #female_fertility_genes=female_fertility_df['feature_symbol']
    ## both
    cond_m_f=((fdf['rgr']<cut_less) & (mdf['rgr']<cut_less))
    both_df=mdf[cond_m_f]
    #both_genes=both_df['feature_symbol']
    cond_neutral=((fdf['rgr']>=cut_gr) & (mdf['rgr']>=cut_gr))
    nopheno_df=mdf[cond_neutral]

    ## get nan index
    from matplotlib.font_manager import FontProperties

    fontP = FontProperties()
    fontP.set_size('xx-small')
    df_nan=mdf[mdf['rgr'].isna()| (fdf['rgr'].isna())]
    # df_extreme=df[df['rgr']<pheno]
    # df_under=df[df['rgr']>no_pheno]
    df_remain=mdf[~(cond_m |cond_f | cond_m_f | cond_neutral |(mdf['rgr'].isna()) | (fdf['rgr'].isna())) ]

    ### now we will plot all points
    f, ax = plt.subplots(figsize=[6.4*1.2,4.8*1.2])
    # p1=ax.scatter(df_nan["knn_graph_X"], df_nan["knn_graph_Y"], s=6,color='#f5f5f5', label='No data')
    p1=ax.scatter(df_nan["knn_graph_X"], df_nan["knn_graph_Y"], s=6,color='#f5f5f5')
    # p1=ax.scatter(df_nan["knn_graph_X"], df_nan["knn_graph_Y"], s=6,color='#A0B1BA', label='No data')
    p5=ax.scatter(nopheno_df["knn_graph_X"], nopheno_df["knn_graph_Y"], s=6,color='#f5f5f5')
    p4=ax.scatter(both_df["knn_graph_X"], both_df["knn_graph_Y"], s=6,color='#f5f5f5' )
    # p5=ax.scatter(nopheno_df["knn_graph_X"], nopheno_df["knn_graph_Y"], s=6,color='#00CD6C', label='MF and FF ( '+'f > '+str(cut_gr)+',m > '+str(cut_gr)+' )' )
    # p6=ax.scatter(df_remain["knn_graph_X"], df_remain["knn_graph_Y"], s=6,color='#A0B1BA', label='others (slight reduce in a sex)' )
    if sex=='F':
        p3=ax.scatter(male_fertility_df["knn_graph_X"], male_fertility_df["knn_graph_Y"], s=6,color='#f5f5f5' )
        p2=ax.scatter(female_fertility_df["knn_graph_X"], female_fertility_df["knn_graph_Y"], s=6,color='#AF58BA', label='FR( '+'f < '+str(cut_less)+',m > '+str(cut_gr)+' )' )
    else:
        p2=ax.scatter(female_fertility_df["knn_graph_X"], female_fertility_df["knn_graph_Y"], s=6,color='#f5f5f5')
        p3=ax.scatter(male_fertility_df["knn_graph_X"], male_fertility_df["knn_graph_Y"], s=6,color='#009ADE', label='MR( '+'f > '+str(cut_gr)+',m < '+str(cut_less)+' )' )
    # if sex=='M':
    #     p3=ax.scatter(male_fertility_df["knn_graph_X"], male_fertility_df["knn_graph_Y"], s=6,color='#009ADE', label='MR( '+'f > '+str(cut_gr)+',m < '+str(cut_less)+' )' )
    # else:
    #     p3=ax.scatter(male_fertility_df["knn_graph_X"], male_fertility_df["knn_graph_Y"], s=6,color='#f5f5f5' )




    # p2=ax.scatter(df_extreme["knn_graph_X"], df_extreme["knn_graph_Y"], s=6,color='r', label='RGR < '+str(pheno))
    # p3=ax.scatter(df_under["knn_graph_X"], df_under["knn_graph_Y"], s=6,color='g', label='RGR > '+str(no_pheno))
    # p3=ax.scatter(df_under["knn_graph_X"], df_under["knn_graph_Y"], s=6,color='#EDF155', label='RGR > '+str(no_pheno))
    # points = ax.scatter(x=df_remain['knn_graph_X'], y=df_remain['knn_graph_Y'], c=df_remain['rgr'], s=6, cmap='plasma')
    # # f.colorbar(points)
    # p2=ax.scatter(df_extreme["knn_graph_X"], df_extreme["knn_graph_Y"], s=6,color='#0B1283', label='RGR < '+str(pheno))

    # plt.legend(handles=[p1, p2,p3], title='Data labels', bbox_to_anchor=(1.1, 1), loc='upper left', prop=fontP)
    # plt.legend(handles=[p1, p2,p3], title='Data labels', bbox_to_anchor=(1.1, 1), loc='upper left')
    if sex=='F':
        plt.legend(handles=[p2], title='RGR data', loc='best',frameon=True,prop=fontP)
    else:
        plt.legend(handles=[p3], title='RGR data', loc='best',frameon=True,prop=fontP)
    # points = ax.scatter(x= male_mca_df['knn_graph_X'], y=male_mca_df['knn_graph_Y'], c=colors, s=6, cmap='Spectral')
    #points = ax.scatter(x=df_remain['knn_graph_X'], y=df_remain['knn_graph_Y'], c=df_remain['rgr'], s=6, cmap='RdYlGn')

    # plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
    plt.grid(False)
    plt.savefig(out_folder+"/"+sex+'_mca_cluster_gene'+str(abs(cut_gr))+str(abs(cut_less))+'.pdf')

    plt.close()



def male_top50_cluster(mdf,fdf):
    '''  '''
    top50=pd.read_excel(open(data_folder+'/Functional gene lists for MCA overlay.xlsx', 'rb'),sheet_name='top50male')
    male_fertility_df=mdf[mdf['feature_symbol'].isin(top50['pbanka_id'].to_list())]
    clust=[i+1 for i in range(20)]
    ind=np.arange(20)
    df=pd.DataFrame(columns=['prob'])
    for c in clust:
        tmp=mca_df[mca_df['Cluster']==c]
        # clusters.append(str(c))
        clust_genes=set(tmp['feature_symbol'].to_list())
        male_genes=set(male_fertility_df['feature_symbol'].to_list())
        comm=clust_genes&male_genes
        df.loc[c,'prob']=len(comm)/len(clust_genes)

    plt.bar(ind, df['prob'].to_list())
    plt.plot(ind, [0.05]*20,'r--')
    plt.xticks(ind, df.index.to_list())
    plt.grid(False)
    plt.savefig(out_folder+'/top50 in cluster.pdf' )

    plt.close()

def male_female_top50(mdf,fdf,sex,cut_gr=-2,cut_less=-2):
    top50=pd.read_excel(open(data_folder+'/Functional gene lists for MCA overlay.xlsx', 'rb'),sheet_name='top50male')

    print('hi')
    top50=pd.read_excel(open(data_folder+'/Functional gene lists for MCA overlay.xlsx', 'rb'),sheet_name='top50male')

    male_fertility_df=mdf[mdf['feature_symbol'].isin(top50['pbanka_id'].to_list())]

    # # get male fertility genes
    # cond_m=((fdf['rgr']>cut_gr) & (mdf['rgr']<cut_less))
    # male_fertility_df=mdf[cond_m]
    # #male_fertility_genes=male_fertility_df['feature_symbol']
    # # female fertility genes
    # cond_f=((fdf['rgr']<cut_less) & (mdf['rgr']>cut_gr))
    # female_fertility_df=mdf[cond_f]
    # #female_fertility_genes=female_fertility_df['feature_symbol']
    # ## both
    # cond_m_f=((fdf['rgr']<cut_less) & (mdf['rgr']<cut_less))
    # both_df=mdf[cond_m_f]
    # #both_genes=both_df['feature_symbol']
    # cond_neutral=((fdf['rgr']>=cut_gr) & (mdf['rgr']>=cut_gr))
    # nopheno_df=mdf[cond_neutral]

    ## get nan index
    from matplotlib.font_manager import FontProperties

    fontP = FontProperties()
    fontP.set_size('xx-small')
    # df_nan=mdf[mdf['rgr'].isna()| (fdf['rgr'].isna())]
    # # df_extreme=df[df['rgr']<pheno]
    # # df_under=df[df['rgr']>no_pheno]
    # df_remain=mdf[~(cond_m |cond_f | cond_m_f | cond_neutral |(mdf['rgr'].isna()) | (fdf['rgr'].isna())) ]
    df_nan=mdf[~mdf['feature_symbol'].isin(top50['pbanka_id'].to_list())]
    ### now we will plot all points
    f, ax = plt.subplots(figsize=[6.4*1.2,4.8*1.2])
    # p1=ax.scatter(df_nan["knn_graph_X"], df_nan["knn_graph_Y"], s=6,color='#f5f5f5', label='No data')
    p1=ax.scatter(df_nan["knn_graph_X"], df_nan["knn_graph_Y"], s=6,color='#f5f5f5')
    # p1=ax.scatter(df_nan["knn_graph_X"], df_nan["knn_graph_Y"], s=6,color='#A0B1BA', label='No data')
    # p5=ax.scatter(nopheno_df["knn_graph_X"], nopheno_df["knn_graph_Y"], s=6,color='#f5f5f5')
    # p4=ax.scatter(both_df["knn_graph_X"], both_df["knn_graph_Y"], s=6,color='#f5f5f5' )
    # p5=ax.scatter(nopheno_df["knn_graph_X"], nopheno_df["knn_graph_Y"], s=6,color='#00CD6C', label='MF and FF ( '+'f > '+str(cut_gr)+',m > '+str(cut_gr)+' )' )
    # p6=ax.scatter(df_remain["knn_graph_X"], df_remain["knn_graph_Y"], s=6,color='#A0B1BA', label='others (slight reduce in a sex)' )
    if sex=='F':
        p3=ax.scatter(male_fertility_df["knn_graph_X"], male_fertility_df["knn_graph_Y"], s=6,color='#f5f5f5' )
        p2=ax.scatter(female_fertility_df["knn_graph_X"], female_fertility_df["knn_graph_Y"], s=6,color='#AF58BA', label='FR( '+'f < '+str(cut_less)+',m > '+str(cut_gr)+' )' )
    else:
        # p2=ax.scatter(female_fertility_df["knn_graph_X"], female_fertility_df["knn_graph_Y"], s=6,color='#f5f5f5')
        p3=ax.scatter(male_fertility_df["knn_graph_X"], male_fertility_df["knn_graph_Y"], s=6,color='#009ADE', label='MR top50' )
    # if sex=='M':
    #     p3=ax.scatter(male_fertility_df["knn_graph_X"], male_fertility_df["knn_graph_Y"], s=6,color='#009ADE', label='MR( '+'f > '+str(cut_gr)+',m < '+str(cut_less)+' )' )
    # else:
    #     p3=ax.scatter(male_fertility_df["knn_graph_X"], male_fertility_df["knn_graph_Y"], s=6,color='#f5f5f5' )




    # p2=ax.scatter(df_extreme["knn_graph_X"], df_extreme["knn_graph_Y"], s=6,color='r', label='RGR < '+str(pheno))
    # p3=ax.scatter(df_under["knn_graph_X"], df_under["knn_graph_Y"], s=6,color='g', label='RGR > '+str(no_pheno))
    # p3=ax.scatter(df_under["knn_graph_X"], df_under["knn_graph_Y"], s=6,color='#EDF155', label='RGR > '+str(no_pheno))
    # points = ax.scatter(x=df_remain['knn_graph_X'], y=df_remain['knn_graph_Y'], c=df_remain['rgr'], s=6, cmap='plasma')
    # # f.colorbar(points)
    # p2=ax.scatter(df_extreme["knn_graph_X"], df_extreme["knn_graph_Y"], s=6,color='#0B1283', label='RGR < '+str(pheno))

    # plt.legend(handles=[p1, p2,p3], title='Data labels', bbox_to_anchor=(1.1, 1), loc='upper left', prop=fontP)
    # plt.legend(handles=[p1, p2,p3], title='Data labels', bbox_to_anchor=(1.1, 1), loc='upper left')
    if sex=='F':
        plt.legend(handles=[p2], title='RGR data', loc='best',frameon=True,prop=fontP)
    else:
        plt.legend(handles=[p3], title='RGR data', loc='best',frameon=True,prop=fontP)
    # points = ax.scatter(x= male_mca_df['knn_graph_X'], y=male_mca_df['knn_graph_Y'], c=colors, s=6, cmap='Spectral')
    #points = ax.scatter(x=df_remain['knn_graph_X'], y=df_remain['knn_graph_Y'], c=df_remain['rgr'], s=6, cmap='RdYlGn')

    # plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
    plt.grid(False)
    plt.savefig(out_folder+"/"+sex+'_mca_cluster_genr_top50'+str(abs(cut_gr))+str(abs(cut_less))+'.pdf')

    plt.close()




def plot_female_male_on_mca(male_mca_df,female_mca_df):

    ## we will plot rgr

    ## cutoff
    cut_grs=[-2]
    cut_lesss=[-2]
    for cut_gr in cut_grs:
        for cut_less in cut_lesss:
            define_much_greater(male_mca_df,female_mca_df,cut_gr,cut_less)


def get_data():
    # get x, y co-ordinates for the and mca gene cluster
    mca_df=pd.read_excel(open(data_folder+'/ginny_gene_mca.xlsx','rb'), sheet_name='mca')
    fertility_df=pd.read_excel(open(data_folder+'/Phenotype_call_plasmogem_id_Vikash_271020.xlsx','rb'),sheet_name='fertility')
    female_df=fertility_df[~(fertility_df['GCKO2_pheno']=='No Data')]
    male_df=fertility_df[~(fertility_df['g145480_pheno']=='No Data')]
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
            male_mca_df.loc[idx,'pheno']=tmp['g145480_pheno'].to_list()[0]
            male_mca_df.loc[idx,'rgr']=tmp['g145480_RGR'].to_list()[0]
            male_mca_df.loc[idx,'sd']=tmp['g145480_sd'].to_list()[0]
        tmp=female_df[female_df['pbanka_id']==female_mca_df.loc[idx,'feature_symbol']]
        if not tmp.empty:
            female_mca_df.loc[idx,'pheno']=tmp['GCKO2_pheno'].to_list()[0]
            female_mca_df.loc[idx,'rgr']=tmp['GCKO2_RGR'].to_list()[0]
            female_mca_df.loc[idx,'sd']=tmp['GCKO2_sd'].to_list()[0]

    return mca_df,male_df,female_df,male_mca_df,female_mca_df



def plot_male_female_scatter(mdf,fdf,cut_gr=-2,cut_less=-2):
    na_df=(mdf['rgr'].isna())| (fdf['rgr'].isna())
    # get male fertility genes
    cond_m=((fdf['rgr']>cut_gr) & (mdf['rgr']<cut_less))
    male_fertility_df=mdf[cond_m]
    #male_fertility_genes=male_fertility_df['feature_symbol']
    # female fertility genes
    cond_f=((fdf['rgr']<cut_less) & (mdf['rgr']>cut_gr))
    female_fertility_df=mdf[cond_f]
    #female_fertility_genes=female_fertility_df['feature_symbol']
    ## both
    cond_m_f=((fdf['rgr']<cut_less) & (mdf['rgr']<cut_less))
    both_df=mdf[cond_m_f]
    #both_genes=both_df['feature_symbol']
    cond_neutral=((fdf['rgr']>=cut_gr) & (mdf['rgr']>=cut_gr))
    nopheno_df=mdf[cond_neutral]

    print(male_fertility_df.shape[0]+female_fertility_df.shape[0]+both_df.shape[0]+nopheno_df.shape[0])

    ## get nan index
    from matplotlib.font_manager import FontProperties

    fontP = FontProperties()
    fontP.set_size('xx-small')

    # df_extreme=df[df['rgr']<pheno]
    # df_under=df[df['rgr']>no_pheno]
    df_remain=mdf[~(cond_m |cond_f | cond_m_f | cond_neutral |(mdf['rgr'].isna()) | (fdf['rgr'].isna())) ]

    ### now we will plot all points
    size=6
    f, ax = plt.subplots(figsize=[6.4,6.4])
    # p1=ax.scatter(df_nan["knn_graph_X"], df_nan["knn_graph_Y"], s=6,color='#f5f5f5', label='No data')
    # p1=ax.scatter(df_nan["knn_graph_X"], df_nan["knn_graph_Y"], s=6,color='#A0B1BA', label='No data')
    p5=ax.scatter(mdf[cond_neutral]['rgr'], fdf[cond_neutral]['rgr'], s=size,color='#00CD6C', label='MF and FF ( '+'f > '+str(cut_gr)+',m > '+str(cut_gr)+' )' )
    # p6=ax.scatter(df_remain["knn_graph_X"], df_remain["knn_graph_Y"], s=6,color='#A0B1BA', label='others (slight reduce in a sex)' )
    p2=ax.scatter(mdf[cond_f]['rgr'], fdf[cond_f]['rgr'], s=size,color='#AF58BA', label='FR( '+'f < '+str(cut_less)+',m > '+str(cut_gr)+' )' )
    p3=ax.scatter(mdf[cond_m]['rgr'], fdf[cond_m]['rgr'], s=size,color='#009ADE', label='MR( '+'f > '+str(cut_gr)+',m < '+str(cut_less)+' )' )
    # p4=ax.scatter(mdf[cond_m_f]['rgr'], fdf[cond_m_f]['rgr'], s=6,color='#FFC61E', label='MR and FR( '+'f < '+str(cut_less)+',m < '+str(cut_less)+' )' )
    p4=ax.scatter(mdf[cond_m_f]['rgr'], fdf[cond_m_f]['rgr'], s=size,color='#F28522', label='MR and FR( '+'f < '+str(cut_less)+',m < '+str(cut_less)+' )' )
    plt.legend(handles=[p5,p2,p3,p4], title='RGR data', loc='best',frameon=True,prop=fontP)
    # points = ax.scatter(x= male_mca_df['knn_graph_X'], y=male_mca_df['knn_graph_Y'], c=colors, s=6, cmap='Spectral')
    #points = ax.scatter(x=df_remain['knn_graph_X'], y=df_remain['knn_graph_Y'], c=df_remain['rgr'], s=6, cmap='RdYlGn')

    plt.xlabel('Relative male fertility')
    plt.ylabel('Relative female fertility')
    # plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
    plt.grid(False)
    plt.savefig(out_folder+"/"+'male_female_scatter'+str(abs(cut_gr))+str(abs(cut_less))+'.pdf')

    plt.close()




def change_in_male_female(mdf,fdf,cut_gr=-2,cut_less=-5):


    ## get nan index
    from matplotlib.font_manager import FontProperties

    fontP = FontProperties()
    fontP.set_size('xx-small')

    df_nan=mdf[mdf['rgr'].isna()| (fdf['rgr'].isna())]

    df_remain_m=mdf[~((mdf['rgr'].isna()) | (fdf['rgr'].isna())) ]
    df_remain_f=fdf[~((mdf['rgr'].isna()) | (fdf['rgr'].isna())) ]

    df_change=df_remain_m.copy()
    df_change['rgr']=df_remain_m['rgr']-df_remain_f['rgr']

    ### now we will plot all points
    f, ax = plt.subplots(figsize=[6.4*1.2,4.8*1.2])
    p1=ax.scatter(df_nan["knn_graph_X"], df_nan["knn_graph_Y"], s=6,color='#f5f5f5', label='No data')
    points = ax.scatter(x=df_change['knn_graph_X'], y=df_change['knn_graph_Y'], c=df_change['rgr'], s=6, cmap='plasma')
    f.colorbar(points)
    plt.legend(handles=[p1], title='RGR data', loc='best',frameon=True,prop=fontP)

    #p2=ax.scatter(df_extreme["knn_graph_X"], df_extreme["knn_graph_Y"], s=6,color='r', label='RGR < '+str(pheno))
    # # p3=ax.scatter(df_under["knn_graph_X"], df_under["knn_graph_Y"], s=6,color='g', label='RGR > '+str(no_pheno))
    # # p3=ax.scatter(df_under["knn_graph_X"], df_under["knn_graph_Y"], s=6,color='#EDF155', label='RGR > '+str(no_pheno))
    # # points = ax.scatter(x=df_remain['knn_graph_X'], y=df_remain['knn_graph_Y'], c=df_remain['rgr'], s=6, cmap='plasma')
    # # # f.colorbar(points)
    # # p2=ax.scatter(df_extreme["knn_graph_X"], df_extreme["knn_graph_Y"], s=6,color='#0B1283', label='RGR < '+str(pheno))
    #
    # # plt.legend(handles=[p1, p2,p3], title='Data labels', bbox_to_anchor=(1.1, 1), loc='upper left', prop=fontP)
    # # plt.legend(handles=[p1, p2,p3], title='Data labels', bbox_to_anchor=(1.1, 1), loc='upper left')
    # plt.legend(handles=[p1, p5,p2,p3,p4], title='RGR data', loc='best',frameon=True,prop=fontP)
    # # points = ax.scatter(x= male_mca_df['knn_graph_X'], y=male_mca_df['knn_graph_Y'], c=colors, s=6, cmap='Spectral')
    #points = ax.scatter(x=df_remain['knn_graph_X'], y=df_remain['knn_graph_Y'], c=df_remain['rgr'], s=6, cmap='RdYlGn')

    # plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
    plt.grid(False)
    plt.savefig(out_folder+"/"+'mca_cluster_difference_male_female.pdf')
    plt.close()


def color_subset_genes(mdf,fdf,metabolic_genes,cut_gr=-2,cut_less=-2,group_name='metabolic genes'):
    cond_m=((fdf['rgr']>cut_gr) & (mdf['rgr']<cut_less) & (mdf['feature_symbol'].isin(metabolic_genes)))
    male_fertility_df=mdf[cond_m]

    cond_f=((fdf['rgr']<cut_less) & (mdf['rgr']>cut_gr) & (mdf['feature_symbol'].isin(metabolic_genes)))
    female_fertility_df=mdf[cond_f]
    cond_m_f=((fdf['rgr']<cut_less) & (mdf['rgr']<cut_less) & (mdf['feature_symbol'].isin(metabolic_genes)))
    both_df=mdf[cond_m_f]
    #both_genes=both_df['feature_symbol']
    cond_neutral=((fdf['rgr']>=cut_gr) & (mdf['rgr']>=cut_gr) & (mdf['feature_symbol'].isin(metabolic_genes)))
    nopheno_df=mdf[cond_neutral]

    ## get nan index
    from matplotlib.font_manager import FontProperties

    fontP = FontProperties()
    fontP.set_size('xx-small')
    df_nan=mdf[mdf['rgr'].isna()| (fdf['rgr'].isna()) | ~(mdf['feature_symbol'].isin(metabolic_genes))]
    ##

    # df_nan=mdf[mdf['rgr'].isna()| (fdf['rgr'].isna()) ]
    tdf=mdf[(mdf['feature_symbol'].isin(metabolic_genes))]
    custom_df=tdf[tdf['rgr'].isna()]

    if group_name=='Flagellar proteins':
        import pdb; pdb.set_trace()
        pointdf=mdf[mdf['feature_symbol']=='PBANKA_0106200']
    # df_extreme=df[df['rgr']<pheno]
    # df_under=df[df['rgr']>no_pheno]
    # df_remain=mdf[~(cond_m |cond_f | cond_m_f | cond_neutral |(mdf['rgr'].isna()) | (fdf['rgr'].isna())) ]

    ### now we will plot all points
    f, ax = plt.subplots(figsize=[6.4*1.2,4.8*1.2])

    p1=ax.scatter(df_nan["knn_graph_X"], df_nan["knn_graph_Y"], s=6,color='#f5f5f5', label='No data')
    # if group_name=='Flagellar proteins':
    p6=ax.scatter(custom_df["knn_graph_X"], custom_df["knn_graph_Y"], s=6,color='#A0B1BA', label='Not in list')
    p5=ax.scatter(nopheno_df["knn_graph_X"], nopheno_df["knn_graph_Y"], s=6,color='#00CD6C', label='MF and FF ( '+'f > '+str(cut_gr)+',m > '+str(cut_gr)+' )' )
    # p6=ax.scatter(df_remain["knn_graph_X"], df_remain["knn_graph_Y"], s=6,color='#A0B1BA', label='others (slight reduce in a sex)' )
    p2=ax.scatter(female_fertility_df["knn_graph_X"], female_fertility_df["knn_graph_Y"], s=6,color='#AF58BA', label='FR( '+'f < '+str(cut_less)+',m > '+str(cut_gr)+' )' )
    p3=ax.scatter(male_fertility_df["knn_graph_X"], male_fertility_df["knn_graph_Y"], s=6,color='#009ADE', label='MR( '+'f > '+str(cut_gr)+',m < '+str(cut_less)+' )' )
    p4=ax.scatter(both_df["knn_graph_X"], both_df["knn_graph_Y"], s=6,color='#F28522', label='MR and FR( '+'f < '+str(cut_less)+',m < '+str(cut_less)+' )' )
    if group_name=='Flagellar proteins':

        ax.scatter(pointdf["knn_graph_X"], pointdf["knn_graph_Y"], s=10,color='#F28522',facecolors='none', edgecolors='r' )

    plt.legend(handles=[p1, p5,p2,p3,p4], title='RGR data', loc='best',frameon=True,prop=fontP)
    plt.title('RGR data for '+ group_name)

    plt.grid(False)
    plt.savefig(out_folder+"/"+'mca_cluster_gene_'+group_name+str(abs(cut_gr))+str(abs(cut_less))+'.pdf')

    plt.close()



def color_group_genes(mdf,fdf,cut_gr=-2,cut_less=-2):
    metabolic_genes,df_subsystem=read_metabolic_genes_subsystem(genelist=None)
    #[apico_df,mito_df,ribo_df,proteinex_df,egress_df,gam_df]
    ldfs=read_group_genes()
    groups=['Apicoplast','Mitochondrion','Ribosomal & DA','Protein export','Egressome','Gametocyte','Microgamete','Flagellar proteins']
    # get new pbanka_id
    for i,gp_name in enumerate(groups):
        genes=ldfs[i]['pbanka_id'].to_list()
        new_ids=[]
        for g in genes:
            if g in prev_to_new.keys():
                new_ids.append(prev_to_new[g])
            else:
                new_ids.append(g)


        color_subset_genes(mdf,fdf,new_ids,cut_gr=-2,cut_less=-2,group_name=gp_name)


def read_group_genes():
    ''' read data from excel files  '''
    # get excel data
    apico_df=pd.read_excel(open(data_folder+'/Functional gene lists for MCA overlay.xlsx', 'rb'),sheet_name='Apicoplast_Bushell et al. 2017')
    mito_df=pd.read_excel(open(data_folder+'/Functional gene lists for MCA overlay.xlsx', 'rb'),sheet_name='Mitochondrion_Bushell et al.')
    ribo_df=pd.read_excel(open(data_folder+'/Functional gene lists for MCA overlay.xlsx', 'rb'),sheet_name='Ribosomal & DA_Bushell et al.')
    proteinex_df=pd.read_excel(open(data_folder+'/Functional gene lists for MCA overlay.xlsx', 'rb'),sheet_name='Protein export_Bushell et al.')
    egress_df=pd.read_excel(open(data_folder+'/Functional gene lists for MCA overlay.xlsx', 'rb'),sheet_name='Egressome_Kehrer et al. 2016')
    gam_df=pd.read_excel(open(data_folder+'/Functional gene lists for MCA overlay.xlsx', 'rb'),sheet_name='Gametocyte_SMC_Pandey et al.')
    mgam_df=pd.read_excel(open(data_folder+'/Functional gene lists for MCA overlay.xlsx', 'rb'),sheet_name='Microgamete_Talman et al. 2014')
    flagellar_df=pd.read_excel(open(data_folder+'/Functional gene lists for MCA overlay.xlsx', 'rb'),sheet_name='Flagellar proteins_Sinden et al')
    return [apico_df,mito_df,ribo_df,proteinex_df,egress_df,gam_df,mgam_df,flagellar_df]


def color_metabolic_genes(mdf,fdf,cut_gr=-2,cut_less=-2):
    metabolic_genes,df_subsystem=read_metabolic_genes_subsystem(genelist=None)
    color_subset_genes(mdf,fdf,metabolic_genes,cut_gr=-2,cut_less=-2)



def enrich_mca_cluster(mca_df,mdf,fdf,cut_gr=-2,cut_less=-2):
    ''' Pefrorming cluster enrichment bassed on screen    '''



    na_df=(mdf['rgr'].isna())| (fdf['rgr'].isna())
    # get male fertility genes
    cond_m=((fdf['rgr']>cut_gr) & (mdf['rgr']<cut_less))
    male_fertility_df=mdf[cond_m]
    #male_fertility_genes=male_fertility_df['feature_symbol']
    # female fertility genes
    cond_f=((fdf['rgr']<cut_less) & (mdf['rgr']>cut_gr))
    female_fertility_df=mdf[cond_f]
    #female_fertility_genes=female_fertility_df['feature_symbol']
    ## both
    cond_m_f=((fdf['rgr']<cut_less) & (mdf['rgr']<cut_less))
    both_df=mdf[cond_m_f]
    #both_genes=both_df['feature_symbol']
    cond_neutral=((fdf['rgr']>=cut_gr) & (mdf['rgr']>=cut_gr))
    nopheno_df=mdf[cond_neutral]

    print(male_fertility_df.shape[0]+female_fertility_df.shape[0]+both_df.shape[0]+nopheno_df.shape[0])

    screen_genens=male_fertility_df['feature_symbol'].to_list()+female_fertility_df['feature_symbol'].to_list()+both_df['feature_symbol'].to_list()+nopheno_df['feature_symbol'].to_list()


    # df=pd.read_csv('/Users/vikash/Documents/Projects/GeneRegulatory/Final_analysis/io/Phenotype_cluster/phenotype_for_clustering.txt',sep='\t')
    # df=df.fillna('NA')
    # # df['consensus phenotype']=df['consensus phenotype'].str.replace('None','NA')
    # tmp=df.copy() ## we saved the original dataframe
    # tmp.set_index('new gene ID',inplace=True)
    # screen_df=tmp[(tmp['consensus phenotype']== 'Males')| (tmp['consensus phenotype']== 'Females')| (tmp['consensus phenotype']== 'Both')|(tmp['consensus phenotype']== 'None')]
    # pheno_df=tmp[(tmp['consensus phenotype']== 'Males')| (tmp['consensus phenotype']== 'Females')| (tmp['consensus phenotype']== 'Both')]
    # ## groupby along Phenotyepes

    clust=[i+1 for i in range(20)]

    types=['male','female','both', 'nosex']
    listgenes=[male_fertility_df['feature_symbol'].to_list(),
               female_fertility_df['feature_symbol'].to_list(),
               both_df['feature_symbol'].to_list(), nopheno_df['feature_symbol'].to_list()
                ]
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
        enrich_soreted.to_csv(out_folder+'/Cluster_enrichment_with_phenotype_'+ types[k] +'.txt',sep='\t',index=None)

        import pdb; pdb.set_trace()

def get_cutoff(mdf,fdf,cut=-1,shift=2):
    # we are going to identify parameetrs
    mdf['cut1']=mdf['rgr']+shift*mdf['sd']
    mdf['cut2']=mdf['rgr']-shift*mdf['sd']

    # female

    fdf['cut1']=fdf['rgr']+shift*fdf['sd']
    fdf['cut2']=fdf['rgr']-shift*fdf['sd']

    male_df=mdf[(mdf['cut1']<cut) & (fdf['cut2']>cut)]
    female_df=mdf[(mdf['cut2']>cut) & (fdf['cut1']<cut)]
    both_df=mdf[(mdf['cut1']<cut) & (fdf['cut1']<cut)]

    tmp_m=mdf[((mdf['Cluster']==11) | (mdf['Cluster']==12) )]
    tmp_f=mdf[((mdf['Cluster']==13) | (mdf['Cluster']==14) )]

    male_genes=set(tmp_m['feature_symbol'].to_list()) & set(male_df['feature_symbol'].to_list())
    female_genes=set(tmp_f['feature_symbol'].to_list()) & set(female_df['feature_symbol'].to_list())
    both_genes= set(both_df['feature_symbol'].to_list())
    return len(male_genes),len(female_genes),len(both_genes)




def sensitivity_analysis(mdf,fdf):
    cuts=np.linspace(-2, 0, num=20)
    shifts=np.linspace(0, 2, num=20)
    x=[]
    y=[]
    v=[]
    v1=[]
    v2=[]
    vsum=[]
    vsex=[]
    # plist=[]
    # nlist=[]

    for cut in cuts:
        for shift in shifts:
            m,f,b=get_cutoff(mdf,fdf,cut,shift)
            x.append(cut)
            y.append(shift)
            v.append(m)
            v1.append(f)
            v2.append(b)
            vsum.append(25*m+25*f+1*b)
            vsex.append(m+f)
            # plist.append(p)
            # nlist.append(n)
    res_df=pd.DataFrame(columns=['cut','shift','m_value','f_value','b_value'])
    res_df['cut']=x
    res_df['shift']=y
    res_df['m_value']=v
    res_df['f_value']=v1
    res_df['b_value']=v2
    res_df['sumall']=vsum
    res_df['sum']=vsex

    pickle.dump(res_df,open(out_folder+'/'+'optim_n.p','wb'))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # f, ax = plt.subplots(figsize=[6.4*1.2,4.8*1.2])
    p1=ax.scatter(res_df["cut"], res_df["shift"], res_df['m_value'], s=6,color='b')

    ax.set_xlabel('Decision cutoff')
    ax.set_ylabel('Order of SD')
    ax.set_zlabel('male_specific_genes')
    # plt.show()
    plt.savefig(out_folder+'/_m_sensitivity.pdf')

    plt.close()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    p1=ax.scatter(res_df["cut"], res_df["shift"], res_df['f_value'], s=6,color='b')

    ax.set_xlabel('Decision cutoff')
    ax.set_ylabel('Order of SD')
    ax.set_zlabel('female_specific_genes')
    # plt.show()
    plt.savefig(out_folder+ '/_f_sensitivity.pdf')

    plt.close()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    p1=ax.scatter(res_df["cut"], res_df["shift"], res_df['b_value'], s=6,color='b')

    ax.set_xlabel('Decision cutoff')
    ax.set_ylabel('Order of SD')
    ax.set_zlabel('both_specific_genes')
    # plt.show()
    plt.savefig(out_folder+  '/_b_sensitivity.pdf')

    plt.close()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    p1=ax.scatter(res_df["cut"], res_df["shift"], res_df['sum'], s=6,color='b')

    ax.set_xlabel('Decision cutoff')
    ax.set_ylabel('Order of SD')
    ax.set_zlabel('sum')
    # plt.show()
    plt.savefig(out_folder+  '/_sum_sensitivity.pdf')

    plt.close()
    return res_df


def best_m_f_plot(mdf,fdf,cut=-1,shift=2):


    # we are going to identify parameetrs
    mdf['cut1']=mdf['rgr']+shift*mdf['sd']
    mdf['cut2']=mdf['rgr']-shift*mdf['sd']

    # female

    fdf['cut1']=fdf['rgr']+shift*fdf['sd']
    fdf['cut2']=fdf['rgr']-shift*fdf['sd']
    cond_m=(mdf['cut1']<cut) & (fdf['cut2']>cut)
    cond_f=(mdf['cut2']>cut) & (fdf['cut1']<cut)
    cond_m_f=(mdf['cut1']<cut) & (fdf['cut1']<cut)
    cond_neutral=(mdf['cut2']>=cut) & (fdf['cut2']>=cut)
    male_fertility_df=mdf[cond_m]
    female_fertility_df=mdf[cond_f]
    both_df=mdf[cond_m_f]
    nopheno_df=mdf[cond_neutral]


    ## get nan index
    from matplotlib.font_manager import FontProperties

    fontP = FontProperties()
    fontP.set_size('xx-small')
    df_nan=mdf[mdf['rgr'].isna()| (fdf['rgr'].isna())]
    # df_extreme=df[df['rgr']<pheno]
    # df_under=df[df['rgr']>no_pheno]
    remain=~(cond_m |cond_f | cond_m_f | cond_neutral |(mdf['rgr'].isna()) | (fdf['rgr'].isna()))
    df_remain=mdf[remain]

    ### now we will plot all points
    f, ax = plt.subplots(figsize=[6.4*1.2,4.8*1.2])
    p1=ax.scatter(df_nan["knn_graph_X"], df_nan["knn_graph_Y"], s=6,color='#f5f5f5', label='No data')
    # p1=ax.scatter(df_nan["knn_graph_X"], df_nan["knn_graph_Y"], s=6,color='#A0B1BA', label='No data')
    p6=ax.scatter(df_remain["knn_graph_X"], df_remain["knn_graph_Y"], s=6,color='#A0B1BA', label='Cant say' )
    p5=ax.scatter(nopheno_df["knn_graph_X"], nopheno_df["knn_graph_Y"], s=6,color='#00CD6C', label='MF and FF' )

    p2=ax.scatter(female_fertility_df["knn_graph_X"], female_fertility_df["knn_graph_Y"], s=6,color='#AF58BA', label='FR' )
    p3=ax.scatter(male_fertility_df["knn_graph_X"], male_fertility_df["knn_graph_Y"], s=6,color='#009ADE', label='MR')
    p4=ax.scatter(both_df["knn_graph_X"], both_df["knn_graph_Y"], s=6,color='#F28522', label='MR and FR' )

    # p2=ax.scatter(df_extreme["knn_graph_X"], df_extreme["knn_graph_Y"], s=6,color='r', label='RGR < '+str(pheno))
    # p3=ax.scatter(df_under["knn_graph_X"], df_under["knn_graph_Y"], s=6,color='g', label='RGR > '+str(no_pheno))
    # p3=ax.scatter(df_under["knn_graph_X"], df_under["knn_graph_Y"], s=6,color='#EDF155', label='RGR > '+str(no_pheno))
    # points = ax.scatter(x=df_remain['knn_graph_X'], y=df_remain['knn_graph_Y'], c=df_remain['rgr'], s=6, cmap='plasma')
    # # f.colorbar(points)
    # p2=ax.scatter(df_extreme["knn_graph_X"], df_extreme["knn_graph_Y"], s=6,color='#0B1283', label='RGR < '+str(pheno))

    # plt.legend(handles=[p1, p2,p3], title='Data labels', bbox_to_anchor=(1.1, 1), loc='upper left', prop=fontP)
    # plt.legend(handles=[p1, p2,p3], title='Data labels', bbox_to_anchor=(1.1, 1), loc='upper left')
    plt.legend(handles=[p1,p6,p5,p2,p3,p4], title='RGR data', loc='best',frameon=True,prop=fontP)
    # points = ax.scatter(x= male_mca_df['knn_graph_X'], y=male_mca_df['knn_graph_Y'], c=colors, s=6, cmap='Spectral')
    #points = ax.scatter(x=df_remain['knn_graph_X'], y=df_remain['knn_graph_Y'], c=df_remain['rgr'], s=6, cmap='RdYlGn')

    # plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
    plt.grid(False)
    plt.savefig(out_folder+"/"+'mca_cluster_gene_best.pdf')

    plt.close()


def plot_male_female_scatter_best(mdf,fdf,cut=-1,shift=2):
    na_df=(mdf['rgr'].isna())| (fdf['rgr'].isna())

    mdf['cut1']=mdf['rgr']+shift*mdf['sd']
    mdf['cut2']=mdf['rgr']-shift*mdf['sd']

    # female

    fdf['cut1']=fdf['rgr']+shift*fdf['sd']
    fdf['cut2']=fdf['rgr']-shift*fdf['sd']
    cond_m=(mdf['cut1']<cut) & (fdf['cut2']>cut)
    cond_f=(mdf['cut2']>cut) & (fdf['cut1']<cut)
    cond_m_f=(mdf['cut1']<cut) & (fdf['cut1']<cut)
    cond_neutral=(mdf['cut2']>=cut) & (fdf['cut2']>=cut)
    male_fertility_df=mdf[cond_m]
    female_fertility_df=mdf[cond_f]
    both_df=mdf[cond_m_f]
    nopheno_df=mdf[cond_neutral]

    print(male_fertility_df.shape[0]+female_fertility_df.shape[0]+both_df.shape[0]+nopheno_df.shape[0])

    ## get nan index
    from matplotlib.font_manager import FontProperties

    fontP = FontProperties()
    fontP.set_size('xx-small')

    # df_extreme=df[df['rgr']<pheno]
    # df_under=df[df['rgr']>no_pheno]
    remain=~(cond_m |cond_f | cond_m_f | cond_neutral |(mdf['rgr'].isna()) | (fdf['rgr'].isna()))
    df_remain=mdf[remain]

    ### now we will plot all points
    size=6
    f, ax = plt.subplots(figsize=[6.4,6.4])
    # p1=ax.scatter(df_nan["knn_graph_X"], df_nan["knn_graph_Y"], s=6,color='#f5f5f5', label='No data')
    # p1=ax.scatter(df_nan["knn_graph_X"], df_nan["knn_graph_Y"], s=6,color='#A0B1BA', label='No data')
    p5=ax.scatter(mdf[cond_neutral]['rgr'], fdf[cond_neutral]['rgr'], s=size,color='#00CD6C', label='MF and FF' )
    p6=ax.scatter(mdf[remain]['rgr'], fdf[remain]['rgr'], s=size,color='#A0B1BA', label='Cant say')
    # p6=ax.scatter(df_remain["knn_graph_X"], df_remain["knn_graph_Y"], s=6,color='#A0B1BA', label='others (slight reduce in a sex)' )
    p2=ax.scatter(mdf[cond_f]['rgr'], fdf[cond_f]['rgr'], s=size,color='#AF58BA', label='FR' )
    p3=ax.scatter(mdf[cond_m]['rgr'], fdf[cond_m]['rgr'], s=size,color='#009ADE', label='MR' )
    # p4=ax.scatter(mdf[cond_m_f]['rgr'], fdf[cond_m_f]['rgr'], s=6,color='#FFC61E', label='MR and FR( '+'f < '+str(cut_less)+',m < '+str(cut_less)+' )' )
    p4=ax.scatter(mdf[cond_m_f]['rgr'], fdf[cond_m_f]['rgr'], s=size,color='#F28522', label='MR and FR')

    plt.legend(handles=[p5,p6,p2,p3,p4], title='RGR data', loc='best',frameon=True,prop=fontP)
    # points = ax.scatter(x= male_mca_df['knn_graph_X'], y=male_mca_df['knn_graph_Y'], c=colors, s=6, cmap='Spectral')
    #points = ax.scatter(x=df_remain['knn_graph_X'], y=df_remain['knn_graph_Y'], c=df_remain['rgr'], s=6, cmap='RdYlGn')

    plt.xlabel('Relative male fertility')
    plt.ylabel('Relative female fertility')
    # plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
    plt.grid(False)
    plt.savefig(out_folder+"/"+'male_female_scatter_best.pdf')

    plt.close()


### define half circle


def get_best_cross(mdf,fdf):

    cuts=np.linspace(-2, 0, num=20)
    shifts=np.linspace(0, 2, num=20)

    res_df=pd.DataFrame(columns=['cut','shift','m','f','b','wm','wf','sum'])
    idx=0

    for cut in cuts:
        for shift in shifts:
            lnum=train_cross(mdf,fdf,mcut=cut,mshift=shift,fcut=cut,fshift=shift)
            res_df.loc[idx,:]=lnum
            idx+=1

    res_df.to_csv(out_folder+'/cross_para.txt',sep='\t')



def train_cross(mdf,fdf,mcut=-0.42,mshift=1.36,fcut=-0.42,fshift=1.15):
    na_df=(mdf['rgr'].isna())| (fdf['rgr'].isna())
    mdf=mdf[~na_df].copy()
    fdf=fdf[~na_df].copy()


    # for male
    mdf['cut1']=mdf['rgr']+mshift*mdf['sd']
    mdf['cut2']=mdf['rgr']-mshift*mdf['sd']
    m_R=mdf['cut1']<mcut
    m_NR=mdf['cut2']>=mcut
    m_NP=~(m_R | m_NR)

    #
    fdf['cut1']=fdf['rgr']+fshift*mdf['sd']
    fdf['cut2']=fdf['rgr']-fshift*mdf['sd']
    f_R=fdf['cut1']<fcut
    f_NR=fdf['cut2']>=fcut
    f_NP=~(f_R | f_NR)

    # get all 9 categories
    cond1= m_R & f_R # super intersting
    cond2= m_R & f_NR # super
    cond3= m_R & f_NP
    cond4= m_NR & f_R
    cond5= m_NR & f_NR # super intersting
    cond6= m_NR & f_NP
    cond7= m_NP & f_R # super
    cond8= m_NP & f_NR
    cond9= m_NP & f_NP # boring

    ## see common between cross pheno

    m_cross=cross_df[cross_df['Cross']=='Yes (male-specific)']
    f_cross=cross_df[cross_df['Cross']=='Yes (female-specific)']
    b_cross=cross_df[cross_df['Cross']=='Yes (male- and female-specific)']
    male_specific_df=mdf[cond1|cond2|cond3]
    female_specific_df=fdf[cond1|cond4|cond7]
    both_specific_df=mdf[cond1]
    wild_fdf=fdf[cond5 | cond2|cond8]
    wild_mdf=mdf[cond5|cond4|cond6]

    comm_m=set(male_specific_df['feature_symbol'])&set(m_cross['Gene'])
    comm_f=set(female_specific_df['feature_symbol'])&set(f_cross['Gene'])
    comm_b=set(both_specific_df['feature_symbol'])&set(b_cross['Gene'])
    comm_w_male=set(wild_mdf['feature_symbol'])&set(f_cross['Gene'])
    comm_w_female=set(wild_fdf['feature_symbol'])&set(m_cross['Gene'])
    all_sum=len(comm_m)+len(comm_f)+len(comm_b)+len(comm_w_male)+len(comm_w_female)
    print(len(comm_m),len(comm_f),len(comm_b),len(comm_w_male),len(comm_w_female),all_sum)
    return mcut,mshift,len(comm_m),len(comm_f),len(comm_b),len(comm_w_male),len(comm_w_female),all_sum


def male_female_half_circle(mdf,fdf,mcut=-2,mshift=1,fcut=-2,fshift=0.96):
    #  # -0.42105263157894735	1.3684210526315788	120 male
    # -0.42105263157894735	1.1578947368421053	76 female
    # '#AF58BA'
    # '#009ADE'
    # '#A0B1BA'
    #
    # gst male female equal
    # '#d8b365' red
    # '#5ab4ac' not red
    # '#f5f5f5' no power
    na_df=(mdf['rgr'].isna())| (fdf['rgr'].isna())
    mdf=mdf[~na_df].copy()
    fdf=fdf[~na_df].copy()


    # for male
    mdf['cut1']=mdf['rgr']+mshift*mdf['sd']
    mdf['cut2']=mdf['rgr']-mshift*mdf['sd']
    m_R=mdf['cut1']<mcut
    m_NR=mdf['cut2']>=mcut
    m_NP=~(m_R | m_NR)

    #
    fdf['cut1']=fdf['rgr']+fshift*mdf['sd']
    fdf['cut2']=fdf['rgr']-fshift*mdf['sd']
    f_R=fdf['cut1']<fcut
    f_NR=fdf['cut2']>=fcut
    f_NP=~(f_R | f_NR)

    # get all 9 categories
    cond1= m_R & f_R # super intersting
    cond2= m_R & f_NR # super
    cond3= m_R & f_NP
    cond4= m_NR & f_R
    cond5= m_NR & f_NR # super intersting
    cond6= m_NR & f_NP
    cond7= m_NP & f_R # super
    cond8= m_NP & f_NR
    cond9= m_NP & f_NP # boring

    ## see common between cross pheno

    m_cross=cross_df[cross_df['Cross']=='Yes (male-specific)']
    f_cross=cross_df[cross_df['Cross']=='Yes (female-specific)']
    b_cross=cross_df[cross_df['Cross']=='Yes (male- and female-specific)']
    male_specific_df=mdf[cond1|cond2|cond3]
    female_specific_df=fdf[cond1|cond4|cond7]
    both_specific_df=mdf[cond1]
    wild_fdf=fdf[cond5 | cond2|cond8]
    wild_mdf=mdf[cond5|cond4|cond6]

    comm_m=set(male_specific_df['feature_symbol'])&set(m_cross['Gene'])
    comm_f=set(female_specific_df['feature_symbol'])&set(f_cross['Gene'])
    comm_b=set(both_specific_df['feature_symbol'])&set(b_cross['Gene'])
    comm_w_male=set(wild_mdf['feature_symbol'])&set(f_cross['Gene'])
    comm_w_female=set(wild_fdf['feature_symbol'])&set(m_cross['Gene'])
    all_sum=len(comm_m)+len(comm_f)+len(comm_b)+len(comm_w_male)+len(comm_w_female)
    print(len(comm_m),len(comm_f),len(comm_b),len(comm_w_male),len(comm_w_female),all_sum)
    f, ax = plt.subplots(figsize=[6.4*1.2,4.8*1.2])
    ## for cond 9
    marker_style = dict(color='#f5f5f5', linestyle='None', marker='o',markersize=6, markerfacecoloralt='#f5f5f5',markeredgecolor='None')
    p9=ax.plot(mdf[cond9]['rgr'], fdf[cond9]['rgr'], fillstyle='left', **marker_style,label=str(mdf[cond9].shape[0]))
    ## for cond 8
    marker_style = dict(color='#f5f5f5', linestyle='None', marker='o',markersize=6, markerfacecoloralt='#5ab4ac' ,markeredgecolor='None')
    p8=ax.plot(mdf[cond8]['rgr'], fdf[cond8]['rgr'], fillstyle='left', **marker_style,label=str(mdf[cond8].shape[0]))
    ## for cond 6
    marker_style = dict(color='#5ab4ac', linestyle='None', marker='o',markersize=6, markerfacecoloralt='#f5f5f5' ,markeredgecolor='None')
    p6=ax.plot(mdf[cond6]['rgr'], fdf[cond6]['rgr'], fillstyle='left', **marker_style,label=str(mdf[cond6].shape[0]))

    ## for cond7
    marker_style = dict(color='#f5f5f5', linestyle='None', marker='o',markersize=6, markerfacecoloralt='#d8b365' ,markeredgecolor='None')
    p7=ax.plot(mdf[cond7]['rgr'], fdf[cond7]['rgr'], fillstyle='left', **marker_style,label=str(mdf[cond7].shape[0]))

    ## for cond 3
    marker_style = dict(color='#d8b365', linestyle='None', marker='o',markersize=6, markerfacecoloralt='#f5f5f5' ,markeredgecolor='None')
    p3=ax.plot(mdf[cond3]['rgr'], fdf[cond3]['rgr'], fillstyle='left', **marker_style,label=str(mdf[cond3].shape[0]))

    ## for cond 5
    marker_style = dict(color='#5ab4ac', linestyle='None', marker='o',markersize=6, markerfacecoloralt='#5ab4ac' ,markeredgecolor='None')
    p5=ax.plot(mdf[cond5]['rgr'], fdf[cond5]['rgr'], fillstyle='left', **marker_style,label=str(mdf[cond5].shape[0]))

    ## for cond 1
    marker_style = dict(color='#d8b365', linestyle='None', marker='o',markersize=6, markerfacecoloralt='#d8b365' ,markeredgecolor='None')
    p1=ax.plot(mdf[cond1]['rgr'], fdf[cond1]['rgr'], fillstyle='left', **marker_style,label=str(mdf[cond1].shape[0]))

    ## for cond2
    marker_style = dict(color='#d8b365', linestyle='None', marker='o',markersize=6, markerfacecoloralt='#5ab4ac' ,markeredgecolor='None')
    p2=ax.plot(mdf[cond2]['rgr'], fdf[cond2]['rgr'], fillstyle='left', **marker_style,label=str(mdf[cond2].shape[0]))

    ## for cond4
    marker_style = dict(color='#5ab4ac' , linestyle='None', marker='o',markersize=6, markerfacecoloralt='#d8b365',markeredgecolor='None')
    p4=ax.plot(mdf[cond4]['rgr'], fdf[cond4]['rgr'], fillstyle='left', **marker_style,label=str(mdf[cond4].shape[0]))

    # plt.legend(handles=[p2,p4,p1,p5,p3,p7,p6,p8,p9], title='Groups', loc='best')
    handles, labels = ax.get_legend_handles_labels()
    # reorders
    idx=[8,7,6,5,4,3,2,1,0]
    new_hadles=[ handles[i] for i in idx]
    new_labels=[  labels [i] for i in idx]
    # import pdb; pdb.set_trace()
    # ax.legend(new_hadles, new_labels)
    plt.legend(handles=new_hadles,labels=new_labels,title='Groups', bbox_to_anchor=(1.01, 1),borderaxespad=0)

    # now we want to plot half circle
    plt.xlabel('Relative male fertility')
    plt.ylabel('Relative female fertility')
    # plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
    plt.grid(False)
    plt.savefig(out_folder+"/"+'male_female_scatter_half_circle.pdf')

    plt.close()


def associate_fertility_phenotype(mdf,fdf):
    ''' what is fertility phenottpe'''
    mcut=-2
    mshift=1
    fcut=-2
    fshift=0.96

    na_df=(mdf['g145480_RGR'].isna())| (fdf['GCKO2_RGR'].isna())
    mdf=mdf[~na_df].copy()
    fdf=fdf[~na_df].copy()


    # for male
    mdf['cut1']=mdf['g145480_RGR']+mshift*mdf['g145480_sd']
    mdf['cut2']=mdf['g145480_RGR']-mshift*mdf['g145480_sd']
    m_R=mdf['cut1']<mcut
    m_NR=mdf['cut2']>=mcut
    m_NP=~(m_R | m_NR)

    #
    fdf['cut1']=fdf['GCKO2_RGR']+fshift*mdf['GCKO2_sd']
    fdf['cut2']=fdf['GCKO2_RGR']-fshift*mdf['GCKO2_sd']
    f_R=fdf['cut1']<fcut
    f_NR=fdf['cut2']>=fcut
    f_NP=~(f_R | f_NR)

    # get male_female_pheno
    pbankaids=set(mdf['pbanka_id'])|set(fdf['pbanka_id'])

    fertility_df=pd.DataFrame(index=pbankaids, columns=['fertility_male_rgr', 'fertility_male_sd','fertility_male_pheno','fertility_female_rgr',
     'fertility_female_sd','fertility_female_pheno'])

    ## assign male values
    fertility_df.loc[mdf[m_R]['pbanka_id'].to_list(),'fertility_male_rgr']=mdf[m_R]['g145480_RGR'].to_list()
    fertility_df.loc[mdf[m_R]['pbanka_id'].to_list(),'fertility_male_sd']=mdf[m_R]['g145480_sd'].to_list()
    fertility_df.loc[mdf[m_R]['pbanka_id'].to_list(),'fertility_male_pheno']='reduced'
    #
    fertility_df.loc[mdf[m_NR]['pbanka_id'].to_list(),'fertility_male_rgr']=mdf[m_NR]['g145480_RGR'].to_list()
    fertility_df.loc[mdf[m_NR]['pbanka_id'].to_list(),'fertility_male_sd']=mdf[m_NR]['g145480_sd'].to_list()
    fertility_df.loc[mdf[m_NR]['pbanka_id'].to_list(),'fertility_male_pheno']='notreduced'
    #
    fertility_df.loc[mdf[m_NP]['pbanka_id'].to_list(),'fertility_male_rgr']=mdf[m_NP]['g145480_RGR'].to_list()
    fertility_df.loc[mdf[m_NP]['pbanka_id'].to_list(),'fertility_male_sd']=mdf[m_NP]['g145480_sd'].to_list()
    fertility_df.loc[mdf[m_NP]['pbanka_id'].to_list(),'fertility_male_pheno']='nopower'

     # assign female values
    fertility_df.loc[fdf[f_R]['pbanka_id'].to_list(),'fertility_female_rgr']=fdf[f_R]['GCKO2_RGR'].to_list()
    fertility_df.loc[fdf[f_R]['pbanka_id'].to_list(),'fertility_female_sd']=fdf[f_R]['GCKO2_sd'].to_list()
    fertility_df.loc[fdf[f_R]['pbanka_id'].to_list(),'fertility_female_pheno']='reduced'
    #
    fertility_df.loc[fdf[f_NR]['pbanka_id'].to_list(),'fertility_female_rgr']=fdf[f_NR]['GCKO2_RGR'].to_list()
    fertility_df.loc[fdf[f_NR]['pbanka_id'].to_list(),'fertility_female_sd']=fdf[f_NR]['GCKO2_sd'].to_list()
    fertility_df.loc[fdf[f_NR]['pbanka_id'].to_list(),'fertility_female_pheno']='notreduced'
    #
    fertility_df.loc[fdf[f_NP]['pbanka_id'].to_list(),'fertility_female_rgr']=fdf[f_NP]['GCKO2_RGR'].to_list()
    fertility_df.loc[fdf[f_NP]['pbanka_id'].to_list(),'fertility_female_sd']=fdf[f_NP]['GCKO2_sd'].to_list()
    fertility_df.loc[fdf[f_NP]['pbanka_id'].to_list(),'fertility_female_pheno']='nopower'
    fertility_df.to_csv(out_folder+'/fertility_new.txt',sep='\t')

if __name__ == "__main__":
    mca_df,male_df,female_df,male_mca_df,female_mca_df=get_data()
    # male_female_top50(male_mca_df,female_mca_df,'M',cut_gr=-2,cut_less=-2)
    male_top50_cluster(male_mca_df,female_mca_df)
    # get_best_cross(male_mca_df,female_mca_df)
    enrich_mca_cluster(mca_df,male_mca_df,female_mca_df,cut_gr=-2,cut_less=-2)
    import pdb; pdb.set_trace()
    associate_fertility_phenotype(male_df,female_df)
    male_female_half_circle(male_mca_df,female_mca_df)
    # import pdb; pdb.set_trace()

    #res_df=sensitivity_analysis(male_mca_df,female_mca_df)
    color_group_genes(male_mca_df,female_mca_df)

    plot_male_female_scatter(male_mca_df,female_mca_df)

    res_df=pickle.load(open(out_folder+'/'+'optim_n.p','rb'))
    res_sum=res_df.sort_values(by='sum',ascending=False)
    res_male=res_df.sort_values(by='m_value',ascending=False)
    res_female=res_df.sort_values(by='f_value',ascending=False)
    res_both=res_df.sort_values(by='b_value',ascending=False)

    #-0.8421052631578947	1.0526315789473684
    # -0.42105263157894735	1.3684210526315788	120 male
    # -0.42105263157894735	1.1578947368421053	76 female

    best_m_f_plot(male_mca_df,female_mca_df,cut=-0.84,shift=1)
    plot_male_female_scatter_best(male_mca_df,female_mca_df,cut=-0.84,shift=1)



    male_female_only(male_mca_df,female_mca_df,'M',cut_gr=-2,cut_less=-2)
    male_female_only(male_mca_df,female_mca_df,'F',cut_gr=-2,cut_less=-2)

    plot_female_male_on_mca(male_mca_df,female_mca_df)

    # change_in_male_female(male_mca_df,female_mca_df)


    get_category_colors(mca_df)

    color_metabolic_genes(male_mca_df,female_mca_df)
    color_group_genes(male_mca_df,female_mca_df)
    # enrichment analysis

# colors=get_color(male_mca_df['rgr'])
#
# f, ax = plt.subplots()
# # points = ax.scatter(x=male_mca_df['knn_graph_X'], y=male_mca_df['knn_graph_Y'], c=colors, s=6, cmap='Spectral')
# points = ax.scatter(x=male_mca_df['knn_graph_X'], y=male_mca_df['knn_graph_Y'], c=colors, s=6, cmap='Spectral')
# import pdb; pdb.set_trace()
# f.colorbar(points)
# # plt.plot(x*0.1, y, 'o-', color='lightgrey', label='No mask')
# # plt.plot(x2*0.4, y2, 'o-', label='Points removed')
# # plt.plot(x*0.7, y3, 'o-', label='Masked values')
# # plt.plot(x*1.0, y4, 'o-', label='NaN values')
# # plt.legend()
# # plt.title('Masked and NaN data')
# # plt.show()
# # plt.xlabel('knn_graph_X')
# # plt.ylabel('knn_graph_Y')
# #
# plt.savefig(out_folder+"/mca_gene_cluster_male_color_gradient.pdf")
# plt.close()
#
# colors=get_color(female_mca_df['rgr'])
# f, ax = plt.subplots()
# points = ax.scatter(x=female_mca_df['knn_graph_X'], y=female_mca_df['knn_graph_Y'], c=colors, s=6, cmap='Spectral')
# f.colorbar(points)
# plt.xlabel('knn_graph_X')
# plt.ylabel('knn_graph_Y')
#
#
#
# plt.savefig(out_folder+"/mca_gene_cluster_female_color_gradient.pdf")
# plt.close()
#
# male_mca_df.to_csv(out_folder+"/mca_gene_data_frame_male.txt",sep='\t')
# female_mca_df.to_csv(out_folder+"/mca_gene_data_frame_female.txt",sep='\t')





# sns.set()
#
# sns_plot=sns.scatterplot(data=male_mca_df, x="knn_graph_X", y="knn_graph_Y", hue="rgr",s=6)
# norm = plt.Normalize(male_mca_df['rgr'].min(), male_mca_df['rgr'].max())
# sm = plt.cm.ScalarMappable(cmap="RdBu", norm=norm)
# sm.set_array([])
#
# # Remove the legend and add a colorbar
# sns_plot.get_legend().remove()
# sns_plot.figure.colorbar(sm)
#
# # plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
# plt.savefig(out_folder+"/mca_gene_cluster_male_color_gradient.pdf")
# plt.close()
#
# sns_plot=sns.scatterplot(data=female_mca_df, x="knn_graph_X", y="knn_graph_Y", hue="rgr",s=6)
# norm = plt.Normalize(female_mca_df['rgr'].min(), female_mca_df['rgr'].max())
# sm = plt.cm.ScalarMappable(cmap="RdBu", norm=norm)
# sm.set_array([])
#
# # Remove the legend and add a colorbar
# sns_plot.get_legend().remove()
# sns_plot.figure.colorbar(sm)
#
# # plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
# plt.savefig(out_folder+"/mca_gene_cluster_female_color_gradient.pdf")
# plt.close()
