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
import scanpy as sc
import math
code=os.getcwd()
upLevel=code.replace('codes','') ####### we are going to upper level of code directory
sys.path.insert(0,upLevel+'/data')
sys.path.insert(1, upLevel+'/Figures')

import matplotlib.pyplot as plt

import matplotlib
import matplotlib.cm as cm
from copy import copy





# print(sys.path)


# input files which will be needed for analysis of cliare screen data

data_folder=sys.path[0]

## output folder where we will write figures and output files
out_folder=sys.path[1]

# ID conversion: we used plasmoDB to convert ID for P. Berghai
mca_df=pd.read_excel(open(data_folder+'/ginny_gene_mca.xlsx','rb'), sheet_name='mca')
obs=mca_df[['feature_symbol','product_description','Cluster']]
obs.set_index('feature_symbol',inplace=True)

fertility_df=pd.read_excel(open(data_folder+'/Phenotype_call_plasmogem_id_Vikash_271020.xlsx','rb'),sheet_name='fertility')
female_df=fertility_df[~(fertility_df['GCKO2_pheno']=='No Data')]
male_df=fertility_df[~(fertility_df['g145480_pheno']=='No Data')]

var=pd.DataFrame(index=['knn_graph_X','knn_graph_Y'],columns=['co-ordinates'])
var.loc['knn_graph_X','co-ordinates']='x'
var.loc['knn_graph_Y','co-ordinates']='y'
adata = sc.AnnData(X=mca_df[['knn_graph_X','knn_graph_Y']].to_numpy(),obs=obs,var=var)
mca_df['Cluster']=mca_df['Cluster'].astype('category')

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
    tmp=female_df[female_df['pbanka_id']==female_mca_df.loc[idx,'feature_symbol']]
    if not tmp.empty:
        female_mca_df.loc[idx,'pheno']=tmp['GCKO2_pheno'].to_list()[0]
        female_mca_df.loc[idx,'rgr']=tmp['GCKO2_RGR'].to_list()[0]


#import pdb; pdb.set_trace()
# color gradient




sns_plot=sns.scatterplot(data=mca_df, x="knn_graph_X", y="knn_graph_Y", hue="Cluster",s=6)
plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
plt.savefig(out_folder+"/mca_gene_cluster.pdf")
plt.close()

# Create an array with the colors you want to use
colorsIdx = {'No Power': '#cccccc', 'Not Reduced': '#b2e2e2','Reduced':'#d8b365'}
# colors = ["#f7f7f7", "#b2e2e2",'#d8b365','#cccccc']
colors = ["#f7f7f7", "#99d594",'#fc8d59','#ffffbf']
# Set your custom color palette
customPalette = sns.set_palette(sns.color_palette(colors))


sns_plot=sns.scatterplot(data=male_mca_df, x="knn_graph_X", y="knn_graph_Y", hue="pheno",s=6)
plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
plt.savefig(out_folder+"/mca_gene_cluster_male.pdf")
plt.close()

sns_plot=sns.scatterplot(data=female_mca_df, x="knn_graph_X", y="knn_graph_Y", hue="pheno",s=6)
plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
plt.savefig(out_folder+"/mca_gene_cluster_female.pdf")
plt.close()



def mask_dataframe_table(df):
    ## get nan index
    import pdb; pdb.set_trace()










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

mask_dataframe_table(male_mca_df)

colors=get_color(male_mca_df['rgr'])

f, ax = plt.subplots()
# points = ax.scatter(x=male_mca_df['knn_graph_X'], y=male_mca_df['knn_graph_Y'], c=colors, s=6, cmap='Spectral')
points = ax.scatter(x=male_mca_df['knn_graph_X'], y=male_mca_df['knn_graph_Y'], c=colors, s=6, cmap='Spectral')
import pdb; pdb.set_trace()
f.colorbar(points)
# plt.plot(x*0.1, y, 'o-', color='lightgrey', label='No mask')
# plt.plot(x2*0.4, y2, 'o-', label='Points removed')
# plt.plot(x*0.7, y3, 'o-', label='Masked values')
# plt.plot(x*1.0, y4, 'o-', label='NaN values')
# plt.legend()
# plt.title('Masked and NaN data')
# plt.show()
# plt.xlabel('knn_graph_X')
# plt.ylabel('knn_graph_Y')
#
plt.savefig(out_folder+"/mca_gene_cluster_male_color_gradient.pdf")
plt.close()

colors=get_color(female_mca_df['rgr'])
f, ax = plt.subplots()
points = ax.scatter(x=female_mca_df['knn_graph_X'], y=female_mca_df['knn_graph_Y'], c=colors, s=6, cmap='Spectral')
f.colorbar(points)
plt.xlabel('knn_graph_X')
plt.ylabel('knn_graph_Y')



plt.savefig(out_folder+"/mca_gene_cluster_female_color_gradient.pdf")
plt.close()

male_mca_df.to_csv(out_folder+"/mca_gene_data_frame_male.txt",sep='\t')
female_mca_df.to_csv(out_folder+"/mca_gene_data_frame_female.txt",sep='\t')





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
