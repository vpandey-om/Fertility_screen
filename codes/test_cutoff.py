import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import numpy as np
import matplotlib.backends.backend_pdf
from bioinfokit import analys, visuz
df=pd.read_csv('/Users/vpandey/projects/githubs/Fertility_screen/data/wt21.txt',sep='\t')
tmp=df[['GCKO2_RGR','GCKO2_diff_min','g145480_RGR','g145480_diff_min']]
# tmp.hist(bins=5)
#plt.show()
malaria_atlas_list= pickle.load(open('/Users/vpandey/projects/gitlabs/evolution_rate/malaria_cell_atlas.pickle', 'rb'))
#ipbe=pickle.load(open('/Users/vpandey/projects/data/ipbeblood_model_info.pickle','wb'))


def get_subset_data(df,col_name,gene_list):
    ## find common genes
    comm=list(set(df[col_name])&set(gene_list))
    list_tmp=[]
    for gene in comm:
        tmp=df[df[col_name]==gene]
        if not tmp.empty:
                list_tmp.append(tmp)
    tmp_df=pd.concat(list_tmp)
    return tmp_df



def plot_diffmin_max(df,sex='GCKO2',file=None):
    ''' plot variance with nod data excluded '''
    plt.clf()
    xx=df[~(df[sex+'_pheno']=='No Data')]
    res_df=xx.copy()
    # yy=res_df[sex+'_diff_min'].values
    # xx=res_df[sex+'_diff_max'].values
    yy=res_df[sex+'_sd'].values
    xx=res_df[sex+'_RGR'].values

    # xx[xx<1e-10]=1e-10;
    # xx=np.log2(xx)
    # df=pd.DataFrame()
    # df['sd']=xx
    # df['RGR']=yy
    plt.scatter(xx, yy, marker='o', s=2)
    # sns_plot=sns.scatterplot(x=xx, y=yy)
    plt.xlim([-13, 2])
    # plt.ylim([-18, 8])
    plt.ylim([0, 5])
    plt.legend([sex+'(genes=%s)'%len(xx)])
    # plt.xlabel('diff_max')
    # plt.ylabel('diff_min')
    plt.xlabel('RGR')
    plt.ylabel('SD')
    # figure = plt.get_figure()

    #plt.show()
    if file:
        plt.savefig(file, dpi=200)


def plot_variance(df,sex='GCKO2',percentile=80,file=None):
    ''' plot variance with nod data excluded '''
    plt.clf()
    xx=df[~(df[sex+'_pheno']=='No Data')]
    res_df=xx.copy()
    xx=xx[sex+'_sd'].values
    xx[xx<1e-10]=1e-10;
    xx=np.log2(xx)
    sns_plot=sns.distplot(xx, hist=True, rug=True,color='blue')
    plt.legend([sex+'(genes=%s)'%len(xx)])
    plt.xlabel('log2(sd)')
    figure = sns_plot.get_figure()



    #plt.show()
    if file:
        figure.savefig(file, dpi=200)

    ## get original_df
    cutoff=np.percentile(xx, percentile);
    ## filter noisy data
    filter_df=res_df[xx<cutoff]
    return filter_df,res_df





def getfigMalariaAtlasGene(genes,out_pdf):

    ## this data is for malaria cell atlas visualizations
    umap_df=malaria_atlas_list[0] # umap_df
    pheno_df=malaria_atlas_list[1] ## pheno_df
    atlas_sparse_data=malaria_atlas_list[2] ### sparse numpy array
    atlas_genes=malaria_atlas_list[3]
    ### find index for a gene
    pdf = matplotlib.backends.backend_pdf.PdfPages(out_pdf)
    # umap_df=umap_df.replace({'Male':'Male_gametocyte'})

    for gene in genes:
        fig = plt.figure()
        if len(atlas_genes[atlas_genes==gene].to_list())>0:
            umap_df['expression']=np.nan
            arr=atlas_sparse_data.todense()
            umap_df['expression']=arr[atlas_genes==gene,:].tolist()[0]
            points = plt.scatter(umap_df["umap0"], umap_df["umap1"],
                     c=umap_df["expression"], s=10,cmap="Reds")
            plt.colorbar(points)
            # fig = sns.scatterplot(data=umap_df, x='umap0', y='umap1', hue='expression')


        else:
            umap_df['expression']=1
            #fig = sns.scatterplot(data=umap_df, x='umap0', y='umap1', hue='expression')
            points = plt.scatter(umap_df["umap0"], umap_df["umap1"],
                     c=umap_df["expression"], s=10,cmap="Reds")
            plt.colorbar(points)

        plt.xlabel('umap1')
        plt.ylabel('umap2')
        plt.title("%s"%gene)
        pdf.savefig(fig)
    pdf.close()



def getMetabolicGenes():
    ''' we want to get metabolic genes data '''
    ### we
    print('hi')





### import excel file
male_gene=pd.read_excel('/Users/vpandey/projects/githubs/Fertility_screen/data/Data for phenotype calls.xlsx', sheet_name='Male_Core',index_col=None)
female_gene=pd.read_excel('/Users/vpandey/projects/githubs/Fertility_screen/data/Data for phenotype calls.xlsx', sheet_name='Female_Core',index_col=None)


# male_gam=pd.read_excel('/Users/vpandey/projects/githubs/Fertility_screen/data/Data for phenotype calls.xlsx', sheet_name='Male_Gametocyte',index_col=None)
# female_gam=pd.read_excel('/Users/vpandey/projects/githubs/Fertility_screen/data/Data for phenotype calls.xlsx', sheet_name='Female_Gametocyte',index_col=None)

oocy_not_reduced=pd.read_excel('/Users/vpandey/projects/githubs/Fertility_screen/data/Data for phenotype calls.xlsx', sheet_name='Oocyst_NotReduced',index_col=None)

# wt_oo=pd.read_excel('/Users/vpandey/projects/githubs/Fertility_screen/data/Data for phenotype calls.xlsx', sheet_name='WT_oo',index_col=None) # wild type oocyst
#
#
final_df=pd.read_excel('/Users/vpandey/projects/githubs/Fertility_screen/analysisCodes/Phenotype_call_final_080920.xlsx',sheet_name='data')
#
# ### old_df
# all_screen=pd.read_csv('/Users/vpandey/projects/githubs/Fertility_screen/Figures/all_screening_data_save.txt',sep='\t')

### get data from


############### final analysis of cutoff

# plot male and female variances
folder='/Users/vpandey/projects/githubs/Fertility_screen/analysisCodes/Figures/'
plot_diffmin_max(final_df,sex='GCKO2',file=folder+'GCKO2_RGRvsSD.pdf') ## filter based on percentile
plot_diffmin_max(final_df,sex='g145480',file=folder+'g145480_RGRvsSD.pdf')
import pdb; pdb.set_trace()
## get common set analysis
male_gene_df=get_subset_data(male_filter_df,'pbanka_id', male_gene['Gene'].to_list())
female_gene_df=get_subset_data(female_filter_df,'pbanka_id', female_gene['Gene'].to_list())
male_oocy_not_reduced_df=get_subset_data(male_filter_df,'pbanka_id', oocy_not_reduced['Gene'].to_list())
female_oocy_not_reduced_df=get_subset_data(female_filter_df,'pbanka_id', oocy_not_reduced['Gene'].to_list())
# we are plotting for male

## plot diffminvs diffmax
female_filter_df,res_df_male=plot_variance(final_df,sex='GCKO2',percentile=80,file=folder+'GCKO2_sd_dist.pdf') ## filter based on percentile
male_filter_df,res_df_female=plot_variance(final_df,sex='g145480',percentile=80,file=folder+'g145480_sd_dist.pdf')



plt.clf()
sns_plot=sns.distplot(male_gene_df['g145480_RGR'].values, hist=True, rug=True,color='red')
sns_plot=sns.distplot(male_oocy_not_reduced_df['g145480_RGR'].values, hist=True, rug=True,color='blue')
plt.legend(['Core genes g145480 (genes=%d)'%male_gene_df.shape[0],'Oocyst WT genes g145480(genes=%d)'%male_oocy_not_reduced_df.shape[0]])
plt.xlabel('g145480_RGR')
figure = sns_plot.get_figure()
figure.savefig(folder+'g145480_cutoff.pdf', dpi=200)

## this is for female
plt.clf()
sns_plot=sns.distplot(female_gene_df['GCKO2_RGR'].values, hist=True, rug=True,color='red')
sns_plot=sns.distplot(female_oocy_not_reduced_df['GCKO2_RGR'].values, hist=True, rug=True,color='blue')
plt.legend(['Core genes GCKO2 (genes=%d)'%female_gene_df.shape[0],'Oocyst WT genes GCKO2(genes=%d)'%female_oocy_not_reduced_df.shape[0]])
plt.xlabel('GCKO2_RGR')
figure = sns_plot.get_figure()
figure.savefig(folder+'GCKO2_cutoff.pdf', dpi=200)







#####



# male_gene_df=get_subset_data(final_df,'pbanka_id', male_gene['Gene'].to_list(),sex='g145480')
# female_gene_df=get_subset_data(final_df,'pbanka_id', female_gene['Gene'].to_list(),sex='GCKO2')
# #actual_genes=['PBANKA_1035200','PBANKA_0704900','PBANKA_1436600','PBANKA_1319500','PBANKA_1112700','PBANKA_1225500','PBANKA_1233600',
# #'PBANKA_1444800','PBANKA_1436300','PBANKA_0517600','PBANKA_0105100','PBANKA_1359700','PBANKA_1323800','PBANKA_1335000']
# # female_gene_df=get_subset_data(final_df,'pbanka_id', actual_genes)
# male_gam_df=get_subset_data(all_screen,'pbanka_id', male_gam['Gene'].to_list())
# female_gam_df=get_subset_data(all_screen,'pbanka_id', female_gam['Gene'].to_list())
# oocy_not_reduced_df=get_subset_data(final_df,'pbanka_id', oocy_not_reduced['Gene'].to_list())
# wt_oo_df=get_subset_data(final_df,'pbanka_id', wt_oo['Gene'].to_list())
# # wt_female=df[['GCKO2_diff_min']]
# # wt_male=df[['g145480_diff_min']]
# # male_test=male_gene_df['g145480_diff_min']
# # female_test=female_gene_df['GCKO2_diff_min']
#
#
# wt_female=df[['GCKO2_RGR']]
# wt_male=df[['g145480_RGR']]
# male_test=male_gene_df['g145480_RGR']
# female_test=female_gene_df['GCKO2_RGR']
# # oocy_not_reduced_test_male=oocy_not_reduced_df['g145480_RGR']
# # oocy_not_reduced_test_female=oocy_not_reduced_df['GCKO2_RGR']
# # oocy_not_reduced_test_male=oocy_not_reduced_df['g145480_diff_min']
# # oocy_not_reduced_test_female=oocy_not_reduced_df['GCKO2_diff_min']
#
# oocy_not_reduced_test_male=wt_oo_df['g145480_diff_min']
# oocy_not_reduced_test_female=wt_oo_df['GCKO2_diff_min']
#
#
# sns_plot=sns.distplot(wt_female, hist=True, rug=True,color='blue')
# sns_plot=sns.distplot(female_test, hist=True, rug=True,color='red')
#
#
# plt.legend(['Wild_type_female','core_gene_female'])
# figure = sns_plot.get_figure()
# figure.savefig('female_cutoff.png', dpi=200)
# # plt.show()
#
# plt.clf()
# sns_plot=sns.distplot(wt_male, hist=True, rug=True,color='blue')
# sns_plot=sns.distplot(male_test, hist=True, rug=True,color='red')
#
#
# plt.legend(['Wild_type_male','core_gene_male'])
# figure = sns_plot.get_figure()
# figure.savefig('male_cutoff.png', dpi=200)
#
# ###
# plt.clf()
# sns_plot=sns.distplot(female_test, hist=True, rug=True,color='red')
# sns_plot=sns.distplot(oocy_not_reduced_test_female, hist=True, rug=True,color='blue')
#
#
#
# ### we want to plot female genes
# genes=female_gene_df['pbanka_id'].to_list()
# out_pdf='/Users/vpandey/projects/githubs/Fertility_screen/Fig_pilot/female_core_genes.pdf'
#
# import pdb; pdb.set_trace()
# getfigMalariaAtlasGene(genes,out_pdf)
#
# ###
