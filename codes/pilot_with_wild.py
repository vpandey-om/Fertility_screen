import os
import sys
import plotly.express as px
import seaborn as sns
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np
import statsmodels.stats.multitest as mtest

code=os.getcwd()
upLevel=code.replace('codes','') ####### we are going to upper level of code directory
sys.path.insert(0,upLevel+'/data')
sys.path.insert(1, upLevel+'/Figures')

# print(sys.path)

from special_fun_fertility  import *
# input files which will be needed for analysis of cliare screen data

data_folder=sys.path[0]

## output folder where we will write figures and output files
out_folder=sys.path[1]

# ID conversion: we used plasmoDB to convert ID for P. Berghai
prev_to_new=pickle.load(open(data_folder+'/prevTonew_PBANKA.pickle','rb'))
db_df=pd.read_csv(data_folder+'/PBANKA_id_conversion.txt', sep='\t')
db_df=db_df.fillna('NA')



def pre_process_pilot(manifests_df,count_df):

    ''' We will take average of the two reads '''
    # manifests_df.set_index('NGI Sample ID',inplace=True)
    req_df=count_df.copy()
    manifests_df=manifests_df.fillna('NA')

    manifests_df.set_index('Sample description',inplace=True)
    final_df=pd.DataFrame(index=req_df.index,columns=['Gene',  'Barcodes'])

    for ngi in manifests_df.index:

       ### get two reads
       # sam= manifests_df.loc[ngi,'Sample description'] # sample description
       sam=ngi
       two_reads=req_df.columns[req_df.columns.str.contains(sam)]
       if len(two_reads)==2:
           # final_df[ngi]=np.nan()
           final_df.loc[:,sam]=req_df.loc[:,two_reads].mean(axis=1)
           # final_df_two_read[ngi+'.read1']=np.nan()
           # final_df_two_read[ngi+'.read2']=np.nan()
           # final_df_two_read.loc[:,sam+'.read1']=req_df.loc[:,two_reads[0]]
           # final_df_two_read.loc[:,sam+'.read2']=req_df.loc[:,two_reads[1]]
       elif len(two_reads)==1:
           # final_df[ngi]=np.nan()
           final_df.loc[:,sam]=req_df.loc[:,two_reads].copy()
           print('Please check there is only one reads for the sample: %s'%sam)
       else:
           print('Number of reads are mismatched for the sample: %s'%sam)

    return final_df,manifests_df


def relative_growth_rate_analysis_pilot(df,manfest_df,prev_to_new,db_df,plot_info=None):
    ''' We are going to do relative growth rate analysis'''

    ## first we need to propagate error from PCR to mosquito feed

    rel_df=df.copy()
    rel_df=rel_df.drop(columns=['Gene','Barcodes'])
    rel_df= rel_df + 1
    rel_df=rel_df.div(rel_df.sum(axis=0), axis=1)
    rel_df_log=np.log2(rel_df)


    ### convert old pbanka to new ids
    geneConv,old_to_new_ids,geneConv_new=getNewIdfromPrevID(rel_df.index,prev_to_new,db_df)
    rel_df=rel_df.rename(old_to_new_ids, axis='index')

    # newIds for control genes
    control_genes= plot_info['control_genes']
    ctrl_genes=[]
    for c in control_genes:
        if c in prev_to_new.keys():
            ctrl_genes.append(prev_to_new[c])
        else:
            ctrl_genes.append(c)

    grp_cols=['sex','d','mf','dr','b','t']
    day_pos=grp_cols.index('d')

    mean_df_d0_mf1,var_df_d0_mf1,mean_df_d0_mf2,var_df_d0_mf2=propagate_error_day0_each_mossifeed(rel_df,manfest_df,grp_cols,day_pos)
    grp_cols=['sex','d','mf','dr','b','t']
    day_pos=grp_cols.index('d')
    mean_df_d13_mf1,var_df_d13_mf1,mean_df_d13_mf2,var_df_d13_mf2=propagate_error_day13_each_mossifeed(rel_df,manfest_df,grp_cols,day_pos)

    ##  we are going to compute rleative growth rate  ###

    ### plot propagated relative abundance.
    propagated_relative_baundance_plot(mean_df_d0_mf1,var_df_d0_mf1,mean_df_d0_mf2,var_df_d0_mf2,mean_df_d13_mf1,var_df_d13_mf1,mean_df_d13_mf2,var_df_d13_mf2,plot_info)
    mf1_RGR,mf1_var=calculate_RGR(mean_df_d0_mf1.copy(),var_df_d0_mf1.copy(),mean_df_d13_mf1.copy(),var_df_d13_mf1.copy(),ctrl_genes)
    mf2_RGR,mf2_var=calculate_RGR(mean_df_d0_mf2.copy(),var_df_d0_mf2.copy(),mean_df_d13_mf2.copy(),var_df_d13_mf2.copy(),ctrl_genes)
    ### now combined fitness for mf1 and mf2
    # take mf1 and mf2 in one dtaframe


    cmb_fitness={}
    backgrounds=['GCKO2','g145480','Cl15cy']

    for b in backgrounds:
        rgr_temp=pd.DataFrame(index=mf1_RGR.index,columns=['mf1','mf2'])
        var_temp=pd.DataFrame(index=mf1_RGR.index,columns=['mf1','mf2'])
        col_mf1 = getColumnsFormDF(mf1_RGR, [b])
        col_mf2 = getColumnsFormDF(mf2_RGR, [b])

        rgr_temp.loc[:,'mf1']=mf1_RGR.loc[:,col_mf1[0].to_list()[0]].copy()
        rgr_temp.loc[:,'mf2']=mf2_RGR.loc[:,col_mf2[0].to_list()[0]].copy()
        var_temp.loc[:,'mf1']=mf1_var.loc[:,col_mf1[0].to_list()[0]].copy()
        var_temp.loc[:,'mf2']=mf2_var.loc[:,col_mf2[0].to_list()[0]].copy()


        cmb_fitness[b]=gaussianMeanAndVariance(rgr_temp,var_temp)

    ## calculate combined fitness

    # pheno_call_df=getPvalZscore(cmb_fitness,upcut=1,lowcut=0.4,pval=0.05,pval1=0.05)
    pheno_call_df=applyPhenocall_CI(cmb_fitness,lower_cut=-1,upper_cut=1)

    return pheno_call_df, [mean_df_d0_mf1,mean_df_d0_mf2]


def pilot():
    ''' We are going do analyis of pilot study for Claire data'''

    ### these are the input files

    manifests_df=pd.read_csv(data_folder+"/data_pilot/manifest_all_with_wild_type.txt",sep='\t')

    d28573=pd.read_csv(data_folder+"/data_pilot/counts_28573.csv")
    d28551=pd.read_csv(data_folder+"/data_pilot/counts_28551.csv")
    ### combine two data frames
    dfn = pd.merge(d28573, d28551, on='gene', how='outer')  # combined first and second files
    df_modi = dfn  #
    df_modi = df_modi.drop(df_modi.index[0])
    df_modi = df_modi.drop(['barcode_x', 'barcode_y'], axis=1)  # drop out some unnecessary columns
    df_modi.set_index(['gene'], inplace=True)
    count_df=df_modi.copy()
    input_df=pd.read_csv(data_folder+'/data_pilot/PbSTM168_cross_phenotypes_final.csv', sep=';')

    comm_genes=set(input_df['Gene_ID'])&set(count_df.index)

    filter_count_df=count_df.loc[comm_genes,:].copy()

    ### filter count df for pilot



    #### end of the input section
    # final_count_df: read1 and read2 are added
    # final_count_df_two_read: reads are sperated
    # manfest_df: maifest_df
    final_df,manfest_df=pre_process_pilot(manifests_df,filter_count_df)

    ## remove genes
    remove_genes=['PBANKA_051200']

    filtered_count_df = final_df.drop(remove_genes)
    filtered_count_df.to_csv(out_folder+"/filterd_count_matrix_pilot_wild.txt",sep='\t')

    ### we are going to perform relative abundance analysis
    ## prev_to_new this is the pickle information which is used when we change old to new ID
    ## db_df: this is the dataframe contains name and description

    ## if you we do not want to plot then plot_info=None
    #plot_info={'pdf':out_folder+"/relative_abundance_of_pool1.pdf",'d':['d0','d13'],'mf':['mf1','mf2'],'sex':['GCKO2','g145480']}
    #plot_info=None
    # relative_abundance_analysis(filtered_count_df,manfest_df,prev_to_new,db_df,plot_info)

    ## we will do diffrent kind of error analysis

    #error_analysis(filtered_count_df,manfest_df,prev_to_new,db_df)
    geneConv,old_to_new_ids,geneConv_new=getNewIdfromPrevID(filtered_count_df.index,prev_to_new,db_df)
    plot_info={'pool':'pilot','file':out_folder+'/pilot_repeat.xlsx','rel_file':out_folder+'/pilot_propagated_error_relative_abundance.pdf','d':['d0','d13'],
    'mf':['mf1','mf2'],'sex':['GCKO2','g145480','Cl15cy13'],'geneConv':geneConv_new,
    'control_genes':['PBANKA_102460','PBANKA_100210','PBANKA_100220','PBANKA_050120']}

    pheno_call_df,input_df=relative_growth_rate_analysis_pilot(filtered_count_df,manfest_df,prev_to_new,db_df,plot_info)

    return pheno_call_df,input_df


if __name__ == '__main__':
    pub_cross=pd.read_excel(data_folder+"/data_pilot/Published_crosses_07122021.xlsx")
    pheno_call_df,input_df=pilot()

    cross=[]
    for gene in pheno_call_df.index.to_list():
        tmp=pub_cross[pub_cross['Gene_ID_new']==gene]

        if tmp.empty:
            cross.append('NA')
        else:

            if tmp['Female'].to_list()[0]=='F':
                cross.append('F')
            elif tmp['Male'].to_list()[0]=='M':
                cross.append('M')
            elif tmp['Female_male'].to_list()[0]=='FM':
                cross.append('FM')
            else:
                print('see')





    RGR_columns=['GCKO2_RGR','g145480_RGR','Cl15cy_RGR']
    RGRs=[]
    Sexes=[]
    publishedcross=[]
    tmpRGRs=[]
    tmpSexes=[]
    tmpCross=[]
    genes=[]
    for item in RGR_columns:
        tmpRGRs.append(pheno_call_df[item].to_list())
        tmpSexes.append([item.replace('_RGR','')]*len(pheno_call_df[item].to_list()))
        tmpCross.append(cross)
        genes.append(pheno_call_df.index.to_list())
    ## flatten the list
    RGRs= [item for sublist in tmpRGRs for item in sublist]
    Sexes=[item for sublist in tmpSexes for item in sublist]
    Cross=[item for sublist in tmpCross for item in sublist]
    allGenes=[item for sublist in genes for item in sublist]
    df_res=pd.DataFrame()
    df_res['Genes']=allGenes
    df_res['Relative fertility']=RGRs
    df_res['Sex']=Sexes
    df_res['published']=Cross
    sex_rename={'GCKO2':'Female','g145480':'Male','Cl15cy':'WT'}
    df_res['Sex']=df_res['Sex'].replace(sex_rename)
    df_res['shape']="o"
    sns.set(font="Arial",font_scale =2)
    sns.set_style("white")
    colors=['#FFFFFF', '#FFFFFF', '#FFFFFF']
    colors1=['#AF58BA', '#009ADE', '#FFC61E','#636363']
    ### shape wise

    plt.figure()
    tmp=df_res[df_res['published']=='NA']
    sns.stripplot(x="Sex", y="Relative fertility", data=tmp,color='#636363',marker="o")

    tmp=df_res[df_res['published']=='F']
    sns.stripplot(x="Sex", y="Relative fertility", data=tmp,color='#636363',marker="^")

    tmp=df_res[df_res['published']=='M']
    sns.stripplot(x="Sex", y="Relative fertility", data=tmp,color='#636363',marker="s")

    tmp=df_res[df_res['published']=='FM']
    sns.stripplot(x="Sex", y="Relative fertility", data=tmp,color='#636363',marker="h")

    sns.violinplot(x="Sex", y='Relative fertility', data=df_res,inner='box',
                         order=['Female','Male','WT'],palette=colors) ## 'quartile' 'box'
    plt.legend().remove()
    plt.tick_params(axis='x', rotation=0)
    # plt.xlabel('Gametocyte', fontsize=20,labelpad=10)
    plt.xlabel('')
    plt.ylabel('')
    plt.ylim((-18,8))

    plt.savefig(out_folder+'/pilot_violin_shape.pdf',
                format='pdf',dpi=300)
    plt.close()



    # sns.set_style("white")
    ## combine data
    # pheno_call_df[pheno_call_df.index.isin(pub_cross['Gene_ID_new'].to_list())]
    s=5
    plt.figure()
    tmp=df_res[df_res['published']=='NA']
    sns.stripplot(x="Sex", y="Relative fertility", color='#636363', data=tmp,size=s)
    tmp=df_res[df_res['published']=='F']
    sns.stripplot(x="Sex", y="Relative fertility", color='#AF58BA', data=tmp,size=s)
    tmp=df_res[df_res['published']=='M']
    sns.stripplot(x="Sex", y="Relative fertility", color='#009ADE', data=tmp,size=s)
    tmp=df_res[df_res['published']=='FM']
    sns.stripplot(x="Sex", y="Relative fertility", color='#FFC61E', data=tmp,size=s)
    # sns.stripplot(x="Sex", y='Fertility rate', data=df_res,color="black", edgecolor="gray",hue="published")

    sns.violinplot(x="Sex", y='Relative fertility', data=df_res,inner='box',
                         order=['Female','Male','WT'],palette=colors) ## 'quartile' 'box'
    #sns.stripplot(x="Sex", y="Relative fertility", hue="published", data=tmp,hue_order=['F','M','FM','NA'],palette=colors1)
    plt.legend().remove()
    plt.tick_params(axis='x', rotation=0)
    # plt.xlabel('Gametocyte', fontsize=20,labelpad=10)
    plt.xlabel('')
    plt.ylabel('')
    plt.ylim((-18,8))

    plt.savefig(out_folder+'/pilot_violin.pdf',
                format='pdf',dpi=300)
    plt.close()
