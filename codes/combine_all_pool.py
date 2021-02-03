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
## end of databse information
colorsIdx = {'No Power': '#cccccc', 'Not Reduced': '#b2e2e2','Reduced':'#d8b365', 'Increased':'#238b45'}

from fertility_pool5_2 import stepwiseAnalysis  as pool5_2
from fertility_pool7_2 import stepwiseAnalysis  as pool7_2


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
    backgrounds=['GCKO2','g145480']

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
    manifests_df=pd.read_csv(data_folder+"/data_pilot/manifest_pilot_small.txt",sep='\t')

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
    filtered_count_df.to_csv(out_folder+"/filterd_count_matrix_pilot1.txt",sep='\t')

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
    'mf':['mf1','mf2'],'sex':['GCKO2','g145480'],'geneConv':geneConv_new,
    'control_genes':['PBANKA_102460','PBANKA_100210','PBANKA_100220','PBANKA_050120']}

    pheno_call_df,input_df=relative_growth_rate_analysis_pilot(filtered_count_df,manfest_df,prev_to_new,db_df,plot_info)

    return pheno_call_df,input_df


def pool1():
    ''' We are going do analyis of pool1 data of Claire '''

    ### these are the input files
    manifests_df=pd.read_csv(data_folder+"/manifests_pool1.txt",sep='\t')
    count_df=pd.read_csv(data_folder+ "/result_240120_barcode_counts_table.txt",sep='\t')
    input_df=pd.read_csv(data_folder+'/input_vector.txt', sep='\t')


    #### end of the input section
    # final_count_df: read1 and read2 are added
    # final_count_df_two_read: reads are sperated
    # manfest_df: maifest_df
    final_count_df,final_count_df_des,final_count_df_two_read,manfest_df=preprocessing(manifests_df,count_df)

    ####  Dropouts and input check
    # input_df: these are genes which was used for pool phenotypes
    percent=0.9 ## this parameters is used to test whether count is too small for 90 % of input samples. Those will be deleted.
    filtered_count_df,filtered_df_read,filtered_count_df_des=filter_input_dropout(final_count_df,final_count_df_des,final_count_df_two_read,input_df,manfest_df,percent)

    ######  write filtered and unfiltered files
    # final_count_df_two_read.to_csv(out_folder+"/unfilterd_count_matrix_pool1.txt",sep='\t')
    # filtered_count_df_des.to_csv(out_folder+"/filterd_count_matrix_pool1.txt",sep='\t')

    ### we are going to perform relative abundance analysis
    ## prev_to_new this is the pickle information which is used when we change old to new ID
    ## db_df: this is the dataframe contains name and description

    # if you we do not want to plot then plot_info=None
    # plot_info={'pdf':out_folder+"/relative_abundance_of_pool1.pdf",'d':['d0','d13'],'mf':['mf1','mf2'],'sex':['GCKO2','g145480']}
    # # plot_info=None
    # relative_abundance_analysis(filtered_count_df,manfest_df,prev_to_new,db_df,plot_info)

    ## we will do diffrent kind of error analysis

    #error_analysis(filtered_count_df,manfest_df,prev_to_new,db_df)

    geneConv,old_to_new_ids,geneConv_new=getNewIdfromPrevID(filtered_count_df.index,prev_to_new,db_df)
    plot_info={'pool':'pool1','file':out_folder+'/pool1_repeat.xlsx','rel_file':out_folder+'/pool1_propagated_error_relative_abundance.pdf','d':['d0','d13'],
    'mf':['mf1','mf2'],'sex':['GCKO2','g145480'],'geneConv':geneConv_new,
    'control_genes':['PBANKA_102460' , 'PBANKA_050120' , 'PBANKA_010110' , 'PBANKA_142210']}
    pheno_call_df,input_df=relative_growth_rate_analysis(filtered_count_df,manfest_df,prev_to_new,db_df,plot_info)
    # import pdb; pdb.set_trace()

    return pheno_call_df,input_df

    ### we are going to start analysing
    # ### now we are going to combine two dataframe with categorical data sets
    # cmd_df=final_df.T.join(manifests_df)


def pool2():
    ''' We are going do analyis of pool1 data of Claire '''

    ### these are the input files
    manifests_df=pd.read_csv(data_folder+"/manifest_pool2.txt",sep='\t')
    count_df=pd.read_csv(data_folder+ "/barcode_counts_table_170620_pool2.txt",sep='\t')
    input_df=pd.read_csv(data_folder+'/input_pool2.txt', sep='\t')


    #### end of the input section
    # final_count_df: read1 and read2 are added
    # final_count_df_two_read: reads are sperated
    # manfest_df: maifest_df
    final_count_df,final_count_df_des,final_count_df_two_read,manfest_df=preprocessing(manifests_df,count_df)

    ####  Dropouts and input check
    # input_df: these are genes which was used for pool phenotypes
    percent=0.9 ## this parameters is used to test whether count is too small for 90 % of input samples. Those will be deleted.
    filtered_count_df,filtered_df_read,filtered_count_df_des=filter_input_dropout(final_count_df,final_count_df_des,final_count_df_two_read,input_df,manfest_df,percent)

    ######  write filtered and unfiltered files
    # final_count_df_two_read.to_csv(out_folder+"/unfilterd_count_matrix_pool2.txt",sep='\t')
    # filtered_count_df_des.to_csv(out_folder+"/filterd_count_matrix_pool2.txt",sep='\t')

    ### we are going to perform relative abundance analysis
    ## prev_to_new this is the pickle information which is used when we change old to new ID
    ## db_df: this is the dataframe contains name and description

    ## if you we do not want to plot then plot_info=None
    # plot_info={'pdf':out_folder+"/relative_abundance_of_pool2.pdf",'d':['d0','d13'],'mf':['mf1','mf2'],'sex':['GCKO2','g145480']}
    # # plot_info=None
    # relative_abundance_analysis(filtered_count_df,manfest_df,prev_to_new,db_df,plot_info)

    ## we will do diffrent kind of error analysis
    # import pdb;pdb.set_trace()
    # error_analysis(filtered_count_df,manfest_df,prev_to_new,db_df)
    geneConv,old_to_new_ids,geneConv_new=getNewIdfromPrevID(filtered_count_df.index,prev_to_new,db_df)
    plot_info={'pool':'pool2','file':out_folder+'/pool2_repeat.xlsx','rel_file':out_folder+'/pool2_propagated_error_relative_abundance.pdf','d':['d0','d13'],
    'mf':['mf1','mf2'],'sex':['GCKO2','g145480'],'geneConv':geneConv_new,
    'control_genes':['PBANKA_102460' , 'PBANKA_050120' , 'PBANKA_010110' , 'PBANKA_142210']}
    pheno_call_df,input_df=relative_growth_rate_analysis(filtered_count_df,manfest_df,prev_to_new,db_df,plot_info)




    return pheno_call_df,input_df


    ### we are going to start analysing
    # ### now we are going to combine two dataframe with categorical data sets
    # cmd_df=final_df.T.join(manifests_df)

def pool3():
    ''' We are going do analyis of pool1 data of Claire '''

    ### these are the input files
    manifests_df=pd.read_csv(data_folder+"/manifest_pool3.txt",sep='\t')
    count_df=pd.read_csv(data_folder+ "/barcode_counts_table_200720_pool3.txt",sep='\t')
    input_df=pd.read_csv(data_folder+'/input_pool3.txt', sep='\t')


    #### end of the input section
    # final_count_df: read1 and read2 are added
    # final_count_df_two_read: reads are sperated
    # manfest_df: maifest_df
    final_count_df,final_count_df_des,final_count_df_two_read,manfest_df=preprocessing(manifests_df,count_df)

    ####  Dropouts and input check
    # input_df: these are genes which was used for pool phenotypes
    percent=0.9 ## this parameters is used to test whether count is too small for 90 % of input samples. Those will be deleted.
    filtered_count_df,filtered_df_read,filtered_count_df_des=filter_input_dropout(final_count_df,final_count_df_des,final_count_df_two_read,input_df,manfest_df,percent)

    ######  write filtered and unfiltered files
    # final_count_df_two_read.to_csv(out_folder+"/unfilterd_count_matrix_pool2.txt",sep='\t')
    # filtered_count_df_des.to_csv(out_folder+"/filterd_count_matrix_pool3.txt",sep='\t')

    ### we are going to perform relative abundance analysis
    ## prev_to_new this is the pickle information which is used when we change old to new ID
    ## db_df: this is the dataframe contains name and description

    ## if you we do not want to plot then plot_info=None
    # plot_info={'pdf':out_folder+"/relative_abundance_of_pool3.pdf",'d':['d0','d13'],'mf':['mf1','mf2'],'sex':['GCKO2','g145480']}
    # # plot_info=None
    # relative_abundance_analysis(filtered_count_df,manfest_df,prev_to_new,db_df,plot_info)

    ## we will do diffrent kind of error analysis

    # error_analysis(filtered_count_df,manfest_df,prev_to_new,db_df)
    geneConv,old_to_new_ids,geneConv_new=getNewIdfromPrevID(filtered_count_df.index,prev_to_new,db_df)
    geneConv,old_to_new_ids,geneConv_new=getNewIdfromPrevID(filtered_count_df.index,prev_to_new,db_df)
    plot_info={'pool':'pool3','file':out_folder+'/pool3_repeat.xlsx','rel_file':out_folder+'/pool3_propagated_error_relative_abundance.pdf','d':['d0','d13'],
    'mf':['mf1','mf2'],'sex':['GCKO2','g145480'],'geneConv':geneConv_new,
    'control_genes':['PBANKA_0804300' , 'PBANKA_1427100' , 'PBANKA_010110' , 'PBANKA_1143400']}
    pheno_call_df,input_df=relative_growth_rate_analysis(filtered_count_df,manfest_df,prev_to_new,db_df,plot_info)

    return pheno_call_df,input_df




    ### we are going to start analysing
    # ### now we are going to combine two dataframe with categorical data sets
    # cmd_df=final_df.T.join(manifests_df)

def pool4():
    ''' We are going do analyis of pool1 data of Claire '''

    ### these are the input files
    manifests_df=pd.read_csv(data_folder+"/manifest_pool4.txt",sep='\t')
    count_df=pd.read_csv(data_folder+ "/barcode_counts_table_170620_pool4.txt",sep='\t')
    input_df=pd.read_csv(data_folder+'/input_pool4.txt', sep='\t')


    #### end of the input section
    # final_count_df: read1 and read2 are added
    # final_count_df_two_read: reads are sperated
    # manfest_df: maifest_df
    final_count_df,final_count_df_des,final_count_df_two_read,manfest_df=preprocessing(manifests_df,count_df)

    ####  Dropouts and input check
    # input_df: these are genes which was used for pool phenotypes
    percent=0.9 ## this parameters is used to test whether count is too small for 90 % of input samples. Those will be deleted.
    filtered_count_df,filtered_df_read,filtered_count_df_des=filter_input_dropout(final_count_df,final_count_df_des,final_count_df_two_read,input_df,manfest_df,percent)

    ######  write filtered and unfiltered files
    # final_count_df_two_read.to_csv(out_folder+"/unfilterd_count_matrix_pool2.txt",sep='\t')
    # filtered_count_df_des.to_csv(out_folder+"/filterd_count_matrix_pool4.txt",sep='\t')
    #
    # ### we are going to perform relative abundance analysis
    # ## prev_to_new this is the pickle information which is used when we change old to new ID
    # ## db_df: this is the dataframe contains name and description
    #
    # ## if you we do not want to plot then plot_info=None
    # plot_info={'pdf':out_folder+"/relative_abundance_of_pool4.pdf",'d':['d0','d13'],'mf':['mf1','mf2'],'sex':['GCKO2','g145480']}
    # # plot_info=None
    # relative_abundance_analysis(filtered_count_df,manfest_df,prev_to_new,db_df,plot_info)

    ## we will do diffrent kind of error analysis

    # error_analysis(filtered_count_df,manfest_df,prev_to_new,db_df)

    geneConv,old_to_new_ids,geneConv_new=getNewIdfromPrevID(filtered_count_df.index,prev_to_new,db_df)
    plot_info={'pool':'pool4','file':out_folder+'/pool4_repeat.xlsx','rel_file':out_folder+'/pool4_propagated_error_relative_abundance.pdf','d':['d0','d13'],
    'mf':['mf1','mf2'],'sex':['GCKO2','g145480'],'geneConv':geneConv_new,
    'control_genes':['PBANKA_0511000' , 'PBANKA_050120' , 'PBANKA_1425200' , 'PBANKA_0404600']}
    pheno_call_df,input_df=relative_growth_rate_analysis(filtered_count_df,manfest_df,prev_to_new,db_df,plot_info)
    return pheno_call_df,input_df


    ### we are going to start analysing
    # ### now we are going to combine two dataframe with categorical data sets
    # cmd_df=final_df.T.join(manifests_df)
def pool6():
    ''' We are going do analyis of pool1 data of Claire '''

    ### these are the input files
    manifests_df=pd.read_csv(data_folder+"/manifest_pool6.txt",sep='\t')
    count_df=pd.read_csv(data_folder+ "/barcode_counts_table_280720_pool6.txt",sep='\t')
    input_df=pd.read_csv(data_folder+'/input_pool6.txt', sep='\t')


    #### end of the input section
    # final_count_df: read1 and read2 are added
    # final_count_df_two_read: reads are sperated
    # manfest_df: maifest_df
    final_count_df,final_count_df_des,final_count_df_two_read,manfest_df=preprocessing(manifests_df,count_df)

    ####  Dropouts and input check
    # input_df: these are genes which was used for pool phenotypes
    percent=0.9 ## this parameters is used to test whether count is too small for 90 % of input samples. Those will be deleted.
    filtered_count_df,filtered_df_read,filtered_count_df_des=filter_input_dropout(final_count_df,final_count_df_des,final_count_df_two_read,input_df,manfest_df,percent)

    ######  write filtered and unfiltered files
    # final_count_df_two_read.to_csv(out_folder+"/unfilterd_count_matrix_pool2.txt",sep='\t')
    # filtered_count_df_des.to_csv(out_folder+"/filterd_count_matrix_pool6.txt",sep='\t')

    ### we are going to perform relative abundance analysis
    ## prev_to_new this is the pickle information which is used when we change old to new ID
    ## db_df: this is the dataframe contains name and description

    # if you we do not want to plot then plot_info=None
    # plot_info={'pdf':out_folder+"/relative_abundance_of_pool6.pdf",'d':['d0','d13'],'mf':['mf1','mf2'],'sex':['GCKO2','g145480']}
    # # plot_info=None
    # relative_abundance_analysis(filtered_count_df,manfest_df,prev_to_new,db_df,plot_info)
    #
    # # ## we will do diffrent kind of error analysis

    # # drop P17509_1003 from manifest and final_count_matrix
    filtered_count_df=filtered_count_df.drop(columns='P17509_1003').copy()
    manfest_df=manfest_df.drop('P17509_1003').copy()
    #
    #
    # # error_analysis(filtered_count_df,manfest_df,prev_to_new,db_df)
    geneConv,old_to_new_ids,geneConv_new=getNewIdfromPrevID(filtered_count_df.index,prev_to_new,db_df)
    plot_info={'pool':'pool6','file':out_folder+'/pool6_repeat.xlsx','rel_file':out_folder+'/pool6_propagated_error_relative_abundance.pdf','d':['d0','d13'],
    'mf':['mf1','mf2'],'sex':['GCKO2','g145480'],'geneConv':geneConv_new,
    'control_genes':['PBANKA_0812800' , 'PBANKA_1031700' , 'PBANKA_010110' , 'PBANKA_1101400']}
    pheno_call_df,input_df=relative_growth_rate_analysis(filtered_count_df,manfest_df,prev_to_new,db_df,plot_info)

    return pheno_call_df,input_df


##### pool 5 to 7

def pool5():
    ''' We are going do analyis of pool1 data of Claire '''

    ### these are the input files
    manifests_df=pd.read_csv(data_folder+"/manifest_pool5_1.txt",sep='\t')
    count_df=pd.read_csv(data_folder+ "/barcode_counts_table_20920_pool5.txt",sep='\t')
    input_df=pd.read_csv(data_folder+'/input_pool5.txt', sep='\t')


    #### end of the input section
    # final_count_df: read1 and read2 are added
    # final_count_df_two_read: reads are sperated
    # manfest_df: maifest_df
    final_count_df,final_count_df_des,final_count_df_two_read,manfest_df=preprocessing(manifests_df,count_df)

    ####  Dropouts and input check
    # input_df: these are genes which was used for pool phenotypes
    percent=0.9 ## this parameters is used to test whether count is too small for 90 % of input samples. Those will be deleted.
    filtered_count_df,filtered_df_read,filtered_count_df_des=filter_input_dropout(final_count_df,final_count_df_des,final_count_df_two_read,input_df,manfest_df,percent)

    ######  write filtered and unfiltered files
    # final_count_df_two_read.to_csv(out_folder+"/unfilterd_count_matrix_pool2.txt",sep='\t')
    filtered_count_df_des.to_csv(out_folder+"/filterd_count_matrix_pool5_1.txt",sep='\t')

    ### we are going to perform relative abundance analysis
    ## prev_to_new this is the pickle information which is used when we change old to new ID
    ## db_df: this is the dataframe contains name and description

    # if you we do not want to plot then plot_info=None
    # plot_info={'pdf':out_folder+"/relative_abundance_of_pool6.pdf",'d':['d0','d13'],'mf':['mf1','mf2'],'sex':['GCKO2','g145480']}
    # # plot_info=None
    # relative_abundance_analysis(filtered_count_df,manfest_df,prev_to_new,db_df,plot_info)
    #
    # # ## we will do diffrent kind of error analysis

    # drop P17509_1003 from manifest and final_count_matrix
    # filtered_count_df=filtered_count_df.drop(columns='P17509_1003').copy()
    # manfest_df=manfest_df.drop('P17509_1003').copy()
    #
    #
    geneConv,old_to_new_ids,geneConv_new=getNewIdfromPrevID(filtered_count_df.index,prev_to_new,db_df)
    geneConv,old_to_new_ids,geneConv_new=getNewIdfromPrevID(filtered_count_df.index,prev_to_new,db_df)
    plot_info={'pool':'pool5_1','file':out_folder+'/pool5_1_repeat.xlsx','rel_file':out_folder+'/pool5_1_propagated_error_relative_abundance.pdf','d':['d0','d13'],
    'mf':['mf1','mf2'],'sex':['GCKO2','g145480'],'geneConv':geneConv_new,
    'control_genes':['PBANKA_1404100' , 'PBANKA_050120' , 'PBANKA_010110' , 'PBANKA_142210']}
    pheno_call_df,input_df=relative_growth_rate_analysis(filtered_count_df,manfest_df,prev_to_new,db_df,plot_info)

    return pheno_call_df,input_df


### pool7


def pool7():
    ''' We are going do analyis of pool1 data of Claire '''

    ### these are the input files
    manifests_df=pd.read_csv(data_folder+"/manifest_pool7_1.txt",sep='\t')
    count_df=pd.read_csv(data_folder+ "/barcode_counts_table_20920_pool7.txt",sep='\t')
    input_df=pd.read_csv(data_folder+'/input_pool7.txt', sep='\t')


    #### end of the input section
    # final_count_df: read1 and read2 are added
    # final_count_df_two_read: reads are sperated
    # manfest_df: maifest_df
    final_count_df,final_count_df_des,final_count_df_two_read,manfest_df=preprocessing(manifests_df,count_df)

    ####  Dropouts and input check
    # input_df: these are genes which was used for pool phenotypes
    percent=0.9 ## this parameters is used to test whether count is too small for 90 % of input samples. Those will be deleted.
    filtered_count_df,filtered_df_read,filtered_count_df_des=filter_input_dropout(final_count_df,final_count_df_des,final_count_df_two_read,input_df,manfest_df,percent)

    ######  write filtered and unfiltered files
    # final_count_df_two_read.to_csv(out_folder+"/unfilterd_count_matrix_pool2.txt",sep='\t')
    filtered_count_df_des.to_csv(out_folder+"/filterd_count_matrix_pool7_1.txt",sep='\t')

    ### we are going to perform relative abundance analysis
    ## prev_to_new this is the pickle information which is used when we change old to new ID
    ## db_df: this is the dataframe contains name and description

    # if you we do not want to plot then plot_info=None
    # plot_info={'pdf':out_folder+"/relative_abundance_of_pool6.pdf",'d':['d0','d13'],'mf':['mf1','mf2'],'sex':['GCKO2','g145480']}
    # # plot_info=None
    # relative_abundance_analysis(filtered_count_df,manfest_df,prev_to_new,db_df,plot_info)
    #
    # # ## we will do diffrent kind of error analysis

    # drop P17509_1003 from manifest and final_count_matrix
    # filtered_count_df=filtered_count_df.drop(columns='P17509_1003').copy()
    # manfest_df=manfest_df.drop('P17509_1003').copy()
    #
    #
    geneConv,old_to_new_ids,geneConv_new=getNewIdfromPrevID(filtered_count_df.index,prev_to_new,db_df)
    geneConv,old_to_new_ids,geneConv_new=getNewIdfromPrevID(filtered_count_df.index,prev_to_new,db_df)
    plot_info={'pool':'pool7_1','file':out_folder+'/pool7_1_repeat.xlsx','rel_file':out_folder+'/pool7_1_propagated_error_relative_abundance.pdf','d':['d0','d13'],
    'mf':['mf1','mf2'],'sex':['GCKO2','g145480'],'geneConv':geneConv_new,
    'control_genes':['PBANKA_102460' , 'PBANKA_050120' , 'PBANKA_010110' , 'PBANKA_142210']}
    pheno_call_df,input_df=relative_growth_rate_analysis(filtered_count_df,manfest_df,prev_to_new,db_df,plot_info)
    return pheno_call_df,input_df




def calPvalQval(pheno_call_df,cutoff=0):
    ''' pheno_call_df'''

    ## select no_power and Both feed noisy
    sex=['GCKO2','g145480']
    for b in sex:
        noisy_df=pheno_call_df[(pheno_call_df[b+'_pheno']=='No Power') & (pheno_call_df[b+'_feed']=='no data')].copy()
        ## we will assign p-value 1 for noisy data
        actual_df=pheno_call_df[~((pheno_call_df[b+'_pheno']=='No Power') & (pheno_call_df[b+'_feed']=='no data'))].copy()
        genes=actual_df.index
        m=actual_df[b+'_RGR'].values
        s=actual_df[b+'_sd'].values
        z=(cutoff-m)/s
        z1=abs(z)
        pvalue=  (1 - st.norm.cdf(z.tolist()))
        pvalue2=  (1 - st.norm.cdf(z1.tolist())) *2
        fdr=mtest.multipletests(pvalue, alpha=0.05, method='fdr_bh')
        fdr2=mtest.multipletests(pvalue2, alpha=0.05, method='fdr_bh')
        pheno_call_df[b+'_pvalue2']=np.nan
        pheno_call_df[b+'_pvalue']=np.nan
        pheno_call_df[b+'_zvalue']=np.nan
        pheno_call_df[b+'_fdr']=np.nan
        pheno_call_df[b+'_fdr2']=np.nan
        pheno_call_df.loc[genes,b+'_pvalue']=pvalue;
        pheno_call_df.loc[genes,b+'_fdr']=fdr[1];
        pheno_call_df.loc[genes,b+'_zvalue']=z;
        pheno_call_df.loc[genes,b+'_pvalue2']=pvalue2;
        pheno_call_df.loc[genes,b+'_fdr2']=fdr2[1];

    return pheno_call_df

def combinepools():
    ''' Combining all pools relative abundance and relative growth rate  '''
    pheno_call_pilot,input_df_pilot=pilot()
    pheno_call_pilot.index.name='pbanka_id'
    pheno_call_pilot['pool']='pilot'
    print ('pilot finished\n\n')

    pheno_call_pool1,input_df_pool1=pool1()
    pheno_call_pool1['pool']='pool1'
    print ('pool1 finished\n\n')

    pheno_call_pool2,input_df_pool2=pool2()
    pheno_call_pool2['pool']='pool2'
    print ('pool2 finished\n\n')
    pheno_call_pool3,input_df_pool3=pool3()
    pheno_call_pool3['pool']='pool3'
    print ('pool3 finished\n\n')
    pheno_call_pool4,input_df_pool4=pool4()
    pheno_call_pool4['pool']='pool4'
    print ('pool4 finished\n\n')
    pheno_call_pool6,input_df_pool6=pool6()
    pheno_call_pool6['pool']='pool6'
    print ('pool6 finished\n\n')
    ###
    pheno_call_pool5,input_df_pool5=pool5()
    pheno_call_pool5['pool']='pool5'
    print ('pool5 finished\n\n')

    pheno_call_pool7,input_df_pool7=pool7()
    pheno_call_pool7['pool']='pool7'
    print ('pool7 finished\n\n')
    pheno_call_pool5_2,input_df_pool5_2=pool5_2()
    pheno_call_pool5_2['pool']='pool5_2'
    print ('pool5_2 finished\n\n')
    pheno_call_pool7_2,input_df_pool7_2=pool7_2()
    pheno_call_pool7_2['pool']='pool7_2'
    print ('pool7_2 finished\n\n')



    pheno_all=pd.concat([pheno_call_pilot,pheno_call_pool1, pheno_call_pool2,pheno_call_pool3,pheno_call_pool4,pheno_call_pool6,pheno_call_pool5,pheno_call_pool7,pheno_call_pool5_2,pheno_call_pool7_2])
    #pheno_all=pd.concat([pheno_call_pool1, pheno_call_pool2,pheno_call_pool3,pheno_call_pool4,pheno_call_pool6])
    # pheno_all=pd.concat([pheno_call_pilot,pheno_call_pool6])
    # pheno_all=pheno_all.replace({'GCKO2_pheno': {'NA': 'No power', 'NE': 'Fertility','E':'Reduced Fertility'}})
    # pheno_all=pheno_all.replace({'g145480_pheno': {'NA': 'No power', 'NE': 'Fertility','E':'Reduced Fertility'}})
    # ### plot female
    ####### Combine one gene for all pheno types
    ### test input cutoff
    #test_for_input_cutoff(input_df_pilot,input_df_pool1, input_df_pool2,input_df_pool3,input_df_pool4,input_df_pool6)

    pheno_all2=apply_weighted_mean(pheno_all)
    pheno_all2=calPvalQval(pheno_all2,cutoff=0)

    # import pdb; pdb.set_trace()

    ### set color

    # colorsIdx = {'No power': '#f5f5f5', 'Not': '#5ab4ac','Reduced Fertility':'#d8b365'}
    ### filter no data entries
    write_df=pheno_all2.copy()
    feeds=['g145480_feed','GCKO2_feed']
    sex_df=[]
    for f in feeds:
        xx=pheno_all2[f].str.split(',').apply(lambda x: set(x))
        bool_ind=[False if item==set(['no data']) else True  for item in xx ]
        tmp=pheno_all2[bool_ind].copy()
        sex_df.append(tmp)
    ### plot female
    xx=sex_df[1].sort_values(by=['GCKO2_RGR']).copy()
    write_df=write_df.sort_values(by=['GCKO2_RGR']).copy()
    write_df=write_df.reset_index()
    write_df['Gene IDs']=write_df['pbanka_id']+'_'+write_df['pool']
    xx=xx.reset_index()

    xx['Gene IDs']=xx['pbanka_id']+'_'+xx['pool']
    fig_GCKO2 = px.scatter(xx,x=xx.index,y='GCKO2_RGR',color="GCKO2_pheno",color_discrete_map=colorsIdx,error_y='GCKO2_sd',hover_name='Gene IDs',render_mode="svg")
    #fig_GCKO2 = px.scatter(xx,x=xx.index,y='GCKO2_RGR',color="GCKO2_pheno",color_discrete_map=colorsIdx,error_y='GCKO2_sd')
    fig_GCKO2.update_layout(
    autosize=False,
    width=2500,
    height=1000,
    margin=dict(
        l=50,
        r=50,
        b=100,
        t=100,
        pad=4
    ),
    paper_bgcolor="LightSteelBlue",
    )



    fig_GCKO2.show()
    ## plot male
    xx=sex_df[0].sort_values(by=['g145480_RGR']).copy()
    xx=xx.reset_index()
    xx['Gene IDs']=xx['pbanka_id']+'_'+xx['pool']
    fig_g145480= px.scatter(xx,x=xx.index,y='g145480_RGR',color="g145480_pheno",color_discrete_map=colorsIdx,error_y='g145480_sd',hover_name='Gene IDs',render_mode="svg")
    fig_g145480.update_layout(
    autosize=False,
    width=2500,
    height=1000,
    margin=dict(
        l=50,
        r=50,
        b=100,
        t=100,
        pad=4
    ),
    paper_bgcolor="LightSteelBlue",
    )
    fig_g145480.show()
    # fig, ax = plt.subplots(figsize=(16,12))
    # plt.scatter(xx.index.to_list(), xx['g145480_RGR'].values, marker='o',c=xx['g145480_pheno'].map(colorsIdx).values);
    # plt.errorbar(xx.index.to_list(), xx['g145480_RGR'].values, yerr=xx['g145480_sd'].values,ecolor=xx['g145480_pheno'].map(colorsIdx).values);
    # # legend1 = ax.legend(*scatter.legend_elements(),
    # #                 loc="lower left", title="Classes")
    # plt.show()
    # import pdb; pdb.set_trace()
    ### add gene anotation file
    #xx=sex_df[0].copy()
    xx=write_df.copy()
    xx['gene_product']='NA'
    xx['gene_symbol']='NA'
    for idx,pbanka_id in enumerate(xx['pbanka_id']):
        tmp=db_df[db_df['Gene ID']== pbanka_id]
        if not tmp.empty:
            xx.loc[idx,'gene_product']=tmp['Product Description'].to_list()[0]
            xx.loc[idx,'gene_symbol']=tmp['Gene Name or Symbol'].to_list()[0]

    ### write in order
    order_col=['pbanka_id','Gene IDs',	'gene_product',	'gene_symbol','pool']
    sex=['GCKO2','g145480']
    for b in sex:
        for item in xx.columns[xx.columns.str.contains(b)].to_list():
            order_col.append(item)
    yy=xx[order_col]
    yy.to_csv(out_folder+'/Phenotype_call_final.txt',sep='\t')



def test_for_input_cutoff(input_df_pilot,input_df_pool1, input_df_pool2,input_df_pool3,input_df_pool4,input_df_pool6):

    input_df_all_mf1=pd.concat([input_df_pilot[0],input_df_pool1[0], input_df_pool2[0],input_df_pool3[0],input_df_pool4[0],input_df_pool6[0]])
    input_df_all_mf2=pd.concat([input_df_pilot[1],input_df_pool1[1], input_df_pool2[1],input_df_pool3[1],input_df_pool4[1],input_df_pool6[1]])

    ### distributions
    x=input_df_all_mf1['GCKO2_d0_NA_NA'].values
    x_grid=np.linspace(x.min(), x.max(), 500)
    trace_d0_mf1_dist_female = go.Histogram( x=x,
    name='# of genes='+str(len(x)),
    marker_color='#bdbdbd',
    opacity=0.75,
    histnorm='probability')

    ### distributions
    x=input_df_all_mf1['g145480_d0_NA_NA'].values
    x_grid=np.linspace(x.min(), x.max(), 500)
    trace_d0_mf1_dist_male = go.Histogram( x=x,
    name='# of genes='+str(len(x)),
    marker_color='#bdbdbd',
    opacity=0.75,
    histnorm='probability')

    ### distributions
    x=input_df_all_mf2['GCKO2_d0_NA_NA'].values
    x_grid=np.linspace(x.min(), x.max(), 500)
    trace_d0_mf2_dist_female = go.Histogram( x=x,
    name='# of genes='+str(len(x)),
    marker_color='#bdbdbd',
    opacity=0.75,
    histnorm='probability')

    ### distributions
    x=input_df_all_mf2['g145480_d0_NA_NA'].values
    x_grid=np.linspace(x.min(), x.max(), 500)
    trace_d0_mf2_dist_male = go.Histogram( x=x,
    name='# of genes='+str(len(x)),
    marker_color='#bdbdbd',
    opacity=0.75,
    histnorm='probability')

    fig = make_subplots(rows=2, cols=2,subplot_titles=("Distribution GCKO2(day0) mf1", "Distribution 145480(day0) mf1","Distribution GCKO2(day0) mf2", "Distribution 145480(day0) mf2"))

    fig.append_trace(trace_d0_mf1_dist_female,row=1, col=1)
    fig.append_trace(trace_d0_mf1_dist_male,row=1, col=2)
    fig.append_trace(trace_d0_mf2_dist_female,row=2, col=1)
    fig.append_trace(trace_d0_mf2_dist_male,row=2, col=2)

    # Update xaxis properties
    # fig.update_yaxes(title_text="relative error",row=1, col=1,range=[0,0.8])
    # for pool6
    fig.update_yaxes(title_text="Distribution",row=1, col=1,range=[0,0.25])
    fig.update_yaxes(title_text="Distribution",row=1, col=2,range=[0,0.25])
    fig.update_yaxes(title_text="Distribution",row=2, col=1,range=[0,0.25])
    fig.update_yaxes(title_text="Distribution",row=2, col=2,range=[0,0.25])

    fig.update_xaxes(title_text="log2 (relative abundance or count)", row=1, col=1,range=[-18,0])
    fig.update_xaxes(title_text="log2 (relative abundance or count)", row=1, col=2,range=[-18,0])
    fig.update_xaxes(title_text="log2 (relative abundance or count)", row=2, col=1,range=[-18,0])
    fig.update_xaxes(title_text="log2 (relative abundance or count)", row=2, col=2,range=[-18,0])

    fig.show()

def pilot_call():
    pheno_call_pilot=pilot()
    pheno_call_pilot.index.name='pbanka_id'
    pheno_call_pilot['pool']='pilot'
    print ('pilot finished\n\n')
    pheno_all=pheno_call_pilot
    # pheno_all=pheno_all.replace({'GCKO2_pheno': {'NA': 'No power', 'NE': 'Fertility','E':'Reduced Fertility'}})
    # pheno_all=pheno_all.replace({'g145480_pheno': {'NA': 'No power', 'NE': 'Fertility','E':'Reduced Fertility'}})
    ### plot female



    ### set color

    # colorsIdx = {'No power': '#f5f5f5', 'Fertility': '#5ab4ac','Reduced Fertility':'#d8b365'}

    ### plot female
    xx=pheno_all.sort_values(by=['GCKO2_RGR']).copy()
    xx=xx.reset_index()
    xx['Gene IDs']=xx['pbanka_id']+'_'+xx['pool']
    #fig_GCKO2 = px.scatter(xx,x=xx.index,y='GCKO2_RGR',color="GCKO2_pheno",color_discrete_map=colorsIdx,error_y='GCKO2_sd',hover_name='Gene IDs')
    fig_GCKO2 = px.scatter(xx,x=xx.index,y='GCKO2_RGR',color="GCKO2_pheno",color_discrete_map=colorsIdx,error_y='GCKO2_sd',hover_name='Gene IDs')
    fig_GCKO2.show()
    ## plot male
    xx=pheno_all.sort_values(by=['g145480_RGR']).copy()
    xx=xx.reset_index()
    xx['Gene IDs']=xx['pbanka_id']+'_'+xx['pool']
    fig_g145480= px.scatter(xx,x=xx.index,y='g145480_RGR',color="g145480_pheno",color_discrete_map=colorsIdx,error_y='g145480_sd',hover_name='Gene IDs')
    fig_g145480.show()

def applyMeanVar_nodata(b,combine_pheno_call,tmp_not_applied,cmb_fitness):
    rgr_temp=pd.DataFrame(index=tmp_not_applied.index,columns=[b])
    var_temp=pd.DataFrame(index=tmp_not_applied.index,columns=[b])

    rgr_temp.loc[:,b]=tmp_not_applied.loc[:,b+'_RGR'].copy()

    var_temp.loc[:,b]=tmp_not_applied.loc[:,b+'_var'].copy()

    cmb_fitness[b]=gaussianMeanAndVariance(rgr_temp.T.copy(),var_temp.T.copy())

    ## find call
    m=cmb_fitness[b][0][b]
    s=cmb_fitness[b][1][b]
    v=cmb_fitness[b][2][b]

    ##
    diff_max=m+2*s
    diff_min=m-2*s
    pheno_type='No Power'

    return m,s,v,diff_max,diff_min,pheno_type

def applyMeanVar_final(b,combine_pheno_call,tmp_not_applied,cmb_fitness):
    rgr_temp=pd.DataFrame(index=tmp_not_applied.index,columns=[b])
    var_temp=pd.DataFrame(index=tmp_not_applied.index,columns=[b])

    rgr_temp.loc[:,b]=tmp_not_applied.loc[:,b+'_RGR'].copy()

    var_temp.loc[:,b]=tmp_not_applied.loc[:,b+'_var'].copy()

    cmb_fitness[b]=gaussianMeanAndVariance(rgr_temp.T.copy(),var_temp.T.copy())

    ## find call
    m=cmb_fitness[b][0][b]
    s=cmb_fitness[b][1][b]
    v=cmb_fitness[b][2][b]

    ##
    if b=='g145480':
        lower_cut=-1#np.log2(0.45)
    else:
        lower_cut=-1.6#np.log2(0.45)
    upper_cut=np.log2(2)
    diff_max=m+2*s
    diff_min=m-2*s

    if (diff_max<lower_cut) and  (diff_min<lower_cut):
        pheno_type='Reduced'
    # elif (diff_max< upper_cut) and  (diff_min >lower_cut) :
    #     pheno_type='Not Reduced'
    elif (diff_min >lower_cut) :
        pheno_type='Not Reduced'
    # elif (diff_max>upper_cut) and  (diff_min>upper_cut) :
    #     pheno_type='Increased'
    else:
        pheno_type='No Power'

    return m,s,v,diff_max,diff_min,pheno_type


def add_experiment(uq_gene,tmp,b,combine_pheno_call,cmb_fitness):

    combine_pheno_call_copy=combine_pheno_call.copy()

    tmp_not_applied=tmp[tmp[b+'_relative_filter']=='Not applied']
    if len(tmp_not_applied.index)>1:
        m,s,v,diff_max,diff_min,pheno_type=applyMeanVar_final(b,combine_pheno_call_copy,tmp_not_applied,cmb_fitness)
        combine_pheno_call_copy.loc[uq_gene,b+'_RGR']=m
        combine_pheno_call_copy.loc[uq_gene,b+'_sd']=s
        combine_pheno_call_copy.loc[uq_gene,b+'_var']=v
        combine_pheno_call_copy.loc[uq_gene,b+'_diff_max']=diff_max
        combine_pheno_call_copy.loc[uq_gene,b+'_diff_min']=diff_min
        combine_pheno_call_copy.loc[uq_gene,b+'_pheno']=pheno_type
        combine_pheno_call_copy.loc[uq_gene,b+'_relative_filter']=tmp[b+'_relative_filter'].str.cat(sep=',')
        combine_pheno_call_copy.loc[uq_gene,b+'_feed']=tmp[b+'_feed'].str.cat(sep=',')
        combine_pheno_call_copy.loc[uq_gene,'pool']=tmp['pool'].str.cat(sep=',')
    elif len(tmp_not_applied.index)==1:
        combine_pheno_call_copy.loc[uq_gene,b+'_RGR']=tmp_not_applied.loc[uq_gene,b+'_RGR'].copy()
        combine_pheno_call_copy.loc[uq_gene,b+'_sd']=tmp_not_applied.loc[uq_gene,b+'_sd'].copy()
        combine_pheno_call_copy.loc[uq_gene,b+'_var']=tmp_not_applied.loc[uq_gene,b+'_var'].copy()
        combine_pheno_call_copy.loc[uq_gene,b+'_diff_max']=tmp_not_applied.loc[uq_gene,b+'_diff_max'].copy()
        combine_pheno_call_copy.loc[uq_gene,b+'_diff_min']=tmp_not_applied.loc[uq_gene,b+'_diff_min'].copy()
        combine_pheno_call_copy.loc[uq_gene,b+'_pheno']=tmp_not_applied.loc[uq_gene,b+'_pheno']
        combine_pheno_call_copy.loc[uq_gene,b+'_relative_filter']=tmp[b+'_relative_filter'].str.cat(sep=',')
        combine_pheno_call_copy.loc[uq_gene,b+'_feed']=tmp[b+'_feed'].str.cat(sep=',')
        combine_pheno_call_copy.loc[uq_gene,'pool']=tmp['pool'].str.cat(sep=',')
    else:
        if len(tmp.index)>1:
            ## apply mean var final
            # if not (b=='g145480'):
            #     print(b,uq_gene)
            # else:
            #     import pdb; pdb.set_trace()
            m,s,v,diff_max,diff_min,pheno_type=applyMeanVar_nodata(b,combine_pheno_call_copy,tmp,cmb_fitness)
            combine_pheno_call_copy.loc[uq_gene,b+'_RGR']=m
            combine_pheno_call_copy.loc[uq_gene,b+'_sd']=s
            combine_pheno_call_copy.loc[uq_gene,b+'_var']=v
            combine_pheno_call_copy.loc[uq_gene,b+'_diff_max']=diff_max
            combine_pheno_call_copy.loc[uq_gene,b+'_diff_min']=diff_min
            combine_pheno_call_copy.loc[uq_gene,b+'_pheno']=pheno_type
            combine_pheno_call_copy.loc[uq_gene,b+'_relative_filter']=tmp[b+'_relative_filter'].str.cat(sep=',')
            combine_pheno_call_copy.loc[uq_gene,b+'_feed']=tmp[b+'_feed'].str.cat(sep=',')
            combine_pheno_call_copy.loc[uq_gene,'pool']=tmp['pool'].str.cat(sep=',')
        else:
            print('all pools have bad data for gene = %s'%uq_gene)
    return combine_pheno_call_copy

def apply_weighted_mean(pheno_all):
    ''' Combine repeated analysis'''
    unique_genes=pheno_all.index.unique()
    backgrounds=['GCKO2','g145480']
    combine_pheno_call=pd.DataFrame(index=unique_genes,columns=pheno_all.columns)
    combine_pheno_save=combine_pheno_call.copy()
    for uq_gene in unique_genes:
        tmp=pheno_all[pheno_all.index==uq_gene].copy()
        if len(tmp.index)>1:

            ## apply weighted sum and variance
            ## filter only not_applied for combining data
            # if uq_gene == 'PBANKA_1109600':
            #      import pdb; pdb.set_trace()
            cmb_fitness={}
            for b in backgrounds:
                    combine_pheno_call=add_experiment(uq_gene,tmp,b,combine_pheno_call,cmb_fitness)
                    combine_pheno_save=combine_pheno_call.copy()
                    ## still we can combine based on tmp
                ##
        else:
            combine_pheno_call.loc[uq_gene,:]=pheno_all.loc[uq_gene,:].copy()
            combine_pheno_save=combine_pheno_call.copy()
    return combine_pheno_save

if __name__ == '__main__':
    combinepools()
    #pilot_call()
