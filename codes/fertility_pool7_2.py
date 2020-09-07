import os
import sys

code=os.getcwd()
upLevel=code.replace('codes','') ####### we are going to upper level of code directory
sys.path.insert(0,upLevel+'/data')
sys.path.insert(1, upLevel+'/Figures')

# print(sys.path)

from special_fun_fertility  import *
from special_fun_repeat  import relative_growth_rate_analysis_repeat
# input files which will be needed for analysis of cliare screen data

data_folder=sys.path[0]

## output folder where we will write figures and output files
out_folder=sys.path[1]

# ID conversion: we used plasmoDB to convert ID for P. Berghai
prev_to_new=pickle.load(open(data_folder+'/prevTonew_PBANKA.pickle','rb'))
db_df=pd.read_csv(data_folder+'/PBANKA_id_conversion.txt', sep='\t')
db_df=db_df.fillna('NA')
## end of databse information
relative_growth_rate_analysis_repeat



def stepwiseAnalysis():
    ''' We are going do analyis of pool1 data of Claire '''

    ### these are the input files
    manifests_df=pd.read_csv(data_folder+"/manifest_pool7_2.txt",sep='\t')
    count_df=pd.read_csv(data_folder+ "/barcode_counts_table_20920_pool7.txt",sep='\t')
    input_df=pd.read_csv(data_folder+'/input_pool7_2.txt', sep='\t')


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
    filtered_count_df_des.to_csv(out_folder+"/filterd_count_matrix_pool7_2.txt",sep='\t')

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
    plot_info={'pool':'pool7_2','file':out_folder+'/pool7_2_repeat.xlsx','rel_file':out_folder+'/pool7_2_propagated_error_relative_abundance.pdf','d':['d0','d13'],
    'mf':['mf1','mf2'],'sex':['GCKO2'],'geneConv':geneConv_new,
    'control_genes':['PBANKA_102460' , 'PBANKA_050120' , 'PBANKA_010110' , 'PBANKA_142210']}
    pheno_call_df,input_df=relative_growth_rate_analysis_repeat(filtered_count_df,manfest_df,prev_to_new,db_df,plot_info)
    return pheno_call_df,input_df







    ### we are going to start analysing
    # ### now we are going to combine two dataframe with categorical data sets
    # cmd_df=final_df.T.join(manifests_df)




if __name__ == '__main__':
    stepwiseAnalysis()
