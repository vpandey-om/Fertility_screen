import os
import sys

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
db_df=pd.read_csv(data_folder+'/GenesByGeneModelChars_Summary.txt', sep='\t')
db_df=db_df.fillna('NA')
new_to_prev= dict((v,k) for k,v in prev_to_new.items())
## end of databse information



def filter_input_dropout_male_pool(df,df_s,df_read,input_df,manfest_df,percent=0.9,rel_cut=1e-5,day0='d0'):
    ''' We are going to find dropouts as well as how many we can map to inputs'''
    ########## input parameters
    # df: original count dataframe (read1+ read2)
    # df_read is dtaframe when reads are sperated
    # input_df: This is the DATAFRAME  which one used as a input for pools
    # manfest_df: manifests dtaframe
    ########

    ########## output parameters
    ## filtered_count_df
    ## filtered_df_read
    ########

    ## we will drop sum string columns and work on numerical values
    tmp=df.copy()
    tmp=tmp.drop(columns=['Gene','Barcodes','pbanka_id'])
    tmp=tmp.div(tmp.sum(axis=0))

    cols=[]

    for k,v in manfest_df.groupby('t').indices.items():
        if not (k =='NA'):
            for id in v:
                cols.append(manfest_df.index[id])
    ### these are the columns

    test1=tmp[cols].copy()
    filtered_count_df=df[(test1<rel_cut).mean(axis=1)<percent]
    filtered_df_read=df_read[(test1<rel_cut).mean(axis=1)<percent]

    drop_out_genes=set(input_df['gene_id'])-set(filtered_count_df['pbanka_id'])

    print('Number of Dropout genes at vitro level : %d' %len(drop_out_genes))
    print('Genes are :',drop_out_genes)


    # now we are going to check dropout gene at day 0 level

    day_cols=[]

    for k,v in manfest_df.groupby('b').indices.items():
        if k=='b1':
            for id in v:
                day_cols.append(manfest_df.index[id])

    #####
    tmp=filtered_count_df.copy()
    tmp=tmp.drop(columns=['Gene','Barcodes','pbanka_id'])
    tmp=tmp.div(tmp.sum(axis=0))
    test1=tmp[day_cols].copy()

    filtered_count_df=filtered_count_df[(test1<rel_cut).mean(axis=1)<percent]
    filtered_df_read=filtered_df_read[(test1<rel_cut).mean(axis=1)<percent]

    drop_out_genes=set(input_df['gene_id'])-set(filtered_count_df['pbanka_id'])

    print('Number of Dropout genes at b1 level : %d' %len(drop_out_genes))
    print('Genes are :',drop_out_genes)

    print('Number of final genes : %d' %len(set(input_df['gene_id'])&set(filtered_count_df['pbanka_id'])))
    contaminated_genes=list(set(filtered_count_df['pbanka_id'])-set(input_df['gene_id']))
    print('contaminated genes  :' ,contaminated_genes)

    ### we are going to drop contaminated genes
    filtered_count_df.set_index('pbanka_id',inplace=True)
    filtered_df_read.set_index('pbanka_id',inplace=True)

    filtered_count_df=filtered_count_df.drop(contaminated_genes)
    filtered_df_read=filtered_df_read.drop(contaminated_genes)

    df_s.set_index('pbanka_id',inplace=True)
    filltered_df_s=df_s.loc[filtered_count_df.index,:].copy()


    return filtered_count_df,filtered_df_read,filltered_df_s



def old_plot_prop_rel_abun_male_pool(mean_df_d0,var_df_d0,mean_df_d13,var_df_d13,geneConv,plot_info):
    ''' We will plot genes sex-specific wise '''

    #### input ####
    # sex_list: sex_list[0]=mean and sex_list[1]=SD
    # geneConv: when we want to give name of gene
    # out_pdf: will generated pdf file
    ########
    out_pdf=plot_info['rel_file']
    mf=plot_info['mf']
    day=plot_info['d']
    sex=plot_info['N']

    pdf = matplotlib.backends.backend_pdf.PdfPages(out_pdf)
    n=2
    genes=mean_df_d0.index;

    labels=[]
    for gene in genes:
        if gene in geneConv.keys():
            labels.append(geneConv[gene])
        else:
            labels.append(gene)

    ### get number of subplots

    num_subplots=2*len(genes) ### number of subplots

    l=len(genes)
    loops=int(l/(n*n))
    rem=l % (n*n)
    k=0

    per_page=num_subplots/16

    for loop in range(loops):
        fig = plt.figure(figsize=(15,15))
        for i in range(1,(n*n)+1):

            plt.subplot(n, n, i)
            x=[1,2]
            ### mf1 female
            y=[mean_df_d0.loc[genes[k],sex[0]+'_'+day[0]+'_NA'],mean_df_d13.loc[genes[k],sex[0]+'_'+day[1]+'_NA']]
            yerr=[var_df_d0.loc[genes[k],sex[0]+'_'+day[0]+'_NA'],var_df_d13.loc[genes[k],sex[0]+'_'+day[1]+'_NA']]
            plt.errorbar(x,y, yerr=yerr, fmt='r--')


            # ### mf1 male
            y=[mean_df_d0.loc[genes[k],sex[1]+'_'+day[0]+'_NA'],mean_df_d13.loc[genes[k],sex[1]+'_'+day[1]+'_NA']]
            yerr=[var_df_d0.loc[genes[k],sex[1]+'_'+day[0]+'_NA'],var_df_d13.loc[genes[k],sex[1]+'_'+day[1]+'_NA']]
            plt.errorbar(x,y, yerr=yerr, fmt='b--')


            plt.ylabel('log2 relative fitness')
            plt.title(labels[k])
            plt.legend((sex[0], sex[1]))
            plt.xticks([1, 1.5, 2],['day0', '', 'day14'],fontsize=15)

            plt.ylim(-18, 1)
            plt.grid(False)
            k=k+1
        pdf.savefig(fig)

    ## for the remaing one



    fig = plt.figure(figsize=(15,15))
    for i in range(1,rem+1):
        plt.subplot(n, n, i)
        x=[1,2]


        ### mf1 female
        y=[mean_df_d0.loc[genes[k],sex[0]+'_'+day[0]+'_NA'],mean_df_d13.loc[genes[k],sex[0]+'_'+day[1]+'_NA']]
        yerr=[var_df_d0.loc[genes[k],sex[0]+'_'+day[0]+'_NA'],var_df_d13.loc[genes[k],sex[0]+'_'+day[1]+'_NA']]
        plt.errorbar(x,y, yerr=yerr, fmt='r--')


        # ### mf1 male
        y=[mean_df_d0.loc[genes[k],sex[1]+'_'+day[0]+'_NA'],mean_df_d13.loc[genes[k],sex[1]+'_'+day[1]+'_NA']]
        yerr=[var_df_d0.loc[genes[k],sex[1]+'_'+day[0]+'_NA'],var_df_d13.loc[genes[k],sex[1]+'_'+day[1]+'_NA']]
        plt.errorbar(x,y, yerr=yerr, fmt='b--')




        plt.title(labels[k])
        plt.legend((sex[0], sex[1]))
        plt.xticks([1, 1.5, 2],['day0', '', 'day14'],fontsize=15)

        plt.ylim(-18, 1)
        plt.grid(False)
        k=k+1
        pdf.savefig(fig)
    pdf.close()




def plot_prop_rel_abun_male_pool(list_mean_df,list_var_df,geneConv,plot_info):
    ''' We will plot genes sex-specific wise '''

    #### input ####
    # sex_list: sex_list[0]=mean and sex_list[1]=SD
    # geneConv: when we want to give name of gene
    # out_pdf: will generated pdf file
    ########
    out_pdf=plot_info['rel_file']
    # # mf=plot_info['mf']
    # day=plot_info['d']
    # sex=plot_info['N']

    pdf = matplotlib.backends.backend_pdf.PdfPages(out_pdf)

    genes=list_mean_df[0].index;
    mini=2
    maxi=-18

    for df in list_mean_df:
        m=df.min()
        n=df.max()
        if mini>m:
            mini=m;
        if maxi<n:
            maxi=n;

    mini=int(mini-1)
    maxi=int(maxi+1)

    n=2
    labels=[]
    for gene in genes:
        if gene in geneConv.keys():
            labels.append(geneConv[gene])
        else:
            labels.append(gene)

    ### get number of subplots

    num_subplots=2*len(genes) ### number of subplots

    l=len(genes)
    loops=int(l/(n*n))
    rem=l % (n*n)
    k=0

    per_page=num_subplots/16

    for loop in range(loops):
        fig = plt.figure(figsize=(15,15))
        for i in range(1,(n*n)+1):

            plt.subplot(n, n, i)
            ## get relative abundance values

            y=[df.loc[genes[k]] for df in list_mean_df]
            yerr=[df.loc[genes[k]] for df in list_var_df]
            x=[ i+1 for i in range(len(y))]
            plt.errorbar(x,y, yerr=yerr, fmt='r*-')

            plt.ylabel('log2 relative fitness')
            plt.title(labels[k])
            # plt.legend((sex[0], sex[1]))
            plt.xticks([1, 1.5, 2,2.5,3,3.5,4],['b1', '', 'b2','','pel', '', 'sup'],fontsize=15)

            plt.ylim(mini, maxi)
            plt.grid(False)
            k=k+1
        pdf.savefig(fig)

    ## for the remaing one



    fig = plt.figure(figsize=(15,15))
    for i in range(1,rem+1):
        plt.subplot(n, n, i)
        x=[1,2]


        y=[df.loc[genes[k]] for df in list_mean_df]
        yerr=[df.loc[genes[k]] for df in list_var_df]
        x=[ i+1 for i in range(len(y))]
        plt.errorbar(x,y, yerr=yerr, fmt='r*-')

        plt.ylabel('log2 relative fitness')
        plt.title(labels[k])
        # plt.legend((sex[0], sex[1]))
        plt.xticks([1, 1.5, 2,2.5,3,3.5,4],['b1', '', 'b2','','pel', '', 'sup'],fontsize=15)
        plt.ylim(mini, maxi)
        # plt.ylim(-18, 2)
        plt.grid(False)
        k=k+1
        pdf.savefig(fig)
    pdf.close()



def calculate_RGR_male_pool(m_df_d0,var_df_d0,m_df_d13,var_df_d13,control_genes):
    ''' computing relative growth rate by divideing day13/day0'''
    # all relative abundance are on log scale

    rel_fitness=m_df_d13-m_df_d0
    rel_var=var_df_d13+var_df_d0

    control_gene_info={}

    control_fitness=rel_fitness[control_genes].copy().to_frame(name='values')
    control_var=rel_var[control_genes].copy().to_frame(name='values')
    l=gaussianMeanAndVariance(control_fitness.T,control_var.T)

    col='values'
    control_gene_info[col]=l # l[0] mean l[1]  SD   l[2] variance



    print(control_gene_info)

    #  we want to normalize by control genes
    normalized_fit=rel_fitness.copy()
    normalized_var=rel_var.copy()

    ctr_mean=control_gene_info[col][0][col]  #  0 mean
    ctr_sd=control_gene_info[col][2][col]  # 2 variance
    ctr_var=control_gene_info[col][2][col]  # 2 variance

    # this is the relative mean
    normalized_fit=rel_fitness-ctr_mean
    # relative variance on log scale
    normalized_var=rel_var+ctr_var

    return normalized_fit,normalized_var



def applyPhenocall_CI_male_pool(mdf,vdf, key='sup_vs_b1', lower_cut=-1):
    ''' This function is used for define the pheno_call '''
    #
    print('phenocall_test')

    m=mdf.copy()
    s=np.sqrt(vdf.copy())
    diff_max=m+2*s
    diff_min=m-2*s
    ##
    viz_df=pd.DataFrame(index=mdf.index)
    viz_df[key+'_RGR'] = m
    viz_df[key+'_sd'] = s
    viz_df[key+'_var'] = vdf.copy()
    viz_df[key+'_diff_max'] = m+2*s
    viz_df[key+'_diff_min'] = m-2*s

    viz_df[key+'_pheno'] = 'Not Reduced'


    red_df=viz_df[(viz_df[key+'_diff_max'] < lower_cut) & (viz_df[key+'_diff_min'] < lower_cut)]

    # # not_red_df=viz_df[(viz_df[key+'_diff_max'] < upper_cut) & (viz_df[key+'_diff_min'] > lower_cut)]
    # not_red_df=viz_df[(viz_df[key+'_diff_min'] > lower_cut)]
    # # increased_df=viz_df[(viz_df[key+'_diff_max'] > upper_cut) & (viz_df[key+'_diff_min'] > upper_cut)]

    viz_df.loc[red_df.index,key+'_pheno'] = 'Reduced'


    return viz_df


def plot_error_scatter(xseries,yseries,xlab,ylab,filename):

    plt.figure()
    plt.scatter(xseries.values,yseries.values, marker='o', c='b')
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.savefig(filename,
                format='pdf',dpi=300)
    plt.close()

def plot_error_hist(xseries,xlab,filename):

    plt.figure()
    plt.hist(xseries.values)
    plt.xlabel(xlab)
    plt.ylabel('Frequency')
    plt.savefig(filename,
                format='pdf',dpi=300)

    plt.close()




def relative_growth_rate_analysis_male_pool(df,manfest_df,prev_to_new,db_df,plot_info=None):
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



    ## now we are going to propagate error

    mean_df_b1,var_df_b1=propagate_error_male_pool(rel_df,manfest_df,grp_col=['b'],values=['b1'])
    mean_df_b2,var_df_b2=propagate_error_male_pool(rel_df,manfest_df,grp_col=['b'],values=['b2'])
    mean_df_pel,var_df_pel=propagate_error_male_pool(rel_df,manfest_df,grp_col=['pellet'],values=['pel'])
    mean_df_sup,var_df_sup=propagate_error_male_pool(rel_df,manfest_df,grp_col=['sup'],values=['sup'])
    list_mean_df=[mean_df_b1,mean_df_b2,mean_df_pel,mean_df_sup]
    list_var_df=[var_df_b1,var_df_b2,var_df_pel,var_df_sup]

    xlab='Relative abundance'
    ylab='Relative error'
    folder='/Users/vpandey/projects/githubs/Fertility_screen_2/preFinals/motility/'

    plot_error_scatter(mean_df_b1,var_df_b1,xlab,ylab,folder+'b1_scatter.pdf')
    plot_error_hist(mean_df_b1,xlab,folder+'b1_hist.pdf')
    plot_error_scatter(mean_df_b2,var_df_b2,xlab,ylab,folder+'b2_scatter.pdf')
    plot_error_hist(mean_df_b2,xlab,folder+'b2_hist.pdf')
    plot_error_scatter(mean_df_pel,var_df_pel,xlab,ylab,folder+'pel_scatter.pdf')
    plot_error_hist(mean_df_pel,xlab,folder+'pel_hist.pdf')
    plot_error_scatter(mean_df_sup,var_df_sup,xlab,ylab,folder+'sup_scatter.pdf')
    plot_error_hist(mean_df_sup,xlab,folder+'sup_hist.pdf')
    

    if 'rel_file' in plot_info.keys():
        print('plotting propagated relative abundance')
        geneConv=plot_info['geneConv']
        plot_prop_rel_abun_male_pool(list_mean_df,list_var_df,geneConv,plot_info)



    pel_vs_b2_RGR,pel_vs_b2_var=calculate_RGR_male_pool(mean_df_b2.copy(),var_df_b2.copy(),mean_df_pel.copy(),var_df_pel.copy(),ctrl_genes)
    viz_df1=applyPhenocall_CI_male_pool(pel_vs_b2_RGR,pel_vs_b2_var, key='pel_vs_b2', lower_cut=-1)

    sup_vs_pel_RGR,sup_vs_pel_var=calculate_RGR_male_pool(mean_df_pel.copy(),var_df_pel.copy(),mean_df_sup.copy(),var_df_sup.copy(),ctrl_genes)
    viz_df2=applyPhenocall_CI_male_pool(sup_vs_pel_RGR,sup_vs_pel_var, key='sup_vs_pel', lower_cut=-1)

    ## combine data
    viz_df=pd.concat([viz_df2,viz_df1],axis=1)



    ## anotate gene name and symbol
    symbols=[]
    gene_names=[]
    for idx in viz_df.index:
        tmp=db_df[db_df['Gene ID']==idx]

        if not tmp.empty:
            symbols.append(tmp['Gene Name or Symbol'].to_list()[0])
            gene_names.append(tmp['Product Description'].to_list()[0])
        else:
            symbols.append('NA')
            gene_names.append('NA')

    viz_df['Gene symbol']=symbols
    viz_df['Gene description']=gene_names

    if 'rgr_file' in plot_info.keys():
        print('writing rgr file---')
        #viz_df.to_csv(plot_info['rgr_file'],sep='\t')
        viz_df.to_excel(plot_info['rgr_file'],sheet_name='RGR_analysis')
    else:
        print('give file name of rgr file---')
    import pdb; pdb.set_trace()
    # cmb_fitness={}
    # backgrounds=['sup_vs_b1']
    #
    # for b in backgrounds:
    #     cmb_fitness[b]=gaussianMeanAndVariance(rgr_temp,var_temp)

    ## calculate combined fitness
    pheno_call_df=applyPhenocall_CI(cmb_fitness,lower_cut=np.log2(0.45),upper_cut=np.log2(2.05))
    # pheno_call_df=getPvalZscore(cmb_fitness,upcut=1,lowcut=0.4,pval=0.05,pval1=0.05)

    ### aplly filter of relative input
    pheno_call_df=apply_filter_testInput(pheno_call_df,mean_df_d0_mf1,mean_df_d0_mf2,rel_cut=-12)

    ### once again apply filter on feeds
    pheno_call_df2=apply_filter_on_feeds(pheno_call_df,mf1_RGR,mf1_var,mf2_RGR,mf2_var,mean_df_d0_mf1,mean_df_d0_mf2)

    # test variance of comined and variance with each step using cutoff
    # rel_cut=-12
    # testforRepeat(pheno_call_df,mean_df_d0_mf1,var_df_d0_mf1,mean_df_d0_mf2,var_df_d0_mf2,mean_df_d13_mf1,var_df_d13_mf1,mean_df_d13_mf2,var_df_d13_mf2,cuvette_mean_df,cuvette_var_df,rel_cut,plot_info)

    ## for pool2 and pool4

    if (plot_info['pool']=='pool2') or (plot_info['pool']=='pool4'):

        grp_cols=['sex','d','mf','dr','e','t']
        day_pos=grp_cols.index('d')
        gDNA_mean_df, gDNA_mean_df_var=propagate_error_gDNA_extraction_method(rel_df,manfest_df,grp_cols,day_pos)
        plot_gDNA_error(gDNA_mean_df, gDNA_mean_df_var)
    ##

    plot_each_step_mean_var(mean_df_d0_mf1,var_df_d0_mf1,mean_df_d0_mf2,var_df_d0_mf2,mean_df_d13_mf1,var_df_d13_mf1,mean_df_d13_mf2,var_df_d13_mf2)

    ## compute relative ratio betwen day13 and day0
    trace_d0_female = go.Scatter(
        x = mean_df_d0_log['GCKO2_d0_NA_NA'],
        y = var_df_d0_log['GCKO2_d0_NA_NA'],
        mode = 'markers',
        marker=dict(size=5,color='red'),
        name='(# of markers=%d)'%mean_df_d0_log.shape[0],
        opacity=0.7,
        text=mean_df_d0_log.index)

    trace_d0_male = go.Scatter(
        x = mean_df_d0_log['g145480_d0_NA_NA'],
        y = var_df_d0_log['g145480_d0_NA_NA'],
        mode = 'markers',
        marker=dict(size=5,color='red'),
        name='(# of markers=%d)'%mean_df_d0_log.shape[0],
        opacity=0.7,
        text=mean_df_d0_log.index)

    trace_d13_female = go.Scatter(
        x = mean_df_d13_log['GCKO2_d13_NA_NA'],
        y = var_df_d13_log['GCKO2_d13_NA_NA'],
        mode = 'markers',
        marker=dict(size=5,color='red'),
        name='(# of markers=%d)'%mean_df_d13_log.shape[0],
        opacity=0.7,
        text=mean_df_d13_log.index)

    trace_d13_male = go.Scatter(
        x = mean_df_d13_log['g145480_d13_NA_NA'],
        y = var_df_d13_log['g145480_d13_NA_NA'],
        mode = 'markers',
        marker=dict(size=5,color='red'),
        name='(# of markers=%d)'%mean_df_d13_log.shape[0],
        opacity=0.7,
        text=mean_df_d13_log.index)

    fig = make_subplots(rows=2, cols=2,subplot_titles=("Input error GCKO2(day0)", "Input error 145480(day0)","Output error GCKO2(day13)","Output error 145480(day13)"))

    fig.append_trace(trace_d0_female,row=1, col=1)
    fig.append_trace(trace_d0_male,row=1, col=2)
    fig.append_trace(trace_d13_female,row=2, col=1)

    fig.append_trace(trace_d13_male,row=2, col=2)

    # fig.show()

    return pheno_call_df2,[mean_df_d0_mf1,mean_df_d0_mf2]




def propagate_error_male_pool(df,manfest_df,grp_col=['b'],values=['b1']):
    ''' We are going to calculate mean and SD for combined analyis from PCR to mosquitofeed '''
    df_log=np.log2(df)
    tmp_mn=pd.DataFrame(index=df.index)
    tmp_var=pd.DataFrame(index=df.index)
    tmp_mn_log=pd.DataFrame(index=df.index)
    ###
    final_mean=pd.DataFrame(index=df.index)
    final_var=pd.DataFrame(index=df.index)
    df_sample=pd.DataFrame(index=['num'])

    for k,v in manfest_df.groupby(grp_col).indices.items():
        if k==values[0]:
            tmp_mani=manfest_df.iloc[v].copy()
            for k1,v1 in tmp_mani.groupby(['pr']).indices.items():

                key=k+'_'+k1
                tmp_mani2=tmp_mani.iloc[v1].copy()
                tmp_mn[key]=df[tmp_mani2.index].mean(axis=1).copy()
                tmp_var[key]=df[tmp_mani2.index].var(axis=1).copy()


    [mean_df,sd_max,var_max]=getCombined_mean_variance(tmp_mn,tmp_var,df_sample)

            # [mean_df,sd_max,var_max]=weighted_mean_variance(tmp_mn_log,tmp_var)
    final_mean=mean_df.copy()
    final_var=var_max.copy()

    final_mean_log=np.log2(final_mean)
    final_var_log=(final_var)/((final_mean*final_mean)*((np.log(10)**2)*np.log10(2)))
    return final_mean_log,final_var_log









def stepwiseAnalysis():
    ''' We are going do analyis of pool1 data of Claire '''

    ### these are the input files
    manifests_df=pd.read_csv(data_folder+"/manifest_male_pool_dis.txt",sep='\t')
    count_df=pd.read_csv(data_folder+ "/barcode_counts_table_050221_male_pool.txt",sep='\t')
    input_df=pd.read_csv(data_folder+'/input_male_pool_dis.txt', sep='\t')

    oldIds=[]
    for item in input_df['gene_id'].to_list():
        if item in new_to_prev.keys():
            oldIds.append(new_to_prev[item])
        else:
            print (item)
            oldIds.append(item)

    input_df['gene_id']=oldIds
    ## change back to old ids



    #### end of the input section
    # final_count_df: read1 and read2 are added
    # final_count_df_two_read: reads are sperated
    # manfest_df: maifest_df
    final_count_df,final_count_df_des,final_count_df_two_read,manfest_df=preprocessing(manifests_df,count_df)

    ####  Dropouts and input check
    # input_df: these are genes which was used for pool phenotypes
    percent=0.9 ## this parameters is used to test whether count is too small for 90 % of input samples. Those will be deleted.
    filtered_count_df,filtered_df_read,filtered_count_df_des=filter_input_dropout_male_pool(final_count_df,final_count_df_des,final_count_df_two_read,input_df,manfest_df,percent)

    ######  write filtered and unfiltered files
    # final_count_df_two_read.to_csv(out_folder+"/unfilterd_count_matrix_pool2.txt",sep='\t')
    filtered_count_df_des.to_csv(out_folder+"/filterd_count_matrix_male_pool_dis.txt",sep='\t')

    geneConv,old_to_new_ids,geneConv_new=getNewIdfromPrevID(filtered_count_df.index,prev_to_new,db_df)

    plot_info={'pool':'male_pool_des','rel_file':out_folder+'/male_pool_des_propagated_error_relative_abundance.pdf','geneConv':geneConv_new,
    'control_genes':['PBANKA_1024600' , 'PBANKA_0501200' , 'PBANKA_0101100' , 'PBANKA_1422100','PBANKA_0616700','PBANKA_1359700'],
    'rgr_file':out_folder+'/male_pool_des_RGR.xlsx'}
    relative_growth_rate_analysis_male_pool(filtered_count_df,manfest_df,prev_to_new,db_df,plot_info)

def stepwiseAnalysis_slowpool():
    ''' We are going do analyis of pool1 data of Claire '''

    ### these are the input files
    manifests_df=pd.read_csv(data_folder+"/manifest_male_pool_slow.txt",sep='\t')
    count_df=pd.read_csv(data_folder+ "/barcode_counts_table_050221_male_pool.txt",sep='\t')
    input_df=pd.read_csv(data_folder+'/input_male_pool_slow.txt', sep='\t')

    oldIds=[]
    for item in input_df['gene_id'].to_list():
        if item in new_to_prev.keys():
            oldIds.append(new_to_prev[item])
        else:
            print (item)
            oldIds.append(item)

    input_df['gene_id']=oldIds
    ## change back to old ids



    #### end of the input section
    # final_count_df: read1 and read2 are added
    # final_count_df_two_read: reads are sperated
    # manfest_df: maifest_df
    final_count_df,final_count_df_des,final_count_df_two_read,manfest_df=preprocessing(manifests_df,count_df)

    ####  Dropouts and input check
    # input_df: these are genes which was used for pool phenotypes
    percent=0.9 ## this parameters is used to test whether count is too small for 90 % of input samples. Those will be deleted.
    filtered_count_df,filtered_df_read,filtered_count_df_des=filter_input_dropout_male_pool(final_count_df,final_count_df_des,final_count_df_two_read,input_df,manfest_df,percent)

    ######  write filtered and unfiltered files
    # final_count_df_two_read.to_csv(out_folder+"/unfilterd_count_matrix_pool2.txt",sep='\t')
    filtered_count_df_des.to_csv(out_folder+"/filterd_count_matrix_male_pool_dis.txt",sep='\t')

    geneConv,old_to_new_ids,geneConv_new=getNewIdfromPrevID(filtered_count_df.index,prev_to_new,db_df)

    plot_info={'pool':'male_pool_slow','rel_file':out_folder+'/male_pool_slow_propagated_error_relative_abundance.pdf','geneConv':geneConv_new,
    'control_genes':['PBANKA_1024600' , 'PBANKA_0501200' , 'PBANKA_0101100' , 'PBANKA_1422100','PBANKA_0616700','PBANKA_1359700'],
    'rgr_file':out_folder+'/male_pool_slow_RGR.xlsx'}
    relative_growth_rate_analysis_male_pool(filtered_count_df,manfest_df,prev_to_new,db_df,plot_info)


if __name__ == '__main__':
    stepwiseAnalysis()
    #stepwiseAnalysis_slowpool()
