import pandas as pd
import numpy as np
import pickle
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
import statsmodels.api as sm
from scipy import stats
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from sklearn.neighbors import KernelDensity
import scipy.stats as st
import scipy.special as sp
from special_fun_fertility import (getNewIdfromPrevID,propagate_error_day0,propagate_error_day13,
propagate_error_day0_each_mossifeed,propagate_error_day13_each_mossifeed,plot_each_step_mean_var,
calculate_RGR,getColumnsFormDF,applyPhenocall_CI,apply_filter_testInput,apply_filter_on_feeds,gaussianMeanAndVariance,kde_sklearn)


def relative_growth_rate_analysis_repeat(df,manfest_df,prev_to_new,db_df,plot_info=None):
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
    grp_cols=['sex','d','mf','dr','b','t']
    day_pos=grp_cols.index('d')

    mean_df_d0,var_df_d0=propagate_error_day0(rel_df,manfest_df,grp_cols,day_pos)
    grp_cols=['sex','d','mf','dr','b','t']
    day_pos=grp_cols.index('d')
    mean_df_d13,var_df_d13=propagate_error_day13(rel_df,manfest_df,grp_cols,day_pos)

    #### now we need to tarnsform at log2 scale
    mean_df_d0_log=np.log2(mean_df_d0)
    var_df_d0_log=(var_df_d0)/((mean_df_d0*mean_df_d0)*((np.log(10)**2)*np.log10(2)))

    mean_df_d13_log=np.log2(mean_df_d13)
    var_df_d13_log=(var_df_d13)/((mean_df_d13*mean_df_d13)*((np.log(10)**2)*np.log10(2)))

    ### we are going to compute each mosquito feed

    grp_cols=['sex','d','mf','dr','b','t']
    day_pos=grp_cols.index('d')
    mean_df_d0_mf1,var_df_d0_mf1,mean_df_d0_mf2,var_df_d0_mf2=propagate_error_day0_each_mossifeed(rel_df,manfest_df,grp_cols,day_pos)
    grp_cols=['sex','d','mf','dr','b','t']
    day_pos=grp_cols.index('d')
    mean_df_d13_mf1,var_df_d13_mf1,mean_df_d13_mf2,var_df_d13_mf2=propagate_error_day13_each_mossifeed(rel_df,manfest_df,grp_cols,day_pos)

    saveres=[[mean_df_d0_mf1,var_df_d0_mf1,mean_df_d0_mf2,var_df_d0_mf2],[mean_df_d13_mf1,var_df_d13_mf1,mean_df_d13_mf2,var_df_d13_mf2]]

    pickle.dump(saveres,open(plot_info['pool']+'.p','wb'))


    ##  we are going to compute rleative growth rate  ###

    ### plot propagated relative abundance.

    propagated_relative_abundance_plot_each_sex(mean_df_d0_mf1,var_df_d0_mf1,mean_df_d0_mf2,var_df_d0_mf2,mean_df_d13_mf1,var_df_d13_mf1,mean_df_d13_mf2,var_df_d13_mf2,plot_info)

    ## propgate error for transfection cuvette
    # grp_cols=['sex','d','mf','dr','b','t']
    # day_pos=grp_cols.index('d')
    # days=['NA']
    #cuvette_mean_df,cuvette_var_df=propagate_error_cuvette(rel_df,manfest_df,grp_cols,day_pos,days)
    ##

    mf1_RGR,mf1_var=calculate_RGR(mean_df_d0_mf1.copy(),var_df_d0_mf1.copy(),mean_df_d13_mf1.copy(),var_df_d13_mf1.copy(),ctrl_genes)
    mf2_RGR,mf2_var=calculate_RGR(mean_df_d0_mf2.copy(),var_df_d0_mf2.copy(),mean_df_d13_mf2.copy(),var_df_d13_mf2.copy(),ctrl_genes)
    ### now combined fitness for mf1 and mf2
    # take mf1 and mf2 in one dtaframe


    cmb_fitness={}
    backgrounds=plot_info['sex']

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
    pheno_call_df=applyPhenocall_CI(cmb_fitness,lower_cut=np.log2(0.45),upper_cut=np.log2(2.05))
    # pheno_call_df=getPvalZscore(cmb_fitness,upcut=1,lowcut=0.4,pval=0.05,pval1=0.05)

    ### aplly filter of relative input
    pheno_call_df=apply_filter_testInput(pheno_call_df,mean_df_d0_mf1,mean_df_d0_mf2,rel_cut=-12,sex=plot_info['sex'])

    ### once again apply filter on feeds
    pheno_call_df2=apply_filter_on_feeds(pheno_call_df,mf1_RGR,mf1_var,mf2_RGR,mf2_var,mean_df_d0_mf1,mean_df_d0_mf2,backgrounds=plot_info['sex'])

    # test variance of comined and variance with each step using cutoff
    # rel_cut=-12
    # testforRepeat(pheno_call_df,mean_df_d0_mf1,var_df_d0_mf1,mean_df_d0_mf2,var_df_d0_mf2,mean_df_d13_mf1,var_df_d13_mf1,mean_df_d13_mf2,var_df_d13_mf2,cuvette_mean_df,cuvette_var_df,rel_cut,plot_info)

    ## for pool2 and pool4

    # if (plot_info['pool']=='pool2') or (plot_info['pool']=='pool4'):
    #
    #     grp_cols=['sex','d','mf','dr','e','t']
    #     day_pos=grp_cols.index('d')
    #     gDNA_mean_df, gDNA_mean_df_var=propagate_error_gDNA_extraction_method(rel_df,manfest_df,grp_cols,day_pos)
    #     plot_gDNA_error(gDNA_mean_df, gDNA_mean_df_var)
    ##

    plot_each_step_mean_var_sex(mean_df_d0_mf1,var_df_d0_mf1,mean_df_d0_mf2,var_df_d0_mf2,mean_df_d13_mf1,var_df_d13_mf1,mean_df_d13_mf2,var_df_d13_mf2)

    ## compute relative ratio betwen day13 and day0
    # trace_d0_female = go.Scatter(
    #     x = mean_df_d0_log['GCKO2_d0_NA_NA'],
    #     y = var_df_d0_log['GCKO2_d0_NA_NA'],
    #     mode = 'markers',
    #     marker=dict(size=5,color='red'),
    #     name='(# of markers=%d)'%mean_df_d0_log.shape[0],
    #     opacity=0.7,
    #     text=mean_df_d0_log.index)
    #
    # trace_d0_male = go.Scatter(
    #     x = mean_df_d0_log['g145480_d0_NA_NA'],
    #     y = var_df_d0_log['g145480_d0_NA_NA'],
    #     mode = 'markers',
    #     marker=dict(size=5,color='red'),
    #     name='(# of markers=%d)'%mean_df_d0_log.shape[0],
    #     opacity=0.7,
    #     text=mean_df_d0_log.index)
    #
    # trace_d13_female = go.Scatter(
    #     x = mean_df_d13_log['GCKO2_d13_NA_NA'],
    #     y = var_df_d13_log['GCKO2_d13_NA_NA'],
    #     mode = 'markers',
    #     marker=dict(size=5,color='red'),
    #     name='(# of markers=%d)'%mean_df_d13_log.shape[0],
    #     opacity=0.7,
    #     text=mean_df_d13_log.index)
    #
    # trace_d13_male = go.Scatter(
    #     x = mean_df_d13_log['g145480_d13_NA_NA'],
    #     y = var_df_d13_log['g145480_d13_NA_NA'],
    #     mode = 'markers',
    #     marker=dict(size=5,color='red'),
    #     name='(# of markers=%d)'%mean_df_d13_log.shape[0],
    #     opacity=0.7,
    #     text=mean_df_d13_log.index)
    #
    # fig = make_subplots(rows=2, cols=2,subplot_titles=("Input error GCKO2(day0)", "Input error 145480(day0)","Output error GCKO2(day13)","Output error 145480(day13)"))
    #
    # fig.append_trace(trace_d0_female,row=1, col=1)
    # fig.append_trace(trace_d0_male,row=1, col=2)
    # fig.append_trace(trace_d13_female,row=2, col=1)
    #
    # fig.append_trace(trace_d13_male,row=2, col=2)

    # fig.show()

    return pheno_call_df2,[mean_df_d0_mf1,mean_df_d0_mf2]



def propagated_relative_abundance_plot_each_sex(mean_df_d0_mf1,var_df_d0_mf1,mean_df_d0_mf2,var_df_d0_mf2,mean_df_d13_mf1,var_df_d13_mf1,mean_df_d13_mf2,var_df_d13_mf2,plot_info):
     '''  We are going to plot propagated relative abundance plot  '''

     if 'rel_file' in plot_info.keys():
         print('plotting propagated relative abundance')
         geneConv=plot_info['geneConv']
         #plot_propgated_relative_abunndance_sex(mean_df_d0_mf1,var_df_d0_mf1,mean_df_d0_mf2,var_df_d0_mf2,mean_df_d13_mf1,var_df_d13_mf1,mean_df_d13_mf2,var_df_d13_mf2,geneConv,plot_info)




def plot_propgated_relative_abunndance_sex(mean_df_d0_mf1,var_df_d0_mf1,mean_df_d0_mf2,var_df_d0_mf2,mean_df_d13_mf1,var_df_d13_mf1,mean_df_d13_mf2,var_df_d13_mf2,geneConv,plot_info):
    ''' We will plot genes sex-specific wise '''

    #### input ####
    # sex_list: sex_list[0]=mean and sex_list[1]=SD
    # geneConv: when we want to give name of gene
    # out_pdf: will generated pdf file
    ########
    out_pdf=plot_info['rel_file']
    mf=plot_info['mf']
    day=plot_info['d']
    sex=plot_info['sex']

    pdf = matplotlib.backends.backend_pdf.PdfPages(out_pdf)
    n=2
    genes=mean_df_d0_mf1.index;

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
            y=[mean_df_d0_mf1.loc[genes[k],sex[0]+'_'+day[0]+'_NA_NA'],mean_df_d13_mf1.loc[genes[k],sex[0]+'_'+day[1]+'_NA_NA']]
            yerr=[var_df_d0_mf1.loc[genes[k],sex[0]+'_'+day[0]+'_NA_NA'],var_df_d13_mf1.loc[genes[k],sex[0]+'_'+day[1]+'_NA_NA']]
            plt.errorbar(x,y, yerr=yerr, fmt='r--')

            ### mf2 female
            y=[mean_df_d0_mf2.loc[genes[k],sex[0]+'_'+day[0]+'_NA_NA'],mean_df_d13_mf2.loc[genes[k],sex[0]+'_'+day[1]+'_NA_NA']]
            yerr=[var_df_d0_mf2.loc[genes[k],sex[0]+'_'+day[0]+'_NA_NA'],var_df_d13_mf2.loc[genes[k],sex[0]+'_'+day[1]+'_NA_NA']]
            plt.errorbar(x,y, yerr=yerr, fmt='r-')

            # ### mf1 male
            # y=[mean_df_d0_mf1.loc[genes[k],sex[1]+'_'+day[0]+'_NA_NA'],mean_df_d13_mf1.loc[genes[k],sex[1]+'_'+day[1]+'_NA_NA']]
            # yerr=[var_df_d0_mf1.loc[genes[k],sex[1]+'_'+day[0]+'_NA_NA'],var_df_d13_mf1.loc[genes[k],sex[1]+'_'+day[1]+'_NA_NA']]
            # plt.errorbar(x,y, yerr=yerr, fmt='b--')
            #
            # ### mf2 male
            # y=[mean_df_d0_mf2.loc[genes[k],sex[1]+'_'+day[0]+'_NA_NA'],mean_df_d13_mf2.loc[genes[k],sex[1]+'_'+day[1]+'_NA_NA']]
            # yerr=[var_df_d0_mf2.loc[genes[k],sex[1]+'_'+day[0]+'_NA_NA'],var_df_d13_mf2.loc[genes[k],sex[1]+'_'+day[1]+'_NA_NA']]
            # plt.errorbar(x,y, yerr=yerr, fmt='b-')
            #

            plt.ylabel('log2 relative fitness')
            plt.title(labels[k])
            plt.legend(('mf1_GCKO2', 'mf2_GCKO2'))
            plt.xticks([1, 1.5, 2],['day0', '', 'day13'],fontsize=15)

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
        y=[mean_df_d0_mf1.loc[genes[k],sex[0]+'_'+day[0]+'_NA_NA'],mean_df_d13_mf1.loc[genes[k],sex[0]+'_'+day[1]+'_NA_NA']]
        yerr=[var_df_d0_mf1.loc[genes[k],sex[0]+'_'+day[0]+'_NA_NA'],var_df_d13_mf1.loc[genes[k],sex[0]+'_'+day[1]+'_NA_NA']]
        plt.errorbar(x,y, yerr=yerr, fmt='r--')

        ### mf2 female
        y=[mean_df_d0_mf2.loc[genes[k],sex[0]+'_'+day[0]+'_NA_NA'],mean_df_d13_mf2.loc[genes[k],sex[0]+'_'+day[1]+'_NA_NA']]
        yerr=[var_df_d0_mf2.loc[genes[k],sex[0]+'_'+day[0]+'_NA_NA'],var_df_d13_mf2.loc[genes[k],sex[0]+'_'+day[1]+'_NA_NA']]
        plt.errorbar(x,y, yerr=yerr, fmt='r-')

        # ### mf1 male
        # y=[mean_df_d0_mf1.loc[genes[k],sex[1]+'_'+day[0]+'_NA_NA'],mean_df_d13_mf1.loc[genes[k],sex[1]+'_'+day[1]+'_NA_NA']]
        # yerr=[var_df_d0_mf1.loc[genes[k],sex[1]+'_'+day[0]+'_NA_NA'],var_df_d13_mf1.loc[genes[k],sex[1]+'_'+day[1]+'_NA_NA']]
        # plt.errorbar(x,y, yerr=yerr, fmt='b--')
        #
        # ### mf2 male
        # y=[mean_df_d0_mf2.loc[genes[k],sex[1]+'_'+day[0]+'_NA_NA'],mean_df_d13_mf2.loc[genes[k],sex[1]+'_'+day[1]+'_NA_NA']]
        # yerr=[var_df_d0_mf2.loc[genes[k],sex[1]+'_'+day[0]+'_NA_NA'],var_df_d13_mf2.loc[genes[k],sex[1]+'_'+day[1]+'_NA_NA']]
        # plt.errorbar(x,y, yerr=yerr, fmt='b-')



        plt.ylabel('log2 relative fitness')
        plt.title(labels[k])
        plt.legend(('mf1_GCKO2', 'mf2_GCKO2'))
        plt.xticks([1, 1.5, 2],['day0', '', 'day13'],fontsize=15)

        plt.ylim(-18, 1)
        plt.grid(False)
        k=k+1
        pdf.savefig(fig)
    pdf.close()




def plot_each_step_mean_var_sex(mean_df_d0_mf1,var_df_d0_mf1,mean_df_d0_mf2,var_df_d0_mf2,mean_df_d13_mf1,var_df_d13_mf1,mean_df_d13_mf2,var_df_d13_mf2):
    ## compute relative ratio betwen day13 and day0
    trace_d0_female_mf1 = go.Scatter(
        x = mean_df_d0_mf1['GCKO2_d0_NA_NA'],
        y = var_df_d0_mf1['GCKO2_d0_NA_NA'],
        mode = 'markers',
        marker=dict(size=5,color='red'),
        name='(# of markers=%d)'%mean_df_d0_mf1.shape[0],
        opacity=0.7,
        text=mean_df_d0_mf1.index)
    trace_d0_female_mf2 = go.Scatter(
        x = mean_df_d0_mf2['GCKO2_d0_NA_NA'],
        y = var_df_d0_mf2['GCKO2_d0_NA_NA'],
        mode = 'markers',
        marker=dict(size=5,color='red'),
        name='(# of markers=%d)'%mean_df_d0_mf2.shape[0],
        opacity=0.7,
        text=mean_df_d0_mf2.index)

    # trace_d0_male_mf1 = go.Scatter(
    #     x = mean_df_d0_mf1['g145480_d0_NA_NA'],
    #     y = var_df_d0_mf1['g145480_d0_NA_NA'],
    #     mode = 'markers',
    #     marker=dict(size=5,color='red'),
    #     name='(# of markers=%d)'%mean_df_d0_mf1.shape[0],
    #     opacity=0.7,
    #     text=mean_df_d0_mf1.index)
    #
    # trace_d0_male_mf2 = go.Scatter(
    #     x = mean_df_d0_mf2['g145480_d0_NA_NA'],
    #     y = var_df_d0_mf2['g145480_d0_NA_NA'],
    #     mode = 'markers',
    #     marker=dict(size=5,color='red'),
    #     name='(# of markers=%d)'%mean_df_d0_mf2.shape[0],
    #     opacity=0.7,
    #     text=mean_df_d0_mf2.index)

    trace_d13_female_mf1 = go.Scatter(
        x = mean_df_d13_mf1['GCKO2_d13_NA_NA'],
        y = var_df_d13_mf1['GCKO2_d13_NA_NA'],
        mode = 'markers',
        marker=dict(size=5,color='red'),
        name='(# of markers=%d)'%mean_df_d13_mf1.shape[0],
        opacity=0.7,
        text=mean_df_d13_mf1.index)

    trace_d13_female_mf2 = go.Scatter(
        x = mean_df_d13_mf2['GCKO2_d13_NA_NA'],
        y = var_df_d13_mf2['GCKO2_d13_NA_NA'],
        mode = 'markers',
        marker=dict(size=5,color='red'),
        name='(# of markers=%d)'%mean_df_d13_mf2.shape[0],
        opacity=0.7,
        text=mean_df_d13_mf2.index)

    # trace_d13_male_mf1 = go.Scatter(
    #     x = mean_df_d13_mf1['g145480_d13_NA_NA'],
    #     y = var_df_d13_mf1['g145480_d13_NA_NA'],
    #     mode = 'markers',
    #     marker=dict(size=5,color='red'),
    #     name='(# of markers=%d)'%mean_df_d13_mf1.shape[0],
    #     opacity=0.7,
    #     text=mean_df_d13_mf1.index)
    #
    # trace_d13_male_mf2 = go.Scatter(
    #     x = mean_df_d13_mf2['g145480_d13_NA_NA'],
    #     y = var_df_d13_mf2['g145480_d13_NA_NA'],
    #     mode = 'markers',
    #     marker=dict(size=5,color='red'),
    #     name='(# of markers=%d)'%mean_df_d13_mf2.shape[0],
    #     opacity=0.7,
    #     text=mean_df_d13_mf2.index)

    ### distributions
    x=mean_df_d0_mf1['GCKO2_d0_NA_NA'].values
    x_grid=np.linspace(x.min(), x.max(), 500)
    trace_d0_mf1_dist_female = go.Histogram( x=x,
    name='dist',
    marker_color='#bdbdbd',
    opacity=0.75,
    histnorm='probability')

    trace_d0_mf1_kde_female= go.Scatter(x=x_grid, y=kde_sklearn(x, x_grid, bandwidth=0.8),
                    mode='lines',name='kde')

    ### distributions
    x=mean_df_d0_mf2['GCKO2_d0_NA_NA'].values
    x_grid=np.linspace(x.min(), x.max(), 500)
    trace_d0_mf2_dist_female = go.Histogram( x=x,
    name='dist',
    marker_color='#bdbdbd',
    opacity=0.75,
    histnorm='probability')

    trace_d0_mf2_kde_female= go.Scatter(x=x_grid, y=kde_sklearn(x, x_grid, bandwidth=0.8),
                    mode='lines',name='kde')
    # ### distributions
    # x=mean_df_d0_mf1['g145480_d0_NA_NA'].values
    # x_grid=np.linspace(x.min(), x.max(), 500)
    # trace_d0_mf1_dist_male = go.Histogram( x=x,
    # name='dist',
    # marker_color='#bdbdbd',
    # opacity=0.75,
    # histnorm='probability')
    #
    # trace_d0_mf1_kde_male= go.Scatter(x=x_grid, y=kde_sklearn(x, x_grid, bandwidth=0.8),
    #                 mode='lines',name='kde')
    # ### distributions
    # x=mean_df_d0_mf2['g145480_d0_NA_NA'].values
    # x_grid=np.linspace(x.min(), x.max(), 500)
    # trace_d0_mf2_dist_male = go.Histogram( x=x,
    # name='dist',
    # marker_color='#bdbdbd',
    # opacity=0.75,
    # histnorm='probability')
    #
    # trace_d0_mf2_kde_male= go.Scatter(x=x_grid, y=kde_sklearn(x, x_grid, bandwidth=0.8),
    #                 mode='lines',name='kde')


    ### distributions
    x=mean_df_d13_mf1['GCKO2_d13_NA_NA'].values
    x_grid=np.linspace(x.min(), x.max(), 500)
    trace_d13_mf1_dist_female = go.Histogram( x=x,
    name='dist',
    marker_color='#bdbdbd',
    opacity=0.75,
    histnorm='probability')

    trace_d13_mf1_kde_female= go.Scatter(x=x_grid, y=kde_sklearn(x, x_grid, bandwidth=0.8),
                    mode='lines',name='kde')

    x=mean_df_d13_mf2['GCKO2_d13_NA_NA'].values
    x_grid=np.linspace(x.min(), x.max(), 500)
    trace_d13_mf2_dist_female = go.Histogram( x=x,
    name='dist',
    marker_color='#bdbdbd',
    opacity=0.75,
    histnorm='probability')

    trace_d13_mf2_kde_female= go.Scatter(x=x_grid, y=kde_sklearn(x, x_grid, bandwidth=0.8),
                    mode='lines',name='kde')


    # ### distributions
    # x=mean_df_d13_mf1['g145480_d13_NA_NA'].values
    # x_grid=np.linspace(x.min(), x.max(), 500)
    # trace_d13_mf1_dist_male = go.Histogram( x=x,
    # name='dist',
    # marker_color='#bdbdbd',
    # opacity=0.75,
    # histnorm='probability')
    #
    # trace_d13_mf1_kde_male= go.Scatter(x=x_grid, y=kde_sklearn(x, x_grid, bandwidth=0.8),
    #                 mode='lines',name='kde')
    #
    # x=mean_df_d13_mf2['g145480_d13_NA_NA'].values
    # x_grid=np.linspace(x.min(), x.max(), 500)
    # trace_d13_mf2_dist_male = go.Histogram( x=x,
    # name='dist',
    # marker_color='#bdbdbd',
    # opacity=0.75,
    # histnorm='probability')
    #
    # trace_d13_mf2_kde_male= go.Scatter(x=x_grid, y=kde_sklearn(x, x_grid, bandwidth=0.8),
    #                 mode='lines',name='kde')





    # fig, ax = plt.subplots()
    # for bandwidth in [0.8,0.9,1.0]:
    #     ax.plot(x_grid, kde_sklearn(x, x_grid, bandwidth=bandwidth),
    #         label='bandwidth={0}'.format(bandwidth), linewidth=3, alpha=0.5)
    # ax.hist(x, 50, fc='gray', histtype='stepfilled', alpha=0.3, normed=True)
    # ax.set_xlim(-1, 32)
    # ax.legend(loc='upper left')
    # # plt.ylim(0, 0.25)
    # plt.xlabel('time')
    # plt.ylabel('distributions of genes that have steepest slope')
    # plt.savefig('/Users/vikash/Documents/Projects/GeneRegulatory/Final_analysis/io/SOMap/figures/wave_continuous_sigmoidal_screen.pdf')

    fig = make_subplots(rows=4, cols=2,subplot_titles=("Input error GCKO2(day0) mf1",
    "Distribution GCKO2(day0) mf1","Input error GCKO2(day0) mf2",
    "Distribution GCKO2(day0) mf2","Output error GCKO2(day13) mf1","Distribution GCKO2(day13) mf1",
    "Output error GCKO2(day13) mf2","Distribution GCKO2(day13) mf2"))

    fig.append_trace(trace_d0_female_mf1,row=1, col=1)
    # fig.append_trace(trace_d0_male_mf1,row=1, col=2)
    fig.append_trace(trace_d0_mf1_dist_female,row=1, col=2)
    fig.append_trace(trace_d0_mf1_kde_female,row=1, col=2)
    # fig.append_trace(trace_d0_mf1_dist_male,row=1, col=4)
    # fig.append_trace(trace_d0_mf1_kde_male,row=1, col=4)


    fig.append_trace(trace_d0_female_mf2,row=2, col=1)
    # fig.append_trace(trace_d0_male_mf2,row=2, col=2)
    fig.append_trace(trace_d0_mf2_dist_female,row=2, col=2)
    fig.append_trace(trace_d0_mf2_kde_female,row=2, col=2)
    # fig.append_trace(trace_d0_mf2_dist_male,row=2, col=4)
    # fig.append_trace(trace_d0_mf2_kde_male,row=2, col=4)



    fig.append_trace(trace_d13_female_mf1,row=3, col=1)
    # fig.append_trace(trace_d13_male_mf1,row=3, col=2)
    fig.append_trace(trace_d13_mf1_dist_female,row=3, col=2)
    fig.append_trace(trace_d13_mf1_kde_female,row=3, col=2)
    # fig.append_trace(trace_d13_mf1_dist_male,row=3, col=4)
    # fig.append_trace(trace_d13_mf1_kde_male,row=3, col=4)

    fig.append_trace(trace_d13_female_mf2,row=4, col=1)
    # fig.append_trace(trace_d13_male_mf2,row=4, col=2)
    fig.append_trace(trace_d13_mf2_dist_female,row=4, col=2)
    fig.append_trace(trace_d13_mf2_kde_female,row=4, col=2)
    # fig.append_trace(trace_d13_mf2_dist_male,row=4, col=4)
    # fig.append_trace(trace_d13_mf2_kde_male,row=4, col=4)

    # Update xaxis properties
    # fig.update_yaxes(title_text="relative error",row=1, col=1,range=[0,0.8])
    # for pool6
    fig.update_yaxes(title_text="relative error",row=1, col=1,range=[0,1.25])
    fig.update_yaxes(title_text="relative error ",row=2, col=1,range=[0,0.8])
    fig.update_yaxes(title_text="relative error "    ,row=3, col=1,range=[0,0.8])
    fig.update_yaxes(title_text="relative error "    ,row=4, col=1,range=[0,0.8])
    # fig.update_yaxes(title_text="relative error",row=1, col=2,range=[0,0.8])
    # fig.update_yaxes(title_text="relative error ",row=2, col=2,range=[0,0.8])
    # fig.update_yaxes(title_text="relative error "    ,row=3, col=2,range=[0,0.8])
    # fig.update_yaxes(title_text="relative error "    ,row=4, col=2,range=[0,0.8])

    fig.update_yaxes(title_text="probability"    ,row=1, col=2)
    fig.update_yaxes(title_text="probability"    ,row=2, col=2)
    fig.update_yaxes(title_text="probability"    ,row=3, col=2)
    fig.update_yaxes(title_text="probability"    ,row=4, col=2)
    # fig.update_yaxes(title_text="probability"    ,row=1, col=4)
    # fig.update_yaxes(title_text="probability"    ,row=2, col=4)
    # fig.update_yaxes(title_text="probability"    ,row=3, col=4)
    # fig.update_yaxes(title_text="probability"    ,row=4, col=4)
    # Update yaxis properties
    fig.update_xaxes(title_text="log2 (relative abundance or count)", row=1, col=1,range=[-18,0])
    fig.update_xaxes(title_text="log2 (relative abundance or count)", row=2, col=1,range=[-18,0])
    fig.update_xaxes(title_text="log2 (relative abundance or count)", row=3, col=1,range=[-18,0])
    fig.update_xaxes(title_text="log2 (relative abundance or count)", row=4, col=1,range=[-18,0])
    fig.update_xaxes(title_text="log2 (relative abundance or count)", row=1, col=2,range=[-18,0])
    fig.update_xaxes(title_text="log2 (relative abundance or count)", row=2, col=2,range=[-18,0])
    fig.update_xaxes(title_text="log2 (relative abundance or count)", row=3, col=2,range=[-18,0])
    fig.update_xaxes(title_text="log2 (relative abundance or count)", row=4, col=2,range=[-18,0])

    # fig.update_xaxes(title_text="log2 (relative abundance or count)", row=1, col=3,range=[-18,0])
    # fig.update_xaxes(title_text="log2 (relative abundance or count)", row=2, col=3,range=[-18,0])
    # fig.update_xaxes(title_text="log2 (relative abundance or count)", row=3, col=3,range=[-18,0])
    # fig.update_xaxes(title_text="log2 (relative abundance or count)", row=4, col=3,range=[-18,0])
    # fig.update_xaxes(title_text="log2 (relative abundance or count)", row=1, col=4,range=[-18,0])
    # fig.update_xaxes(title_text="log2 (relative abundance or count)", row=2, col=4,range=[-18,0])
    # fig.update_xaxes(title_text="log2 (relative abundance or count)", row=3, col=4,range=[-18,0])
    # fig.update_xaxes(title_text="log2 (relative abundance or count)", row=4, col=4,range=[-18,0])


    fig.update_layout(height=1000, width=1000, title_text="Diffrent steps of input/output error analysis")
    fig.show()
