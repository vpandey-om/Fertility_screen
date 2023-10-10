import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle
from pylab import *
# plt.style.use('seaborn-whitegrid')
cmap = cm.get_cmap('seismic', 5)
from matplotlib.font_manager import FontProperties
import matplotlib as mpl

# Set the default font family to Arial
mpl.rcParams['font.family'] = 'Arial'
## import final xlsx file then plot s curve
fertility_df=pd.read_excel('/Users/vpandey/projects/githubs/Fertility_screen/preFinals/Phenotype_call_final_100621.xlsx',sheet_name='data')

## choice of color
# arial_font = FontProperties(family='Arial', size=1)

# get male dataframe
female_columns=['GCKO2_RGR','GCKO2_sd','GCKO2_pheno']


def plotScurve(df,ylab,pheno_names=None,pheno_colors=None,pheno_legends=None,plot_df='test.pdf'):
    ## sort dataframe
    if not pheno_names:
        pheno_names=df['pheno'].unique()
        pheno_legends=pheno_names
        pheno_labels=pheno_names
        pheno_colors=[]
        for i in range(len(pheno_names)):
            rgba = cmap(i)
            pheno_colors.append(matplotlib.colors.rgb2hex(rgba))
    sorted_df=df.sort_values(by=['rgr']).copy()
    sorted_df['x']=np.arange(sorted_df.shape[0])
    # get for each phenoype and plot
    fig, ax = plt.subplots(figsize=[8, 2.75])
    # ax = subplot(1,1,1)

    for i,item in enumerate(pheno_names):
        p_df=sorted_df[sorted_df['pheno']==item].copy()
        # h1=ax.errorbar(p_df['x'].values, p_df['rgr'].values, yerr=p_df['sd'].values*2, fmt='o', color=pheno_colors[i],
        #          ecolor=pheno_colors[i], elinewidth=1, capsize=3,label=pheno_legends[i],ms=2)
        h1=ax.errorbar(p_df['x'].values, p_df['rgr'].values, yerr=p_df['sd'].values*2, fmt='o', color=pheno_colors[i],
                 ecolor=pheno_colors[i], elinewidth=1, capsize=3,ms=2)

    # handles, labels = ax.get_legend_handles_labels()
    # ax.legend(handles[::-1], labels[::-1])
    # plt.grid(axis='y', color='0.95')
    # plt.xlabel('Ranked genes [1-%d]'%sorted_df.shape[0], fontsize=16)
    fsize=18
    # plt.xlabel('Ranked genes', fontsize=fsize)
    # plt.ylabel(ylab, fontsize=fsize)
    plt.xticks(fontsize=fsize)
    plt.yticks(fontsize=fsize)
    plt.ylim(-18, 8)
    # plt.xlim(-14, 1150)
    ## reoder legends

    plt.savefig(plot_df, dpi=300)
    plt.close()


def plot_half_circle(df,plot_df='test.df'):
    ''' we would like to fill half circles  '''
    ## get male female redcued
    fig, ax = plt.subplots(figsize=[6.4, 6.4])
    #"#f7f7f7"
    # ax.scatter(df['g145480_RGR'], df['GCKO2_RGR'], c=df['male_color'], s=30,edgecolor='#F5F5F5', alpha=0.6,linewidth=0.2,marker=MarkerStyle("o", fillstyle="right"))
    # ax.scatter(df['g145480_RGR'], df['GCKO2_RGR'], c=df['female_color'], s=30,edgecolor='#F5F5F5', alpha=0.6,linewidth=0.2,marker=MarkerStyle("o", fillstyle="left"))
    #

    for i in df.index:
        marker_style = dict(color=df.loc[i,'female_color'],markersize=6, markerfacecoloralt=df.loc[i,'male_color'],markeredgecolor='#F5F5F5',alpha=0.8,markeredgewidth=0.2)
        ax.plot(df.loc[i,'g145480_RGR'],  df.loc[i,'GCKO2_RGR'], 'o',fillstyle='left', **marker_style)
    fsize=30
    # plt.xlabel('Relative male fertility',fontsize=fsize)
    # plt.ylabel('Relative female fertility',fontsize=fsize)

    ticks=[5, 0, -5, -10, -15]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    plt.xticks(fontsize=fsize)
    plt.yticks(fontsize=fsize)
    # ax.set_xlabel('Male fertility',fontsize=15)
    # ax.set_ylabel('Female fertility',fontsize=15)
    # ax.tick_params(axis='y',labelsize=15)
    # ax.tick_params(axis='x',labelsize=15)
    # ax.set_xlim(-13,3)
    # ax.set_ylim(-13,3)
    plt.xlim(-13,3)
    plt.ylim(-13,3)
    plt.savefig(plot_df, dpi=300)
    plt.close()
    #plt.show()






## preprocessing
male_columns=['g145480_RGR','g145480_sd','g145480_pheno_new']
male_df=fertility_df[male_columns].copy()
male_df=male_df.rename(columns = {'g145480_RGR': 'rgr','g145480_sd':'sd','g145480_pheno_new':'pheno'})


# remove columns with no data
male_df=male_df[~(male_df['pheno']=='No data')].copy()

#plotScurve(male_df) ## bydefault

# #pheno_names=['Reduced','Not Reduced','No Power']
# # pheno_names=['No Power','Not Reduced','Reduced',]
# # pheno_colors=['#cccccc','#b2e2e2','#d8b365',]
# # pheno_legends=['No power','Not reduced','Reduced']
# pheno_labels=['Reduced','Not reduced','No power']

# pheno_names=['Reduced','Not reduced']
# pheno_colors=['#FC8D62','#66C2A5']
# pheno_legends=['Reduced','Not reduced']

pheno_names=['Not reduced','Reduced',]
pheno_colors=['#66C2A5','#FC8D62']
pheno_colors=['#0571b0','#ca0020']
pheno_legends=['Not reduced','Reduced']
# color_gam_dict={'Reduced':'#ca0020','Not reduced':'#0571b0'}
plot_df='/Users/vpandey/projects/githubs/Fertility_screen/preFinals/male_S_2.pdf'

plotScurve(male_df,'Relative male fertility',pheno_names,pheno_colors,pheno_legends,plot_df)


### for female fertility

female_columns=['GCKO2_RGR','GCKO2_sd','GCKO2_pheno_new']
female_df=fertility_df[female_columns].copy()
female_df=female_df.rename(columns = {'GCKO2_RGR': 'rgr','GCKO2_sd':'sd','GCKO2_pheno_new':'pheno'})

# remove columns with no data
female_df=female_df[~(female_df['pheno']=='No data')].copy()

plot_df='/Users/vpandey/projects/githubs/Fertility_screen/preFinals/female_S_2.pdf'

plotScurve(female_df,'Relative female fertility',pheno_names,pheno_colors,pheno_legends,plot_df)

## get common between male and female
common_df=fertility_df[~((fertility_df['g145480_pheno_new']=='No data') | (fertility_df['GCKO2_pheno_new']=='No data'))]

### plot half circle
## colors
common_df['male_color']=common_df['g145480_pheno_new'] ## Not reduced
common_df['female_color']=common_df['GCKO2_pheno_new']  ## Not reduced

common_df=common_df.replace({'male_color': dict(zip(pheno_names, pheno_colors)),
'female_color': dict(zip(pheno_names, pheno_colors))})
plot_df='/Users/vpandey/projects/githubs/Fertility_screen/preFinals/male_female_half_circle_2.pdf'
plot_half_circle(common_df,plot_df)

import pdb; pdb.set_trace()



##### we

# def plotMaleFemaleScatter(df):
#     # fig = plt.figure(figsize=(8,8))
#     fig, ax = plt.subplots(figsize=(9, 9))
#
#     ###
#     # change df_values
#
#     df['145480_d13_pheno']=df['145480_d13_pheno'].replace({'E': 'IM', 'NE': 'FM','NA':'RM'})
#     df['GCKO2_d13_pheno']=df['GCKO2_d13_pheno'].replace({'E': 'IF', 'NE': 'FF', 'NA': 'RF'})
#     # viz_df['Published_cross_phenotype']=viz_df['Published_cross_phenotype'].replace({'N': 'NA'})
#
#     # cmap={'FM':"#66c2a5", 'RM':"#8da0cb", 'IM':"#fc8d62"}
#     cmap_male={'FM':"#1b9e77", 'RM':"#7570b3", 'IM':"#d95f02"}
#     cmap_female={'FF':"#1b9e77", 'RF':"#7570b3", 'IF':"#d95f02"}
#
#
#     for i,item in enumerate(df['GCKO2_d13_pheno'].to_list()):
#         marker_style = dict(color=cmap_female[item],markersize=6, markerfacecoloralt=cmap_male[df['145480_d13_pheno'][i]],markeredgecolor='white',alpha=0.8)
#         ax.plot(df['GCKO2_d13_rel'][i],df['145480_d13_rel'][i], 'o',fillstyle='left', **marker_style)
#
#
#     ax.set_xlabel('Relative growth rate (Female)',fontsize=15)
#     ax.set_ylabel('Relative growth rate (Male)',fontsize=15)
#     ax.tick_params(axis='y',labelsize=12)
#     ax.tick_params(axis='x',labelsize=12)
#     ax.set_xlim(-11,3)
#     ax.set_ylim(-11,3)
#     # handles, labels = ax.get_legend_handles_labels()
#     # ax.legend(handles, labels)
#     infertile= mpatches.Patch(color="#fc8d62", label='Infertility (F/M)')
#     fertile=mpatches.Patch(color="#66c2a5", label='Normal fertility (F/M)')
#     reduced_fertile=mpatches.Patch(color="#8da0cb", label='Reduced fertility (F/M)')
#     plt.legend(handles=[infertile,fertile,reduced_fertile],loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=3, borderaxespad=0, frameon=False, prop={'size': 6})
#
#     fig.savefig(plot_folder + "scatter_plot_male_female_RGR_pool1.pdf")
#
#

####
