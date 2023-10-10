import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle
from pylab import *
# plt.style.use('seaborn-whitegrid')
cmap = cm.get_cmap('seismic', 5)

## import final xlsx file then plot s curve
fertility_df=pd.read_excel('/Users/vpandey/projects/githubs/Fertility_screen/preFinals/Phenotype_call_final_100621.xlsx',sheet_name='data')

## choice of color


# get male dataframe
female_columns=['GCKO2_RGR','GCKO2_sd','GCKO2_pheno']


def plotScurve(df,pheno_names=None,pheno_colors=None,pheno_legends=None,plot_df='test.pdf'):
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
    ax = subplot(1,1,1)

    for i,item in enumerate(pheno_names):
        p_df=sorted_df[sorted_df['pheno']==item].copy()
        h1=ax.errorbar(p_df['x'].values, p_df['rgr'].values, yerr=p_df['sd'].values*2, fmt='o', color=pheno_colors[i],
                 ecolor=pheno_colors[i], elinewidth=1, capsize=3,label=pheno_legends[i],ms=2)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1])
    plt.grid(axis='y', color='0.95')
    plt.xlabel('Ranked genes [1-%d]'%sorted_df.shape[0], fontsize=16)
    plt.ylabel('Motality rate', fontsize=16)
    plt.ylim(-5, 1)
    # plt.xlim(-14, 1150)
    ## reoder legends

    plt.savefig(plot_df, dpi=300)
    plt.close()


def plot_half_circle(df,plot_df='test.df'):
    ''' we would like to fill half circles  '''
    ## get male female redcued
    fig, ax = plt.subplots(figsize=[3, 3])
    #"#f7f7f7"
    # ax.scatter(df['g145480_RGR'], df['GCKO2_RGR'], c=df['male_color'], s=30,edgecolor='#F5F5F5', alpha=0.6,linewidth=0.2,marker=MarkerStyle("o", fillstyle="right"))
    # ax.scatter(df['g145480_RGR'], df['GCKO2_RGR'], c=df['female_color'], s=30,edgecolor='#F5F5F5', alpha=0.6,linewidth=0.2,marker=MarkerStyle("o", fillstyle="left"))
    #

    for i in df.index:
        marker_style = dict(color=df.loc[i,'female_color'],markersize=6, markerfacecoloralt=df.loc[i,'male_color'],markeredgecolor='#F5F5F5',alpha=0.8,markeredgewidth=0.2)
        ax.plot(df.loc[i,'g145480_RGR'],  df.loc[i,'GCKO2_RGR'], 'o',fillstyle='left', **marker_style)

    # plt.xlabel('Male fertility',fontsize=15)
    # plt.ylabel('Female fertility',fontsize=15)
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




## read motality dispensble pool.
motility_df=pd.read_excel('/Users/vpandey/projects/githubs/Fertility_screen_2/preFinals/11.02.21_male_pool_dis_RGR.xlsx')
tmp=motility_df[motility_df['sup_vs_pel_diff_min']>-1]
motility_df['sup_vs_pel_pheno_new']='Reduced'
motility_df.loc[tmp.index,'sup_vs_pel_pheno_new']='Not reduced'

tmpdf=motility_df.copy()
tmpdf=tmpdf.rename(columns = {'sup_vs_pel_RGR': 'rgr','sup_vs_pel_sd':'sd','sup_vs_pel_pheno_new':'pheno'})
pheno_names=['Not reduced','Reduced',]
pheno_colors=['#66C2A5','#FC8D62']
pheno_legends=['Not reduced','Reduced']

plot_df='/Users/vpandey/projects/githubs/Fertility_screen_2/preFinals/Motality_S_Curve.pdf'

plotScurve(tmpdf,pheno_names,pheno_colors,pheno_legends,plot_df)
