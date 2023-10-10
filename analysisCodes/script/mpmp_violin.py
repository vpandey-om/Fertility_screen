import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotnine as pn
import pandas.api.types as pdtypes
import pickle
# read male and female pathway
# male_df=pd.read_excel('/Users/vpandey/projects/githubs/Fertility_screen_2/preFinals/mpmp/male_mpmp_enrich_CS.xlsx')
# female_df=pd.read_excel('/Users/vpandey/projects/githubs/Fertility_screen_2/preFinals/mpmp/female_mpmp_enrich_CS.xlsx')
## load mpmp data
mpmp_data=pickle.load(open('/Users/vpandey/projects/enrichmentData/pbe/mpmp_pbe.pickle','rb'))
pathway_to_genes=mpmp_data[0]
print(len(pathway_to_genes))
## remove some of the categories
mpmp_remove=pd.read_excel('/Users/vpandey/projects/enrichmentData/pbe/MPMP.xlsx',sheet_name='rCS',header=None)

for k in mpmp_remove[0].to_list():
    if k in pathway_to_genes:
       del pathway_to_genes[k]

print(len(pathway_to_genes))

male_df=pd.read_csv('/Users/vpandey/projects/githubs/Fertility_screen_2/preFinals/male_mpmp_enrich.txt',sep='\t')
female_df=pd.read_csv('/Users/vpandey/projects/githubs/Fertility_screen_2/preFinals/female_mpmp_enrich.txt',sep='\t')

top5male=['Kinetochores power chromosome movements in mitosis',
'Regulation of spindle microtubule dynamics',
'Proteins involved in steps during passage through prophase',
'Putative flagellar proteins']
top5female=['Genes encoding respiratory chain proteins',
'F0F1-ATPase',
'Nuclear genes with mitochondrial signal sequences',
'Mitochondrial TCA cycle',
'Double strand break repair and homologous recombination']

orderPath=top5male+top5female
# get geens
top_male=male_df[male_df['Pathways'].isin(top5male)]
top_female=female_df[female_df['Pathways'].isin(top5female)]
# get fertility screen
pvals_male= top_male['Pvals'].apply('{:.1e}'.format).tolist()
pvals_female= top_female['Pvals'].apply('{:.1e}'.format).tolist()
indices = [3,2,1,0]
indices1 = [4,3,2,1,0]

female_mpmp = [top5female[i] for i in indices1]
female_pval_mpmp = [pvals_female[i] for i in indices1]
male_mpmp = [top5male[i] for i in indices]
male_pval_mpmp = [pvals_male[i] for i in indices]
sexlist=['Male']*len(indices)+['Female']*len(indices1)
orderPath=male_mpmp+female_mpmp
orderPval=male_pval_mpmp+female_pval_mpmp
y=[0]*(len(indices)+len(indices1))
order=[i+1 for i in range(len(indices)+len(indices1))]
custom_pvals = pd.DataFrame({
    'Sex': sexlist,
    'x': orderPath,
    'y': y,
    'order':order,
    'label': orderPval
})


fertility_screen=pd.read_excel('/Users/vpandey/projects/githubs/Fertility_screen/preFinals/Phenotype_call_final_100621.xlsx',sheet_name='data')

#### create dataframe with pathway and rgr and genes

def pathway_gene_rgr(top_male,top_female,mpmp_path,sexdf):
    pathways=[]
    rgrs_male=[]
    rgrs_female=[]
    sex=[]
    pheno_male=[]
    pheno_female=[]
    pheno=[]
    rgrs=[]
    paths=[]
    sexes=[]
    gg=[]
    path_genes=[]

    ## for male
    for ind in  top_male.index:
        path=top_male.loc[ind,'Pathways']
        all_genes=mpmp_path[path]
        genes=top_male.loc[ind,'Genes'].split(',')
        genes=[item.strip() for item in genes]
        tmp=sexdf[sexdf['pbanka_id'].isin(all_genes)]
        for item in tmp.index:
            pathways.append(path)
            rgrs_male.append(tmp.loc[item,'g145480_RGR'])
            rgrs_female.append(tmp.loc[item,'GCKO2_RGR'])
            sex.append('Male')
            pheno_male.append(tmp.loc[item,'g145480_pheno_new'])
            pheno_female.append(tmp.loc[item,'GCKO2_pheno_new'])
            rgrs.append(tmp.loc[item,'g145480_RGR'])
            rgrs.append(tmp.loc[item,'GCKO2_RGR'])
            sexes.append('Male')
            sexes.append('Female')
            pheno.append(tmp.loc[item,'g145480_pheno_new'])
            pheno.append(tmp.loc[item,'GCKO2_pheno_new'])
            paths.append(path)
            paths.append(path)
            gg.append(tmp.loc[item,'pbanka_id'])
            gg.append(tmp.loc[item,'pbanka_id'])
            path_genes.append(tmp.loc[item,'pbanka_id'])

    for ind in  top_female.index:
        path=top_female.loc[ind,'Pathways']
        all_genes=mpmp_path[path]
        genes=top_female.loc[ind,'Genes'].split(',')
        genes=[item.strip() for item in genes]
        tmp=sexdf[sexdf['pbanka_id'].isin(all_genes)]
        for item in tmp.index:
            pathways.append(path)
            rgrs_male.append(tmp.loc[item,'g145480_RGR'])
            rgrs_female.append(tmp.loc[item,'GCKO2_RGR'])
            sex.append('Male')
            pheno_male.append(tmp.loc[item,'g145480_pheno_new'])
            pheno_female.append(tmp.loc[item,'GCKO2_pheno_new'])
            rgrs.append(tmp.loc[item,'g145480_RGR'])
            rgrs.append(tmp.loc[item,'GCKO2_RGR'])
            sexes.append('Male')
            sexes.append('Female')
            pheno.append(tmp.loc[item,'g145480_pheno_new'])
            pheno.append(tmp.loc[item,'GCKO2_pheno_new'])
            paths.append(path)
            paths.append(path)
            gg.append(tmp.loc[item,'pbanka_id'])
            gg.append(tmp.loc[item,'pbanka_id'])
            path_genes.append(tmp.loc[item,'pbanka_id'])

    res_df=pd.DataFrame()
    res_df['Pathway']=pathways
    res_df['Fertility_rate_male']=rgrs_male
    res_df['Fertility_rate_female']=rgrs_female
    res_df['Sex']=sex
    res_df['Phenotype_male']=pheno_male
    res_df['Phenotype_female']=pheno_female
    res_df_all=pd.DataFrame()
    res_df_all['Paths']=paths
    res_df_all['Sexes']=sexes
    res_df_all['Phenotype']=pheno
    res_df_all['Fertility rate']=rgrs
    res_df_all['genes']=gg
    return res_df_all

def plot_violin(df,filepdf):
    ''' plot violin all'''
    color_gam_dict={'Reduced':'#ca0020','Not reduced':'#0571b0'}
    df['Pathway'] = df['Paths'].astype(pdtypes.CategoricalDtype(categories=orderPath))
    df['Phenotype'] = df['Phenotype'].astype(pdtypes.CategoricalDtype(categories=['Reduced','Not reduced']))
    df['Sex'] = df['Sexes'].astype(pdtypes.CategoricalDtype(categories=['Female','Male']))

    plot=(pn.ggplot(df, pn.aes(y='Fertility rate', x="Pathway"))
     + pn.coord_flip()
     # + pn.geom_violin(df,style='full',position='dodge')
     + pn.geom_violin()
     + pn.geom_point(pn.aes(color='Phenotype'),size=0.75)
     + pn.scale_color_manual(values=color_gam_dict)
     # Add custom annotations to specific facets
     + pn.geom_text( pn.aes(x='order+0.3', y='y', label='label'), data=custom_pvals,size=6, color='red')

     # + pn.scale_fill_manual(values=color_gam_dict)
     # + pn.annotate('text', y=max(df['Fertility rate'])+0.5, x=df["Pathway"], label='NAAAAAA')
     + pn.facet_wrap('Sex')
     + pn.theme(
            # figure_size=(11, 4.8), ### 4.8
            # legend_direction="vertical",
            # legend_box_spacing=0.4,
            legend_position='none',
            axis_line=pn.element_line(size=1, colour="black"),
            # panel_grid_major=pn.element_line(colour="#d3d3d3"),
            panel_grid_major=pn.element_blank(),
            panel_grid_minor=pn.element_blank(),
            panel_border=pn.element_blank(),
            panel_background=pn.element_blank(),
            # plot_title=pn.element_text(size=15, family="Arial",
            #                         face="bold"),
            text=pn.element_text(family="Arial", size=11),
            axis_text_x=pn.element_text(colour="black", size=10),
            axis_text_y=pn.element_text(colour="black", size=10),
        )
    )


    # Set the figure size


    plot.save(filepdf, dpi=300)





res_df_all=pathway_gene_rgr(top_male,top_female,pathway_to_genes,fertility_screen)
### when we plot male fertility remove No data
female_male_plot_df=res_df_all[~(res_df_all['Phenotype']=='No data')]
# female_plot_df=res_df_all[~(res_df_all['Phenotype_female']=='No data')]
pdffile='/Users/vpandey/projects/githubs/Fertility_screen_2/analysisCodes/Figures/mpmp_image_2.pdf'
plot_violin(female_male_plot_df,pdffile)
import pdb; pdb.set_trace()


## for male
for ind in  top_male.index:
    path=top_male.loc[ind,'Pathways']
    genes=top_male.loc[ind,'Genes'].split(',')
    genes=[item.strip() for item in genes]
    tmp=fertility_screen[fertility_screen['pbanka_id'].isin(genes)]
    for item in tmp.index:
        pathways.append(path)
        total_genes=pathway_to_genes[path]
        rgrs_male.append(tmp.loc[item,'g145480_RGR'])
        rgrs_female.append(tmp.loc[item,'GCKO2_RGR'])
        sex.append('Male')
        pheno_male.append(tmp.loc[item,'g145480_pheno_new'])
        pheno_female.append(tmp.loc[item,'GCKO2_pheno_new'])
        rgrs.append(tmp.loc[item,'g145480_RGR'])
        rgrs.append(tmp.loc[item,'GCKO2_RGR'])
        sexes.append('Male')
        sexes.append('Female')
        pheno.append(tmp.loc[item,'g145480_pheno_new'])
        pheno.append(tmp.loc[item,'GCKO2_pheno_new'])
        paths.append(path)
        paths.append(path)
        gg.append(tmp.loc[item,'pbanka_id'])
        gg.append(tmp.loc[item,'pbanka_id'])
        path_genes.append(tmp.loc[item,'pbanka_id'])
# For female
for ind in  top_female.index:
    path=top_female.loc[ind,'Pathways']
    genes=top_female.loc[ind,'Genes'].split(',')
    genes=[item.strip() for item in genes]
    tmp=fertility_screen[fertility_screen['pbanka_id'].isin(genes)]
    for item in tmp.index:
        pathways.append(path)
        rgrs_female.append(tmp.loc[item,'GCKO2_RGR'])
        rgrs_male.append(tmp.loc[item,'g145480_RGR'])
        sex.append('Female')
        pheno_female.append(tmp.loc[item,'GCKO2_pheno_new'])
        pheno_male.append(tmp.loc[item,'g145480_pheno_new'])
        rgrs.append(tmp.loc[item,'g145480_RGR'])
        rgrs.append(tmp.loc[item,'GCKO2_RGR'])
        sexes.append('Male')
        sexes.append('Female')
        pheno.append(tmp.loc[item,'g145480_pheno_new'])
        pheno.append(tmp.loc[item,'GCKO2_pheno_new'])
        paths.append(path)
        paths.append(path)
        gg.append(tmp.loc[item,'pbanka_id'])
        gg.append(tmp.loc[item,'pbanka_id'])
        path_genes.append(tmp.loc[item,'pbanka_id'])


res_df=pd.DataFrame()
res_df['Pathway']=pathways
res_df['Fertility_rate_male']=rgrs_male
res_df['Fertility_rate_female']=rgrs_female
res_df['Sex']=sex
res_df['Phenotype_male']=pheno_male
res_df['Phenotype_female']=pheno_female
res_df_all=pd.DataFrame()
res_df_all['Paths']=paths
res_df_all['Sexes']=sexes
res_df_all['Phenotype']=pheno
res_df_all['Fertility rate']=rgrs
res_df_all['genes']=gg
path_gene_df=pd.DataFrame()
path_gene_df['Pathway']=pathways
path_gene_df['genes']=path_genes
path_gene_df.to_excel('/Users/vpandey/projects/githubs/Fertility_screen_2/preFinals/Pathgene_for_claire.xlsx')


part_male_vio=res_df_all[res_df_all['Paths'].isin(top5male)]
part_female_vio=res_df_all[res_df_all['Paths'].isin(top5female)]

part_male=res_df_all[res_df_all['Sexes'].isin(['Male'])]
part_female=res_df_all[res_df_all['Sexes'].isin(['Female'])]

tmp=res_df_all[res_df_all['Phenotype']=='No data']


unique_gg=set(gg)


sns.set(font="Arial",font_scale =2)
sns.set_style("white")


color_dict={'Male':'#009ADE','Female':'#AF58BA'}
# color_dict={'Male':'#FC8D62','Female':'#66C2A5'}
# color_dict={'Male':'#FFFFFF','Female':'#FFFFFF'}
color_gam_dict={'Reduced':'#FC8D62','Not reduced':'#66C2A5'}
plt.figure()

res_df_all1=res_df_all.copy()

res_df_all1['Pathways'] = res_df_all1['Paths'].astype(pdtypes.CategoricalDtype(categories=orderPath))
res_df_all1['Sex'] = res_df_all1['Sexes'].astype(pdtypes.CategoricalDtype(categories=['Male','Female']))
res_df_all1['Phenotype'] = res_df_all1['Phenotype'].astype(pdtypes.CategoricalDtype(categories=['Reduced','Not reduced']))

# res_df_all1['sextype']=[ color_dict[item] for item in res_df_all1['Sexes'] ]
# res_df_all1['phenotype']=[ color_gam_dict[item] for item in res_df_all1['Pheno'] ]


import pdb;pdb.set_trace()

# plot=(pn.ggplot(res_df_all1, pn.aes(x="Paths",y="Fertility_rate"))
#  + pn.geom_violin(res_df_all1)
# )


# plot=(pn.ggplot(res_df_all1, pn.aes(y='Fertility rate', x="Pathways",color='Phenotype',fill='Sex'))
#  + pn.coord_flip()
#  + pn.geom_violin(res_df_all1,position='dodge',style='full')
#  + pn.geom_point()
#  + pn.scale_color_manual(values=color_gam_dict)
#  + pn.scale_fill_manual(values=color_dict)
# )


plot=(pn.ggplot(res_df_all1, pn.aes(y='Fertility rate', x="Pathways",fill='Sex'))
 + pn.coord_flip()
 + pn.geom_violin(res_df_all1,style='full',position='dodge')
 + pn.geom_point()
 + pn.scale_color_manual(values=color_gam_dict)
 + pn.scale_fill_manual(values=color_dict)
 + pn.theme(
        legend_direction="vertical",
        legend_box_spacing=0.4,
        axis_line=pn.element_line(size=1, colour="black"),
        panel_grid_major=pn.element_line(colour="#d3d3d3"),
        panel_grid_minor=pn.element_blank(),
        panel_border=pn.element_blank(),
        panel_background=pn.element_blank(),
        plot_title=pn.element_text(size=15, family="Tahoma",
                                face="bold"),
        text=pn.element_text(family="Arial", size=11),
        axis_text_x=pn.element_text(colour="black", size=10),
        axis_text_y=pn.element_text(colour="black", size=10),
    )
)



plot.save('/Users/vpandey/projects/githubs/Fertility_screen_2/analysisCodes/Figures/image.pdf', dpi=300)

import pdb; pdb.set_trace()
# sns.violinplot(x="Fertility_rate_male", y="Pathway", data=res_df,color='#FFFFFF') ## 'quartile' 'box'
#
#
#
# sns.stripplot(x="Fertility_rate_male",
#                 y="Pathway",
#                 data=res_df,
#                    hue='Phenotype_male',edgecolor="gray",size=8,palette=color_gam_dict)

# sns.violinplot(x="Fertility_rate", y="Paths", data=part_male_vio,palette=color_dict,hue='Sexes',order=orderPath,scale="width",dodge=True) ## 'quartile' 'box'
# sns.violinplot(x="Fertility_rate", y="Paths", data=part_female_vio,palette=color_dict,hue='Sexes',hue_order=['Female','Male'],order=orderPath,scale="width",dodge=True) ## 'quartile' 'box

sns.violinplot(x="Fertility_rate", y="Paths", data=res_df_all,palette=color_dict,hue='Sexes',order=orderPath,scale="width",dodge=True) ## 'quartile' 'box'

# sns.violinplot(x="Fertility_rate", y="Paths", data=res_df_all,palette=color_dict,hue='Sexes',order=orderPath,scale="width",dodge=True) ## 'quartile' 'box'

sns.stripplot(x="Fertility_rate",
                y="Paths",
                data=res_df_all,
                   hue='Pheno',edgecolor="gray",size=8,palette=color_gam_dict,hue_order=['Reduced','Not reduced'],order=orderPath,dodge=False)

# sns.stripplot(x="Fertility_rate",
#                 y="Paths",
#                 data=part_male_vio,
#                    hue='Pheno',edgecolor="gray",size=8,palette=color_gam_dict,hue_order=['Reduced','Not reduced'],order=orderPath,dodge=True)
# sns.stripplot(x="Fertility_rate",
#                 y="Paths",
#                 data=part_female_vio,
#                    hue='Pheno',edgecolor="gray",size=8,palette=color_gam_dict,hue_order=['Not reduced','Reduced'],order=orderPath,dodge=True)
#
#

# sns.violinplot(x="Fertility_rate", y="Paths", data=part_male,color='#FFFFFF',order=orderPath,scale="width") ## 'quartile' 'box'
#
#
#
# sns.stripplot(x="Fertility_rate",
#                 y="Paths",
#                 data=part_male,
#                    hue='Pheno',edgecolor="gray",size=8,palette=color_gam_dict,order=orderPath)


plt.yticks([])

plt.legend().remove()
# plt.tick_params(axis='x', rotation=0)
# plt.xlabel('Gametocyte', fontsize=20,labelpad=10)
plt.xlabel('')
plt.ylabel('')
plt.xlim((-18,5))
import pdb; pdb.set_trace()

plt.savefig('/Users/vpandey/projects/githubs/Fertility_screen_2/analysisCodes/Figures/mpmp_violin1.pdf',
            format='pdf',dpi=300)

plt.close()

# sns.violinplot(x="Fertility_rate", y="Paths", data=part_female,color='#FFFFFF',order=orderPath,scale="width") ## 'quartile' 'box'
#
#
#
# sns.stripplot(x="Fertility_rate",
#                 y="Paths",
#                 data=part_female,
#                    hue='Pheno',edgecolor="gray",size=8,palette=color_gam_dict,order=orderPath)
#
#
#
# plt.legend().remove()
# # plt.tick_params(axis='x', rotation=0)
# # plt.xlabel('Gametocyte', fontsize=20,labelpad=10)
# plt.xlabel('')
# plt.ylabel('')
# plt.xlim((-18,8))
#
# plt.savefig('/Users/vpandey/projects/githubs/Fertility_screen_2/analysisCodes/Figures/mpmp_female_violin.pdf',
#             format='pdf',dpi=300)
#
#
#
# plt.close()
