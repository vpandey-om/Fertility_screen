
malaria_atlas_list= pickle.load(open('/Users/vikash/git-hub/pbeDB/ipbedb/external/malaria_cell_atlas.pickle'), 'rb'))

def getfigMalariaAtlasGene(gene):

    ## this data is for malaria cell atlas visualizations
    umap_df=malaria_atlas_list[0] # umap_df
    pheno_df=malaria_atlas_list[1] ## pheno_df
    atlas_sparse_data=malaria_atlas_list[2] ### sparse numpy array
    atlas_genes=malaria_atlas_list[3]
    ### find index for a gene
    umap_df=umap_df.replace({'Male':'Male_gametocyte'})
    if len(atlas_genes[atlas_genes==gene].to_list())>0:
        umap_df['expression']=np.nan
        arr=atlas_sparse_data.todense()
        umap_df['expression']=arr[atlas_genes==gene,:].tolist()[0]
        fig = px.scatter(umap_df, x='umap0', y='umap1', color='expression')
        fig.update_traces(marker=dict(size=4))
        fig.update_layout(title=gene,width=500, height=400,showlegend=True)

    else:
        umap_df['expression']=1
        fig = px.scatter(umap_df, x='umap0', y='umap1', color='expression')
        fig.update_traces(marker=dict(size=4))
        fig.update_layout(title='Not found',width=500, height=400,showlegend=True)
