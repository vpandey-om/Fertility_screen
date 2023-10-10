
## get all phenotype data 
from sqlalchemy import create_engine
import pandas as pd
import numpy as np
# from sklearn.datasets import make_blobs
# from matplotlib import pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# from sklearn.decomposition import PCA

# from sklearn.preprocessing import StandardScaler
# from sklearn import cluster
# from sklearn.neighbors import NearestNeighbors,kneighbors_graph
# import networkx as nx 
# import json
# import py4cytoscape as p4c
# import umap
# reducer = umap.UMAP()
# from fa2 import ForceAtlas2

# import pickle
# import os
# cur_dir = os.path.dirname(__file__)


 
##databse connection
POSTGRES = {
    'user': 'vikash',
    'pw': 'om16042020',
    'db': 'phenodb',
    'host': 'localhost',
    'port': '5432',
}

DATABASE_URI='postgresql+psycopg2://%(user)s:%(pw)s@%(host)s:%(port)s/%(db)s' % POSTGRES
engine = create_engine(DATABASE_URI)

phenodata_df = pd.read_sql_query('''select * from phenotable''',con=engine)

genedf = pd.read_sql_query('''select * from geneanot''',con=engine)
genedf.to_csv('gene_des.txt',sep='\t',index=None)


#### read all phenotype data



from sklearn.cluster import KMeans

def kmeans_missing(X, n_clusters, max_iter=10):
    """Perform K-Means clustering on data with missing values.

    Args:
      X: An [n_samples, n_features] array of data to cluster.
      n_clusters: Number of clusters to form.
      max_iter: Maximum number of EM iterations to perform.

    Returns:
      labels: An [n_samples] vector of integer labels.
      centroids: An [n_clusters, n_features] array of cluster centroids.
      X_hat: Copy of X with the missing values filled in.
    """

    # Initialize missing values to their column means
    missing = ~np.isfinite(X)
    mu = np.nanmean(X, 0, keepdims=1)
    X_hat = np.where(missing, mu, X)

    for i in range(max_iter):
        if i > 0:
            # initialize KMeans with the previous set of centroids. this is much
            # faster and makes it easier to check convergence (since labels
            # won't be permuted on every iteration), but might be more prone to
            # getting stuck in local minima.
            cls = KMeans(n_clusters, init=prev_centroids)
        else:
            # do multiple random initializations in parallel
            cls = KMeans(n_clusters)

        # perform clustering on the filled-in data
        labels = cls.fit_predict(X_hat)
        centroids = cls.cluster_centers_

        # fill in the missing values based on their cluster centroids
        X_hat[missing] = centroids[labels][missing]

        # when the labels have stopped changing then we have converged
        if i > 0 and np.all(labels == prev_labels):
            break

        prev_labels = labels
        prev_centroids = cls.cluster_centers_

    return labels, centroids, X_hat


pheno_col_rename={'B2toMG_RGR':'Oocyst conversion rate', 'B2toMG_SD':'Oocyst SD',
        'B2toMG_pheno': 'Oocyst phenotype', 'MGtoSG_RGR':'Sporozoite conversion rate',
        'MGtoSG_SD':'Sporozoite SD','MGtoSG_pheno':'Sporozoite phenotype', 'SGtoB2_RGR':'Liver conversion rate',
        'SGtoB2_SD':'Liver SD', 'SGtoB2_pheno': 'Liver phenotype', 'blood_RGR':'Blood RGR',
        'blood_SD':'Blood SD', 'blood_pheno':'Blood phenotype', 'female_fertility_RGR': 'Female fertility rate',
        'female_fertility_SD' :'Female fertility SD' , 'female_fertility_pheno':'Female fertility phenotype',
        'male_fertility_RGR' :'Male fertility rate','male_fertility_SD':'Male fertility SD',
        'male_fertility_pheno':'Male fertility phenotype', 'female_gam_RGR':'Female gametocyte conversion rate',
        'female_gam_SD':'Female gametocyte SD', 'female_gam_pheno':'Female gametocyte phenotype',
        'male_gam_RGR':'Male gametocyte conversion rate', 'male_gam_SD':'Male gametocyte SD',
        'male_gam_pheno':'Male gametocyte phenotype','new_pbanka_id':'Gene ID'}

arranged_pheno_col=['Gene ID','Blood RGR','Blood SD','Blood phenotype','Female gametocyte conversion rate',
        'Female gametocyte SD','Female gametocyte phenotype','Male gametocyte conversion rate',
        'Male gametocyte SD','Male gametocyte phenotype','Female fertility rate','Female fertility SD',
        'Female fertility phenotype','Male fertility rate','Male fertility SD','Male fertility phenotype',
        'Oocyst conversion rate','Oocyst SD','Oocyst phenotype','Sporozoite conversion rate','Sporozoite SD',
        'Sporozoite phenotype','Liver conversion rate','Liver SD','Liver phenotype']

stage_pheno=['Blood','Female Gametocyte','Male Gametocyte','Female fertility',
        'Male Fertility','Oocyst','Sporozoite','Liver']

pheno_rate_columns=['Blood RGR','Female gametocyte conversion rate','Male gametocyte conversion rate',
        'Female fertility rate','Male fertility rate','Oocyst conversion rate','Sporozoite conversion rate',
        'Liver conversion rate']
pheno_rate_min=[0.003293,-5.49574,-4.487991,-12.408198,-11.598957,-11.944122,-8.06721,-13.715141]
pheno_rate_max=[1.246831,2.453507,1.552065,2.943137,2.103551, 5.634955,5.183004,12.836909]
pheno_sd_columns=['Blood SD','Female gametocyte SD','Male gametocyte SD','Female fertility SD',
        'Male fertility SD','Oocyst SD','Sporozoite SD','Liver SD']

pheno_call_columns=['Blood phenotype','Female gametocyte phenotype','Male gametocyte phenotype',
        'Female fertility phenotype','Male fertility phenotype','Oocyst phenotype','Sporozoite phenotype',
        'Liver phenotype']

pheno_df=pd.read_csv('all_pheno_data.txt',sep='\t')
motility_df=pd.read_excel('/Users/vpandey/projects/githubs/Fertility_screen_2/preFinals/11.02.21_male_pool_dis_RGR.xlsx')


df=pheno_df.copy()
df=df.rename(columns=pheno_col_rename)
# df=df.fillna('')
### remove non essential genes
# phenotypes='Dispensable', 'Slow', 'Essential', 'Insufficient data', nan,
#        'Fast']
m_arr=df['Blood RGR'].values
s_arr=df['Blood SD'].values
data= np.empty((len(m_arr),100))
data[:] = np.nan
for j in range(len(m_arr)):
    if np.isnan(m_arr[j]) and np.isnan(s_arr[j]):
        tmp=np.empty((1,100,))
        tmp[:] = np.nan
    else:
        tmp = np.random.normal(m_arr[j], s_arr[j], 100)  
    data[j,:]=tmp

import pdb;pdb.set_trace()
phenotypes=['Dispensable', 'Slow', 'Insufficient data','Fast']
df_red=df[(df['Blood phenotype'].isin(phenotypes))&(~(df['Blood RGR'].isna()))]
df_gene_data=df_red[['Gene ID']+pheno_rate_columns+pheno_sd_columns].copy()
df_gene_data.set_index('Gene ID',inplace=True)
df_gene_data.to_csv('genes_phenotype_clutsering.txt',sep='\t')
### count is na column wise 
xx=df_gene_data.isna().sum(axis=1)
remove_genes=xx[xx>10]
final_df_pre=df_gene_data.drop(remove_genes.index.to_list())
final_df=final_df_pre.interpolate()
final_df=final_df[~final_df['Oocyst conversion rate'].isna()]
print(final_df.isna().sum(axis=0))

#### now we can apply clustering on data sets



# plot the inferred points, color-coded according to the true cluster labels



def make_fake_data(fraction_missing, n_clusters=5, n_samples=1500,
                   n_features=3, seed=None):
    # complete data
    gen = np.random.RandomState(seed)
    # X, true_labels = make_blobs(n_samples, n_features, n_clusters,
    #                             random_state=gen)
    X, true_labels=make_blobs(n_samples=1500, centers=5, n_features=3,  random_state=0)
    # with missing values
    missing = gen.rand(*X.shape) < fraction_missing
    Xm = np.where(missing, np.nan, X)
    return X, true_labels, Xm


# X, true_labels, Xm = make_fake_data(fraction_missing=0.3, n_clusters=5, n_samples=1500,seed=0)
Xm=df_gene_data.values
labels, centroids, X_hat = kmeans_missing(Xm, n_clusters=20)
X = StandardScaler().fit_transform(X_hat)
  # normalize dataset for easier parameter selection
# connectivity matrix for structured Ward





spectral = cluster.SpectralClustering(
        n_clusters=10,
        eigen_solver="arpack",
        affinity="nearest_neighbors",
    )
clustering=spectral.fit(X)

nbrs = NearestNeighbors(n_neighbors=5, algorithm='auto').fit(X)

distances, indices = nbrs.kneighbors(X)

G = nx.from_numpy_array(nbrs.kneighbors_graph(X).toarray())
# forceatlas2 = ForceAtlas2(
#                         # Behavior alternatives
#                         outboundAttractionDistribution=True,  # Dissuade hubs
#                         linLogMode=False,  # NOT IMPLEMENTED
#                         adjustSizes=False,  # Prevent overlap (NOT IMPLEMENTED)
#                         edgeWeightInfluence=1.0,

#                         # Performance
#                         jitterTolerance=1.0,  # Tolerance
#                         barnesHutOptimize=True,
#                         barnesHutTheta=1.2,
#                         multiThreaded=False,  # NOT IMPLEMENTED

#                         # Tuning
#                         scalingRatio=2.0,
#                         strongGravityMode=False,
#                         gravity=1.0,

#                         # Log
#                         verbose=True)

# positions = forceatlas2.forceatlas2_networkx_layout(G, pos=None, iterations=2000)
# nx.draw_networkx_nodes(G, positions, node_size=20, with_labels=False, node_color="blue", alpha=0.4)
# nx.draw_networkx_edges(G, positions, edge_color="green", alpha=0.05)
# plt.axis('off')
# plt.show()
p4c.create_network_from_networkx(G)



netcyjs=nx.cytoscape_data(G)  




# Serializing json
json_object = json.dumps(netcyjs, indent=4)
 
# Writing to sample.json
with open("sample.json", "w") as outfile:
    outfile.write(json_object )



###




### apply PCA 
pca = PCA(n_components=3)
pca.fit(X_hat)
scoreX=pca.transform(X_hat)
print(pca.explained_variance_ratio_)


###
# fig=plt.figure()
# ax = fig.add_subplot()
# # ax[0].scatter3D(X[:, 0], X[:, 1], X[:, 2], c=true_labels, cmap='gist_rainbow')
# ax.scatter(scoreX[:, 0], scoreX[:, 1], c=labels,
#                 cmap='gist_rainbow')
# ax.set_title('5 clusters')

# fig.tight_layout()
# plt.show()


###

###
# plot the inferred points, color-coded according to the true cluster labels
# fig, ax = plt.subplots(1, 1, subplot_kw={'projection':'3d', 'aspect':'equal'})
fig=plt.figure()
ax = fig.add_subplot(projection='3d')
# ax[0].scatter3D(X[:, 0], X[:, 1], X[:, 2], c=true_labels, cmap='gist_rainbow')
ax.scatter(scoreX[:, 0], scoreX[:, 1], scoreX[:, 2], c=clustering.labels_,
                cmap='gist_rainbow')
ax.set_title('5 clusters')

fig.tight_layout()
plt.show()
### apply pca














