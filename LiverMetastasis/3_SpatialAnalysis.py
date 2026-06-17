import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import scanpy as sc
import SOAPy_st as sp

sp.pp.make_network(adata,sample_key='dataset',method='knn',
                   cutoff=20,scale=1.0,cluster_key='cluster_L3')

### CCI
sp.tl.neighborhood_analysis(adata,sample_key='dataset',method='include',
                            cluster_key='cluster_L3',n_iters=100,n_jobs=30)


### niche identification
sp.tl.get_c_niche(adata,k_max=30,celltype_key='cluster_L3',sample_key='dataset')

neig_df=adata.uns['SOAPy']['niche']
neig_df=neig_df[np.unique(adata.obs['cluster_L3'])]
cluster_num=10
km=KMeans(n_clusters=cluster_num,max_iter=300)
clusters_res=km.fit_predict(neig_df)
cluster_means=neig_df.groupby(clusters_res).mean()
cluster_means=cluster_means.apply(stats.zscore,axis=0)
