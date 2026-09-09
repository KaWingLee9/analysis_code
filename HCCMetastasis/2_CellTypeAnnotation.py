import os
import re

import numpy as np
import pandas as pd
import scipy
import matplotlib.pyplot as plt

import anndata as ad
import scanpy as sc

import matplotlib
matplotlib.rcParams['pdf.fonttype']='TrueType'

# this function is totally copied from https://github.com/scverse/scanpy/issues/2142
# def scale_by_batch(adata,batch_key,
#                    max_value=None) -> ad.AnnData:
#     return ad.concat(
#         {
#             k: sc.pp.scale(adata[idx],copy=True,max_value=10)
#             for k, idx in adata.obs.groupby(batch_key).indices.items()
#         },
#         merge="first"
#     )

def minmax_scale_per_batch(adata, batch_key='batch',
                           quantile=0.99,layer=None):
    
    from sklearn.preprocessing import minmax_scale
    
    X=adata.layers[layer].copy() if layer else adata.X.copy()
    X[np.isnan(X)]=0
    batches=adata.obs[batch_key].unique()
    
    for batch in batches:
        batch_mask = adata.obs[batch_key] == batch
        X_batch = X[batch_mask]
        
        gene_quantiles = np.quantile(X_batch,quantile,axis=0)

        if scipy.sparse.issparse(X_batch):
            X_batch = X_batch.toarray()

        for i in range(X_batch.shape[1]):
            X_batch[X_batch[:, i] > gene_quantiles[i],i]=gene_quantiles[i]
        
        X_batch_scaled = minmax_scale(X_batch,axis=0)
        
        if scipy.sparse.issparse(X):
            X[batch_mask] = scipy.sparse.csr_matrix(X_batch_scaled)
        else:
            X[batch_mask] = X_batch_scaled
    
    adata.layers['minmax_cell_batch'] = X
    return adata

### cluster_L1
var_for_clustering=[
    # Epithelial cell
    'PanCK','GPC3','HepPar1','AFP','CK19',
    # Stromal cell
    'SMA','Collagen_I','CD31',
    # Immune cell
    'CD45','CD3','CD20','CD68','CD66b','CD117']
adata.uns['var_for_clustering']=var_for_clustering
adata.var['highly_variable']=[True if i in var_for_clustering else False for i in adata.var_names]

adata.X[np.isnan(adata.X)]=0
adata=minmax_scale_per_batch(adata,batch_key='dataset',quantile=0.97)
adata.X=adata.layers['minmax_cell_batch']
sc.pp.scale(adata,max_value=5)

sc.tl.pca(adata,svd_solver='arpack',use_highly_variable=True)
sc.pp.neighbors(adata,n_neighbors=20,n_pcs=15)

for r in [0.5,0.8,1,1.2,1.5,1.8]:
    print(r)
    sc.tl.louvain(adata,resolution=r,key_added='louvain_'+str(r))
    adata.write('./Version4/adata_annotated_louvain_0.97_scale.h5ad')
    
sc.tl.umap(adata)

sc.pl.umap(adata,color=var_for_clustering,cmap='bwr',vmin=0,vmax=1)

### cluster_L2
var_for_clustering=[
    'CD3','CD4','CD8A','FoxP3', # T cell
    'CD57', # NK cell
    'CD20', # B cell
    'CD68','CD206','CD163','HLA-DR', # macrophage
    'CD11c', # DC
    'CD14','CD16', # monocyte
    'CD15','CD66b', # granulocyte/neutrophil
    'CD117' # mast cell
]
adata_immune.uns['var_for_clustering']=var_for_clustering
adata_immune.var['highly_variable']=[True if i in var_for_clustering else False for i in adata_immune.var_names]

sc.tl.pca(adata_immune,svd_solver='arpack',use_highly_variable=True)
sc.pp.neighbors(adata_immune,n_neighbors=20,n_pcs=15)

for r in [0.5,0.8,1,1.2,1.5,1.8]:
    print(r)
    sc.tl.louvain(adata_immune,resolution=r,key_added='louvain_'+str(r))
    adata_immune.write('./Version4/adata_immune.h5ad')
    
sc.tl.umap(adata_immune)


### cluster_L3_Epi
var_for_clustering=[
    'PD-L1','HLA-DR','Ki67','CK19','Vimentin','PanCK',
    'Arg1','AFP','E-cadherin','GPC3','HepPar1'
]

adata_epi.uns['var_for_clustering']=var_for_clustering
adata_epi.var['highly_variable']=[True if i in var_for_clustering else False for i in adata_epi.var_names]

sc.tl.pca(adata_epi,svd_solver='arpack',use_highly_variable=True)
sc.pp.neighbors(adata_epi,n_neighbors=20,n_pcs=5)

for r in [1.5,1.8,2,2.2,2.5]:
    print(r)
    sc.tl.louvain(adata_epi,resolution=r,key_added='louvain_'+str(r))
    
sc.tl.umap(adata_epi)

### cluster_L3_Myeloid
var_for_clustering=[
    'CD68','CD206','CD163', # macrophage
    'CD11c', # DC
    'CD14','CD16', # monocyte
    'HLA-DR',# 'Ki67',
    'PD-L1'
]
adata_myeloid.uns['var_for_clustering']=var_for_clustering
adata_myeloid.var['highly_variable']=[True if i in var_for_clustering else False for i in adata_myeloid.var_names]
sc.tl.pca(adata_myeloid,svd_solver='arpack',use_highly_variable=True)
sc.pp.neighbors(adata_myeloid,n_neighbors=20,n_pcs=5)

for r in [0.5,0.8,1,1.2,1.5,1.8,2,2.2,2.5]:
    print(r)
    sc.tl.louvain(adata_myeloid,resolution=r,key_added='louvain_'+str(r))
    
sc.tl.umap(adata_myeloid)

### cluster_L3_T
var_for_clustering=[
    # 'CD3',
    'CD4','CD8A',
    'FoxP3',# 'CD57',
    'HLA-DR','CD69','GranzymeB',
    'PD-1',# 'CTLA-4',
    'Ki67'
]
adata_T.uns['var_for_clustering']=var_for_clustering
adata_T.var['highly_variable']=[True if i in var_for_clustering else False for i in adata_T.var_names]

sc.tl.pca(adata_T,svd_solver='arpack',use_highly_variable=True)
sc.pp.neighbors(adata_T,n_neighbors=50,n_pcs=5)

for r in [0.5,0.8,1,1.2,1.5]:
    print(r)
    sc.tl.louvain(adata_T,resolution=r,key_added='louvain_'+str(r))
    
sc.tl.umap(adata_T)
