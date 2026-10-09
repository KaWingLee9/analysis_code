import os
import sys
import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

import anndata as ad
import scanpy as sc

import palantir
import cellrank as cr
import scvelo as scv

adata=sc.read('./adata_pseudotime.h5ad') # converted from seurat_obj in 1C_HumanDataFibro_analysis.R

sc.pp.scale(adata,max_value=10)
sc.pp.pca(adata)
sc.pp.neighbors(adata,use_rep="X_pca",n_pcs=30,n_neighbors=20,metric="euclidean")

# build diffusion map
sc.tl.diffmap(adata)
sc.pl.diffmap(adata,color=["cell_type_final"],components=['0,1','2,3','4,5','6,7','8,9'],wspace=0.5,ncols=3,frameon=True)

# 1. PAGA
sc.tl.paga(adata,groups="cell_type_final")
sc.pl.paga(adata)

# 2. find root
sc.pl.scatter(adata,basis="diffmap",color=["cell_type_final"],components=[1,6])
root_ixs=adata.obsm["X_diffmap"][:,1].argmax()
adata.uns["iroot"]=root_ixs
scv.pl.scatter(adata,basis="diffmap",color=["cell_type_final",root_ixs],components=["1,6"],frameon=True)

# 3. DPT
sc.tl.dpt(adata)

# 4. CellRank
adata.obsm['X_diffmapclip']=adata.obsm['X_diffmap'][:,[1,8]]
pk=cr.kernels.PseudotimeKernel(adata,time_key="dpt_pseudotime")
pk.compute_transition_matrix()

adata.obsm['X_diffmapclip']=adata.obsm['X_diffmap'][:,[3,1,8]]
pk.plot_projection(basis="diffmapclip",recompute=True,
                   color="cell_type_final",legend_loc="right margin",
                   save='cellrank.pdf')