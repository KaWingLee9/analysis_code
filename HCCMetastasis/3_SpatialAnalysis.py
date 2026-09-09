import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import scanpy as sc
import SOAPy_st as sp

sp.pp.make_network(adata,sample_key='dataset',method='knn',
                   cutoff=20,scale=1.0,cluster_key='cluster_L3')

### CCI (modified from SOAPy, not shown in the final result)
# exec(open('https://github.com/KaWingLee9/in_house_tools/blob/main/multiplexed_images_pipeline/neighbourhood_pval.py').read())
# sp.tl.neighborhood_analysis(adata,sample_key='dataset',method='include',
#                             cluster_key='cluster_L3',n_iters=100,n_jobs=30)
# df_result=mean_connectivity(adata,method='pval',cluster_key='cluster_L3',sample_key='dataset',
#                             pval_thres=0.01)

### niche identification
sp.tl.get_c_niche(adata,k_max=10,celltype_key='cluster_L3',sample_key='dataset')

neig_df=adata.uns['SOAPy']['niche']
neig_df=neig_df[np.unique(adata.obs['cluster_L3'])]
cluster_num=10
km=KMeans(n_clusters=cluster_num,max_iter=300)
clusters_res=km.fit_predict(neig_df)
cluster_means=neig_df.groupby(clusters_res).mean()
cluster_means=cluster_means.apply(stats.zscore,axis=0)

# in-situ visualization
exec(open('https://github.com/KaWingLee9/in_house_tools/blob/main/multiplexed_images_pipeline/multiplexed_image_processing.py').read())

d='S15M12'
adata_d=sc.read(d+'.h5ad')
pixel_panel_3={'CD20':'blue','CD3':'yellow','Collagen_I':'green',
               'HepPar1':'red'}
plot_cell_pixel(adata_d,pixel_panel_3,min_quantile=0.02,max_quantile=0.98,show_boundary=False,show_legend=False,return_mat=False)
# p2=plot_cell_pixel(adata_d,pixel_panel_3,min_quantile=0.02,max_quantile=0.98,show_boundary=False,show_legend=True,return_mat=True)
# tifffile.imwrite('p2.tiff',p2.astype(np.uint8))
plot_cell_cluster(adata_d,tag='cluster_L3',col=cell_palette,show_boundary=True,show_legend=True,return_mat=False)
plot_cell_cluster(adata_d,tag='niche_cluster',col=niche_palette,show_boundary=False,show_legend=True,return_mat=False)

d='P86M9'
adata_d=sc.read(d+'.h5ad')
pixel_panel_3={'CD15':'yellow','Collagen_I':'green',
               'CD206':'cyan','PanCK':'red'}
plot_cell_pixel(adata_d,pixel_panel_3,min_quantile=0.02,max_quantile=0.98,show_boundary=False,show_legend=False,return_mat=False)
plot_cell_cluster(adata_d,tag='cluster_L3',col=cell_palette,show_boundary=True,show_legend=True,return_mat=False)
plot_cell_cluster(adata_d,tag='niche_cluster',col=niche_palette,show_boundary=False,show_legend=True,return_mat=False)

### the following codes are run in R
# library(dplyr)
# library(ComplexHeatmap)
# source('https://github.com/KaWingLee9/in_house_tools/blob/main/visulization/custom_fun.R')

# niche_result=clusters_res
# niche_result[,'sample_type']=sample_clustering[ niche_result[,'sample'] ,'c1']
# niche_result[,'sample_type']=as.character(niche_result[,'sample_type'])
# niche_result[,'niche_cluster']=as.character(niche_result[,'niche_cluster'])
# niche_result=na.omit(niche_result)

# nihce_result_summed=niche_result %>% group_by(sample,sample_type,niche_cluster) %>% summarise(Num=n()) %>% data.frame() %>% 
#     group_by(sample,sample_type) %>% mutate(Percentage=Num/sum(Num)) %>% data.frame()

# ht1=SumHeatmap(nihce_result_summed,group.col='sample_type',variable.col='niche_cluster',value.col='Percentage',
#            test.mode='ONEvsOTHER',test.method='wilcox.test',permutated=FALSE,transpose=TRUE,show.significance=TRUE,
#            sig.level=c(0.001,0.01),sig.label=c('**','*'),p.adj=TRUE,p.adj.method='bonferroni',scale=TRUE,
#            name='Scaled\npercentage',
#                cluster_rows=FALSE,cluster_columns=FALSE,
#               col=circlize::colorRamp2( c( seq(-1,-0.8,length=4) ,0, seq(0.2,3,length=4) ), 
#                                         rev (c('#82498CFF','#B574C2FF','#D2A9DBFF','#E8D4EDFF','#F1F1F1FF','#C8E1C9FF','#91C392FF','#4E9D4FFF','#256C26FF'))  ),
#             width=unit(6,"cm"))

# ht2=Heatmap(cluster_means %>% ResetOrder(by='col') # %>% ResetOrder(by='row')
#             ,
#         col=circlize::colorRamp2(c( seq(-1.3,-0.3,length=5),0,seq(0.8,4,length=5) ),
#                                  c('#063061','#2266AC','#4392C5','#93C5DF','#D1E6F1','#F7F7F7',
#                                    '#FCDCC7','#F6A581','#D4614D','#B2182C','#68011F')),
#         cluster_rows=FALSE,cluster_columns=FALSE,
#             width=unit(15,"cm"))

# ht1+ht2

