import warnings
warnings.filterwarnings("ignore")
import scanpy as sc
import pandas as pd
import SOAPy_st as sp
import numpy as np
import matplotlib.pyplot as plt

# load edge labels from Loupe Browser
label=pd.read_csv('./Manual_label2.csv',index_col=0,header=0)
adata.obs.loc[label.index.tolist(),'manual_label2']=label.iloc[:,0]
adata=adata[adata.obs['manual_label2'].isin(['Region','Edge']),]

# visualize labels
sc.pl.spatial(adata,color='manual_label2') # crop_coord=[9320,12550,9430,12300]

### PD axis
# define distance of PD axis: distance to the nearest spot of the edge
def find_nearest_points(list1, list2, k=4):
    # Author: Heqi Wang
    # list1: query sites;
    nearest_points_indices = []
    nearest_points_distances = []
    
    for point1 in list1:
        distances = np.linalg.norm(list2 - point1, axis=1)
        nearest_indices = np.argsort(distances)[:k]
        nearest_distances = distances[nearest_indices]
        nearest_points_indices.append(nearest_indices)
        nearest_points_distances.append(nearest_distances)
    
    return np.array(nearest_points_indices), np.array(nearest_points_distances)

list2=adata.obsm['hires_spatial']=adata.obsm['spatial']*adata.uns['spatial']['WSSS_THYst9699525']['scalefactors']['tissue_hires_scalef']
list2=adata.obsm['hires_spatial'][adata.obs['manual_label2']=='Edge',:]
list1=adata.obsm['hires_spatial']
adata.obs['distance']=find_nearest_points(list1, list2, k=1)[1]

# fit polynomial models using SOAPy
mask = sp.tl.get_mask_from_domain(adata, clusters='Edge', KSize=35, cluster_key='manual_label2')
sp.tl.spatial_tendency(adata, mask, radius=1000, location='out', frac=4)
x=adata.uns['SOAPy']['poly']['df_param_poly']
x.loc[:,'spot_exp_ratio']=sc.pp.calculate_qc_metrics(adata)[1]['n_cells_by_counts']/adata.shape[0]
adata.obsm['hires_spatial']=adata.obsm['spatial']*adata.uns['spatial']['WSSS_THYst9699525']['scalefactors']['tissue_hires_scalef']
x.to_csv('./Edge_level.csv')
dis_and_express.to_csv('./exp_dis.csv')

### analysis of AP axis
# define distance of AP axis: distance to the x-axis of the edge
list1=adata.obsm['hires_spatial']
nearest_points_distances = []
for point1 in list1:
    nearest_points_distances.append(-(point1[0]-933.95803668))
adata.obs['distance']=nearest_points_distances

# codes modified from SOAPy
def polynomial_regression(dis_and_express,gene_name,frac,sd=None,drop_zero= False,):

    import statsmodels.api as sm

    mat_only_gene_express=dis_and_express.loc[:, ['distance', gene_name]].values

    if sd != None:
        y = mat_only_gene_express[:, 1]
        mu = np.mean(y)
        std = np.std(y)
        mat_only_gene_express[(mat_only_gene_express[:, 1] < mu - sd * std), 1] = mu - sd * std
        mat_only_gene_express[(mat_only_gene_express[:, 1] > mu + sd * std), 1] = mu + sd * std

    if drop_zero:
        mat_only_gene_express = mat_only_gene_express[~(mat_only_gene_express[:, 1] == 0)]

    x = mat_only_gene_express[:, 0]
    y = mat_only_gene_express[:, 1]
    ind = x.argsort()
    x = x[ind]
    y = y[ind]
    X = x

    if frac > 1:
        for index in range(frac-1):
            X = np.column_stack((X, x ** (index + 2)))
    X = sm.add_constant(X)

    model = sm.OLS(y.astype(float), X)
    res = model.fit()
    y_fitted = res.fittedvalues

    y_min = y_fitted.min()
    y_max = y_fitted.max()
    ran = y_max - y_min

    aic = res.aic
    bic = res.bic

    corr = np.corrcoef(x=x, y=y)[0, 1]

    return res.f_pvalue, x, y_fitted, res.params, ran, corr, len(y), aic, bic

dis_and_express=pd.DataFrame(adata.X.toarray())
dis_and_express.columns=adata.var_names
dis_and_express['distance']=adata.obs['distance']

dic_crd={}
list_name=[]
list_p_value_Poly=[]
list_param_Poly=[]
list_range_Poly=[]
list_corr_gene_distance=[]
list_num_of_spots=[]
list_aic=[]
list_bic=[]

mark=0
for gene_name in adata.var_names:
    f_pvalue,xest,yest,param,ran,corr,spot_num,aic,bic=polynomial_regression(
        dis_and_express=dis_and_express,gene_name=gene_name,frac=4)
    if mark == 0:
        dic_crd['Xest']=xest
        mark=1
    dic_crd[gene_name]=yest
    list_name.append(gene_name)
    list_p_value_Poly.append(f_pvalue)
    list_param_Poly.append(param)
    list_range_Poly.append(ran)
    list_corr_gene_distance.append(corr)
    list_num_of_spots.append(spot_num)
    list_aic.append(aic)
    list_bic.append(bic)

dict_PR = {
    'p_value':list_p_value_Poly,
    'param':list_param_Poly,
    'range':list_range_Poly,
    'correlation':list_corr_gene_distance,
    'Spots number':list_num_of_spots,
    'AIC':list_aic,
    'BIC':list_bic
}

df_param_poly=pd.DataFrame(dict_PR, index=list_name)

x=df_param_poly
x.loc[:,'spot_exp_ratio']=sc.pp.calculate_qc_metrics(adata)[1]['n_cells_by_counts']/adata.shape[0]
adata.obsm['hires_spatial']=adata.obsm['spatial']*adata.uns['spatial']['WSSS_THYst9699525']['scalefactors']['tissue_hires_scalef']
dis_and_express.loc[:,'x']=adata.obsm['hires_spatial'][:,0]
dis_and_express.loc[:,'y']=adata.obsm['hires_spatial'][:,1]
x.to_csv('./Edge_level_1.csv')
dis_and_express.to_csv('./exp_dis_1.csv')