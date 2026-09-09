library(dplyr)
library(ggplot2)
library(Seurat)


# data preprocessing
seurat_list=lapply(seurat_list,function(seurat_obj){
    seurat_obj=subset(seurat_obj,subset= (nFeature_RNA>=300)&(nCount_RNA>=500)&(pct_MT<=0.05) )
    return(seurat_obj)
})

# data integration using RPCA
seurat_list=lapply(seurat_list,SCTransform,method="glmGamPoi")
anchor_features=SelectIntegrationFeatures(seurat_list,nfeatures=4000)
seurat_list=seurat_list %>% 
    PrepSCTIntegration(anchor.features=anchor_features) %>% 
    lapply(FUN=RunPCA,features=anchor_features) 

feature_integrate=lapply(seurat_list,function(seurat_obj){
  seurat_obj@assays$SCT@SCTModel.list$model1@feature.attributes
})
feature_integrate=Reduce(intersect,x=lapply(X=feature_integrate,FUN=rownames))

anchor_set=FindIntegrationAnchors(object.list=seurat_list,
                                  reduction="rpca",
                                  dims=1:30,
                                  normalization.method='SCT',
                                  anchor.features=anchor_features,
                                  k.anchor=5,
                                  k.filter=300,
                                  k.score=30)

seurat_obj=IntegrateData(anchorset=anchor_set,
                         dims=1:30, 
                         normalization.method="SCT",
                         features.to.integrate=feature_integrate # anchor_features
                        )

# SCTransform raw data for downstream analysis
seurat_obj <- SCTransform(seurat_obj,assay = "RNA",new.assay.name = "SCT_all",
                          verbose = FALSE,reference.SCT.model = NULL,
                          return.only.var.genes = FALSE, # return all genes for downstream analysis
                          variable.features.n = 3500, # features number considered to be variable
                          vars.to.regress = c('pct_MT', 'nFeature_RNA') # factors to remove, for example ribosomal genes
                         	)

DefaultAssay(seurat_obj)='integrated'
seurat_obj=RunPCA(seurat_obj,npcs=30)
seurat_obj=FindNeighbors(seurat_obj,dims=1:20)
seurat_obj=FindClusters(seurat_obj,algorithm=1,resolution=0.5)
seurat_obj=RunUMAP(seurat_obj,dims=1:20)
