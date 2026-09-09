library(dplyr)
library(Seurat)

# label transfer using Seurat
transfer_anchor=FindTransferAnchors(reference=ref_sce,query=query_sce,dims=1:30,
                                    normalization.method='SCT',
                                    reference.reduction='pca')

query_prediction=TransferData(anchorset=transfer_anchor,
                              refdata=Idents(ref_sce),
                              dims=1:30)

query_sce=IntegrateEmbeddings(anchorset=transfer_anchor,reference=ref_sce,
                              query=query_sce,new.reduction.name="ref.pca")
query_sce=ProjectUMAP(reference=ref_sce,query=query_sce,
                      reference.reduction="pca",query.reduction="ref.pca",
                      reduction.model="umap")

# merge reference and query
ref_count=ref_sce@assays$RNA@counts
query_count=query_sce@assays$RNA@counts
count_mat=RowMergeSparseMatrices(ref_count,query_count)

ref_meta=ref_sce@meta.data[,c('nCount_RNA','pct_MT','pct_RB','batch','predict_celltype')]
query_meta=query_sce@meta.data[,c('nCount_RNA','pct_MT','pct_RB','batch','predict_celltype')]

seurat_obj=CreateSeuratObject(counts=count_mat,meta.data=rbind(ref_meta,query_meta))

options(future.globals.maxSize=1000*1024^3)
seurat_obj=SCTransform(seurat_obj,assay ="RNA",new.assay.name = "SCT_all",
                          verbose = FALSE,reference.SCT.model = NULL,
                          return.only.var.genes = FALSE, variable.features.n = 3500, vars.to.regress = c('pct_MT', 'nFeature_RNA') )

