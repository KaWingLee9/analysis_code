library(dplyr)
library(ggplot2)
library(aplot)
library(ggplotify)

library(Seurat)
library(monocle)

exp_mat=as(as.matrix(seurat_obj@assays$RNA@counts), 'sparseMatrix')
pd=new('AnnotatedDataFrame',data=seurat_obj@meta.data)
fData=data.frame(gene_short_name=row.names(exp_mat),row.names=row.names(exp_mat))
fd=new('AnnotatedDataFrame',data=fData)
cds=newCellDataSet(exp_mat,phenoData=pd,featureData=fd,lowerDetectionLimit=0.5,expressionFamily=negbinomial.size())
rm(exp_mat,pd,fd,fData)

# precrocessing cds object
cds=estimateSizeFactors(cds)
cds=estimateDispersions(cds)
cds=detectGenes(cds)

# gene selection
expressed_genes=row.names(subset(fData(cds),num_cells_expressed>=50))
interested_cells=rownames(pData(cds))

### select one of the two following methods:
# genes along pseudotime
disp_table=dispersionTable(cds)
ordergene=subset(disp_table,mean_expression>=0.1)
cds=setOrderingFilter(cds,ordergene)
plot_ordering_genes(cds)
plot_pc_variance_explained(cds,return_all=FALSE) 

cds=reduceDimension(cds,max_components=2,method='DDRTree',residualModelFormulaStr='~batch+nCount_RNA')
cds=orderCells(cds)

pData(cds)[,'cell_type_final']=Idents(seurat_obj)[rownames(pData(cds))]

