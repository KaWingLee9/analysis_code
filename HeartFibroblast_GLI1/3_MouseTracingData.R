library(dplyr)
library(ggplot2)
library(aplot)
library(ggplotify)

library(Seurat)
library(monocle)

source('https://github.com/KaWingLee9/analysis_code/blob/main/HeartFibroblast_GLI1/monocle2_function.R')

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

cds=reduceDimension(cds,max_components=2,method='DDRTree')
cds=orderCells(cds)

pData(cds)[,'cell_type_final']=Idents(seurat_obj)[rownames(pData(cds))]

p=plot_cell_trajectory(cds,color_by='Pseudotime')+
    scale_color_gradientn(colors=paletteer::paletteer_c("grDevices::Inferno",30))

# pseudotime ~ cell type, treatment time, gene expression
options(repr.plot.height=5,repr.plot.width=8)
p1=plot_cell_type(cds,grouping_col='cluster_annotated',
               order=rev(c('qFibro','Fos/Jun Fibro','Pcolce2+ Fibro','Apoe+ Fibro','Postn+ Fibro',
                       'Ccl2+ Fibro','Myofibroblast')))+
    scale_fill_manual(values=cell_color)
p2=plot_cell_type(cds,grouping_col='batch',order=rev(c('MID0','MID3','MID7','MID14')))

p3=plot_gene_heatmap(cds,gene_names=c('Scara5','Pi16','Gli1','Thbs4',
                                      'Postn','Col1a1','Fn1','Acta2','Cthrc1'))

p=p2 %>% insert_bottom(p1) %>% insert_bottom(p3)