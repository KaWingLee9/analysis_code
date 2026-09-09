library(dplyr)
library(Seurat)

ref_sce=readRDS('./bigSeu.anno.rds')
query_sce=readRDS('./query_sce_1.rds')

fibro.human.color=c('qFibro'='#E64B35FF','PCOLCE3+ Fibro'='#3C5488FF','POSTN+ Fibro'='#4DBBD5FF','CCL2+ Fibro'='#F39B7FFF',
                        'PLA2G2A+ Fibro'='#91D1C2FF','APOE+ Fibro'='#8491B4FF','ELN+ Myofibroblast'='#00A087FF')

ref_sce=RunUMAP(ref_sce,reduction="pca",dims=1:30,min.dist=0.01,seed.use=202320,n.neighbors=50L,return.model=TRUE)

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

ref_sce=readRDS('./bigSeu.anno.rds')
query_sce=readRDS('./query_sce_2.rds')
query_prediction=read.csv('./query_prediction_all.csv',row.names=1,check.names=FALSE)
query_sce@meta.data[,'batch']=Idents(query_sce)
query_sce@meta.data[,'predict_celltype']=query_prediction[rownames(query_sce@meta.data),'predicted.id']
Idents(query_sce)=query_sce@metat.data[,'predict_celltype']

ref_count=ref_sce@assays$RNA@counts
query_count=query_sce@assays$RNA@counts
count_mat=RowMergeSparseMatrices(ref_count,query_count)

ref_meta=ref_sce@meta.data[,c('nCount_RNA','pct_MT','pct_RB','batch','predict_celltype')]
query_meta=query_sce@meta.data[,c('nCount_RNA','pct_MT','pct_RB','batch','predict_celltype')]

seurat_obj=CreateSeuratObject(counts=count_mat,meta.data=rbind(ref_meta,query_meta))

umap_coord=rbind(ref_sce@reductions$umap@cell.embeddings,query_sce@reductions$ref.umap@cell.embeddings)

seurat_obj@reductions$umap=CreateDimReducObject(embeddings=umap_coord,key='umap')

options(future.globals.maxSize=1000*1024^3)
seurat_obj=SCTransform(seurat_obj,assay ="RNA",new.assay.name = "SCT_all",
                          verbose = FALSE,reference.SCT.model = NULL,
                          return.only.var.genes = FALSE, variable.features.n = 3500, vars.to.regress = c('pct_MT', 'nFeature_RNA') )


#============
# extract fibroblast
query_sce=readRDS('./query_sce_2.rds')
query_sce@meta.data[,'predict_celltype']=query_prediction[rownames(query_sce@meta.data),'predicted.id']
Idents(query_sce)=query_prediction[rownames(query_sce@meta.data),'predicted.id']
query_fibro=subset(query_sce,idents='Fibroblasts')

ref_fibro=readRDS('./Fibro_new_result.rds')
query_fibro=readRDS('./query_fibro.rds')

ref_fibro=RunUMAP(ref_fibro,reduction="pca",dims=1:30,min.dist=0.01,seed.use=202320,n.neighbors=50L,return.model=TRUE)

transfer_anchor=FindTransferAnchors(reference=ref_fibro,query=query_fibro,dims=1:30,
                                    normalization.method='SCT',
                                    reference.reduction='pca')

query_prediction=TransferData(anchorset=transfer_anchor,
                              refdata=Idents(ref_fibro),
                              dims=1:30)

query_fibro=IntegrateEmbeddings(anchorset=transfer_anchor,reference=ref_fibro,
                              query=query_fibro,new.reduction.name="ref.pca")
query_fibro=ProjectUMAP(reference=ref_fibro,query=query_fibro,
                      reference.reduction="pca",query.reduction="ref.pca",
                      reduction.model="umap")

saveRDS(query_fibro,'query_fibro_1.rds') 

write.csv(query_prediction,'./query_prediction_fibro.csv')

library(dplyr)
library(Seurat)

ref_fibro=readRDS('./Fibro_new_result.rds')
query_fibro=readRDS('./query_fibro_1.rds')
query_prediction=read.csv('./query_prediction_fibro.csv')

ref_count=ref_fibro@assays$RNA@counts
query_count=query_fibro@assays$RNA@counts
count_mat=RowMergeSparseMatrices(ref_count,query_count)

ref_meta=ref_fibro@meta.data[,c('nCount_RNA','batch')]
ref_meta[,'final_anno']=Idents(ref_fibro)
query_fibro@meta.data[,'batch']=query_fibro@meta.data[,'orig.ident']
query_meta=query_fibro@meta.data[,c('nCount_RNA','batch')]
query_meta[,'final_anno']=query_prediction[,'predicted.id']

seurat_obj=CreateSeuratObject(counts=count_mat,meta.data=rbind(ref_meta,query_meta))

umap_coord=rbind(ref_fibro@reductions$umap@cell.embeddings,query_fibro@reductions$ref.umap@cell.embeddings)

seurat_obj@reductions$umap=CreateDimReducObject(embeddings=umap_coord,key='umap')


seurat_obj@meta.data[,'atlas']=NA

seurat_obj@meta.data[ seurat_obj@meta.data[,'batch'] %in% 
        c('D5_scRNAseq','D5_snRNAseq','D7_scRNAseq','D7_snRNAseq','patient1_1','patient1_2','patient2') ,'atlas']='core'

seurat_obj@meta.data[ is.na(seurat_obj@meta.data[,'atlas']) ,'atlas']='query'

p=DimPlot(seurat_obj,split.by='atlas')+scale_color_manual(values=fibro.human.color)
ggsave('./fibro_plot.pdf',p,width=7,height=7)

p=DimPlot(seurat_obj)+scale_color_manual(values=fibro.human.color)+NoLegend()
ggsave('./fibro_plot_combined.pdf',p,width=7,height=7)

Idents(seurat_obj)=seurat_obj@meta.data[,'final_anno']

seurat_obj@reductions$umap@cell.embeddings[ seurat_obj@meta.data[,'atlas']=='core' ,1]=-seurat_obj@reductions$umap@cell.embeddings[ seurat_obj@meta.data[,'atlas']=='core'  ,1]



library(dplyr)
library(Seurat)
library(GSEABase)
library(AUCell)

seurat_obj=readRDS('seurat_obj_fibro_mapping.rds')
exp_mat=seurat_obj@assays$RNA@data

gene_ls=GeneSet(c('PI16','COL15A1','CD34','SCA1','DPT','CD55','HIC1','SCARA5','ACKR3',
          'CAV1','S100A4','DPP4','PRX1','THY1','CD44','WNT2','ANXA3'),setName="Fibro_progenitor")

AUCell_result=AUCell_run(exp_mat,gene_ls)

AUCell_result@assays@data@listData$AUC
seurat_obj@meta.data[,'AUCell']=AUCell_result@assays@data@listData$AUC[1,]

p=ggplot(seurat_obj@meta.data,aes(x=reorder(final_anno,AUCell,median),
                                y=AUCell,fill=final_anno))+
    geom_boxplot()+
    theme_classic()+
    theme(axis.text.x=element_text(angle=45,hjust=1))+
    scale_fill_manual(values=fibro.human.color)+
    guides(fill='none')+
    xlab('')+
    ylab('AUCell score')
ggsave('./AUCellscore.pdf',p)

p=ggplot(seurat_obj@meta.data,aes(x=reorder(final_anno,AUCell,median),
                                y=AUCell,fill=final_anno))+
    geom_boxplot()+
    theme_classic()+
    theme(axis.text.x=element_text(angle=45,hjust=1))+
    scale_fill_manual(values=fibro.human.color)+
    guides(fill='none')+
    xlab('')+
    ylab('AUCell score')+
    facet_wrap(.~atlas)
ggsave('./AUCellscore_split.pdf',p,width=14,height=7)

### transcription factor

library(OmnipathR)
library(decoupleR)

exp_mat=as.matrix(seurat_obj@assays$RNA@data)
tf_network=get_dorothea(organism='human',levels=c('A','B','C'))

tf_network=read.csv('./Dorothea_TFnetwork.csv',row.names=1,check.names=FALSE)
tf_network=filter(tf_network,confidence %in% c('A','B','C'))

tf_act=run_wmean(mat=exp_mat,net=tf_network)
tf_act_wide=tf_act %>% filter(statistic=='corr_wmean') %>% 
                                         reshape2::dcast(condition~source,value.var='score') %>% 
                                         data.frame(row.names=1,check.names=FALSE)  # %>% matrix()

colnames(tf_act_wide)=paste0(colnames(tf_act_wide),'(+)')
seurat_obj[['regulon']]=CreateAssayObject(counts=t(tf_act_wide))
DefaultAssay(seurat_obj)='regulon'

Idents(seurat_obj)=seurat_obj@meta.data[,'final_anno']
regulon.markers=FindAllMarkers(seurat_obj,logfc.threshold=0,min.pct=0)
top_regulon=regulon.markers %>% filter(avg_log2FC>=0) %>% group_by(cluster) %>% top_n(n=5,wt=avg_log2FC) %>% .[,'gene',drop=TRUE] %>% unique()

#======== query
seurat_obj_query=subset(seurat_obj, subset=atlas=='query'  )

exp_mat=as.matrix(seurat_obj_query@assays$RNA@data)
# tf_network=get_dorothea(organism='human',levels=c('A','B','C'))

tf_network=read.csv('./Dorothea_TFnetwork.csv',row.names=1,check.names=FALSE)
tf_network=filter(tf_network,confidence %in% c('A','B','C'))

tf_act=run_wmean(mat=exp_mat,net=tf_network)
tf_act_wide=tf_act %>% filter(statistic=='wmean') %>% 
                                         reshape2::dcast(condition~source,value.var='score') %>% 
                                         data.frame(row.names=1,check.names=FALSE)  # %>% matrix()

colnames(tf_act_wide)=paste0(colnames(tf_act_wide),'(+)')
seurat_obj_query[['regulon']]=CreateAssayObject(counts=t(tf_act_wide))
DefaultAssay(seurat_obj_query)='regulon'

Idents(seurat_obj_query)=seurat_obj_query@meta.data[,'final_anno']
regulon.markers=FindAllMarkers(seurat_obj_query,logfc.threshold=0,min.pct=0)
top_regulon=regulon.markers %>% filter(avg_log2FC>=0) %>% group_by(cluster) %>% top_n(n=5,wt=avg_log2FC) %>% .[,'gene',drop=TRUE] %>% unique()

#========



cell_order=c('qFibro','PCOLCE3+ Fibro','PLA2G2A+ Fibro','NAV3+ Fibro','APOE+ Fibro','IFN Fibro','CCL2+ Fibro',
             'POSTN+ Fibro','ELN+ Myofibroblast')
# regulon_order=c('GLI1(+)','SOX13(+)','MAX(+)','MNT(+)','KLF4(+)','TCF12(+)',
#                 'STAT6(+)','ELK4(+)','NR1H3(+)','STAT2(+)','IRF9(+)','IRF1(+)',
#                 'RFX1(+)','SNAI1(+)','SMAD3(+)','TBP(+)',
#                 'REL(+)','SRF(+)','FOXK1(+)','HSF1(+)') # from top_regulon

# "LHX2(+)",   "YBX1(+)",   "EBF1(+)",   "ETV1(+)",   "GLI2(+)",   "SNAI1(+)",
# "TBP(+)",    "RFX1(+)",   "SMAD3(+)",  "ETV6(+)",   "GLI1(+)",   "ETV4(+)",
# "KLF4(+)",   "RARG(+)",   "NANOG(+)",  "MAFB(+)",   "HOXB13(+)", "BCL6(+)",
# "NCOA2(+)",  "STAT6(+)",  "ELK4(+)",   "ELK1(+)",   "SMAD1(+)",  "FOXL2(+)",
# "NR1H3(+)",  "RARB(+)",   "ATF6(+)",   "REL(+)",    "WT1(+)",    "HIC1(+)",
# "E2F3(+)",   "SRF(+)",    "RXRG(+)",   "FOXK1(+)",  "STAT2(+)",  "IRF9(+)",
# "IRF1(+)",   "IRF2(+)",   "STAT4(+)"


regulon_order=c("LHX2(+)",   "YBX1(+)",   "EBF1(+)",   "ETV1(+)",
"GLI1(+)",   "ETV4(+)", "KLF4(+)",   "RARG(+)",
"STAT6(+)",  "ELK4(+)",   "ELK1(+)",   "SMAD1(+)",  "FOXL2(+)",
"NANOG(+)",  "MAFB(+)",   "HOXB13(+)", "BCL6(+)", "NCOA2(+)",
"NR1H3(+)",  "RARB(+)",   "ATF6(+)",
"STAT2(+)",  "IRF9(+)", "IRF1(+)",   "IRF2(+)",   "STAT4(+)",
"REL(+)","E2F3(+)",
"SNAI1(+)", "TBP(+)",    "RFX1(+)",   "SMAD3(+)",  "ETV6(+)",
"SRF(+)",    "RXRG(+)",   "FOXK1(+)") # from top_regulon

p=averageHeatmap(seurat_obj,assays='regulon',
               markerGene=regulon_order,
               # cluster_rows=TRUE,cluster_columns=TRUE,
               cluster.order=cell_order,
               # annoCol=TRUE,myanCol=fibro.human.color,
               htRange=c(-2,0,3))
pdf('./regulon_heatmap.pdf')
print(p)
dev.off()

regulon.markers %>% filter(avg_log2FC>=0) %>% group_by(cluster) %>% top_n(n=5,wt=avg_log2FC) %>% data.frame()

x=intersect(rownames(seurat_obj@meta.data),  rownames(seurat_obj@reductions$umap@cell.embeddings))
seurat_obj=subset(seurat_obj,cells=x)