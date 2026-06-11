library(dplyr)
library(Seurat)
library(OmnipathR)
library(decoupleR)

cell_color=c('qFibro'='#E64B35FF','Postn+ Fibro'='#4DBBD5FF','Myofibroblast'='#00A087FF','Pcolce2+ Fibro'='#3C5488FF',
             'Ccl2+ Fibro'='#F39B7FFF','Apoe+ Fibro'='#8491B4FF','Fos/Jun Fibro'='#91D1C2FF','IFN Fibro'='#7E6148FF')

# fibroblast subtype annotation


seurat_obj=RunUMAP(seurat_obj,dims=1:20,reduction='harmony')

fibro.markers=c('Col1a1','Col1a2','Acta2','Ifit1','Ifit3','Isg15',
                'Ccl2','Pcolce2','Ackr3','Scara5','Cd55',
                'Fos','Jun','Junb','Apoe',
                'Postn','Thbs4','Col15a1','Scn7a')
DefaultAssay(seurat_obj)='RNA'
p=DotPlot(seurat_obj,features=fibro.markers,col.min=0)
p=p+theme(axis.text.x=element_text(angle=45,hjust=1))

# similarity of human and mouse fibroblast
query_sce=seurat_obj

annotation_ref=read.csv('gene_matches_1v1_mouse2human.csv',row.names=1,check.names=FALSE) # gene matching from: https://xingyanliu.github.io/CAME/tutorials.html
exp_mat=as.matrix(query_sce@assays$RNA@counts)
ind=intersect(rownames(exp_mat),rownames(annotation_ref))
exp_mat_transformed=exp_mat[ind,]
rownames(exp_mat_transformed)=annotation_ref[ind,1]

query_sce_human=CreateSeuratObject(counts=exp_mat_transformed,assay='RNA')
query_sce_human@meta.data=query_sce@meta.data
query_sce_human=SCTransform(query_sce_human,assay="RNA",new.assay.name="SCT")

transfer_anchor=FindTransferAnchors(reference=ref_sce,query=query_sce_human,dims=1:30,
                                    normalization.method='SCT',
                                    reference.reduction='pca')
query_prediction=TransferData(anchorset=transfer_anchor,
                              refdata=Idents(ref_sce),
                              dims=1:30)

confusion_matrix=table(query_prediction[,c('predicted.id','original.id')]) %>% as.data.frame.array()
confusion_matrix=confusion_matrix %>% apply(1,function(x) x/sum(x)) %>% t() %>% data.frame()
confusion_matrix=scale(confusion_matrix %>% t()) %>% t()

p=Heatmap(confusion_matrix %>% ResetOrder(by='row') %>% ResetOrder(by='col') %>% scale(),
        cluster_rows=FALSE,
        cluster_columns=FALSE,
        col=circlize::colorRamp2(# c(seq(-1,3,length.out=10)),
                                 c(c(seq(-1,2,length.out=4),seq(1.8,2.5,length.out=6))),
                                 c(viridis::viridis_pal(option='viridis_pal')(10)))
        )

# transcription factor
library(OmnipathR)
library(decoupleR)

tf_network=get_dorothea(organism='mouse',levels=c('A','B','C'))

exp_mat=as.matrix(seurat_obj@assays$RNA@data)
tf_act=run_wmean(mat=exp_mat,net=tf_network)

tf_act_wide=tf_act %>% filter(statistic=='wmean') %>% 
                                         reshape2::dcast(condition~source,value.var='score') %>% 
                                         data.frame(row.names=1,check.names=FALSE)

colnames(tf_act_wide)=paste0(colnames(tf_act_wide),'(+)')
seurat_obj[['regulon']]=CreateAssayObject(counts=t(tf_act_wide))

df=data.frame('Gli1_regulon'=seurat_obj_subset@assays$regulon@data['Gli1(+)',],
              fibro_type=Idents(seurat_obj_subset))
df$binary=as.factor(df[,'fibro_type']=='Pcolce2+ Fibro')
options(repr.plot.width=5,repr.plot.height=5)
p=ggplot(data=df,aes(x=reorder(fibro_type,-Gli1_regulon),y=Gli1_regulon))+
    geom_boxplot()+
    theme_classic()+
    theme(axis.text.x=element_text(angle=45,hjust=1))
