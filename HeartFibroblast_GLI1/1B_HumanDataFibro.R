library(dplyr)
library(ggplot2)
library(Seurat)
library(scRNAtoolVis)

### fibroblast subtype definition
seurat_obj=subset(seurat_obj_all)


### gene markers of fibroblast
future::plan("multicore",workers=10) 
options(future.globals.maxSize=10000*1024^2)
sce.markers=FindAllMarkers(object=seurat_obj,only.pos=FALSE,min.pct=0.25,logfc.threshold=0.25)

fib_feature=c('ACSM3','SCN7A','COL15A1',
              'PCOLCE2','CD55','ACKR3',
              'COL1A1','ACTA2','PDGFRB',
              'POSTN','THBS4','FGF14',
              'NAV3','KAZN','PLA2G2A','FOS','JUNB',
              'APOE','CCL2','IFIT3','IFIT1','ISG15')

cell_order=c('qFibro','PCOLCE3+ Fibro','ELN+ Myofibroblast',
             'POSTN+ Fibro','NAV3+ Fibro','PLA2G2A+ Fibro','APOE+ Fibro',
             'CCL2+ Fibro','IFN Fibro') %>% rev()
p=DotPlot(seurat_obj,features=fib_feature,col.min=0)+
    theme(axis.text.x=element_text(angle=45,hjust=1))+
    scale_y_discrete(limits=cell_order)
p

### transcriptional regulon analysis
library(scRNAtoolVis)
library(OmnipathR)
library(decoupleR)

exp_mat=as.matrix(seurat_obj@assays$RNA@data)
tf_network=get_dorothea(organism='human',levels=c('A','B','C'))
tf_act=run_wmean(mat=exp_mat,net=tf_network)

tf_act_wide=tf_act %>% filter(statistic=='wmean') %>% 
                                         reshape2::dcast(condition~source,value.var='score') %>% 
                                         data.frame(row.names=1,check.names=FALSE)  # %>% matrix()
tf_act_wide[tf_act_wide>=30]=30
tf_act_wide[tf_act_wide<=-5]=-5

colnames(tf_act_wide)=paste0(colnames(tf_act_wide),'(+)')
seurat_obj[['regulon']]=CreateAssayObject(counts=t(tf_act_wide))
DefaultAssay(seurat_obj)='regulon'

regulon.markers=FindAllMarkers(seurat_obj)
top_regulon=regulon.markers %>% filter(avg_log2FC>=0) %>% group_by(cluster) %>% top_n(n=5,wt=avg_log2FC) %>% .[,'gene',drop=TRUE] %>% unique()

cell_order=c('qFibro','PCOLCE3+ Fibro','PLA2G2A+ Fibro','NAV3+ Fibro','APOE+ Fibro','CCL2+ Fibro',
             'IFN Fibro','POSTN+ Fibro','ELN+ Myofibroblast')
regulon_order=c('GLI1(+)','SOX13(+)','MAX(+)','MNT(+)','KLF4(+)','TCF12(+)',
                'STAT6(+)','ELK4(+)','NR1H3(+)','STAT2(+)','IRF9(+)','IRF1(+)',
                'RFX1(+)','SNAI1(+)','SMAD3(+)','TBP(+)',
                'REL(+)','SRF(+)','FOXK1(+)','HSF1(+)') # from top_regulon

p=averageHeatmap(seurat_obj,assays='regulon',
               markerGene=regulon_order,
               # cluster_rows=TRUE,
               # cluster_columns=TRUE,
               cluster.order=cell_order,
               annoCol=TRUE,myanCol=fibro.human.color,
               htRange=c(-2,0,3))
p

p=FeaturePlot(seurat_obj,'GLI1(+)',min.cutoff=0,max.cutoff=1)
p

### fibroblast progenitor score
library(irGSEA)
gene_ls=list('Fibro_progenitor'=c('PI16','COL15A1','CD34','SCA1','DPT',
                                  'CD55','DPP4','HIC1','SCARA5','ACKR3',
                                  'CAV1','S100A4','DPP4','PRX1','THY1',
                                  'CD44','WNT2','ANXA3'))
seurat_obj=irGSEA.score(object=seurat_obj,assay="SCT_all",
                        slot="data",custom=TRUE,geneset=gene_ls,
                        method="AUCell", 
                        aucell.MaxRank=ceiling(0.05*nrow(seurat_obj)),
                        seeds=123,ncores=20)

seurat_obj@meta.data[,'AUCell']=as.matrix(seurat_obj@assays$AUCell@data)['Fibro-progenitor',] 
p=ggplot(seurat_obj@meta.data,aes(x=reorder(cell_type,AUCell,median),
                                y=AUCell,fill=cell_type))+
    geom_boxplot()+
    theme_classic()+
    theme(axis.text.x=element_text(angle=45,hjust=1))+
    scale_fill_manual(values=fibro.human.color)+
    guides(fill='none')+
    xlab('')+
    ylab('AUCell score')

p