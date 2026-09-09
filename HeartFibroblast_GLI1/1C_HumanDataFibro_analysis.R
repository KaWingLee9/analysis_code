library(dplyr)
library(ggplot2)
library(Seurat)
library(scRNAtoolVis)

### data integration of fibroblast from core dataset
seurat_obj=subset(seurat_obj_all,ident='Fibroblast') # seurat_obj_all from 1A_HumanAtlas_core.R
seurat_list=SplitObject(seurat_obj,split.by="batch")
seurat_list=lapply(X=seurat_list, function(x){
	x=SCTransform(x,
		assay = "RNA",
		new.assay.name = "SCT",
		variable.features.n = 3500,
        vars.to.regress = c('pct_MT', 'nFeature_RNA'),
		return.only.var.genes = F,
		seed.use = 202320,
		verbose = FALSE
	)
})

features=SelectIntegrationFeatures(
	object.list=seurat_list,
	nfeatures=2000)
seurat_list=PrepSCTIntegration(
	object.list=seurat_list,
	anchor.features=features)
seurat_list=lapply(X=seurat_list,FUN=RunPCA,features=features)
seurat_obj=FindIntegrationAnchors(object.list=seurat_obj,
	normalization.method='SCT',
	anchor.features=features, dims=1:30,
	reduction="rpca", 
    k.anchor=10
)
seurat_obj=IntegrateData(anchorset=seurat_obj,
	normalization.method='SCT',
	features.to.integrate=features,
	dims=1:30)

seurat_obj=RunPCA(seurat_obj,verbose=FALSE,features=features)
seurat_obj=RunUMAP(seurat_obj,reduction="pca",dims=1:30,min.dist=0.01,seed.use=202320,n.neighbors=50L)
seurat_obj=FindNeighbors(seurat_obj,reduction="pca", dims = 1:30)
seurat_obj=FindClusters(seurat_obj,resolution=c(0.2, 0.5, 1))

seurat_obj=SCTransform(seurat_obj,assay="RNA",new.assay.name="SCT_all",variable.features.n=2000,
                       vars.to.regress=c('pct_MT', 'nFeature_RNA'),
                       return.only.var.genes=F,seed.use=202320,verbose=FALSE)

### label transfer of fibroblast from query dataset
ref_fibro=seurat_obj
query_fibro=subset(query_sce,idents='Fibroblasts') # query_sce from 1B_HumanData_query.R

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

cell_order=c('qFibro','SCARA5+ Fibro','Myofibroblast',
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

cell_order=c('qFibro','SCARA5+ Fibro','PLA2G2A+ Fibro','NAV3+ Fibro','APOE+ Fibro','CCL2+ Fibro',
             'IFN Fibro','POSTN+ Fibro','Myofibroblast')
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