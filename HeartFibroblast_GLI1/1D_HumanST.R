library(dplyr)
library(Seurat)
library(spacexr)

### cell type deconvolution
RCTD_result=create.RCTD(SpatialRNA(data.frame(x=rep(1,dim(stRNA_sce)),y=rep(1,dim(stRNA_sce))),GetAssayData(stRNA_sce,assay='RNA',slot='counts'),use_fake_coords=TRUE), 
            Reference(GetAssayData(scRNA_sce,assay='RNA',slot='counts'),Idents(scRNA_sce)), 
            max_cores=20,test_mode=FALSE) %>%
  run.RCTD(doublet_mode='full')

deconv_result=data.frame(RCTD_result@results$weights,check.names=FALSE)
deconv_result=normalize_weights(deconv_result)

### pie plot for all cell type in ST sample
deconv_result=cbind(deconv_result[,non_fibro],
                    Fibroblasts=deconv_result[,fibro] %>% apply(1,sum))

deconv_result=reshape2::melt(as.matrix(deconv_result),varnames=c('spot','cell_type'),value.name='Proportion')
deconv_result=cbind(deconv_result,stRNA_sce@images$image@coordinates[deconv_result[,'spot'],c('imagerow','imagecol')])

options(repr.plot.width=10,repr.plot.height=10)
p1=ggplot(deconv_result,aes(x=imagecol,y=-imagerow,group=spot))+
    jjPlot::geom_jjPointPie(aes(fill=cell_type,pievar=Proportion,width=0.14),
                            color=NA,add.circle=FALSE)+
    geom_point(size=3.4,shape=21,color='#9C9C9C',fill=NA,stroke=0.1)+
    theme_classic()+
    theme(axis.line=element_blank(),axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())+
    scale_fill_manual(values=cell.type.color)

p1

### proportion of fibroblast subtypes
stRNA_sce[['decvon_result']]=CreateAssayObject(counts=t(deconv_result))
DefaultAssay(stRNA_sce)='decvon_result'
p2=SpatialFeaturePlot(stRNA_sce,features=fibro,ncol=5)
p2

### cell type composition between normal and patient samples
deconv_result_normal=read.csv(file.path(pathRCTD, 'normal_RCTD_result_FBsubtype.csv'), row.names=1,check.names=FALSE) %>% as.matrix()
deconv_result_normal=reshape2::melt(deconv_result_normal,varnames=c('spot','cell_type'),value.name='proportion')
deconv_result_normal[,'condition']='Normal'
deconv_result_patient=read.csv(file.path(pathRCTD, 'patient_RCTD_result_FBsubtype.csv'),row.names=1,check.names=FALSE) %>% as.matrix()
deconv_result_patient=reshape2::melt(deconv_result_patient,varnames=c('spot','cell_type'),value.name='proportion')
deconv_result_patient[,'condition']='Patient'
deconv_result=rbind(deconv_result_normal,deconv_result_patient)

options(repr.plot.width = 15, repr.plot.height = 4)
ggboxplot(data=filter(deconv_result, !cell_type %in% c('APOE+ Fibro','CCL2+ Fibro','ELN+ Myofibroblast','IFN Fibro',
                                                'NAV3+ Fibro','PCOLCE2+ Fibro','PLA2G2A+ Fibro','POSTN+ Fibro','qFibro','other')),
          x='condition', y='proportion', fill='condition', color="condition", width=0.6, outlier.alpha = 0.5,  outlier.size = 0.75)+
    scale_color_manual(values = c('Normal'="#EFAB67", "Patient"="#40220F"))+
scale_fill_manual(values = c('Normal'="#FFD5B3", "Patient"="#7F4528"))+
    facet_grid(.~cell_type)+
    stat_compare_means(method='t.test',paired=FALSE,label='p.signif',label.x=1.5,
                       symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*",  "ns")))+
    theme(strip.background=element_blank(),axis.text.x=element_text(angle=45,hjust=1))+
    guides(fill='none')

deconv_result2 <- deconv_result[which(deconv_result$cell_type %in% c('APOE+ Fibro','CCL2+ Fibro','ELN+ Myofibroblast','IFN Fibro',
                                                'NAV3+ Fibro','PCOLCE2+ Fibro','PLA2G2A+ Fibro','POSTN+ Fibro','qFibro')), ]
deconv_result2$cell_type <- factor(deconv_result2$cell_type,
				   levels = c('APOE+ Fibro','CCL2+ Fibro','ELN+ Myofibroblast','IFN Fibro',
                                                'NAV3+ Fibro','PCOLCE2+ Fibro','PLA2G2A+ Fibro','POSTN+ Fibro','qFibro'))

ggboxplot(data=deconv_result2,
          x='condition',y='proportion',fill='condition', color="condition", width=0.6, outlier.alpha = 0.5,  outlier.size = 0.75)+
    facet_grid(.~cell_type)+
    scale_color_manual(values = c('Normal'="#EFAB67", "Patient"="#40220F"))+
    scale_fill_manual(values = c('Normal'="#FFD5B3", "Patient"="#7F4528"))+
    stat_compare_means(method='t.test',paired=FALSE,label='p.signif',label.x=1.5,
                       symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*",  "ns")))+
    theme(strip.background=element_blank(),axis.text.x=element_text(angle=45,hjust=1))+
    guides(fill='none')# +coord_cartesian(ylim=c(0,0.3))

### fibroblast composition between normal and fibrotic region in patient sample
deconv_pat=read.csv(file.path(pathRCTD, 'patient_RCTD_result_FBsubtype.csv'), row.names=1, check.names=FALSE)
df=reshape2::melt(as.matrix(deconv_pat),varnames=c('spot','cell_type'),value.name='proportion')
df[,'region_anno']=region_anno[df[,'spot'],'region_anno']

df2 <- df[which(df$cell_type %in% c('APOE+ Fibro','CCL2+ Fibro','ELN+ Myofibroblast','IFN Fibro',
                                                'NAV3+ Fibro','PCOLCE2+ Fibro','PLA2G2A+ Fibro','POSTN+ Fibro','qFibro')), ]
df2$cell_type <- factor(df2$cell_type,
				   levels = c('APOE+ Fibro','CCL2+ Fibro','ELN+ Myofibroblast','IFN Fibro',
                                                'NAV3+ Fibro','PCOLCE2+ Fibro','PLA2G2A+ Fibro','POSTN+ Fibro','qFibro'))

ggboxplot(data=df2,
          x='region_anno',y='proportion',fill='region_anno', color='region_anno', outlier.size = 0.75, width = 0.66)+
scale_fill_manual(values = c('Fibro_R'="#aa2176", "NonFibro_R"="#ffd2d2"))+
scale_color_manual(values = c('Fibro_R'="#631140", "NonFibro_R"="#FC9D9D"))+
    facet_grid(.~cell_type)+
    stat_compare_means(method='t.test',paired=FALSE,label='p.signif',label.x=1.5,
                       symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*",  "ns")))+
    theme(strip.background=element_blank(),axis.text.x=element_text(angle=45,hjust=1))+
    guides(fill='none')
