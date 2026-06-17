library(dplyr)
library(ggplot2)
library(ComplexHeatmap)

tendency_result=read.csv('./Edge_level_1.csv',row.names=1,check.names=FALSE)
# tendency_result[,'spot_exp_ratio']=tendency_result[,'spot_exp']/tendency_result[,'Spots number']

gene_anno=filter(tendency_result,spot_exp_ratio>=0.05,p_value<=0.05
                             )
# tendency_result_filter=tendency_result_filter[! grepl('^MT-|^RP[SL]',rownames(tendency_result_filter)),]
# tendency_result_filter=tendency_result

exp_dis=read.csv('./exp_dis_1.csv',row.names=1,check.names=FALSE)
exp_dis=exp_dis[,c('x','y','distance',rownames(tendency_result))]

y=reshape2::melt(exp_dis,id.vars=c('x','y','distance'),variable.name='gene',value.name='exp')
z=reshape2::dcast(y,distance~gene,value.var='exp',fun.aggregate=mean) %>% 
    arrange(distance) %>% .[,-1] %>% scale() %>% t()  %>% data.frame()

z=na.omit(z)
l=rownames(gene_anno)

l=intersect(l,rownames(z))

z_norm=t(z) %>% data.frame()
z_norm=z_norm[,l]
z_norm$'Distance_rank'=1:nrow(z_norm)

z_norm=sapply(setdiff(colnames(z_norm),'Distance_rank'),function(i) 
    loess(as.formula( paste0(i,'~Distance_rank') ),data=z_norm,span=0.5,normalize=FALSE)$fitted) %>% data.frame()
z_norm[,'Distance_rank']=NULL
z_norm=t(z_norm) %>% data.frame()

gene_anno=data.frame(gene_name=l,max_distance=apply(z_norm[l,],1,function(x) which(x==max(x))))

gene_anno[,'p_value']=tendency_result[rownames(gene_anno),'p_value']
# gene_anno=filter(gene_anno,p_value<=0.05)

gene_anno=na.omit(gene_anno)
gene_anno[,'label']=''
gene_anno[which(gene_anno[,'p_value']<=0.05),'label']='*'
gene_anno[which(gene_anno[,'p_value']<=0.01),'label']='**'
gene_anno[which(gene_anno[,'p_value']<=0.001),'label']='***'

gene_anno=gene_anno %>% arrange(max_distance)
# gene_anno=gene_anno %>% filter(gene_name %in% rownames(tendency_result_filter))
row_order=rownames(gene_anno) %>% rev()

z_norm[z_norm<=-1.5]=-1.5
z_norm[z_norm>=1.5]=1.5

# heatmap visualization: 1) LOESS smooth across distance; 2) z-score normalized.
ht=Heatmap(z_norm[row_order,],# row_order=row_order,# column_split=x,
        # row_split=gene_anno[row_order,'group'],
        show_row_names=FALSE,
        cluster_columns=FALSE,show_column_names=FALSE,cluster_rows=FALSE,name='Scaled\nexpression',
        rect_gp = gpar(col="#FFFFFF00",lwd=0),
        # left_annotation=rowAnnotation('Layer'=gene_anno[row_order,'layer_group'],
        #                               'Phalange'=gene_anno[row_order,'phalange_group'],
        #                               col=list('Phalange'=c('DP'='#3969AC','MP'='#E9C39B','PP'='#CC6677',
        #                                                     'Layer2'='#7F3C8D','Layer-shared'='#b85292',
        #                                                     'Layer2'='#7F3C8D',
        #                                                     'Layer3-shared'='#11A579'),
        #                                        'Layer'=c('Layer2'='#7F3C8D','Layer3'='#11A579',
        #                                                  'Layer3-shared'='#11A579','Layer-shared'='#b85292')
        #                                       ),
        #                               show_legend=TRUE),
        right_annotation=rowAnnotation('  '=anno_text(row_order)) )+
        rowAnnotation('Significance'=anno_text(gene_anno[row_order,'label']))

# LOESS-smoothed regression curve
