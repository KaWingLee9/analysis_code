library(ggplot2)
library(ComplexHeatmap)

### data preprocessing
library(dplyr)
library(biomaRt)
library(GenomicFeatures)
source('https://github.com/KaWingLee9/in_house_tools/blob/main/BulkRNASeq/RNANormalization.R')

count_mat=dplyr::left_join( count_mat , df[,c('ensembl_gene_id','external_gene_name')],by='ensembl_gene_id')
tpm_mat=NormalizeCount(count_mat,species='mouse')

### QC
group1=c('D3','D3','D3','D3','D3','D3','D3',
        'D3','D7','D7','D7','D7','D7','D7',
        'D7','D7','D7','D7','D7','D7','D14',
        'D14','D14','D14','D14','D14','D14',
        'D14','D14','D14','D14','D14','D0',
        'D0','D0','D0','D0','D0','D0','D0',
        'D0','D0','D0','D0','D3','D3','D3','D3')
group1=setNames(group1,colnames(count_mat))
group1=factor(group1,levels=c('D0','D3','D7','D14'))

p=SampleQC_PCA(tpm_mat,group1,PC=c('PC1','PC2'))

### temporally related gens
library(maSigpro2)

time=gsub('D','',group1) %>% as.numeric()
design=data.frame(time=time,replicate=1:length(time),group=rep(1,length.out=length(time)))
rownames(design)=colnames(count_mat)
design=make.design.matrix(design,degree=2)
tc.p=p.vector(tpm_mat,design,counts=FALSE,MT.adjust='fdr')
tc.tstep=T.fit(tc.p)
sigs=get.siggenes(tc.tstep,rsq=0.15,vars="all")

j=4
cluster=see.genes(sigs$sig.genes,show.fit=TRUE,k=j,cluster.method="hclust",dis=design$dis,cluster.data=1,newX11=F)

genes_cluster=cluster$cut %>% data.frame()

genes_cluster=genes_cluster[rownames(tpm_mat),,drop=FALSE] %>% na.omit()
genes=rownames(genes_cluster)

df=t(scale(t(log(tpm_mat+1)[genes,])))
df_anno_x=data.frame(time=group1)
time=df_anno_x[colnames(df),'time']
time=factor(time,levels=c('D0','D3','D7','D14'))
df_anno_y=data.frame(cluster=genes_cluster)
gcluster=paste0('C',df_anno_y[rownames(df),'cluster'])

ht=Heatmap(df,
           cluster_rows=FALSE,cluster_columns=FALSE,
           show_column_names=FALSE,show_row_names=FALSE,
           column_split=time,
           top_annotation=HeatmapAnnotation(time=time),
           row_split=genes_cluster[rownames(df),],
           left_annotation=rowAnnotation(gcluster=genes_cluster[rownames(df),]))

### transcription factor analysis with decoupleR
library(decoupleR)

tf_network=get_dorothea(organism='human',levels=c('A','B','C'))

exp_mat=as.matrix(tpm_mat)
exp_mat=log(exp_mat+1)
tf_act=run_wmean(mat=exp_mat,net=tf_network)

tf_act_wide=tf_act %>% filter(statistic=='wmean') %>% 
                                         reshape2::dcast(condition~source,value.var='score') %>% 
                                         data.frame(row.names=1,check.names=FALSE)
colnames(tf_act_wide)=paste0(colnames(tf_act_wide),'(+)')
tf_act_wide=tf_act_wide[paste0('Gli1MI.UMI_counts.sc',1:48),]
tf_act_wide[,'time']=c('D3','D3','D3','D3','D3','D3','D3','D3','D7',
'D7','D7','D7','D7','D7','D7','D7','D7','D7',
'D7','D7','D14','D14','D14','D14','D14','D14',
'D14','D14','D14','D14','D14','D14','D0','D0',
'D0','D0','D0','D0','D0','D0','D0','D0','D0',
'D0','D3','D3','D3','D3')

regulon_selected=c('Gli1(+)','Pou4f2(+)','D3YY78(+)','Cebpe(+)','Nr2f2(+)',
                                'E2f2(+)','Nr2c2(+)','Tfdp1(+)','E2f4(+)','Q6P1H7(+)',
                                'Batf(+)','Lyl1(+)','G5E8H3(+)','Spi1(+)','Q3U1Q3(+)',
                                'E2f6(+)','Smad4(+)','Zbtb7a(+)','Maff(+)','Sox11(+)')

df=reshape2::melt(tf_act_wide,id.vars='time',variable.name='regulon',value.name='regulon_score')
df=df %>% filter(regulon %in% regulon_selected)
df[,'time']=factor(df[,'time'],levels=c('D0','D3','D7','D14'))
df[,'regulon']=factor(df[,'regulon'],levels=regulon_selected)

ht=SumHeatmap(df,group.col='time',variable.col='regulon',value.col='regulon_score',
              test.mode='ONEvsOTHER',permutated=FALSE,test.method='t.test',p.adj=FALSE,
              cluster_rows=FALSE,cluster_columns=FALSE,name='Regulon\nactivity')

### gene co-expression
library(igraph)
library(tidygraph)
library(ggraph)
source('https://github.com/KaWingLee9/in_house_tools/blob/main/BulkRNASeq/network.R')
source('https://github.com/KaWingLee9/in_house_tools/tree/main/visulization/custom_fun.R')


regulators=msigdbr::msigdbr(species='Homo sapiens',category='C3',subcategory='TFT:GTRD') %>% data.frame %>% .[,'gs_name'] %>% unique %>% gsub('_TARGET_GENES','',.)
test_result=BuildNetwork(tpm_mat,method='pearson')

df=test_result[['r']][c('Gli1','Atoh8'),] %>% reshape2::melt(varnames=c('gene1','gene2'),value.name='r')
df[,'P']=df %>% apply(1,function(x) {test_result[['P']] [x[1],x[2]] } )

gs_1=df %>% filter(gene1=='Gli1',r>=0.25,P<=0.05) %>% .[,'gene2']
gs_2=df %>% filter(gene1=='Atoh8',r>=0.25,P<=0.05) %>% .[,'gene2']
gs_3=intersect(gs_1,gs_2)

gs_1=setdiff(gs_1,gs_3)
gs_2=setdiff(gs_2,gs_3)

x=c(setNames(rep(1,length.out=length(gs_1)),gs_1),
  setNames(rep(2,length.out=length(gs_2)),gs_2),
      setNames(rep(NA,length.out=2),c('Atoh8','Gli1')) )

df_1=df %>% filter( ((r>0.25)|(r<=-0.25)),P<=0.05)
df_1[ df_1[,'r']<=0 ,'sign']=-1
df_1[ df_1[,'r']>=0 ,'sign']=1

df_1=df_1 %>% filter(gene2 %in% names(x))
g=graph_from_data_frame(df_1)
edge_attr(g,'sign')=as.character(edge_attr(g,'sign'))

vertex_attr(g,'community')[vertex_attr(g,'name')=='Gli1']=1
vertex_attr(g,'community')[vertex_attr(g,'name')=='Atoh8']=2
node_coords=layout_circular_community(vertex_attr(g,'community'),R=15,k=0.3)
rownames(node_coords)=vertex_attr(g,'name')
node_coords[,'community']=vertex_attr(g,'community')
node_coords[,'gene']=rownames(node_coords)
node_coords['Atoh8',1:2]=c(-15,40)
node_coords['Atoh8',3]=3
node_coords['Gli1',1:2]=c(15,40)
node_coords['Gli1',3]=3

p=ggraph(g,layout='manual',x=node_coords[,'x'],y=node_coords[,'y'])+
    geom_edge_fan2(aes(color=sign),width=0.05)+
    geom_node_point(size=0.001)+
    scale_edge_color_manual(values=c('#35978F50','#BFB12D50'))+
    theme_void()