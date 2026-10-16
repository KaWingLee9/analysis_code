# functions for finding the branch containing most cells of root cells, by Zhixuan Tang/Ping Lin
# group.by: cell type information; group.value: the cell type in "group.by"
# usage: cds=orderCells(cds,root_state=FindRoot(cds,group.by='cluster_annotated',group.value='qFibro'))
FindRoot=function(cds,group.by,group.value){
  if (length(unique(pData(cds)$State))>1){
    counts=table(pData(cds)$State, pData(cds)[[group.by]])[,group.value]
    return(as.numeric(names(counts)[which
          (counts==max(counts))]))
  } else {
    return (1)
  }
}

# plot cell type distribution across the trajectory
# required packages: ggridges
plot_cell_type=function(cds,grouping_col='seurat_clusters',pseudotime_col='Pseudotime',
                        order=NULL){
    library(ggridges)
    
    df=pData(cds)[,c(grouping_col,pseudotime_col)]

    if (!is.null(order)){
        df[,grouping_col]=factor(df[,grouping_col],levels=order)
    }

    df=df %>% na.omit()
    
    p=ggplot(df,aes_string(x=pseudotime_col,y=grouping_col,fill=grouping_col))+
        geom_density_ridges(scale=2,alpha=0.7,color='white',na.rm=TRUE)+
        theme_bw()

    return(p)
}

# plot smoothed gene expression in heatmap
plot_gene_heatmap=function(cds,gene_names){
    new_data=data.frame(Pseudotime=seq(min(pData(cds)$Pseudotime),max(pData(cds)$Pseudotime),length.out=100))
    df=genSmoothCurves(cds[gene_names,],new_data,
                       trend_formula="~sm.ns(Pseudotime,df=3)" )
    df=df %>% plyr::aaply(1,scale)
    df[df>=3]=3
    df[df<=-3]=-3
    df=reshape2::melt(df,varnames=c('gene','Pseudotime'),value.name='expression')
    df[,'gene']=factor(df[,'gene'],levels=rev(gene_names))

    old_min=min(df[,'Pseudotime'])
    old_max=max(df[,'Pseudotime'])
    new_min=min( pData(cds)[,'Pseudotime'] )
    new_max=max( pData(cds)[,'Pseudotime'] )
    df[,'old_Pseudotime']=df[,'Pseudotime']
    df[,'Pseudotime']=(df[,'old_Pseudotime']-old_min)/(old_max-old_min)*(new_max - new_min)+new_min

    p3=ggplot(df,aes_string(x='Pseudotime',y='gene',fill='expression'))+
        geom_tile()+
        theme_classic()+
        theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.title=element_blank(),
              axis.text.x=element_blank())+
        scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0)

    return(p3)
}