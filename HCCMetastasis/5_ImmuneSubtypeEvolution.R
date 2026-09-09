library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
source('https://github.com/KaWingLee9/in_house_tools/blob/main/visulization/custom_fun.R')

### TME distance
# comparison of Primary-Primary, Primary-Metastasis, Metastasis-Metastasis
x=dist(cell_type_comp) %>% as.matrix()
x=x %>% reshape2::melt(varnames=c('sample1','sample2'),value.name='distance')
x[,'sample1_type']=x[,'sample1'] %>% gsub('[PS][0-9]+','',.) %>% gsub('[0-9]+','',.)
x[,'sample2_type']=x[,'sample2'] %>% gsub('[PS][0-9]+','',.) %>% gsub('[0-9]+','',.)
x[,'pair_type']=paste0(x[,'sample1_type'],'-',x[,'sample2_type'])
x[,'sample1_patient']=x[,'sample1'] %>% gsub('[MT][0-9]+','',.) # %>% gsub('[0-9]+','',.)
x[,'sample2_patient']=x[,'sample2'] %>% gsub('[MT][0-9]+','',.) # %>% gsub('[0-9]+','',.)
x=x %>% filter(sample1_patient==sample2_patient,pair_type!='M-T')

x_patient=x %>% group_by(sample1_patient,pair_type) %>% summarize(distance=mean(distance)) %>% data.frame(check.names=FALSE)
df_1=x_patient %>% filter(pair_type=='T-T') # %>% data.frame(row.names=1,check.names=FALSE)
colnames(df_1)=c('patient','p1','d1')
df_2=x_patient %>% filter(pair_type=='T-M') # %>% data.frame(row.names=1,check.names=FALSE)
colnames(df_2)=c('patient','p2','d2')
df_3=x_patient %>% filter(pair_type=='M-M') # %>% data.frame(row.names=1,check.names=FALSE)
colnames(df_3)=c('patient','p3','d3')
df=df_1 %>% dplyr::full_join(df_2,by='patient') %>% dplyr::full_join(df_3,by='patient')
df=na.omit(df)

p=ggplot(x_patient,aes(x=forcats::fct_relevel(pair_type,'T-T','T-M','M-M'),y=distance))+
    geom_boxplot()+
    geom_point()+
    theme_bw()+
    xlab('')+
    geom_segment(data=df,aes(x=p1,xend=p2,y=d1,yend=d2),color="#DDDDDD",alpha=0.5,linewidth=0.4)+
    geom_segment(data=df,aes(x=p2,xend=p3,y=d2,yend=d3),color="#DDDDDD",alpha=0.5,linewidth=0.4)+
    theme_classic()

df=x_patient %>% filter(pair_type %in% c('T-T','T-M')) %>% reshape2::dcast(sample1_patient~pair_type,value.var='distance')
t.test(x=df[,'T-T'],y=df[,'T-M'],paired=TRUE)
df=x_patient %>% filter(pair_type %in% c('M-M','T-M')) %>% reshape2::dcast(sample1_patient~pair_type,value.var='distance')
t.test(x=df[,'M-M'],y=df[,'T-M'],paired=TRUE)

# comparison of TME distance and ITH


### consistency between primary and metastatic tumor ~ dissemination mode of patients
df=sample_info %>% group_by(patient,sample_type,TME_type_number) %>% summarise(Num=n()) %>% data.frame()

transfer_df=matrix(0,nrow=6,ncol=6)
rownames(transfer_df)=as.character(1:6)
colnames(transfer_df)=as.character(1:6)
patients=unique(df[,'patient'])

for (pa in patients){
    
    df_tmp=df %>% filter(patient==pa)
    
    if ( ( ('primary' %in% df_tmp[,'sample_type']) & ('metastatic' %in% df_tmp[,'sample_type']) ) ){
        p=df_tmp[ df_tmp[,'sample_type']=='primary' , 'TME_type_number' ]
        m=df_tmp[ df_tmp[,'sample_type']=='metastatic' , 'TME_type_number' ]
        
        for (i in p){
            for (j in m){
                transfer_df[i,j]=transfer_df[i,j]+1
            }
        }
    }
}

# Fisher's exact test: immune subtypes composition between primary and metastatic tumor
# retain the patients with both of primary and metastatic tumor sample number larger than 2
pvals=sapply(unique(df[,'patient']),function(p){
    
    df_patient=df %>% filter(patient==p)
    df_patient=df_patient %>% reshape2::dcast(sample_type~TME_type_number,value.var='Num') %>% data.frame(row.names=1,check.names=FALSE)
    df_patient[is.na(df_patient)]=0

    if (nrow(df_patient)>=2){ # with both primary and metastatic samples

        # remove samples with less than 2 samples in primary or metastatic sites
        if ( sum( (df_patient %>% apply(1,sum))==1 ) !=0 ){
            pval=NA
        }

        # composition in primary and metastatic sites are totally sme
        else if (ncol(df_patient)==1){
            pval=1
        }

        # TME subtypes totally different
        else if (sum(apply(df_patient,2,function(x) {0 %in% x}))==ncol(df_patient)){
            pval=0
        }

        else {
            test_result=fisher.test(df_patient,alternative='two.sided')
            # test_result=chisq.test(df_patient)
            pval=test_result[['p.value']]
        }
        

    } else {
        pval=NA
    }
    
    return(pval)

})

clinical_info=read.csv('./clinical_data.csv',row.names=1,check.names=FALSE)
clinical_info[,'pval']=pvals[clinical_info[,'patient']]
clinical_info=clinical_info[,c('patient','TimeDivergence','ColonalType','PhyleticType','pval')]
clinical_info=na.omit(clinical_info)

# ~TimeDivergence
clinical_info=clinical_info %>% group_by(TimeDivergence) %>% mutate(Rate=1/n())
ggplot(clinical_info %>% filter(TimeDivergence %in% c('EarlyDivergence','LateDivergence')),
       aes(x=TimeDivergence,y=Rate,fill=sig))+
    geom_bar(stat='identity')

m=clinical_info %>% filter(TimeDivergence %in% c('EarlyDivergence','LateDivergence'))
a=m %>% filter(TimeDivergence=='EarlyDivergence',sig=='0') %>% nrow()
b=m %>% filter(TimeDivergence=='EarlyDivergence',sig=='1') %>% nrow()
c=m %>% filter(TimeDivergence=='LateDivergence',sig=='0') %>% nrow()
d=m %>% filter(TimeDivergence=='LateDivergence',sig=='1') %>% nrow()
fisher.test( matrix(c(a,b,c,d),nrow=2),alternative='less' )

# ~ColonalType
clinical_info=clinical_info %>% group_by(ColonalType) %>% mutate(Rate=1/n())
ggplot(clinical_info %>% filter(ColonalType %in% c('monoclonal','polyclonal')),
       aes(x=ColonalType,y=Rate,fill=sig))+
    geom_bar(stat='identity')

m=clinical_info %>% filter(ColonalType %in% c('monoclonal','polyclonal'))
a=m %>% filter(ColonalType=='monoclonal',sig=='0') %>% nrow()
b=m %>% filter(ColonalType=='monoclonal',sig=='1') %>% nrow()
c=m %>% filter(ColonalType=='polyclonal',sig=='0') %>% nrow()
d=m %>% filter(ColonalType=='polyclonal',sig=='1') %>% nrow()
fisher.test( matrix(c(a,b,c,d),nrow=2),alternative='greater' )

# ~PhyleticType
clinical_info=clinical_info %>% group_by(PhyleticType) %>% mutate(Rate=1/n())
ggplot(clinical_info %>% filter(PhyleticType %in% c('monophyletic','polyphyletic')),
       aes(x=PhyleticType,y=Rate,fill=sig))+
    geom_bar(stat='identity')

m=clinical_info %>% filter(PhyleticType %in% c('monophyletic','polyphyletic'))
a=m %>% filter(PhyleticType=='monophyletic',sig=='0') %>% nrow()
b=m %>% filter(PhyleticType=='monophyletic',sig=='1') %>% nrow()
c=m %>% filter(PhyleticType=='polyphyletic',sig=='0') %>% nrow()
d=m %>% filter(PhyleticType=='polyphyletic',sig=='1') %>% nrow()
fisher.test( matrix(c(a,b,c,d),nrow=2),alternative='less' )

# Immune subtype composition between primary and metastatic tumor
df=sample_info %>% group_by(patient,sample_type) %>% mutate(Rate=1/n()) %>% data.frame()
df[,'sample_type']=factor(df[,'sample_type'],levels=c('primary','metastatic'))
df[,'sig']= as.character(clinical_info[ df[,'patient'] , 'sig' ])
df=na.omit(df)

ggplot(df,aes(x=sample_type,y=Rate,fill=TME_type_number))+
    geom_bar(stat='identity',color='black')+
    # facet_grid(sig~patient)+
    facet_wrap(sig~patient)+
    theme_bw()+
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),
          strip.background=element_blank(),# strip.text=element_blank()
          )

### transition probability matrix of immune subtypes
comp_df=df %>%
  group_by(patient, sample_type) %>%
  mutate(total = sum(Num)) %>%
  ungroup() %>%
  mutate(prob = Num / total)
comp_wide=comp_df %>%
  select(patient, sample_type, TME_type_number, prob) %>%
  pivot_wider(
    names_from = TME_type_number,
    values_from = prob,
    values_fill = 0
  )

P_mat=comp_wide %>%
  filter(sample_type == "primary") %>%
  arrange(patient)
M_mat=comp_wide %>%
  filter(sample_type == "metastatic") %>%
  arrange(patient)
P_mat=data.frame(P_mat)
M_mat=data.frame(M_mat)
rownames(P_mat)=P_mat[,'patient']
rownames(M_mat)=M_mat[,'patient']

P_mat=P_mat[intersect(rownames(P_mat),rownames(M_mat)),]
M_mat=M_mat[intersect(rownames(P_mat),rownames(M_mat)),]
P=as.matrix(P_mat[,-c(1,2)])
M=as.matrix(M_mat[,-c(1,2)])

# M=PT, resolve T with OLS
T_est=solve(t(P) %*% P) %*% t(P) %*% M
Heatmap(T_est,cluster_rows=FALSE,cluster_columns=FALSE,row_names_side='left',
        col=circlize::colorRamp2(seq(0,0.6,length=7),c('#F4FAED','#D6EFD0','#B2E1B9','#77CAC5','#42A6CB','#1373B2','#084384')
                                ))

### Immune subtypes in primary tumor ~ tumor dissemination mode
sample_clustering=read.csv('../cell_type_comp_scale_clustering_pearson_6.csv',row.names=1,check.names=FALSE)
cell_type_comp[,'TME_type_number']=sample_clustering[cell_type_comp[,'sample'],'c1']

cell_type_comp_tmp=cell_type_comp %>% filter(metastatic_site=='primary',TimeDivergence %in% c('EarlyDivergence','LateDivergence'))
p1=ContingencyPlot(x=cell_type_comp_tmp[ , 'TimeDivergence'],
               y=cell_type_comp_tmp[ , 'TME_type_number'],return_plot=FALSE)
p1['group']='Divergence'

cell_type_comp_tmp=cell_type_comp %>% filter(metastatic_site=='primary',ColonalType %in% c('monoclonal','polyclonal'))
p2=ContingencyPlot(x=cell_type_comp_tmp[ , 'ColonalType'],
               y=cell_type_comp_tmp[ , 'TME_type_number'],return_plot=FALSE)
p2['group']='Clonal'

cell_type_comp_tmp=cell_type_comp %>% filter(metastatic_site=='primary',PhyleticType %in% c('monophyletic','polyphyletic'))
p3=ContingencyPlot(x=cell_type_comp_tmp[ , 'PhyleticType'],
               y=cell_type_comp_tmp[ , 'TME_type_number'],return_plot=FALSE)
p3['group']='Phyletic'

test_result=Reduce(rbind,list(p1,p2,p3))
p=ggplot(data=test_result,aes(x=group_x,y=group_y,fill=OR,size=pval,color=sig))+
    geom_point(shape=21,stroke=1)+
    # scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=1)+
    scale_fill_gradientn(colors=c('blue','white','red'),
                         values=scales::rescale( c( min(test_result[,'OR']),1,max(test_result[,'OR']) ),to=c(0,1) ) )+
    scale_size(trans='reverse')+
    scale_size_continuous(range=c(6,0.1))+
    theme_bw()+
    theme(axis.title=element_blank())+
    scale_color_manual(values=c('TRUE'='black','FALSE'='#FFFFFF00'),na.value='#FFFFFF00',limits=c(TRUE,FALSE))+
    facet_grid(.~group,scale='free')
p=p+scale_size(range=c(8,0.001),breaks=c(0.05,0.2,0.5,1))+theme(axis.text.x=element_text(angle=45,hjust=1))
