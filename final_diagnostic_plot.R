df1=read.csv('LIHC_for_python.csv')
df1 |> head()
df2=read.csv('Ps_for_python.csv')
df2 |> head()
library(glmnet)
fit_function=function(type='liver'){
  if (type=='liver'){
    df1=read.csv('LIHC_for_python.csv')
  }
  else{
    df1=read.csv('Ps_for_python.csv')
  }
  input=df1[,-c(1,ncol(df1))]
  input1=input |> as.matrix()
  label1=df1[,ncol(df1)]
  cvfit=cv.glmnet(input1,label1,family='binomial',
                  type.measure = 'deviance',nfolds = 10)
  return(cvfit)
}
plot_figure1=function(cvfits,type='liver'){
  ### grep the fit 
  x=coef(cvfits$glmnet.fit)
  tmp=as.data.frame(as.matrix(x))
  ## remove the intercept columns
  tmp=tmp[-1,]
  ### creat the gene column
  tmp$coef=rownames(tmp)
  ### wid to long data
  tmp=reshape::melt(tmp,id='coef')
  ### s1 to1
  tmp$variable=as.numeric(gsub('s','',tmp$variable))
  tmp$coef <- gsub('_','-',tmp$coef) 
  
  ### grep the lambda acrording to index
  tmp$lambda=cvfits$lambda[tmp$variable+1]
  nozero=coef(cvfits,s=cvfits$lambda.min)
  index=as.matrix(nozero)
  index=index[-1,]
  filtered_genes=index[index!=0] |> names()
  if(type=='liver'){
  filtered_genes=c('KPNA2','FAM83D','AURKA',
          'PTTG1', 'CCNB1','BIRC5',
          'CDC6','CDKN3')}
  else{
  filtered_genes=c('KPNA2','FAM83D','AURKA',
           'PTTG1','BUB1', 'CENPF',
           'TOP2A','RRM2','KIF23')}
  tmp$Selected=tmp$coef %in% filtered_genes
  
  p1=ggplot(tmp,aes(log(lambda),value,color = coef,linetype=Selected)) + 
    geom_vline(xintercept = log(cvfits$lambda.min),
               size=0.8,color='grey60',
               alpha=0.8,linetype=2)+
    geom_line(size=0.8)+
    xlab("Log lambda") + 
    ylab('Coefficients')+ 
    theme_bw()+ 
    scale_color_manual(values= c(pal_npg()(10),
                                 pal_d3()(10),
                                 pal_jco()(10),
                                 pal_lancet()(10),
                                 pal_aaas()(10),
                                 pal_simpsons()(10),
                                 pal_gsea()(10),
                                 pal_jama()(10)),
                       breaks = filtered_genes,
                       name='Selected_genes')+
    scale_x_continuous(expand = c(0.01,0.01))+ 
    scale_y_continuous(expand = c(0.01,0.01))+ 
    theme(panel.grid = element_blank(), 
          axis.title = element_text(size=16), 
          axis.text = element_text(size=15), 
          legend.text = element_text(size=13), 
          legend.position = 'right')+ 
    guides(col=guide_legend(ncol = 1))
  res=list(x=p1,y=filtered_genes)
  return(res)
}
plot_figure2=function(cvfits){
  xx <- data.frame(lambda=cvfits[["lambda"]],
                   cvm=cvfits[["cvm"]],
                   cvsd=cvfits[["cvsd"]], 
                   cvup=cvfits[["cvup"]],
                   cvlo=cvfits[["cvlo"]],
                   nozezo=cvfits[["nzero"]]) 
  xx$ll<- log(xx$lambda) 
  xx$NZERO<- paste0(xx$nozezo,' vars')
  p2=ggplot(xx,aes(ll,cvm,color=NZERO))+ 
    geom_errorbar(aes(x=ll,ymin=cvlo,ymax=cvup),
                  width=0.05,size=1,color='grey')+ 
    geom_vline(xintercept = xx$ll[which.min(xx$cvm)],
               size=0.8,color='grey60',alpha=0.8,
               linetype=2)+ 
    scale_x_continuous(labels = function(x) sprintf("%.1f", x)
    )+ 
    geom_point(size=2,color='#DB3124')+ 
    xlab("Log Lambda")+
    ylab('Binomial Deviance')+ 
    theme_bw()+ 
    scale_color_manual(values= c(pal_npg()(10),
                                 pal_d3()(10),
                                 pal_lancet()(10),
                                 pal_aaas()(10)))+
    theme(panel.grid = element_blank(), 
          axis.title = element_text(size=16), 
          axis.text = element_text(size=15), 
          legend.title = element_blank(), 
          legend.text = element_text(size=13), 
          legend.position = 'bottom')+
    guides(col=guide_legend(ncol = 1))
  return(p2)
}
plot_figure3=function(index_genes,type='liver'){
  ### remove the first column
  if (type=='liver'){
    df=read.csv('LIHC_for_python.csv')
  }
  else{
    df=read.csv('Ps_for_python.csv')
  }
  df=df[,-1]
  df1=df[,-ncol(df)]
  label1=df[,ncol(df)]
  
  #### then plot the ridges_for the other

  target_df=df1
  df_select=target_df[,index_genes]
  
  ### rename
  label_select=label1
  plot_density=cbind(df_select,label_select)
  plot_density=plot_density |> as.data.frame()
  ### rename label_select to sample_type
  plot_density=plot_density |>
    dplyr::rename(sample_type=label_select)
  #### discrete
  if(type=='liver'){
  plot_density$sample_type1=if_else(plot_density$sample_type=='1','Tumor','Normal')
  }
  else{
    plot_density$sample_type1=if_else(plot_density$sample_type=='1','Ps','Normal')
  }
  ### remove the non-numeric col
  plot_density=plot_density |>
    dplyr::select(-c(sample_type))
  ### extract the gene cols
  select_cols=colnames(plot_density)[1:ncol(plot_density)-1]
  ### wider to longer data format
  plot_density1=plot_density %>%
    gather(key='Gene',value = 'Value',all_of(select_cols))

  library(ggridges)
  ### relevl the factor
  plot_density1$sample_type1=plot_density1$sample_type1 |> as.factor()
  if (type=='liver'){
  plot_density1$sample_type1=relevel(plot_density1$sample_type1,ref = 'Tumor')
  }else{
    plot_density1$sample_type1=relevel(plot_density1$sample_type1,ref = 'Ps')
  }
  index1=rev(index_genes)
  plot_density1$Gene=factor(plot_density1$Gene,levels=index1)
  
  p3=ggplot(plot_density1,aes(x=Value,y=Gene,fill=sample_type1))+
    geom_density_ridges(alpha=0.7)+
    theme_ridges()+
    scale_fill_manual(values=c('#DB3124','#90BEE0'))+
    theme(legend.position = "top",
          legend.title =element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size=16), 
          axis.text = element_text(size=15),
          legend.text = element_text(size=13))+
    labs(x='Expression:log2(TPM+1)')
  return(p3)
}
plot_figure4=function(cvfits,type='liver'){
  if (type=='liver'){
  jp_exp=read.csv('validation_icgc_exp.csv')
  rownames(jp_exp)=jp_exp$X
  jp_exp=jp_exp[,-1]
  y_label_jp=jp_exp |> 
    colnames() %>%
    grepl('normal',.)%>%
    ifelse(.,0,1)
  test_set2=jp_exp |> t()
  test_set2_log=log2(test_set2+1)

  predict_test=predict(cvfits,newx = test_set2_log,s='lambda.min',
                       type='response')
  
  label=y_label_jp |> 
    as.matrix() 
  }
  else{
    ps_val=read.csv('validate_ps.csv')
    ps_val=ps_val[,-1]
    input_data=ps_val[,-ncol(ps_val)]
    input_datax=as.matrix(input_data)
    label=ps_val[,ncol(ps_val)]
    predict_test=predict(cvfits,newx = input_datax,s='lambda.min',
                         type='response')
  }
    
  
  library(pROC)
  roc_obj=roc(label,predict_test)
  p4<- ggroc(roc_obj, legacy.axes = TRUE,
             color='#DB3124',
             size=1)+
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype=4)+
    ggtitle('ROC')+ theme_bw()+
    annotate("text",x=0.75,y=0.375,label=paste("AUC = ", round(roc_obj$auc,3)),size=7)+
    theme(panel.grid = element_blank(), 
          axis.title = element_text(size=16), 
          axis.text = element_text(size=15), 
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size=13))
  return(p4)
}


cvfit=fit_function()
p_return=plot_figure1(cvfit)
p1=p_return$x
index_genes=p_return$y
p2=plot_figure2(cvfit)
p3=plot_figure3(index_genes)
p4=plot_figure4(cvfit)



cvfit1=fit_function(type = 'Ps')
p_return1=plot_figure1(cvfit1,type='Ps')
p5=p_return1$x
index_genes1=p_return1$y
p6=plot_figure2(cvfit1)
p7=plot_figure3(index_genes1,type = 'Ps')
p8=plot_figure4(cvfit1,type='Ps')


p1/p2/p3/p4/p5/p6/p7/p8 + 
  plot_layout(ncol=4)+
  plot_annotation(tag_levels = 'A')
ggsave('result/figure/comb.pdf',height = 8,width = 14)



df1=df2
input=df1[,-c(1,ncol(df1))]
input1=input |> as.matrix()
label1=df1[,ncol(df1)]

cvfit=cv.glmnet(input1,label1,family='binomial',
                type.measure = 'deviance',nfolds = 10)
cvfit
plot(cvfit,col='#DB3124')
plot

pdf('result/figure/cv_lasso_liver.pdf')
par(cex.axis = 1.8, cex.lab = 2.1, cex.main = 3.5)
plot(cvfit,col='#DB3124')
dev.off()
pdf('result/figure/cv_lasso_ps1.pdf')
par(cex.axis = 1.8, cex.lab = 2.1, cex.main = 3.5)
plot(cvfit$glmnet.fit)
dev.off()
plot(cvfit$glmnet.fit)
cvfit$glmnet.fit
cvfit
### test
x=coef(cvfit$glmnet.fit)
x
tmp=as.data.frame(as.matrix(x))
## remove the intercept columns
tmp=tmp[-1,]

tmp$coef=rownames(tmp)

tmp=reshape::melt(tmp,id='coef')

tmp$variable=as.numeric(gsub('s','',tmp$variable))
tmp$variable
tmp$coef <- gsub('_','-',tmp$coef) 

tmp$lambda=cvfit$lambda[tmp$variable+1]
tmp$lambda
tmp$norm=apply(abs(x[-1,]),2,sum)[tmp$variable+1]
tmp$norm

tmp$Selected=tmp$coef %in% c('KPNA2','FAM83D',
                'AURKA','PTTG1',
                'CCNB1','BIRC5',
                'CDC6','CDKN3')

line_types=ifelse(tmp$selected,'dotted','solid')
tmp$line_size=ifelse(tmp$selected,1,0.5)
geom_line(aes(size=Selected)) 
  scale_size_discrete(range = c(0.6, 1))

p1=ggplot(tmp,aes(log(lambda),value,color = coef,linetype=Selected)) + 
  geom_vline(xintercept = log(cvfit$lambda.min),
             size=0.8,color='grey60',
             alpha=0.8,linetype=2)+
  geom_line(size=0.8)+
  xlab("Log lambda") + 
  ylab('Coefficients')+ 
  theme_bw()+ 
  scale_color_manual(values= c(pal_npg()(10),
                               pal_d3()(10),
                               pal_jco()(10),
                               pal_lancet()(10),
                               pal_aaas()(10),
                               pal_simpsons()(10),
                               pal_gsea()(10),
                               pal_jama()(10)),
                     breaks = c('KPNA2','FAM83D',
                                'AURKA','PTTG1',
                                'CCNB1','BIRC5',
                                'CDC6','CDKN3'),
                     name='Selected_genes')+
  scale_x_continuous(expand = c(0.01,0.01))+ 
  scale_y_continuous(expand = c(0.01,0.01))+ 
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size=16), 
        axis.text = element_text(size=15), 
        legend.text = element_text(size=13), 
        legend.position = 'right')+ 
  guides(col=guide_legend(ncol = 1))
ggsave('result/figure/1_22_test_lasso.pdf',height = 4,width =5.5)  

### plot next
xx <- data.frame(lambda=cvfit[["lambda"]],
                 cvm=cvfit[["cvm"]],
                 cvsd=cvfit[["cvsd"]], 
                 cvup=cvfit[["cvup"]],
                 cvlo=cvfit[["cvlo"]],
                 nozezo=cvfit[["nzero"]]) 
xx$ll<- log(xx$lambda) 
xx$NZERO<- paste0(xx$nozezo,' vars')


p2=ggplot(xx,aes(ll,cvm,color=NZERO))+ 
  geom_errorbar(aes(x=ll,ymin=cvlo,ymax=cvup),
                width=0.05,size=1,color='grey')+ 
  geom_vline(xintercept = xx$ll[which.min(xx$cvm)],
             size=0.8,color='grey60',alpha=0.8,
             linetype=2)+ 
  scale_x_continuous(labels = function(x) sprintf("%.1f", x)
                     )+ 
  geom_point(size=2,color='#DB3124')+ 
  xlab("Log Lambda")+
  ylab('Binomial Deviance')+ 
  theme_bw()+ 
  scale_color_manual(values= c(pal_npg()(10),
                               pal_d3()(10),
                               pal_lancet()(10),
                               pal_aaas()(10)))+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size=16
                                  ), 
        axis.text = element_text(size=15
                                 ), 
        legend.title = element_blank(), 
        legend.text = element_text(size=13), 
        legend.position = 'bottom')+
  guides(col=guide_legend(ncol = 1))
p1+p2
ggsave('result/figure/test_loglambda.pdf',height = 4,width = 6)

coef(cvfit,s=cvfit$lambda.min)

###
### plot ridge_plot
getwd()

df=read.csv('LIHC_for_python.csv')
### remove the first column
df=df[,-1]
log_df1=df[,-ncol(df)]
label1=df[,ncol(df)]

#### then plot the ridges_for the other
index=c('KPNA2','FAM83D','AURKA',
        'PTTG1', 'CCNB1','BIRC5',
        'CDC6','CDKN3')
index1=c('KPNA2','FAM83D','AURKA','PTTG1','BUB1',
         'CENPF','TOP2A','RRM2','KIF23')


log_df1=log2(log_df1+1)

target_df=log_df1
df_select=target_df[,index]

label_select=label1
plot_density=cbind(df_select,label_select)
plot_density=plot_density |> as.data.frame()
plot_density |>head()
plot_density=plot_density |>
  dplyr::rename(sample_type=label_select)

plot_density$sample_type1=if_else(plot_density$sample_type=='1','Tumor','Normal')

plot_density=plot_density |>
  dplyr::select(-c(sample_type))
colnames(plot_density)
select_cols=colnames(plot_density)[1:8]


plot_density1=plot_density %>%
  gather(key='Gene',value = 'Value',all_of(select_cols))
plot_density1
library(ggplot2)
library(ggridges)

plot_density1$sample_type1=plot_density1$sample_type1 |> as.factor()
plot_density1$sample_type1=relevel(plot_density1$sample_type1,ref = 'Tumor')

index1=rev(index)
plot_density1$Gene=factor(plot_density1$Gene,levels=index1)
p3=ggplot(plot_density1,aes(x=Value,y=Gene,fill=sample_type1))+
  geom_density_ridges(alpha=0.7)+
  theme_ridges()+
  scale_fill_manual(values=c('#DB3124','#90BEE0'))+
  theme(legend.position = "top",
        legend.title =element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size=16), 
        axis.text = element_text(size=15),
        legend.text = element_text(size=13))+
  labs(x='Expression:log2(TPM+1)')

p1+p2+p3+p4
p1/p2/p3/p4 + 
  plot_layout(ncol=4)+ 
  plot_annotation(tag_levels = 'A')+
  theme(plot.tag = element_text(size = 8))
ggsave('result/figure/ridges_hcc.pdf',height = 3,width = 3)
ggsave('result/figure/test_comb.pdf',height = 4,width = 14)



### the code for p4
ICGC_clinic=read.csv('survival/donor.tsv/donor.tsv',sep = '\t')

sample_barcode=read.csv('survival/specimen.tsv/specimen.tsv',sep='\t')
exp=read.csv('survival/exp_seq.tsv/exp_seq.tsv',sep ='\t')

exp |> head()
exp$analysis_id |> factor() |> levels()
index=grepl('Cancer',exp$analysis_id)
index1=grepl('Liver',exp$analysis_id)

exp_tumor=exp[index,]
exp_normal=exp[index1,]
library(dplyr)
library(tidyverse)
exp |> colnames()
exp1_tumor=exp_tumor |> dplyr::select(c('icgc_donor_id','gene_id','raw_read_count'))
exp1_normal=exp_normal |> dplyr::select(c('icgc_donor_id','gene_id','raw_read_count'))
exp1_tumor |> head()

#spread(exp1,gene_id,raw_read_count)
#exp1 |> head()

library(reshape2)
test_tumor=exp1_tumor |> 
  group_by(icgc_donor_id,gene_id) |>
  summarise(meanvalues=mean(raw_read_count))

test_tumor |> head()
raw_count_tumor=pivot_wider(test_tumor,names_from = gene_id,
                            values_from = meanvalues,
                            values_fill = 0)
raw_count_tumor=as.data.frame(raw_count_tumor)
#raw_count_tumor=dcast(exp1_tumor,icgc_donor_id~exp1_tumor$gene_id,
#               value.var = 'raw_read_count')

rownames(raw_count_tumor)=raw_count_tumor[,1]
raw_count_tumor=raw_count_tumor[,-1]
raw_count_tumor |> dim()
raw_count_tumor=raw_count_tumor |> 
  t() |>
  as.data.frame()


#raw_count_normal=dcast(exp1_normal,icgc_donor_id~exp1_normal$gene_id,
#                      value.var = 'raw_read_count')
exp1_normal=exp1_normal |> 
  group_by(icgc_donor_id,gene_id) |>
  summarise(meanvalues=mean(raw_read_count))

raw_count_normal=pivot_wider(exp1_normal,names_from = gene_id,
                             values_from = meanvalues,
                             values_fill = 0)
raw_count_normal=as.data.frame(raw_count_normal)
rownames(raw_count_normal)=raw_count_normal[,1]
raw_count_normal=raw_count_normal[,-1]
raw_count_normal |> dim()
raw_count_normal=raw_count_normal |> 
  t() |>
  as.data.frame()
colnames(raw_count_normal)=raw_count_normal |>
  colnames() %>%
  paste0(.,'normal')
raw_count_normal |> dim()
raw_count_tumor |> dim()
index=intersect(rownames(raw_count_normal),rownames(raw_count_tumor))
raw_count_normal=raw_count_normal[index,]
raw_count_tumor=raw_count_tumor[index,]

raw_count_tumor$gene_symbol=rownames(raw_count_tumor)
raw_count_normal$gene_symbol=rownames(raw_count_normal)
total_counts=left_join(raw_count_tumor,raw_count_normal,by="gene_symbol")
rownames(total_counts)=total_counts$gene_symbol
total_counts=total_counts |> 
  dplyr::select(-c('gene_symbol'))
total_counts |>
  class()
total_counts[1:4,1:4]
#total_counts[c('KPNA2','AURKA','PTTG1'),]
#index_order=order(rowMeans(total_counts),decreasing = T)
#total_counts[index_order,]



tpm_count=count_to_tpm(total_counts)
tpm_count=tpm_count |>
  as.data.frame()
tpm_count
genelists=read.csv('intertested_genes_23.csv')
jp_exp=tpm_count[genelists$x,]
jp_exp
write.csv(jp_exp,'validation_icgc_exp.csv')
y_label_jp=jp_exp |> 
  colnames() %>%
  grepl('normal',.)%>%
  ifelse(.,0,1)
test_set2=jp_exp |> t()
test_set2_log=log2(test_set2+1)

test_set2
predict_test=predict(cvfit,newx = test_set2_log,s='lambda.min',
                     type='response')
predict_test |> dim()
predict_test |> class()
label=y_label_jp
label=label |> 
  as.matrix() 
table(label)
library(pROC)
library(ggplot2)
roc_obj=roc(label,predict_test)
roc_obj

p4<- ggroc(roc_obj, legacy.axes = TRUE,
           color='#DB3124',
           size=1)+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype=4)+
  ggtitle('ROC')+ theme_bw()+
  annotate("text",x=0.75,y=0.375,label=paste("AUC = ", round(roc_obj$auc,3)),size=7)+
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size=16), 
        axis.text = element_text(size=15), 
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=13))
p1+p2+p3+p4
p1/p2/p3/p4 + 
  plot_layout(ncol=4)+ 
  plot_annotation(tag_levels = 'A')
ggsave('result/figure/test_comb.pdf',height = 4,width = 14)
