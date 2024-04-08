### final file for survival and forest plot
library(glmnet)
library(survival)


###TCGA_LIHC####
LIHC_clinical=read.csv('survival/LIHC_survival.txt',sep='\t')
LIHC_clinical |> dim()
LIHC_clinical |> head()
source('function_defined.R')


### get the exp
exp_LIHC=get_TCGA_exp('TCGA-LIHC')
dataExpr1=as.data.frame(exp_LIHC)
dataExpr1 |> dim()
dataExpr1=dataExpr1 |> 
  colnames() |>
  substr(1,16) %>% 
  grepl('TCGA-..-....-0..',x=.) %>%
  dataExpr1[,.]
dataExpr1 |> dim()
library(dplyr)
tpm_LIHC=count_to_tpm(dataExpr1)
### filter low exp
tpm_LIHC[1:5,1:5]
index=filter_low_exp(dataExpr1)
tpm_LIHC=as.data.frame(tpm_LIHC)
tpm_LIHC=tpm_LIHC[index,]


### read the 23 genes
genelists=read.csv('intertested_genes_23.csv')
final_exp=tpm_LIHC[genelists$x,]
final_exp |> dim()
final_exp

### check the IDs if duplicated
IDS=final_exp |> 
  colnames() %>%
  substr(.,1,15) 
IDS |> duplicated() |> any()
final_exp=final_exp[,!duplicated(IDS)]
final_exp |> dim()
IDS %in% LIHC_clinical$sample

LIHC_clinical1=LIHC_clinical |>
  dplyr::filter(sample %in% IDS)

### remove the NA
index=LIHC_clinical1$OS.time |>
  is.na() |>
  which()
LIHC_clinical1=LIHC_clinical1[-index,]
LIHC_clinical1=LIHC_clinical1[LIHC_clinical1$OS.time!=0,]
#### match the order of the samples
final_exp=t(final_exp) |> as.data.frame()
library(magrittr)
rownames(final_exp) %<>%
  substr(.,1,15)
final_exp=final_exp[LIHC_clinical1$sample,]
final_exp |> head()
### check the sample order
all(rownames(final_exp)==LIHC_clinical1$sample) 
y=data.matrix(Surv(time=LIHC_clinical1$OS.time,
                   event = LIHC_clinical1$OS))
#y=data.matrix(Surv(time=LIHC_clinical1$PFI.time,
 #                    event = LIHC_clinical1$PFI))
y |> head()
log_final_exp=log2(final_exp+1)
#final_exp1=final_exp[,c('KPNA2','AURKA','PTTG1','FAM83D')]
fit=cv.glmnet(x=as.matrix(log_final_exp),
              y,
              family = 'cox',
              alpha = 1,
              nfolds = 10)

fit
pdf('result/figure/lasso_cox_1_111.pdf')
par(cex.axis = 2.0, cex.lab = 2.3, cex.main = 3.7)
plot(fit)

fit

plot_figure1=function(cvfits){
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
  print(length(filtered_genes))
  filtered_genes=c('CDKN3','KIF2C','KPNA2',
                   'FAM83D','DLGAP5','NCAPH',
                   'STIL','TOP2A','CENPE')
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
    ylab('Partial Likelihood Deviance')+ 
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
p1=plot_figure1(fit)
p11=p1$x
p2=plot_figure2(fit)
p11+p2
ggsave('result/figure/diag_1_22.pdf',height = 4,width = 10)
dev.off()
pdf('result/figure/lasso_cox1_11.pdf')
par(cex.axis = 2.0, cex.lab = 2.3, cex.main = 3.7)
plot(fit$glmnet.fit)
dev.off()

coefs=coef(fit,s=fit$lambda.min)
coefs
index=which(as.numeric(coefs)!=0)
index
coefs_retain=as.numeric(coefs)[index]
sig_gene_mult_cox=rownames(coefs)[index]

lasso_cox=log_final_exp[,sig_gene_mult_cox]
lasso_cox |> head()
filter=LIHC_clinical1 |> 
  dplyr::select(c('OS','OS.time'))
#filter=LIHC_clinical1 |> 
  #dplyr::select(c('PFI','PFI.time'))

data_cox=cbind(filter,lasso_cox)
data_cox |> head()
multicox=coxph(Surv(OS.time,OS)~.,data=data_cox)
#multicox=coxph(Surv(PFI.time,PFI)~.,data=data_cox)
multicox$means

library(broom)
multi_tidy=tidy(multicox)
index=p.adjust(multi_tidy$p.value)<0.05
new_lasso_cox=lasso_cox[,index] 
new_lasso_cox |> head()
new_data_cox=cbind(filter,new_lasso_cox)

new_multicox=coxph(Surv(OS.time,OS)~.,data=new_data_cox)
new_multicox

new_data_cox=new_data_cox[,-ncol(new_data_cox)]
new_multicox=coxph(Surv(OS.time,OS)~.,data=new_data_cox)
new_multicox
final_multi=tidy(new_multicox)
library(survminer)
ggforest(new_multicox, #直接用前面多因素cox回归分析的结果
         main = "Hazard ratio",
         cpositions = c(0.02,-0.15, 0.25), #前三列的位置，第二列是样品数，设了个负值，相当于隐藏了
         fontsize = 0.8, #字体大小
         refLabel = "reference", 
         noDigits = 2) 
write.csv(final_multi,'result/figure/for_foreast.csv')

### plot forestplt
mul_cox1=summary(new_multicox)
colnames(mul_cox1$conf.int)
multi1=as.data.frame(round(mul_cox1$conf.int[,c(1,3,4)],2))
multi1
library(forestplot)
library(stringr)
library(survival)
library(tableone)
multi2=ShowRegTable(new_multicox,
             exp=T,
             digits=2,
             pDigits = 3,
             printToggle = T,
             quote = F,
             ciFun = confint)

multi1
result=cbind(multi1,multi2)
result
result<-tibble::rownames_to_column(result, var = "Genes")
result[,c(1,5,6)]

fig1<- forestplot(result[,c(1,5,6)], #告诉函数，合成的表格result的第1，5，6列还是显示数字
                  mean=result[,2],   #告诉函数，表格第2列为HR，它要变成森林图的小方块
                  lower=result[,3],  #告诉函数表格第3列为5%CI，
                  upper=result[,4],  #表格第5列为95%CI，它俩要化作线段，穿过方块
                  zero=1,            #告诉函数，零线或参考线为HR=1即x轴的垂直线
                  boxsize=0.3,       #设置小黑块的大小
                  graph.pos=2,
                  title='Multivariate analysis')       #森林图应插在图形第2列
fig1


riskscore=predict(new_multicox,type = 'risk',
                  newdata = new_data_cox)
riskscores=as.data.frame(riskscore)
rownames(riskscores)=rownames(final_exp)
rownames(new_data_cox)=rownames(final_exp)



rownames(LIHC_clinical1)=LIHC_clinical1$sample
riskScore_cli=cbind(LIHC_clinical1,riskscores)
riskScore_cli
#按照中位数分为高低风险两组
riskScore_cli$riskScore2 <- ifelse(riskScore_cli$riskscore > median(riskScore_cli$riskscore),
                                   "High","Low")

riskScore_cli 

data.pca=prcomp(new_data_cox[3:5],scale. = T)
summary(data.pca)
riskScore_cli
pca.res=data.pca$x
all(rownames(pca.res)==rownames(riskScore_cli))

pca.res=cbind(pca.res,as.data.frame(riskScore_cli$riskScore2))
colnames(pca.res)[4]='group'
pca.var=data.pca$sdev^2 %>% as.data.frame()
pca.var$var=round(pca.var$./sum(pca.var)*100,2)
pca.var$pc=colnames(pca.res)[1:(ncol(pca.res)-1)]
pca.var
library(ggplot2)
ggplot(pca.res, aes(PC1, PC2, color = group))+ 
  # 选择X轴Y轴并映射颜色和形状
  geom_point(size = 4)+ # 画散点图并设置大小
  theme_bw() + # 加上边框
  # 自动提取主成分解释度进行绘图
  labs(x = paste('PC1(', pca.var$var[1],'%)', sep = ''),
       y = paste('PC2(', pca.var$var[2],'%)', sep = '')) +
  theme(
    panel.grid = element_blank(),
    legend.position = 'right')+ggtitle('TCGA-LIHC')+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,size=20),
        axis.text=element_text(size=19),
        axis.title =element_text(size=21),
        legend.text =element_text(size=19) ,
        legend.title = element_text(size=20))+
  stat_ellipse(level = 0.95,data=pca.res,
               aes(fill=group),
               alpha=0.1,
               geom='polygon',
               size=0.8)
ggsave('result/figure/pca_LIHC_1_9.pdf',height =  5,width = 6)
library('Rtsne')
data_cox_tsne=Rtsne(new_data_cox[3:5])
data_cox_tsne_res=as.data.frame(data_cox_tsne$Y)
colnames(data_cox_tsne_res)=c('tSNE1','tSNE2')
data_cox_tsne_res$group=paste0(pca.res$group,'_risk')

ggplot(data_cox_tsne_res,aes(tSNE1,tSNE2,color=group))+
  geom_point()+theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,size=22),
        axis.text=element_text(size=20),
        axis.title =element_text(size=22),
        legend.text =element_text(size=20) ,
        legend.title = element_text(size=21))+
  stat_ellipse(level = 0.95,data=data_cox_tsne_res,
               aes(fill=group),
               alpha=0.1,
               geom='polygon',
               size=0.8)+
  ggtitle('TCGA-LIHC')
ggsave('result/figure/tsne_LIHC.pdf',height = 6,width = 8)
#KM分析
riskScore_cli$OS.time1=riskScore_cli$OS.time/30
riskScore_cli$RiskScore=riskScore_cli$riskScore2

unicox=riskScore_cli |> dplyr::select(c(OS.time,OS,riskscore))
coxph(Surv(OS.time,OS)~.,data=unicox)
LIHC_clinical1 |> head()
all(rownames(unicox)==rownames(LIHC_clinical1))
LIHC_clinical1 |> head()

index=riskScore_cli$OS.time1<108
riskScore_cli_modified=riskScore_cli[index,]
fit <- survfit(Surv(OS.time1, as.numeric(OS)) ~ RiskScore, data=riskScore_cli)

library("survminer")
lasso_KM <- ggsurvplot(fit, data = riskScore_cli,
                       pval = T,
                       linetype = 1,
                       risk.table = F,
                       surv.median.line = "hv",
                       legend.title=" ",
                       title="Overall survival(TCGA_LIHC)", #标题
                       ylab="Survival probability",xlab = " Time (Months)" ,
                       break.x.by = 12,
                       censor.size = 6,
                       font.size = 50,
                       font.x=55,
                       font.y=55,
                       font.legend=38,
                       font.main=c(55),
                       font.tickslab=c(45),
                       conf.int = F,
                       Legend.labs=NULL)

lasso_KM$plot+theme(legend.position = 'none')+
  theme(plot.title = element_text(hjust = 0.5))
ggsave('result/figure/survival_TCGA_LIHC_1_9.pdf',
       height =11,width = 11)  

library(timeROC)
with(riskScore_cli,
    ROC_riskscore <<- timeROC(T = OS.time,
                              delta = OS,
                              marker = riskscore,
                              cause = 1,
                              weighting = "marginal",
                              times = c(365,730,1080),
                              ROC = TRUE,
                              iid = TRUE)
)
ROC_riskscore

dat = data.frame(fpr = as.numeric(ROC_riskscore$FP),
                tpr = as.numeric(ROC_riskscore$TP),
                time = rep(as.factor(c(365,730,1080)),each = nrow(ROC_riskscore$TP)))

library(ggplot2)
ggplot() + 
 geom_line(data = dat,aes(x = fpr, y = tpr,color = time),size = 1) + 
 scale_color_manual(name = NULL,values = c("#92C5DE", "#F4A582", "#66C2A5"),
                    labels = paste0("AUC of ",c(1,2,3),"-y survival: ",
                                    format(round(ROC_riskscore$AUC,3),nsmall = 2)))+
 geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey",linetype='dotted')+
 theme_bw()+
 theme(panel.grid = element_blank(),
       legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
       legend.position = c(0.765,0.125))+
 scale_x_continuous(expand = c(0.005,0.005))+
 scale_y_continuous(expand = c(0.005,0.005))+
 labs(x = "1 - Specificity",
      y = "Sensitivity")+
 theme(axis.text = element_text(size = 34),
       axis.title = element_text(size = 36),
       legend.text =element_text(size = 30),
       title = element_text(size=30),
       panel.grid = element_blank())+
 coord_fixed()+
 ggtitle('timeROC')


ggsave('result/figure/time_ROC_TCGA_LIHC_1_9.pdf',height = 8,width = 8)










###ICGC_LIRI_JP####
### extract the exp and clinical info
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
jp_exp=tpm_count[genelists$x,]
jp_exp
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
dev.new()



pc<- ggroc(roc_obj, legacy.axes = TRUE,
          color='red',
          size=1)+
 geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype=4)+
 ggtitle('ROC')+ theme_bw()+
 annotate("text",x=0.75,y=0.375,label=paste("AUC = ", round(roc_obj$auc,3)),size=7)+
 theme(axis.text = element_text(size = 24),
       axis.title = element_text(size = 26),
       legend.text =element_text(size = 22),
       title = element_text(size=20),
       panel.grid = element_blank())
pc
ggsave('result/figure/ICTC_jp_all_12_31.pdf',height = 5,width = 5)
pc
### get y variables
sample_barcode |> colnames()

tumor_barcode=sample_barcode |> 
 dplyr::select(c('icgc_donor_id','specimen_type'))

### select hte primary-solid tissues
index=tumor_barcode$specimen_type=="Primary tumour - solid tissue"
tumor_barcode=tumor_barcode[index,]
jp_exp |> colnames()
jp_exp[,tumor_barcode$icgc_donor_id]
colnames(jp_exp) %in% tumor_barcode$icgc_donor_id

#### normal run
index=jp_exp |> 
 colnames() %>%
 grepl('normal',.)
index
jp_exp1=jp_exp[!index]
jp_exp1 |> dim()
ICGC_clinic |> colnames()
ICGC_clinic1=ICGC_clinic |> 
 dplyr::select(c("icgc_donor_id","donor_vital_status","donor_survival_time"))
rownames(ICGC_clinic1)=ICGC_clinic1$icgc_donor_id
ICGC_clinic1=ICGC_clinic1[colnames(jp_exp1),] 
all(rownames(ICGC_clinic1)==colnames(jp_exp1))
ICGC_clinic1$OS_event=if_else(ICGC_clinic1$donor_vital_status=='alive',0,1)
### convert the data formt for validation
ICGC_clinic1

y_validate=data.matrix(Surv(time=ICGC_clinic1$donor_survival_time,
                           event = ICGC_clinic1$OS_event))

jp_exp1=jp_exp1 |> t()
jp_exp1 |> dim()
sig_gene_mult_cox
jp_exp1=jp_exp1[,c('KIF2C','KPNA2','CDKN3')]
jp_exp1_log=log2(jp_exp1+1)
new_multicox

ICGC_clinic2=ICGC_clinic1 |> 
 dplyr::select("OS_event","donor_survival_time") |>
 dplyr::rename(OS=OS_event) |>
 dplyr::rename(OS.time=donor_survival_time)




data_cox1=cbind(ICGC_clinic2,jp_exp1_log)
data_cox1 |> dim()



sink('result/survival_univariate_info.txt',append = T)
for (i in 3:ncol(data_cox)){
 index_vector_LIHC=c(1,2)
 index_vector_LIHC=c(index_vector_LIHC,i)
 final_data22=data_cox[,index_vector_LIHC]
 glm_logistic=coxph(Surv(OS.time,OS)~.,data=final_data22)
 print(summary(glm_logistic))
}
sink()

ggforest(multicox,
        data=data_cox,
        main = "Hazard ratio",        # 标题
        cpositions = c(0.01, 0.1, 0.3), #左边三列的相对距离
        fontsize = 0.8, # 字体大小
        refLabel = "reference", #显示因子的参考水平
        noDigits = 3)
riskscore=predict(new_multicox,type = 'risk',
                 newdata = data_cox1)

riskscores=as.data.frame(riskscore)
rownames(riskscores)=rownames(jp_exp1)
riskscores


riskScore_cli=cbind(ICGC_clinic2,riskscores)
riskScore_cli
riskScore_cli$riskScore2 <- ifelse(riskScore_cli$riskscore > median(riskScore_cli$riskscore),
                                   "High","Low")


riskScore_cli$OS.time1=riskScore_cli$OS.time/30
library(ggplot2)
install.packages('ProComp')
data_cox1

data.pca=prcomp(data_cox1[3:5],scale. = T)
summary(data.pca)

pca.res=data.pca$x
all(rownames(pca.res)==rownames(riskScore_cli))
pca.res=cbind(pca.res,as.data.frame(riskScore_cli$riskScore2))
colnames(pca.res)[4]='group'
pca.var=data.pca$sdev^2 %>% as.data.frame()
pca.var$var=round(pca.var$./sum(pca.var)*100,2)
pca.var$pc=colnames(pca.res)[1:(ncol(pca.res)-1)]
pca.var
ggplot(pca.res, aes(PC1, PC2, color = group))+ 
 # 选择X轴Y轴并映射颜色和形状
 geom_point(size = 4)+ # 画散点图并设置大小
 theme_bw() + # 加上边框
 # 自动提取主成分解释度进行绘图
 labs(x = paste('PC1(', pca.var$var[1],'%)', sep = ''),
      y = paste('PC2(', pca.var$var[2],'%)', sep = '')) +
  ggtitle('ICGC_LIRI-JP')+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5,size=20),
        axis.text=element_text(size=19),
        axis.title =element_text(size=21),
        legend.text =element_text(size=19) ,
        legend.title = element_text(size=20))+
  stat_ellipse(level = 0.95,data=pca.res,
               aes(fill=group),
               alpha=0.1,
               geom='polygon',
               size=0.8)
# 设置图例位置，此处为相对位置
ggsave('result/figure/pca_validation.pdf',height = 5,width = 6)



#按照中位数分为高低风险两组
riskScore_cli$riskScore2 <- ifelse(riskScore_cli$riskscore > median(riskScore_cli$riskscore),
                                  "High","Low")


riskScore_cli$OS.time1=riskScore_cli$OS.time/30



fit <- survfit(Surv(OS.time1, as.numeric(OS)) ~ riskScore2, data=riskScore_cli)
lasso_KM <- ggsurvplot(fit, data = riskScore_cli,
                      pval = T,
                      linetype = 1,
                      risk.table = F,
                      surv.median.line = "hv",
                      legend.title=" ",
                      title="Overall survival(ICGC_LIRI-JP)", #标题
                      ylab="Survival probability",xlab = " Time (Months)" ,
                      break.x.by = 12,
                      censor.size = 5,
                      font.size = 18,
                      font.x=32,
                      font.y=32,
                      font.legend=30,
                      font.main=c(32),
                      font.tickslab=c(32),
                      conf.int = F,
                      Legend.labs=NULL)

lasso_KM$plot+theme(legend.position = 'none')+
 theme(plot.title = element_text(hjust = 0.5))
ggsave('result/figure/validate_1_10.pdf',height = 8,width = 8)




library(timeROC)
with(riskScore_cli,
    ROC_riskscore <<- timeROC(T = OS.time,
                              delta = OS,
                              marker = riskscore,
                              cause = 1,
                              weighting = "marginal",
                              times = c(365,730,1080),
                              ROC = TRUE,
                              iid = TRUE))
ROC_riskscore
dev.new()

dat = data.frame(fpr = as.numeric(ROC_riskscore$FP),
                tpr = as.numeric(ROC_riskscore$TP),
                time = rep(as.factor(c(365,730,1080)),each = nrow(ROC_riskscore$TP)))

library(ggplot2)
ggplot() + 
 geom_line(data = dat,aes(x = fpr, y = tpr,color = time),size = 1) + 
 scale_color_manual(name = NULL,values = c("#92C5DE", "#F4A582", "#66C2A5"),
                    labels = paste0("AUC of ",c(1,2,3),"-y survival: ",
                                    format(round(ROC_riskscore$AUC,3),nsmall = 2)))+
 geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey",linetype='dotted')+
 theme_bw()+
 theme(panel.grid = element_blank(),
       legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
       legend.position = c(0.765,0.125))+
 scale_x_continuous(expand = c(0.005,0.005))+
 scale_y_continuous(expand = c(0.005,0.005))+
 labs(x = "1 - Specificity",
      y = "Sensitivity")+
  theme(axis.text = element_text(size = 34),
        axis.title = element_text(size = 36),
        legend.text =element_text(size = 30),
        title = element_text(size=30),
        panel.grid = element_blank())+
 coord_fixed()+
 ggtitle('timeROC(ICGC_LIRI-JP)')


ggsave('result/figure/time_ROC_1_10.pdf',height = 8,width = 8)
ROC_riskscore


