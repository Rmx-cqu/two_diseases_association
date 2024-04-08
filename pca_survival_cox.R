library(glmnet)
library(survival)


###TCGA_LIHC####
LIHC_clinical=read.csv('survival/LIHC_survival.txt',sep='\t')
LIHC_clinical |> dim()
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

### check the IDs 
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
y |> head()

#final_exp1=final_exp[,c('KPNA2','AURKA','PTTG1','FAM83D')]
fit=cv.glmnet(x=as.matrix(final_exp),
       y,
       family = 'cox',
       alpha = 1,
       nfolds = 20)

plot(fit)
coefs=coef(fit,s=fit$lambda.min)
coefs
index=which(as.numeric(coefs)!=0)
index
coefs_retain=as.numeric(coefs)[index]
sig_gene_mult_cox=rownames(coefs)[index]

lasso_cox=final_exp[,sig_gene_mult_cox]
lasso_cox |> head()
filter=LIHC_clinical1 |> 
  dplyr::select(c('OS','OS.time'))

data_cox=cbind(filter,lasso_cox)
data_cox |> head()
multicox=coxph(Surv(OS.time,OS)~.,data=data_cox)
multicox
riskscore=predict(multicox,type = 'risk',
        newdata = data_cox)
riskscores=as.data.frame(riskscore)
rownames(riskscores)=rownames(final_exp)
rownames(data_cox)=rownames(final_exp)



rownames(LIHC_clinical1)=LIHC_clinical1$sample
riskScore_cli=cbind(LIHC_clinical1,riskscores)

#按照中位数分为高低风险两组
riskScore_cli$riskScore2 <- ifelse(riskScore_cli$riskscore > median(riskScore_cli$riskscore),
                                   "High","Low")

data_cox
data.pca=prcomp(data_cox[3:6],scale. = T)
summary(data.pca)
riskScore_cli
pca.res=data.pca$x
all(rownames(pca.res)==rownames(riskScore_cli))

pca.res=cbind(pca.res,as.data.frame(riskScore_cli$riskScore2))
colnames(pca.res)[5]='group'
pca.var=data.pca$sdev^2 %>% as.data.frame()
pca.var$var=round(pca.var$./sum(pca.var)*100,2)
pca.var$pc=colnames(pca.res)[1:(ncol(pca.res)-1)]
pca.var
ggplot(pca.res, aes(PC1, PC2, color = group, shape = group))+ 
  # 选择X轴Y轴并映射颜色和形状
  geom_point(size = 4)+ # 画散点图并设置大小
  theme_bw() + # 加上边框
  # 自动提取主成分解释度进行绘图
  labs(x = paste('PC1(', pca.var$var[1],'%)', sep = ''),
       y = paste('PC2(', pca.var$var[2],'%)', sep = '')) +
  theme(
    panel.grid = element_blank(),
    legend.position = 'top')+ggtitle('TCGA-LIHC')+stat_ellipse(level = 0.95)
ggsave('result/figure/pca_LIHC.pdf',height =  4,width = 4)
library('Rtsne')
data_cox_tsne=Rtsne(data_cox[3:6])
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
fit <- survfit(Surv(OS.time1, as.numeric(OS)) ~ RiskScore, data=riskScore_cli)

BiocManager::install('survminer')
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
ggsave('result/figure/survival_TCGA_LIHC_12_28_1.pdf',
       height =9,width = 9)  
lasso_KM <- ggsurvplot(fit,data = riskScore_cli,
                       pval = T,
                       risk.table = F,
                       surv.median.line = "hv", #添加中位生存曲线
                       #palette=c("red", "blue"),  #更改线的颜色
                       #legend.labs=c("High risk","Low risk"), #标签
                       legend.title="RiskScore",
                       title="Overall survival", #标题
                       ylab="Survival probability",xlab = " Time (Days)", #更改横纵坐标
                       censor.shape = 124,censor.size = 2,conf.int = FALSE, #删失点的形状和大小
                       break.x.by = 720#横坐标间隔)
lasso_KM
BiocManager::install('timeROC')
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
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 32),
        legend.text =element_text(size = 26),
        title = element_text(size=26),
        panel.grid = element_blank())+
  coord_fixed()+
  ggtitle('timeROC')


ggsave('result/figure/time_ROC_TCGA_LIHC_12_28.pdf',height = 8,width = 8)







plot(ROC_riskscore, time = 365, col = "red", add = F,title = "")
plot(ROC_riskscore, time = 730, col = "blue", add = T)
plot(ROC_riskscore, time = 1080, col = "purple", add = T)
legend("bottomright",c("1-Year","3-Year","5-Year"),col=c("red","blue","purple"),lty=1,lwd=2)
text(0.5,0.2,paste("1-Year AUC = ",round(ROC_riskscore$AUC[1],3)))
text(0.5,0.15,paste("3-Year AUC = ",round(ROC_riskscore$AUC[2],3)))
text(0.5,0.1,paste("5-Year AUC = ",round(ROC_riskscore$AUC[3],3)))


###ICGC_LIRI_JP####
### extract the exp and clinical info
ICGC_clinic=read.csv('survival/donor.tsv/donor.tsv',sep = '\t')
sample_barcode=read.csv('survival/specimen.tsv/specimen.tsv',sep='\t')
exp=read.csv('survival/exp_seq.tsv/exp_seq.tsv',sep ='\t')
'Cancer' %in% exp$analysis_id
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

index=jp_exp |> 
  colnames() %>%
  grepl('normal',.)
index
jp_exp1=jp_exp[!index]
jp_exp1 |> dim()
ICGC_clinic |> colnames()
ICGC_clinic1=ICGC_clinic |> 
  dplyr::select(c("icgc_donor_id","donor_vital_status","donor_survival_time" ))
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

jp_exp1=jp_exp1[,sig_gene_mult_cox]
jp_exp1 |> dim()


ICGC_clinic2=ICGC_clinic1 |> 
  dplyr::select("OS_event","donor_survival_time") |>
  dplyr::rename(OS=OS_event) |>
  dplyr::rename(OS.time=donor_survival_time)

select_ICGC=ICGC_clinic2[tumorbarcode,] 
for_survival |> dim()
for_survival
jp_exp1=jp_exp1[,3]
for_survival1=for_survival[,sig_gene_mult_cox]
data_cox1=cbind(select_ICGC,for_survival1)

data_cox1=cbind(ICGC_clinic2,jp_exp1)
data_cox1 |> dim()
multicox1=coxph(Surv(OS.time,OS)~.,data=data_cox1)
multicox
data_cox |> head()
summary(multicox)
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
riskscore=predict(multicox,type = 'risk',
                  newdata = data_cox1)

riskscore
riskscores=as.data.frame(riskscore)
rownames(riskscores)=rownames(jp_exp1)
riskscores


riskScore_cli=cbind(ICGC_clinic2,riskscores)
riskScore_cli
library(ggplot2)
install.packages('ProComp')
data_cox1[3:6]

data.pca=prcomp(data_cox1[3:6],scale. = T)
summary(data.pca)
riskScore_cli
pca.res=data.pca$x
all(rownames(pca.res)==rownames(riskScore_cli))
pca.res=cbind(pca.res,as.data.frame(riskScore_cli$riskScore2))
colnames(pca.res)[5]='group'
pca.var=data.pca$sdev^2 %>% as.data.frame()
pca.var$var=round(pca.var$./sum(pca.var)*100,2)
pca.var$pc=colnames(pca.res)[1:(ncol(pca.res)-1)]
pca.var
ggplot(pca.res, aes(PC1, PC2, color = group, shape = group))+ 
  # 选择X轴Y轴并映射颜色和形状
  geom_point(size = 4)+ # 画散点图并设置大小
  theme_bw() + # 加上边框
  # 自动提取主成分解释度进行绘图
  labs(x = paste('PC1(', pca.var$var[1],'%)', sep = ''),
       y = paste('PC2(', pca.var$var[2],'%)', sep = '')) +
  theme(
        panel.grid = element_blank(),
        legend.position = c(0.85,0.85))+ggtitle('ICGC_LIRI-JP')
  # 设置图例位置，此处为相对位置
ggsave('result/figure/pca_validation.pdf',height = 4,width = 4)



#按照中位数分为高低风险两组
riskScore_cli$riskScore2 <- ifelse(riskScore_cli$riskscore > median(riskScore_cli$riskscore),
                                   "High","Low")


riskScore_cli$OS.time1=riskScore_cli$OS.time/30

data_cox_tsne=Rtsne(data_cox1[3:6])
data_cox_tsne_res=as.data.frame(data_cox_tsne$Y)
colnames(data_cox_tsne_res)=c('tSNE1','tSNE2')
data_cox_tsne_res$group=paste0(riskScore_cli$riskScore2,'_risk')

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
  ggtitle('ICGC_LIRI-JP')
ggsave('result/figure/tsne_JP.pdf',height = 6,width = 8)


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
ggsave('result/figure/validate_12_28_1.pdf',height = 8,width = 8)




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
  theme(axis.text = element_text(size = 30),
        axis.title = element_text(size = 32),
        legend.text =element_text(size = 26),
        title = element_text(size=26),
        panel.grid = element_blank())+
  coord_fixed()+
  ggtitle('timeROC(ICGC_LIRI-JP)')


ggsave('result/figure/time_ROC_12_28.pdf',height = 8,width = 8)
ROC_riskscore

pdf('result/figure/test_time_roc.pdf',height = 7,width = 7)
plot(ROC_riskscore, time = 365, col = "red", add = F,title = "")
plot(ROC_riskscore, time = 730, col = "blue", add = T)
plot(ROC_riskscore, time = 1080, col = "purple", add = T)

axis(1, at = seq(0, 1, by = 0.125), tcl = 0)
axis(1, at = seq(0, 1, by = 0.25), labels = FALSE, tcl = -0.5, lwd = 2)
axis(2, at = seq(0, 1, by = 0.125), tcl = 0)
axis(2, at = seq(0, 1, by = 0.25), labels = FALSE, tcl = -0.5, lwd = 2)

legend("bottomright",c("1-Year","2-Year","3-Year"),col=c("red","blue","purple"),lty=1,lwd=2)
text(0.5,0.2,paste("1-Year AUC = ",round(ROC_riskscore$AUC[1],3)))
text(0.5,0.15,paste("3-Year AUC = ",round(ROC_riskscore$AUC[2],3)))
text(0.5,0.1,paste("5-Year AUC = ",round(ROC_riskscore$AUC[3],3)))
#axis(2, at = seq(0, 1, by = 0.25))
#abline(h = seq(0, 1, by = 0.125), col = "gray", lty = 3)
#grid(nx = 8, ny = 8, col = "gray",lty='solid',lwd=0.8)
#grid(nx = 4, ny = 4, col = "gray",lwd=1.5,lty='solid')
xticks <- axTicks(1)
yticks <- axTicks(2)
abline(v = xticks, col = "gray", lty = "solid")
# 添加水平网格线
abline(h = yticks, col = "gray", lty = "solid")
dev.off()

