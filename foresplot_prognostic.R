library(TCGAbiolinks)
clinical_data <- GDCquery_clinic(project = "TCGA-LIHC", type = "clinical")

clinical_data |> head()
clinical_data |> colnames()
test_clinical=clinical_data |> dplyr::select(c(bcr_patient_barcode,gender,age_at_index,
                                 age_at_diagnosis,ajcc_pathologic_stage))
 

LIHC_clinical1 |> head()
test_clinical |> head()

total_clinical=left_join(riskScore_cli,test_clinical,by=c('X_PATIENT'='bcr_patient_barcode'))
total_clinical |> dim()
filter_cli=total_clinical |> dplyr::select(c(sample,X_PATIENT,
                                  OS,OS.time,riskscore,
                                  riskScore2,OS.time1,
                                  gender,age_at_index,
                                  ajcc_pathologic_stage))

filter_cli=filter_cli |> na.omit() 
filter_cli$ajcc_pathologic_stage |> factor()
filter_cli$age=filter_cli$age_at_index
filter_cli$stage=if_else(filter_cli$ajcc_pathologic_stage %in% c('Stage I','Stage II'),'Stage_I_II','stage_III_Iv')
filter_cli$riskScore2=factor(filter_cli$riskScore2) |> relevel(,ref='Low') 

LIHC_fit=coxph(Surv(OS.time,OS)~riskScore2+age+gender+stage,data=filter_cli)

covariates=c('riskscore','age','gender','stage')
filter_final=filter_cli |> dplyr::select(c(OS.time,OS,riskscore,age,gender,stage))
library(broom)
for (i in c(3:6)){
  index_col=c(1,2,i)
  filter_test=filter_final[,index_col]
  res_cox=coxph(Surv(OS.time,OS)~.,data=filter_test)
  tydy1=tidy(res_cox)
  print(tydy1)
}



### Then for ICGC
ICGC_clinic=read.csv('survival/donor.tsv/donor.tsv',sep = '\t')
ICGC_clinic$donor_tumour_stage_at_diagnosis
ICGC_clinic$donor_sex
ICGC_clinic$donor_age_at_diagnosis
ICGC_clinic$icgc_donor_id
log_jp_exp1=log2(jp_exp1+1)
ICGC_clinic1
riskScore_cli |> head()
riskScore_cli$icgc_donor_id=rownames(riskScore_cli)
ICGC_clinic_filter=ICGC_clinic |> dplyr::select(c(icgc_donor_id,
                                                  donor_age_at_diagnosis,
                                                  donor_sex,
                                                  donor_tumour_stage_at_diagnosis))
ICGC_clinic_filter$age=ICGC_clinic_filter$donor_age_at_diagnosis
ICGC_clinic_filter$sex=ICGC_clinic_filter$donor_sex
ICGC_clinic_filter$stage=if_else(ICGC_clinic_filter$donor_tumour_stage_at_diagnosis %in% c(1,2),'stage_I_II','stage_III_IV')
ICGC_clinic_filter1=ICGC_clinic_filter |> dplyr::select(c(age,sex,
                                                          stage,
                                                       icgc_donor_id))
total_clinical=left_join(riskScore_cli,ICGC_clinic_filter1)
rownames(total_clinical)=total_clinical$icgc_donor_id
total_clinical=total_clinical |>dplyr::select(-c(icgc_donor_id))
total_clinical |> head()
total_clinical$risk_group=if_else(total_clinical$riskscore>median(total_clinical$riskscore),'High','Low')
total_clinical1=total_clinical |> dplyr::select(-c(riskscore,riskScore2,OS.time1))
total_clinical1$risk_group=factor(total_clinical1$risk_group)
total_clinical1$risk_group=relevel(total_clinical1$risk_group,ref='Low')
test_cox=coxph(Surv(OS.time,OS)~.,data=total_clinical1)
test_cox
for (i in c(3:6)){
  index_col=c(1,2,i)
  filter_final=total_clinical
  filter_test=filter_final[,index_col]
  res_cox=coxph(Surv(OS.time,OS)~.,data=filter_test)
  tydy1=tidy(res_cox)
  print(res_cox)
}









### plot forestplt
LIHC_fit
mul_cox1=summary(LIHC_fit)
mul_cox1
colnames(mul_cox1$conf.int)
multi1=as.data.frame(round(mul_cox1$conf.int[,c(1,3,4)],2))
multi1
library(forestplot)
library(stringr)
library(survival)
library(tableone)
multi2=ShowRegTable(LIHC_fit,
                    exp=T,
                    digits=2,
                    pDigits = 3,
                    printToggle = T,
                    quote = F,
                    ciFun = confint)

multi1
result=cbind(multi1,multi2)

result<-tibble::rownames_to_column(result, var = "Characteristics")
result=rbind(c('characteristics',NA,NA,NA,"HR (95%CI)", "P.value"),result)
result[, 2:4] <- lapply(result[, 2:4], as.numeric)

result=result[c(1,3,4,5,2),]
result
result$Characteristics[1]="Multivariate analysis"
result$Characteristics=c("Multivariate analysis",
                         "Age",
                         'Gender(Male vs Female)',
                         'Tumor stage(III/IV vs I/II)',
                         'Risk score(High vs Low)')
pdf('result/figure/multi_figure_lihc_forest.pdf',,height = 9,width = 7)
P1 <-   forestplot(result[,c(1,5,6)], #告诉函数，合成的表格result的第1，5，6列还是显示数字
                   mean=result[,2],   #告诉函数，表格第2列为HR，它要变成森林图的小方块
                   lower=result[,3],  #告诉函数表格第3列为5%CI，
                   upper=result[,4],
                   # 
                   is.summary = c(T,rep(F,nrow(result)-1)),
                   ##定义x轴
                   # xlab="log2(OR)",
                   xlab="HR",
                   
                   #定义标题
                   # title=paste0(id,":",method),
                   hrzl_lines=list( "2" = gpar(lwd=2, col="black")), # 画水平线
                   
                   zero = 1, # log2转换后，纵参考线改为0
                   boxsize = 0.3,
                   lineheight = unit(8,"mm"),
                   colgap = unit(4,"mm"),
                   lwd.zero = 2,
                   lwd.ci = 2,
                   # clip=c(-2,2), # 选择特定的区域
                   # coord_trans(x='log2'),
                   col = fpColors(box=  "#2166AC",
                                  # lines=  "#2166AC",
                                  zero = "gray50",
                                  # box = "#458B00",
                                  # summary = "8B008B",
                                  lines = "black",
                   ),
                   graphwidth=unit(55,"mm"),
                   # 文字大小
                   txt_gp=fpTxtGp(label=gpar(cex=3.0),
                                  ticks=gpar(cex=3.2),
                                  xlab=gpar(cex = 3.0),
                                  title=gpar(cex = 3.2)),
                   lwd.xaxis = 2,
                   lty.ci = "solid",
                   graph.pos = 2)
P1
dev.off()
