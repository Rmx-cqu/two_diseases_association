## load counts table from GEO
# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE83645", "file=GSE83645_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

tbl |> dim()
tbl[1:4,1:4]
source('function_defined.R')
index=filter_low_exp(tbl)
tbl=tbl[index,]
tbl |> dim()
tbl |> head()

library(clusterProfiler)
library('org.Hs.eg.db')
trans_df=bitr(rownames(tbl),
              fromType = 'ENTREZID',
              toType = 'SYMBOL',
              OrgDb = org.Hs.eg.db)

trans_df |> dim()
trans_df |> head()
tbl |> dim()
tbl2=tbl |>  as.data.frame()
tbl2$gene_name=rownames(tbl2)
#tbl3=tbl2[tbl2$gene_name %in% trans_df$ENTREZID,]
#all(tbl3$gene_name==trans_df$ENTREZID)
#tbl3$gene_name
#trans_df$ENTREZID
tbl2=tbl2 |>
  dplyr::filter(gene_name %in% trans_df$ENTREZID) |>
  dplyr:: select(-c(gene_name))
tbl2 |> dim()
trans_df |> dim()
### resign the rownames
rownames(tbl2)=trans_df$SYMBOL

tbl2 |> dim()
trans_df |> dim()
library(dplyr)
tbl3=count_to_tpm(tbl2)
tbl3 |> colSums()

tbl3 |> colnames()






library(GEOquery)
GSE_name='GSE83645'
options('download.file.method.GEOquery'='libcurl')
get=getGEO(GSE_name,getGPL = F)
Get=get[[1]]
info=pData(Get)

info
info=info |> 
  dplyr::select(c('tissue type:ch1')) %>%
  dplyr::mutate(id=rownames(.))
info |> dim()
info
tbl3 |> dim()
tbl31=t(tbl3)
tbl31=tbl31 |> as.data.frame()
tbl31 |> dim()
tbl31$id=rownames(tbl31)
new_cbind=left_join(tbl31,info,by='id')
new_cbind |> dim()

label_psoaris=new_cbind[,c(ncol(new_cbind)-1,ncol(new_cbind))]
label_psoaris
label_psoaris$label=if_else(label_psoaris$`tissue type:ch1`=='psoriasis',1,0)
label_psoaris$label
tbl31 |> dim()

rownames(tbl31)==label_psoaris$id
tbl32=tbl31[,genelists$x]

tbl32=tbl32 |> as.matrix()
tbl32 |> dim()
cbind(tbl32,)
tbl32=log2(tbl32,label_psoaris$label)

predict_test1=predict(cvfit1,newx = tbl32,s='lambda.min',
                      type='response')


predict_test1 |> dim()
library(pROC)
roc_obj=roc(label_psoaris$label,predict_test1)
table(label_psoaris$label)
table(predict_test1)
psorisis_label=psorisis_label |> as.matrix()
roc_obj
library(ggplot2)
ggr
roc_data=data.frame(specificity=roc_obj$specificities,
           sensitivity=roc_obj$sensitivities)
ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "red")

pc<- ggroc(roc_obj, legacy.axes = TRUE,
           color='red',
           size=1)+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype=4)+
  theme_bw() +ggtitle('ROC')+ 
  annotate("text",x=0.75,y=0.375,label=paste("AUC = ", round(roc_obj$auc,3)),size=7)+
  theme(axis.text = element_text(size = 24),
        axis.title = element_text(size = 26),
        legend.text =element_text(size = 22),
        title = element_text(size=20))
pc
ggsave('result/figure/psoriasis_gse83645.pdf',height = 5,width = 5)
label_psoaris
predict_test1
