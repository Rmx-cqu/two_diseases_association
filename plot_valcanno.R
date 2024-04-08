#### plot valcanlo plot
source('function_defined.R') 
setwd('D:/two_diseases/')
## Get TCGA_LIHC exp non-compared ----
library(dplyr)
### get the exp
exp=get_TCGA_exp('TCGA-LIHC')
dataExpr1=as.data.frame(exp)
dataExpr1 |> dim()
### filter the low exp
index=filter_low_exp(dataExpr1)
dataExpr1=dataExpr1[index,]
dataExpr1 |> dim()


### get the trait info
library(dplyr)
traits = dataExpr1 |> 
  colnames() |>
  substr(1,16) %>% 
  grepl('TCGA-..-....-0..',x=.) |>
  ifelse('liver_tumor','normal') %>%
  data.frame(row.names =colnames(dataExpr1),status=.)
traits |> head()


###
datTraits1=traits
traits$sample_name=rownames(traits)
all(rownames(traits)==colnames(dataExpr1))
traits$status=factor(traits$status)
traits$status=relevel(traits$status,ref = 'normal')

DEG=get_DEGS(dataExpr1,traits)
DEG11=get_DEG_not_filter(dataExpr1,traits)



#####
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE54456", "file=GSE54456_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)
#write.csv(tbl,'GSE54456_counts.csv')
tbl |> dim()
index1=filter_low_exp(tbl)
tbl1=tbl[index1,]
tbl1 |> dim()
### write a function to filter the low cpm values

###Psoriasis for WGCNA
library(GEOquery)
GSE_name='GSE54456'
options('download.file.method.GEOquery'='libcurl')
get=getGEO(GSE_name,getGPL = F)
Get=get[[1]]
info=pData(Get)
info=info |> dplyr::select(c(1,37))
info |> nrow()
info
### test the order
info1=info[match(colnames(tbl1),rownames(info)),]
match(colnames(tbl1),rownames(info1))
### the new info
traits2=info1
traits2$sample_name=rownames(traits2)
traits2=traits2 |>
  dplyr::rename(status=`tissue type:ch1`) |>
  dplyr::select(c(status,sample_name))
traits2=traits2 %>%
  mutate(
    status=case_when(
      status=='normal skin' ~'normal_skin',
      status=='lesional psoriatic skin' ~ 'lesional_psoriatic_skin'
    )
  )

traits2$status=factor(traits2$status)
traits2$status=relevel(traits2$status,ref = 'normal_skin')
traits2$status
match(traits2$sample_name,colnames(tbl1))
all(traits2$sample_name==colnames(tbl1))
tbl1 |> dim()
library(org.Hs.eg.db)
library(clusterProfiler)
trans_df=bitr(rownames(tbl1),
              fromType = 'ENTREZID',
              toType = 'SYMBOL',
              OrgDb = org.Hs.eg.db)
trans_df |> dim()
tbl1 |> dim()
tbl2=tbl1 |>  as.data.frame()
tbl2$gene_name=rownames(tbl2)
tbl2=tbl2 |>
  dplyr::filter(gene_name %in% trans_df$ENTREZID) |>
  dplyr:: select(-c(gene_name))
tbl2 |> dim()
trans_df |> dim()
rownames(tbl2)=trans_df$SYMBOL
DEG3=get_DEGS(tbl2,traits2)
DEG33=get_DEG_not_filter(tbl2,traits2)

final_intersect=intersect(rownames(DEG3),rownames(DEG))
final_intersect |> length()


entrez_genelists1=bitr(final_intersect,fromType = 'SYMBOL',
                       toType = 'ENTREZID',
                       OrgDb = 'org.Hs.eg.db')
entrez_genelists1
kk1 <- enrichKEGG(gene = entrez_genelists1$ENTREZID,
                  organism = 'hsa', #KEGG可以用organism = 'hsa'
                  pvalueCutoff = 1,
                  qvalueCutoff = 0.05)
index=kk1@result$ID=='hsa04110'
index1=kk1@result[index,]$geneID %>% 
  strsplit(.,'/') |> 
  unlist() 
index1
cell_cycle=entrez_genelists1 |> 
  filter(ENTREZID %in% index1)
cell_cycle



plot_valcano=function(input_df,DEG,cell_cycle,file_name){
  ### assign the symbol column
  input_df=input_df |> as.data.frame()
  input_df$symbol=rownames(input_df)
  input_df=input_df |>
    mutate(change=as.factor(if_else(padj<0.01 & abs(log2FoldChange)>1,
                                    if_else(log2FoldChange>1,'Up','Down'),'Not Sig')))
  ### assign the symbol column
  DEG_DF=as.data.frame(DEG) 
  DEG_DF$symbol=rownames(DEG_DF)
  ### get the cell cycle related DEGS
  DEG_DF_cell_cycle=DEG_DF |>
    filter(symbol %in% cell_cycle$SYMBOL)
  ### color  the cell cycle related DEGS
  input_df=input_df |> 
    mutate(annotation=symbol %in% DEG_DF_cell_cycle$symbol)
  index=input_df$annotation==T
  
  ### select the top5 genes
  cell_cycle_anno=input_df %>%
    filter(annotation==T) |>
    top_n(5,log2FoldChange)
  ###annotate the top5 genes
  input_df=input_df |> 
    mutate(annotation_gene=symbol %in% cell_cycle_anno$symbol)
  
  ### plot
  ggplot(input_df,aes(log2FoldChange,-log10(padj)))+
    geom_hline(yintercept = -log10(0.01),linetype='dashed')+
    geom_vline(xintercept = c(-1,1),linetype='dashed')+
    geom_point(aes(color=change),
               size=5,
               alpha=0.8)+
    theme_bw(base_size = 15)+
    scale_color_manual(values=c("blue","grey","black"))+
    theme(
          legend.position = 'right',
          axis.text = element_text(size = 21),
          axis.title = element_text(size = 23),
          legend.text =element_text(size = 19) )+
    xlab('Log2FC')+
    ylab('-Log10(p.adj)')+
    geom_point(data=subset(input_df,annotation==T),
               aes(x=log2FoldChange,
                   y=-log10(padj)),
               color='red',
               size=5)+
    ggrepel::geom_label_repel(
      data=subset(input_df,annotation_gene==T),
      aes(log2FoldChange,-log10(padj),label=symbol),
      size=8,
      box.padding = unit(1, "lines"),
      point.padding = unit(0.5, "lines"),
      segment.color = "black")
  file_name=file.path('result/figure',file_name)
  ggsave(file_name,height = 5,width = 7)
}

plot_valcano(DEG11,DEG,cell_cycle,'test.pdf')
plot_valcano(DEG33,DEG3,cell_cycle,'test1.pdf')






