### Figure2A'
setwd('D:/two_diseases/')
source('function_defined.R')
save_dir='result/figure'

write.csv(DEG,'result/figure/LIHC_DEGS.csv')
write.csv(DEG3,'result/figure/Ps_DEGS.csv')

venn_plot=read.csv('result/figure/LIHC_DEGS.csv')
venn_plot1=read.csv('result/figure/Ps_DEGS.csv')

plot_list=list('TCGA_LIHC'=venn_plot$X,'PS'=venn_plot1$X)
save_dir='result/figure/figure2A.pdf'
plot_venn_diagram(plot_list,save_dir)


### figure 2B,C

final_intersect=intersect(venn_plot$X,
                          venn_plot1$X)
final_intersect
library(clusterProfiler)
library(org.Hs.eg.db)
entrez_genelists1=bitr(final_intersect,fromType = 'SYMBOL',
                       toType = 'ENTREZID',
                       OrgDb = 'org.Hs.eg.db')
entrez_genelists1
kk <- enrichKEGG(gene = entrez_genelists1$ENTREZID,
                  organism = 'hsa', #KEGG可以用organism = 'hsa'
                  pvalueCutoff = 1,
                  qvalueCutoff = 0.05)
ego_ALL <- enrichGO(gene = entrez_genelists1$ENTREZID, 
                     OrgDb = org.Hs.eg.db, 
                     #keytype = 'ENSEMBL',
                     ont = "ALL", #也可以是 CC  BP  MF中的一种
                     pAdjustMethod = "BH", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                     pvalueCutoff = 0.05, #P值会过滤掉很多，可以全部输出
                     qvalueCutoff = 0.05,
                     readable = TRUE)
dotplot(ego_ALL1,font.size=15)
kk@result |> head()
ego_ALL@result |> head()

Generatio=c(31/237,65/388,55/388,42/388,53/388,
            7/25,)

kk1@result |> head()
go=ego_ALL@result |> as.data.frame()
go_1=ego_ALL1@result |> as.data.frame()

go=go |> dplyr::select(c(Description,GeneRatio))
go_1=go_1 |> dplyr::select(c(Description,GeneRatio))
go=slice(go,1:5)

go_1
go
y=c(31/237,65/388,55/388,42/388,53/388,43/388,
    7/25,19/38,17/38,13/38,16/38,14/38)
y
left_join(go,go_1,by='Description')

getwd()
setwd('D:/two_diseases/HCC/GSE233422_RAW (2)/')
### set a base dataframe
test='GSM6134590_HCC_1_N_Genes_ReadCounts.txt.gz'
test=read.table(gzfile(test))
total_count=data.frame(gene_symbol=test$V1)
total_count |> head()
filenames=list.files()
#### combine the gz file
library(dplyr)
library(tidyverse)
for (filename in filenames){
  res=str_match(filename,'(GSM[0-9]{7}_)(.*)(_Genes_ReadCounts.txt.gz)')
  samplename=res[,3]
  count=read.table(gzfile(filename))
  colnames(count)=c('gene_symbol',samplename)
  total_count=left_join(total_count,count,by='gene_symbol')
}

rownames(total_count)=total_count[,1]
total_count=total_count[,-1]
index=filter_low_exp(total_count)

total_count=total_count[index,]
total_count |> dim()
traits3=total_count |> 
  colnames() %>%
  endsWith('N') |>
  if_else('normal','tumor') %>%
  data.frame(sample_name=colnames(total_count),status=.)

### set the rownames
rownames(traits3)=traits3$sample_name
all(rownames(traits3)==colnames(total_count))
traits3$status=factor(traits3$status)
traits3$status
traits3
### get the DEGS
validated_HCC=get_DEGS(total_count,traits3)
validated_HCC
setwd('D:/two_diseases/')
####
### DEGS for Psoriasis
validate=read.csv('D:/two_diseases/Psoriasis/GSE66511_Psoriasis_counts.txt/GSE66511_Psoriasis_counts.txt',
                  sep='\t')

library(GEOquery)
GSE_name='GSE66511'
options('download.file.method.GEOquery'='libcurl')
get=getGEO(GSE_name,getGPL = F)
Get=get[[1]]
info=pData(Get)
info=info |> dplyr::select(c(1,45))
### LP VS NLP
info$`disease:ch1`
info=info[info$`disease:ch1`!='C',]
### LP VS C
info=info[info$`disease:ch1`!='NLP',]
validate[(duplicated(validate$Symbol)),] |> head()

index=order(rowMeans(validate[,-1]),decreasing = TRUE)
expr_ordered=validate[index,]
keep=!duplicated(expr_ordered$Symbol)
expr_max=expr_ordered[keep,]

rownames(expr_max)=expr_max$Symbol
expr_max=expr_max[,-1]
index=filter_low_exp(expr_max)
expr_max=expr_max[index,]
expr_max=expr_max[,info$title]
expr_max[1:4,1:4]
info$`disease:ch1`=factor(info$`disease:ch1`)

info$status=info$`disease:ch1`
rownames(info)=info$title

info$status=relevel(info$status,ref = 'NLP')
info$status=relevel(info$status,ref='C')
info$status
validate_PS=get_DEGS(expr_max,info)
validate_PS1=get_DEGS(expr_max,info)
intersect2=intersect(rownames(validate_PS),rownames(validated_HCC))
intersect2
intersect3=intersect(rownames(validate_PS1),rownames(validated_HCC))
intersect4=intersect(intersect2,intersect3)
intersect2=intersect(intersect4,final_intersect)

entrez_genelists1=bitr(intersect2,fromType = 'SYMBOL',
                       toType = 'ENTREZID',
                       OrgDb = 'org.Hs.eg.db')
entrez_genelists1
kk1 <- enrichKEGG(gene = entrez_genelists1$ENTREZID,
                  organism = 'hsa', #KEGG可以用organism = 'hsa'
                  pvalueCutoff = 1,
                  qvalueCutoff = 0.05)
ego_ALL1 <- enrichGO(gene = entrez_genelists1$ENTREZID, 
                     OrgDb = org.Hs.eg.db, 
                     #keytype = 'ENSEMBL',
                     ont = "ALL", #也可以是 CC  BP  MF中的一种
                     pAdjustMethod = "BH", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                     pvalueCutoff = 0.05, #P值会过滤掉很多，可以全部输出
                     qvalueCutoff = 0.05,
                     readable = TRUE)
kk1@result |>head()
ego_ALL1@result |> head()

Generatio=c(31/237,65/388,7/25,19/38)
x=(Generatio[3]+Generatio[4])/2
y=(Generatio[1]+Generatio[2])/2
x-y


combn(Generatio,2)
permutation_res=c(0.24083625000000003,0.12836204999999995,
  0.09163795,0.09163795,0.12836204999999995,
  0.24083625000000003)
permutation_res
plot_df1=data.frame(values=permutation_res)
ggplot(plot_df1, aes(x = values)) +
  geom_histogram(stat='bin',binwidth = 0.01)
plot_df=data.frame(GeneRatio=c(31/237,65/388,7/25,19/38),
                   Pathways=rep(c('KEGG:cell_cycle','Go:nuclear_divsion'),2),
                   Gene_sets=rep(c('407(genes)','39(genes)'),each=2))
plot_df=data.frame(GeneRatio=y,
                   Gene_sets=rep(c('407(genes)','39(genes)'),each=6),
                   Pathways=rep(c('cell_cylce',go$Description),2),
                   source=rep(c('KEGG',rep('Go',5)),2))

plot_df$Pathways=factor(plot_df$Pathways,levels=rev(c('cell_cylce',go$Description)))


library(ggplot2)
library(stringr)
p=ggplot(plot_df,aes(x=Pathways,weight=GeneRatio))+
  geom_bar(aes(fill=Gene_sets),width = 0.7)+
  scale_fill_manual(values = c('#EA7369','#90BEE0'))+
  theme_bw()+
  ylab('GeneRatio')+
  annotate('text',x=6,y=0.6,label='p<0.001',size=4)+
  coord_flip()+
  scale_x_discrete(labels=function(x) str_wrap(x, width=23))+
  theme(axis.text = element_text(size = 28),
         axis.title = element_text(size = 28),
        legend.text =element_text(size = 27),
        legend.title = element_text(size = 26))

p
ggsave('result/1_18.pdf',height = 7,width = 14)

# red 219 49 36
# blue 144 190 224
red_code=rgb(144,190,224)

pathway_per=read.csv('pathways_permutate.csv')
pathway_per$different_values=pathway_per$test
p2=ggplot(pathway_per,aes(x=different_values))+
  geom_histogram(aes(y=..density..),
                 binwidth = 0.01,fill='#90BEE0',
                 color='black',
                 alpha=0.4)+
  geom_density(color='#EA7369')+
  theme_bw()+
  scale_x_continuous(breaks=seq(0, 0.3, 0.05))+
  annotate(
    "segment", x=0.26, xend=0.26, y=2.5,yend=0.5, 
    color="#EA7369",arrow = arrow(length = unit(0.5, "cm"), angle = 20))+
  annotate('text',x=0.26,y=2.7,label='OV')+
  theme(axis.text = element_text(size = 28),
        axis.title = element_text(size = 28),
        legend.text =element_text(size = 27),
        legend.title = element_text(size = 26))
p2
ggsave('per_pathways_histogram.pdf',height = 7,width = 7)
