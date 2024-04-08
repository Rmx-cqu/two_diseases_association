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
        plot.title = element_text(hjust = 0.5))+
  labs(x='Expression:log2(TPM+1)')

p1+p2+p3
ggsave('result/figure/ridges_hcc.pdf',height = 3,width = 3)
ggsave('result/figure/test_comb.pdf',height = 3,width = 10)

### the code for p4
