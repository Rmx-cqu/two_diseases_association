### plot ridge_plot
getwd()
LIHC_all_info=cbind(log_df1,label1)
write.csv(LIHC_all_info,'LIHC_for_python.csv')
Ps_all_info=cbind(df2,label2)
write.csv(Ps_all_info,'Ps_for_python.csv')

index=c('KPNA2','FAM83D','AURKA','PTTG1','BUB1',
        'CENPF','TOP2A','RRM2','KIF23')
label2,df2
df2 |> colnames()
target_df=df2
df_select=target_df[,index]


label_select=label2
plot_density=cbind(df_select,label_select)
plot_density=plot_density |> as.data.frame()
plot_density |>head()
plot_density=plot_density |>
  dplyr::rename(sample_type=label_select)

plot_density$sample_type1=if_else(plot_density$sample_type=='1','Ps','Normal')

plot_density=plot_density |>
  dplyr::select(-c(sample_type))
colnames(plot_density)
select_cols=colnames(plot_density)[1:9]


plot_density1=plot_density %>%
  gather(key='Gene',value = 'Value',all_of(select_cols))
library(ggplot2)
library(ggridges)

plot_density1$sample_type1=plot_density1$sample_type1 |> as.factor()
plot_density1$sample_type1=relevel(plot_density1$sample_type1,ref = 'Ps')

index1=rev(index)
plot_density1$Gene=factor(plot_density1$Gene,levels=index1)
ggplot(plot_density1,aes(x=Value,y=Gene,fill=sample_type1))+
  geom_density_ridges(alpha=0.7)+
  theme_ridges()+
  theme(legend.position = "top",
        legend.title =element_blank(),
        panel.grid = element_blank())+
  labs(x='Expression:log2(TPM+1)')

plot_density |> colnames()
test_t=plot_density[,c(1,10)]
index=test_t$sample_type1=='Ps'
x1=test_t[index,][,1]
x2=test_t[!index,][,1]
t_result=t.test(x1,x2)
t_result$p.value
t.test(x1,x2)
empty_vector=vector()
for(i in 1:(ncol(plot_density)-1)){
  index_vector=c(i,ncol(plot_density))
  test_t=plot_density[,index_vector]
  index=test_t$sample_type1=='Ps'
  ### extract
  x1=test_t[index,][,1]
  x2=test_t[!index,][,1]
  t_result=wilcox.test(x1,x2)
  empty_vector=c(empty_vector,t_result$p.value)
}
empty_vector
p.adjust(empty_vector, method = "BH")
ggsave('result/figure/ridges_ps.pdf',height = 3,width = 3)





#### then plot the ridges_for the other
index=c('KPNA2','FAM83D','AURKA',
        'PTTG1', 'CCNB1','BIRC5',
        'CDC6','CDKN3')
index1=c('KPNA2','FAM83D','AURKA','PTTG1','BUB1',
        'CENPF','TOP2A','RRM2','KIF23')

df1
log_df1=log2(df1+1)

df1 |> colnames()
log_df1[,2] |> quantile(probs=c(0.25,0.75,0.95))

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
library(ggplot2)
library(ggridges)

plot_density1$sample_type1=plot_density1$sample_type1 |> as.factor()
plot_density1$sample_type1=relevel(plot_density1$sample_type1,ref = 'Tumor')

index1=rev(index)
plot_density1$Gene=factor(plot_density1$Gene,levels=index1)
ggplot(plot_density1,aes(x=Value,y=Gene,fill=sample_type1))+
  geom_density_ridges(alpha=0.7)+
  theme_ridges()+
  theme(legend.position = "top",
        legend.title =element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  labs(x='Expression:log2(TPM+1)')


plot_density |> colnames()

empty_vector=vector()
for(i in 1:(ncol(plot_density)-1)){
  index_vector=c(i,ncol(plot_density))
  test_t=plot_density[,index_vector]
  index=test_t$sample_type1=='Tumor'
  ### extract
  x1=test_t[index,][,1]
  x2=test_t[!index,][,1]
  t_result=wilcox.test(x1,x2)
  empty_vector=c(empty_vector,t_result$p.value)
}
empty_vector
p.adjust(empty_vector, method = "BH")
ggsave('result/figure/ridges_hcc.pdf',height = 3,width = 3)
