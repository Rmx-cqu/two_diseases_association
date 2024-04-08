### plot the rank conversion
plot_res=function(test=FALSE){
  library(patchwork)
  library(ggplot2)
if(test){
  df=data.frame(group=c('HCC','Ps'),
                'Rank_Convervation_index'=c(0.870,0.934))
  df1=data.frame(group=c('HCC_N_VS_Ps_N','HCC_vs_Ps'),
                 'Cosine_similarity'=c(0.857,0.892))
  permu=read.csv('permutation.csv')
  #p_value_location=0.0213
  save_file='result/figure/first_per.pdf'
}else{
  df=data.frame(group=c('HCC','Ps'),
                'Rank_Convervation_index'=c(0.856,0.948))
  df1=data.frame(group=c('HCC_NVSPs_N','HCCvsPs'),
                 'Cosine_similarity'=c(0.850,0.921))
  permu=read.csv('permutation_validate.csv')
  save_file='result/figure/validate_per.pdf'
}

  permu$different_values=permu$X0
  p0=ggplot(df,aes(x=group,y=Rank_Convervation_index))+
    geom_bar(stat='identity',width = 0.5,fill='#90BEE0')+
    coord_cartesian(ylim=c(0.5,1))+
    geom_text(aes(label=Rank_Convervation_index), vjust=1.6, color="white", size=3.5)+
    theme_bw()
  p1=ggplot(df1,aes(x=group,y=Cosine_similarity))+
    geom_bar(stat='identity',width = 0.5,fill='#90BEE0')+
    coord_cartesian(ylim=c(0.5,1))+
    geom_text(aes(label=Cosine_similarity), vjust=1.6, color="white", size=3.5)+
    theme_bw()+
    geom_signif(y_position = c(0.95),
                xmin = c(1),xmax=c(2),
                annotation=c('p<0.001'),tip_length = 0)
  p2=ggplot(permu,aes(x=different_values))+
    geom_histogram(aes(y=..density..),
                   binwidth = 0.001,fill='#90BEE0',
                   color='black',
                   alpha=0.4)+
    geom_density(color='#EA7369')+
    theme_bw()
  p=p0+p1+p2+plot_annotation(tag_levels = 'A')
  ggsave(save_file,height = 3,width = 8)
  return(p)
  }

plot_return=plot_res(test = T)
plot_return1=plot_res(test=F)
plot_return/plot_return1
ggsave('result/figure/per_comb_1_19.pdf',height = 5,width = 8)

diff=0.892-0.857
diff
permu=read.csv('permutation.csv')
permu$different_values=permu$X0
permu_res=c(permu$different_values)
count(permu_res<diff)
permu=read.csv('permutation_validate.csv')
diff=0.921-0.850
diff
permu$different_values=permu$X0
permu_res=c(permu$different_values)
count(permu_res>diff)


df=data.frame(group=c('HCC','Ps'),
              'Rank_Convervation_index'=c(0.870,0.934))
df1=data.frame(group=c('HCC_NVSPs_N','HCCvsPs'),
               'Cosine_similarity'=c(0.857,0.892))
permu=read.csv('permutation.csv')
#p_value_location=0.0213
save_dir='result/figure/'
permu$different_values=permu$X0
p0=ggplot(df,aes(x=group,y=Rank_Convervation_index))+
  geom_bar(stat='identity',width = 0.5,fill='#90BEE0')+
  coord_cartesian(ylim=c(0.5,1))+
  geom_text(aes(label=Rank_Convervation_index), vjust=1.6, color="white", size=3.5)+
  theme_bw()
p0
library(ggpubr)
p1=ggplot(df1,aes(x=group,y=Cosine_similarity))+
  geom_bar(stat='identity',width = 0.5,fill='#90BEE0')+
  coord_cartesian(ylim=c(0.5,1))+
  geom_text(aes(label=Cosine_similarity), vjust=1.6, color="white", size=3.5)+
  theme_bw()+
  geom_signif(y_position = c(0.95),
              xmin = c(1),xmax=c(2),
              annotation=c('p<0.001'),tip_length = 0)

p2=ggplot(permu,aes(x=different_values))+
  geom_histogram(aes(y=..density..),
                 binwidth = 0.001,fill='#90BEE0',
                 color='black',
                 alpha=0.4)+
  geom_density(color='#EA7369')+
  theme_bw()
p2
  geom_segment(aes(x = p_value_location ,xend = 1,p_value_location, y = 0, yend = 10), linetype = "dashed", color = "red")+
  annotate('text',x=p_value_location,y=14,label='p=0.001',size=4)
library(patchwork)
p0+p1+p2+plot_annotation(tag_levels = 'A')
save_file=file.path(save_dir,'permu1.pdf')
ggsave(save_file,height = 3,width = 8)

