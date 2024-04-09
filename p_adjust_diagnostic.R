library(table1)
#install.packages('table1')


multi_tidy=tidy(multicox)
multi_tidy$p.adj=p.adjust(multi_tidy$p.value)
index=order(multi_tidy$p.adj)
multi_tidy1=multi_tidy[index,]

multi_tidy2=tidy(new_multicox)
multi_tidy2$p.adj=p.adjust(multi_tidy2$p.value)
index=order(multi_tidy2$p.adj)
multi_tidy2=multi_tidy2[index,]
multi_tidy2

multi_tidy3=tidy(new_multicox)
multi_tidy3$p.adj=p.adjust(multi_tidy3$p.value)
index=order(multi_tidy3$p.adj)
multi_tidy3=multi_tidy3[index,]
multi_tidy3

multi_tidy11=multi_tidy1 |> dplyr::select(c(term,p.adj))
multi_tidy22=multi_tidy2 |> dplyr::select(c(term,p.adj))
multi_tidy33=multi_tidy3 |> dplyr::select(c(term,p.adj))

multi_tidy11
multi_tidy22
test1=full_join(multi_tidy11,multi_tidy22,by='term')
final_table=full_join(test1,multi_tidy33)
final_table |> colnames()
colnames(final_table)=c('genes','model1','model2','final_model')
table1(~. ,data=final_table)
install.packages("sjPlot")
library(sjPlot)
summary_table <- tab_df(final_table,digits = 3)
tab_df
summary_table
ggsave('result/figure/test_table.pdf',plot=summary_table)
install.packages("pander")
library(pander)
final_table |> colnames()
summary_table <- pander(final_table)
final_table1=final_table[,c(2,3,4)] |> round(digits = 3)
rownames(final_table1)=final_table$genes
write.csv(final_table1,'result/figure/summary.csv')
