forest_plot <- function(data_input){
  
  # 数据预处理
  for (cols in c("or","or_lci95","or_uci95","egger_intercept","se_pleio","pval_pleio","Q_MR_Egger","Q_pval_MR_Egger",'Q_pval_ivw')) {
    data_input[,cols] <- as.numeric(sprintf("%.3f",  data_input[,cols]))
    
  }
  
  data_input$or_95ci <- paste0(data_input$or,"(",data_input$or_lci95,"-",data_input$or_uci95,")")
  data_input$pval_mr <-  scales::scientific(  data_input$pval_mr,digits = 3)
  data_input$Q_pval <-  scales::scientific(  data_input$Q_pval_MR_Egger,digits = 3)
  data_input$nsnp <- ifelse(!duplicated(data_input$exposure),data_input$nsnp,"")
  data_input$egger_intercept <- ifelse(!duplicated(data_input$exposure),data_input$egger_intercept,"")
  data_input$pval_pleio <- ifelse(!duplicated(data_input$exposure),data_input$pval_pleio,"")
  #data_input$Q_MR_Egger <- ifelse(!duplicated(data_input$exposure),data_input$Q_MR_Egger,"")
  #data_input$Q_pval_MR_Egger <- ifelse(!duplicated(data_input$exposure),data_input$Q_pval_MR_Egger,"")
  data_input$heterogeneity_ivw <- ifelse(!duplicated(data_input$exposure),data_input$Q_pval_ivw,"")
  data_input$exposure <- ifelse(!duplicated(data_input$exposure),data_input$exposure,"")
  
  
  # 准备森林图的第一行
  if (T) {
    data_input <- rbind(colnames(data_input),data_input)
    data_input[1,c("exposure","outcome","nsnp","method","pval_mr","or_95ci","egger_intercept","pval_pleio",'heterogeneity_ivw')] <-
      c("Exposure","Outcome","Nsnp","Method","Pval_mr","OR(95%CI)","Egger_intercept","Pval_pleio",'Heterogeneity_IVW')
    
    
    
    data_input[is.na(data_input$pval_mr),"or"] <- 1
    data_input[is.na(data_input$pval_mr),"or_lci95"] <- 1
    data_input[is.na(data_input$pval_mr),"or_uci95"] <- 1
    data_input[is.na(data_input$pval_mr),"pval_mr"] <- 1
    data_input[data_input$method=="Inverse variance weighted","method"] <- "IVW*"
  }
  
  
  
  
  
  # 作图
  library(forestplot)
  P1 <-   forestplot(labeltext = as.matrix(data_input[,c("exposure","outcome","nsnp","method","pval_mr","or_95ci",
                                                         'heterogeneity_ivw', "pval_pleio")]),
                     
                     lower = as.numeric(data_input[,"or_lci95"]),
                     upper = as.numeric(data_input[,"or_uci95"]),
                     mean = as.numeric(data_input[,"or"]),
                     # 
                     is.summary = c(T,rep(F,nrow(data_input)-1)),
                     ##定义x轴
                     # xlab="log2(OR)",
                     xlab="OR",
                     
                     #定义标题
                     # title=paste0(id,":",method),
                     hrzl_lines=list( "2" = gpar(lwd=2, col="black")), # 画水平线
                     
                     zero = 1, # log2转换后，纵参考线改为0
                     boxsize = 0.4,
                     lineheight = unit(8,"mm"),
                     colgap = unit(2,"mm"),
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
                     graphwidth=unit(40,"mm"),
                     # 文字大小
                     txt_gp=fpTxtGp(label=gpar(cex=1.2),
                                    ticks=gpar(cex=1.4),
                                    xlab=gpar(cex = 1.2),
                                    title=gpar(cex = 1.4)),
                     lwd.xaxis = 2,
                     lty.ci = "solid",
                     graph.pos = 7)
  return(P1)
  
}