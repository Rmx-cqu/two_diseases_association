## define the function(DEGS,WCGNA,TCGA)----
### define the function for DEGS
get_DEGS=function(data,info){
  library(DESeq2)
  dds <-DESeqDataSetFromMatrix(countData=data, 
                               colData=info, 
                               design=~status,
                               tidy=FALSE)
  dds <- DESeq(dds)
  resLFc1=results(dds)
  resLFc1=subset(resLFc1,!is.na(padj))
  resLFc1=resLFc1[abs(resLFc1$log2FoldChange)>1,]
  resLFc1=resLFc1[order(resLFc1$log2FoldChange,decreasing = TRUE),]
  filter=resLFc1[resLFc1$padj<0.01,]
  return(filter)
  
}

get_DEG_not_filter=function(data,info){
  library(DESeq2)
  dds <-DESeqDataSetFromMatrix(countData=data, 
                               colData=info, 
                               design=~status,
                               tidy=FALSE)
  dds <- DESeq(dds)
  resLFc1=results(dds)
  resLFc1=subset(resLFc1,!is.na(padj))
  #resLFc1=resLFc1[abs(resLFc1$log2FoldChange)>1,]
  resLFc1=resLFc1[order(resLFc1$log2FoldChange,decreasing = TRUE),]
  #filter=resLFc1[resLFc1$padj<0.01,]
  return(resLFc1)
  
}
### define function to acquire TCGA-LIHC counts
get_TCGA_exp=function(TCGA_ID,data_format='count'){
  filename=paste0(TCGA_ID,'_mRNA.Rdata')
  file_path=file.path('D:/gsva',filename)
  load(file_path)
  ### get TCGA-LIHC exp
  se=data
  library(SummarizedExperiment)
  ### get rowData: gene information
  row=rowData(se)
  ### keep the protein_coding genes
  se_mrna=se[row$gene_type=='protein_coding',]
  ### using assay to extract the expression matrix count(TPM)
  if (data_format=='count'){
    tpm=assay(se_mrna,'unstranded')}
  else{
    tpm=assay(se_mrna,'tpm_unstrand')
  }
  
  ### get the symbol_name
  symbol_name=rowData(se_mrna)$gene_name
  ### cbind
  tpm_symbol=cbind(data.frame(symbol_name),as.data.frame(tpm))
  ### find the duplicated genes
  tpm_symbol[duplicated(tpm_symbol$symbol_name ),] |> nrow()
  ### use mean statistic for the duplicated genes
  #exp=aggregate(.~symbol_name,mean,data=tpm_symbol)
  
  ### use the higher expression for the duplicated names
  index=order(rowMeans(tpm_symbol[,-1]),decreasing = TRUE)
  expr_ordered=tpm_symbol[index,]
  keep=!duplicated(expr_ordered$symbol_name)
  exp=expr_ordered[keep,]
  
  ### pre-process the exp-matrix
  rownames(exp)=exp$symbol_name
  exp=exp[,-1]
  ### convert the count to integer
  #exp=round(exp)
  return(exp)
}

### load packages
get_module_color=function(dataExpr,design,project_name,power='auto_selection'){
  library(WGCNA)
  powers=c(c(1:10),seq(from = 12, to=30, by=2))
  sft=pickSoftThreshold(dataExpr,powerVector = powers,
                          verbose = 5,
                          networkType = "signed",
                          corFnc = 'bicor')
  
  sizeGrWindow(9,5)
  par(mfrow=c(1,2))
  cex1=0.9
  
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.90,col="red")
  
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  print(sft$powerEstimate)
  saveTOM_filename=paste0(project_name,'_TOM')
  ### build the networks
  if (power=='auto_selection'){
    sft_power=sft$powerEstimate
  }
  else{
    sft_power=7
  }
  net = blockwiseModules(dataExpr,
                         power = sft_power,
                         corType = "bicor",
                         TOMType = "signed", 
                         networkType = "signed",
                         minModuleSize = 200,
                         reassignThreshold = 0, 
                         mergeCutHeight = 0.15,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         pamRespectsDendro = FALSE,
                         saveTOMFileBase = saveTOM_filename,
                         nThreads=31,
                         verbose = 3)
  
  mergedColors=labels2colors(net$colors)
  plotDendroAndColors(net$dendrograms[[1]],
                      mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  
  moduleColors <- labels2colors(net$colors)
  
  MEs0 = moduleEigengenes(dataExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩阵(样本vs模块)
  ### get correlation
  MEs= MEs[,-ncol(MEs)]
  moduleTraitCor = cor(MEs, design , use = "p");
  ### 
  nSamples= nrow(dataExpr)
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  
  #### save the moduletraitPvalue
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  
  dim(textMatrix) = dim(moduleTraitCor)
  save_name=paste0(project_name,'-Module-trait-relationships_2_withoutlog.pdf')
  pdf(save_name,width = 7,height = 7,pointsize = 12)
  
  par(mar = c(6, 8.5, 3, 3))
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
  return(list(moduleColors,MEs0))
}

### return the module_genes
get_module_genes=function(dataExpr,moduleColors,color){
  module_genes=dataExpr |> 
    colnames() %>% 
    .[moduleColors==color]
  return(module_genes)
  
}

### load gene_list for GSVA
get_genelist=function(ID='hsa04110'){
  library('KEGGREST')
  GS=keggGet(ID)
  genes=unlist(lapply(GS[[1]]$GENE,function(x) strsplit(x,';')))
  genelist <- genes[1:length(genes)%%3 ==2]    
  genelist
  class(genelist)
  genelist=list(genelist)
  names(genelist)='Cell_cycle'
  return(genelist)
}

### return the TPM value for a row counts
count_to_tpm=function(raw_counts){
  if (file.exists('efflen_10_17.csv')){
    geneid_efflen=read.csv('efflen_10_17.csv')
  }
  else{
    library(rtracklayer)
    
    gtf=rtracklayer::import('gencode.v44.annotation.gtf.gz')
    gtf=as.data.frame(gtf)
    
    exon=gtf[gtf$type=='exon',c('start','end','gene_name')]
    
    exon_bygeneid=split(exon,exon$gene_name)
    library(parallel)
    cl <- makeCluster(0.75*detectCores()) 
    efflen=parLapply(cl,exon_bygeneid,function(x){
      tmp=apply(x,1,function(y){y[1]:y[2]})
      length(unique(unlist(tmp)))
    })
    geneid_efflen=data.frame(geneid=names(efflen),
                             efflen=as.numeric(efflen))
    write.csv(geneid_efflen,'efflen_10_17.csv')
  }
  
  gene_df=data.frame(geneid=rownames(raw_counts))
  gene_df=left_join(gene_df,geneid_efflen,by='geneid')
  gene_df=gene_df |> na.omit()
  raw_counts=raw_counts[gene_df$geneid,]
  
  kb=gene_df$efflen/1000
  df=apply(raw_counts,2,function(x){
    (x/kb)/sum(x/kb)*1000000
  })
  return(df)
}

### plot the wcgna heatmaps
plot_heatmaps=function(expression,WCGNA_res,target_DEGS,project_name){
  
  ### access the color labels of each gene
  color_labels=WCGNA_res[[1]]
  ### each gene  row: samples; col:genes
  gene_HC=expression |> colnames()
  ### re-assign the gene label to color
  names(color_labels)=gene_HC
  
  ### count the number of DEGS of the modules
  x=color_labels[rownames(target_DEGS)] |> table()
  ### count the number of genes belong to modules
  y=color_labels |> table()
  x1=as.data.frame(x)
  y1=as.data.frame(y)
  ## join 
  x_y=left_join(y1,x1,by=c('color_labels'='Var1'))
  x_y[is.na(x_y)]=0
  ### calculate the ratio
  x_y$ratio=x_y$Freq.y/x_y$Freq.x
  
  ### transpose(row:genes,col:samples)
  gsva_tpm=expression |> t()
  
  library(GSVA)
  library(parallel)
  genelist=get_genelist()
  scores=gsva(expr=gsva_tpm, 
              gset.idx.list=genelist, 
              kcdf="Gaussian" ,#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
              verbose=T, 
              method='gsva',
              parallel.sz = parallel::detectCores())
  
  scores=scores |> t()
  ### extract the egen-gene
  ege_gene=WCGNA_res[[2]]
  moduleTraitCor=cor(ege_gene,scores,method = c('spearman'))
  nSamples= ncol(gsva_tpm)
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "")
  ### cbind the ratio 
  textMatrix=cbind(textMatrix,round(x_y$ratio,2))
  moduleTraitCor1=cbind(moduleTraitCor,round(x_y$ratio,2))
  #print(ege_gene)
  ### remove the grey module
  index=which(rownames(moduleTraitCor1)=='MEgrey')
  moduleTraitCor1=moduleTraitCor1[-index,]
  textMatrix=textMatrix[-index,]
  ege_gene=ege_gene[,-which(names(ege_gene) == "MEgrey")]
  ###
  dim(textMatrix) = dim(moduleTraitCor1)
  save_name=paste('Cell_cycle',project_name,'without_normal.pdf',sep='_')
  
  pdf(save_name,width = 7,height = 7,pointsize = 12)
  xlabs=c('Cell_cycle','DEGS_ratio')
  length(xlabs)
  par(mar = c(6, 8.5, 3, 3))
  labeledHeatmap(Matrix = moduleTraitCor1,
                 xLabels = xlabs,
                 yLabels = names(ege_gene),
                 ySymbols = names(ege_gene),
                 colorLabels = FALSE,
                 colors = greenWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 1.8,
                 zlim = c(-1,1),
                 cex.lab.x = 1.7,
                 cex.lab.y = 1.7,
                 cex.lab = 1.7,
                 cex.legendLabel = 1.7,
                 xLabelsAngle = 0,
                 main = paste("Module-trait relationships"))
  dev.off()
  
}

### filter the low exp according to cpm(half of the samples>1)
filter_low_exp=function(exp){
  ## row:gene; col:sample
  tbl=as.data.frame(exp)
  cpm=apply(tbl, 2, function(x){
    x/sum(x)*1000000
  })
  index=rowSums( cpm >= 1 ) >= ncol(tbl)/2
  return(index)
}

### plot the venn diagram
plot_venn_diagram=function(res_list,save_name){
  library(VennDiagram)
  if(length(res_list)==2){
    fill_color=c( "blue",'yellow')
  } else if(length(res_list)==3){
    fill_color=c( "blue",'yellow',"red")
  } else{
    fill_color=c( "blue",'yellow',"red", "green")
  }
  venn.plot=venn.diagram(
    x = res_list,
    filename = NULL,##韦恩图的名字
    lty = 1,
    lwd = 1,
    col = "black",  ##圈的颜色
    fill = fill_color,##对应每个圈的颜色，有几个数据集，就需要有相应数量的颜色
    alpha = 0.60,
    cat.col = "black",##此处设置每个数据集的名称颜色，也可以使用c（）函数输入三种颜色
    cat.cex = 1.4,
    cat.fontface = "bold",
    margin = 0.07,
    cex = 2
  )
  pdf(save_name)
  grid.draw(venn.plot)
  dev.off()
}
