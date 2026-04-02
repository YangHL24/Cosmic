###################################################乳腺癌亚型特异性ceRNA模块推断###########################################

## Input data
load("lncRNA_se.RData")
load("mRNA_se.RData")
load("miRNA_se.RData")
load("miRTar_se.RData")
load("BRCA_subtype_index.RData")
load('survival_result.RData')
load("lncREXP_final.RData")
load("mREXP_final.RData")
load('subtype_col.RData')

# Install and load the required packages
library(org.Hs.eg.db)
library(reactome.db)
library(SummarizedExperiment)
library(GSEABase)
# Details of the miRSM R package can be found at http://bioconductor.org/packages/miRSM/ and https://github.com/zhangjunpeng411/miRSM
library(miRSM)
library(igraph)
# Details of the miRspongeR R package can be found at http://bioconductor.org/packages/miRspongeR/ and https://github.com/zhangjunpeng411/miRspongeR.
library(miRspongeR)
library(ggplot2)

############################## 所有亚型的ceRNA模块推断 ###############################
set.seed(1234)
# lncRNA_se: SummarizedExperiment object of lncRNA: rows are sample, columns are lncRNA
# mRNA_se: SummarizedExperiment object of mRNA: rows are sample, columns are mRNA
# miRNA_se: SummarizedExperiment object of miRNA: rows are sample, columns are miRNA
# miRTar_se: SummarizedExperiment object of putative miRNA-target interactions
# Output: modulegenes_WGCNA_BRCA is a list of gene co-expression modules
modulegenes_WGCNA_BRCA <- module_WGCNA(lncRNA_se,mRNA_se)
# Output:miRSM_WGCNA_SDC is a list of ceRNA modules
miRSM_WGCNA_SDC <- miRSM(miRNA_se, lncRNA_se, mRNA_se, miRTar_se,
                         modulegenes_WGCNA_BRCA,
                         num_shared_miRNAs = 3, pvalue.cutoff = 0.05,
                         method = "SDC", MC.cutoff = 0.8,
                         SMC.cutoff = 0.1, RV_method = "RV")

################################ 乳腺癌亚型特异性ceRNA模块的推断 ################################
library(foreach)
library(doParallel)
## 设置并行计算的核数
num.cores <- 3
cl <- makeCluster(num.cores)
registerDoParallel(cl)
# nsamples: 乳腺癌亚型的数量
nsamples <- 5
modulegenes_all_WGCNA<-modulegenes_WGCNA_BRCA

## Infer the gene co-expression module after knocking out the sample set
# lncRNA_se: SummarizedExperiment object of lncRNA: rows are sample, columns are lncRNA
# mRNA_se: SummarizedExperiment object of mRNA: rows are sample, columns are mRNA
# miRNA_se: SummarizedExperiment object of miRNA: rows are sample, columns are miRNA
# miRTar_se: SummarizedExperiment object of putative miRNA-target interactions
# index: A list of sample indices for each breast cancer subtype in the gene expression matrix
# Output: results_modulegenes_exceptk_WGCNA results_modulegenes_exceptk_WGCNA is a list of gene co-expression modules identified after removing subtype k
results_modulegenes_exceptk_WGCNA <- foreach(i = seq(nsamples), .packages = "miRSM") %dopar% {
  module_WGCNA(lncRNA_se[-index[[i]]$index, ], mRNA_se[-index[[i]]$index, ])
}
## Infer the ceRNA module after knocking out the sample set
# Output: miRSM_SDC_exceptk_WGCNA is a list of ceRNA modules identified after removing subtype k
miRSM_SDC_exceptk_WGCNA<- foreach(i = seq(nsamples), .packages = "miRSM") %dopar% { miRSM(miRNA_se[-index[[i]]$index, ],
                                                                                          lncRNA_se[-index[[i]]$index,], mRNA_se[-index[[i]]$index, ],
                                                                                          miRTar_se, results_modulegenes_exceptk_WGCNA[[i]],
                                                                                          method = "SDC",
                                                                                          SMC.cutoff = 0.1)}

stopCluster(cl)
stopImplicitCluster()

miRSM_SDC_all_WGCNA <- miRSM_WGCNA_SDC
Modulegenes_all_SDC_WGCNA <- miRSM_SDC_all_WGCNA[[2]]
Modulegenes_exceptk_SDC_WGCNA <- lapply(seq(nsamples), function(i) miRSM_SDC_exceptk_WGCNA[[i]][[2]])
# Output: Modules_SS_SDC_WGCNA is the list of ceRNA modules of subtype k
Modules_SS_SDC_WGCNA <- miRSM_SS(Modulegenes_all_SDC_WGCNA, Modulegenes_exceptk_SDC_WGCNA)
Modules_SS_SDC_WGCNA


####Visualization of the number of breast cancer-specific ceRNA modules
#Modules_SS_SDC_WGCNA: List of ceRNA modules identified for each breast cancer subtype
##Check whether the input length matches
subtype_names = c("Normal", "LumA", "Her2+", "LumB", "Basal")
if (length(subtype_names) != length(Modules_SS_SDC_WGCNA)) {
  stop("subtype_names length must match modules_list length")
}

ceRNA_modules <- data.frame(
  Subtype = subtype_names,
  "Number of ceRNA modules" = sapply(Modules_SS_SDC_WGCNA, length),
  row.names = NULL,
  check.names = FALSE
)
except_ceRNA_Modulegenes_plot <- ggplot(data = ceRNA_modules, mapping = aes(x = Subtype, y =`Number of ceRNA modules` , fill = Subtype, group = factor(1))) +
  ggtitle("BRCA Subtype-specific CeRNA Modules") +
  scale_fill_manual(values = c("#80FFFF", "#BFFFFF", "#FFD5FF", "#FFBFFF", "#FF80FF"))+
  theme_bw() +
  geom_bar(stat="identity",width=0.4)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +  # 向下扩展y轴范围
  theme(text = element_text(family = "serif"),
        axis.text = element_text(size = 17),
        axis.title = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 21),  # 居中并增大标题字体大小
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13))+  # 调整图例文本大小为12
  geom_text(aes(label = `Number of ceRNA modules`), vjust = -0.5, color = "black", size = 7,family = "serif")
except_ceRNA_Modulegenes_plot

##################################亚型特异性ceRNA模块功能解析 ###############################################
#### 异质性分析
subtype_mapping <- list("Normal" = "Sample 1",
                        "LumA" = "Sample 2",
                        "Her2+" = "Sample 3",
                        "LumB" = "Sample 4",
                        "Basal" = "Sample 5")

specific_ceRNA_module_list <- list()
for (subtype in names(subtype_mapping))
{
  ceRNA_data <- Modules_SS_SDC_WGCNA[[subtype_mapping[[subtype]]]]
  # Process each module and combine ceRNA and mRNA
  merged_data <- sapply(ceRNA_data, function(x) paste(c(x$ceRNA, x$mRNA), collapse = ""))
  specific_ceRNA_module_list[[subtype]] <- merged_data
}
library(UpSetR)
except_ceRM_upset<-upset(fromList(specific_ceRNA_module_list),
                         sets=names(specific_ceRNA_module_list),
                         point.size=3.5,
                         line.size=1.5,
                         mainbar.y.label="Number of shared modules",
                         main.bar.color = 'black',
                         matrix.color="black",
                         sets.x.label="Set size",
                         sets.bar.color=c('red', 'orange','green','blue','purple'),
                         order.by = "freq",
                         text.scale=c(3.2,3,2,2.5,2.5,3)
)
except_ceRM_upset

### 基于亚型特异性ceRNA模块计算亚型的相似度
# Normal_ceRNA: Normal亚型的ceRNA模块列表
# LumA_ceRNA: LumA亚型的ceRNA模块列表
# Her2_ceRNA: Her2亚型的ceRNA模块列表
# LumB_ceRNA: LumB亚型的ceRNA模块列
# Basal_ceRNA: Basal亚型的ceRNA模块列
Normal_ceRNA<-Modules_SS_SDC_WGCNA$`Sample 1`
LumA_ceRNA<-Modules_SS_SDC_WGCNA$`Sample 2`
Her2_ceRNA<-Modules_SS_SDC_WGCNA$`Sample 3`
LumB_ceRNA<-Modules_SS_SDC_WGCNA$`Sample 4`
Basal_ceRNA<-Modules_SS_SDC_WGCNA$`Sample 5`

## 分别计算Normal亚型与LumA、Her2、LumB、Basal亚型的相似度
sim_Normal_LumA<-module_group_sim(Normal_ceRNA, LumA_ceRNA, sim.method = "Simpson")
sim_Normal_Her2<-module_group_sim(Normal_ceRNA, Her2_ceRNA, sim.method = "Simpson")
sim_Normal_LumB<-module_group_sim(Normal_ceRNA, LumB_ceRNA, sim.method = "Simpson")
sim_Normal_Basal<-module_group_sim(Normal_ceRNA, Basal_ceRNA, sim.method = "Simpson")

## 分别计算LumA亚型与Her2、LumB、Basal亚型的相似度
sim_LumA_Her2<-module_group_sim(LumA_ceRNA, Her2_ceRNA,sim.method = "Simpson")
sim_LumA_LumB<-module_group_sim(LumA_ceRNA, LumB_ceRNA,sim.method = "Simpson")
sim_LumA_Basal<-module_group_sim(LumA_ceRNA, Basal_ceRNA,sim.method = "Simpson")

## 分别计算Her2亚型与LumB和Basal亚型的相似度
sim_Her2_LumB<-module_group_sim(Her2_ceRNA, LumB_ceRNA,sim.method = "Simpson")
sim_Her2_Basal<-module_group_sim(Her2_ceRNA, Basal_ceRNA,sim.method = "Simpson")

## 计算LumB亚型和Basal亚型的相似度
sim_LumB_Basal<-module_group_sim(LumB_ceRNA, Basal_ceRNA,sim.method = "Simpson")

#### 富集分析
##乳腺癌亚型特异性ceRNA模块功能富集分析
library(openxlsx)
miRSM_Normal_FEA <- module_FA(Normal_ceRNA, Analysis.type = "FEA")
miRSM_LumA_FEA <- module_FA(LumA_ceRNA, Analysis.type = "FEA")
miRSM_Her_FEA <- module_FA(Her2_ceRNA, Analysis.type = "FEA")
miRSM_LumB_FEA <- module_FA(LumB_ceRNA, Analysis.type = "FEA")
miRSM_Basal_FEA <- module_FA(Basal_ceRNA, Analysis.type = "FEA")
# Normal
Nor_FEA_BP<-miRSM_Normal_FEA[[1]]
for (i in 1:length(Nor_FEA_BP)) {
  Nor_FEA_BP_data <- data.frame(Nor_FEA_BP[[i]])
  
  file_name <- paste0("Normal_", i, "FEA_BP_result.xlsx")
  write.xlsx(Nor_FEA_BP_data, file = file_name, rowNames = FALSE)
}
Nor_FEA_KEGG<-miRSM_Normal_FEA[[2]]
for (i in 1:length(Nor_FEA_KEGG)) {
  Nor_FEA_KEGG_data <- data.frame(Nor_FEA_KEGG[[i]])
  
  file_name <- paste0("Normal_", i, "FEA_KEGG_result.xlsx")
  write.xlsx(Nor_FEA_KEGG_data, file = file_name, rowNames = FALSE)
}
Nor_FEA_Reactome<-miRSM_Normal_FEA[[3]]
for (i in 1:length(Nor_FEA_Reactome)) {
  Nor_FEA_Reactome_data <- data.frame(Nor_FEA_Reactome[[i]])
  
  file_name <- paste0("Normal_", i, "FEA_Reactome_result.xlsx")
  write.xlsx(Nor_FEA_Reactome_data, file = file_name, rowNames = FALSE)
}
# Her2
Her_FEA_BP<-miRSM_Her_FEA[[1]]
for (i in 1:length(Her_FEA_BP)) {
  Her_FEA_BP_data <- data.frame(Her_FEA_BP[[i]])
  
  file_name <- paste0("Her_", i, "FEA_BP_result.xlsx")
  write.xlsx(Her_FEA_BP_data, file = file_name, rowNames = FALSE)
}
Her_FEA_KEGG<-miRSM_Her_FEA[[2]]
for (i in 1:length(Her_FEA_KEGG)) {
  Her_FEA_KEGG_data <- data.frame(Her_FEA_KEGG[[i]])
  
  file_name <- paste0("Her_", i, "FEA_KEGG_result.xlsx")
  write.xlsx(Her_FEA_KEGG_data, file = file_name, rowNames = FALSE)
}
Her_FEA_Reactome<-miRSM_Her_FEA[[3]]
for (i in 1:length(Her_FEA_Reactome)) {
  Her_FEA_Reactome_data <- data.frame(Her_FEA_Reactome[[i]])
  
  file_name <- paste0("Her_", i, "FEA_Reactome_result.xlsx")
  write.xlsx(Her_FEA_Reactome_data, file = file_name, rowNames = FALSE)
}
# LumB
LB_FEA_BP<-miRSM_LumB_FEA[[1]]
for (i in 1:length(LB_FEA_BP)) {
  LumB_FEA_BP_data <- data.frame(LB_FEA_BP[[i]])
  
  file_name <- paste0("LumB_", i, "FEA_BP_result.xlsx")
  write.xlsx(LumB_FEA_BP_data, file = file_name, rowNames = FALSE)
}
LB_FEA_KEGG<-miRSM_LumB_FEA[[2]]
for (i in 1:length(LB_FEA_KEGG)) {
  LumB_FEA_KEGG_data <- data.frame(LB_FEA_KEGG[[i]])
  
  file_name <- paste0("LumB_", i, "FEA_KEGG_result.xlsx")
  write.xlsx(LumB_FEA_KEGG_data, file = file_name, rowNames = FALSE)
}
LB_FEA_Reactome<-miRSM_LumB_FEA[[3]]
for (i in 1:length(LB_FEA_Reactome)) {
  LumB_FEA_Reactome_data <- data.frame(LB_FEA_Reactome[[i]])
  
  file_name <- paste0("LumB_", i, "FEA_Reactome_result.xlsx")
  write.xlsx(LumB_FEA_Reactome_data, file = file_name, rowNames = FALSE)
}
#Basal
Ba_FEA_BP<-miRSM_Basal_FEA[[1]]
for (i in 1:length(Ba_FEA_BP)) {
  Basal_FEA_BP_data <- data.frame(Ba_FEA_BP[[i]])
  
  file_name <- paste0("Basal_", i, "FEA_BP_result.xlsx")
  write.xlsx(Basal_FEA_BP_data, file = file_name, rowNames = FALSE)
}
Ba_FEA_KEGG<-miRSM_Basal_FEA[[2]]
for (i in 1:length(Ba_FEA_KEGG)) {
  Basal_FEA_KEGG_data <- data.frame(Ba_FEA_KEGG[[i]])
  
  file_name <- paste0("Basal_", i, "FEA_KEGG_result.xlsx")
  write.xlsx(Basal_FEA_KEGG_data, file = file_name, rowNames = FALSE)
}
Ba_FEA_Reactome<-miRSM_Basal_FEA[[3]]
for (i in 1:length(Ba_FEA_Reactome)) {
  Basal_FEA_Reactome_data <- data.frame(Ba_FEA_Reactome[[i]])
  
  file_name <- paste0("Basal_", i, "FEA_Reactome_result.xlsx")
  write.xlsx(Basal_FEA_Reactome_data, file = file_name, rowNames = FALSE)
}

## 乳腺癌亚型特异性ceRNA模块疾病富集分析
miRSM_Normal_DEA <- module_FA(Normal_ceRNA, Analysis.type = "DEA")
miRSM_LumA_DEA <- module_FA(LumA_ceRNA, Analysis.type = "DEA")
miRSM_Her_DEA <- module_FA(Her2_ceRNA, Analysis.type = "DEA")
miRSM_LumB_DEA <- module_FA(LumB_ceRNA, Analysis.type = "DEA")
miRSM_Basal_DEA <- module_FA(Basal_ceRNA, Analysis.type = "DEA")
###Normal
Nor_DEA_HDO<-miRSM_Normal_DEA[[1]]
for (i in 1:length(Nor_DEA_HDO)) {
  Nor_DEA_HDO_data <- data.frame(Nor_DEA_HDO[[i]])
  
  file_name <- paste0("Normal_", i, "DEA_HDO_result.xlsx")
  write.xlsx(Nor_DEA_HDO_data, file = file_name, rowNames = FALSE)
}
Nor_DEA_DisGeNET<-miRSM_Normal_DEA[[2]]
for (i in 1:length(Nor_DEA_DisGeNET)) {
  Nor_DEA_DisGeNET_data <- data.frame(Nor_DEA_DisGeNET[[i]])
  
  file_name <- paste0("Normal_", i, "DEA_DisGeNET_result.xlsx")
  write.xlsx(Nor_DEA_DisGeNET_data, file = file_name, rowNames = FALSE)
}
Nor_DEA_NCG<-miRSM_Normal_DEA[[3]]
for (i in 1:length(Nor_DEA_NCG)) {
  Nor_DEA_NCG_data <- data.frame(Nor_DEA_NCG[[i]])
  
  file_name <- paste0("Normal_", i, "DEA_NCG_result.xlsx")
  write.xlsx(Nor_DEA_NCG_data, file = file_name, rowNames = FALSE)
}
# LumA
LA_DEA_HDO<-miRSM_LumA_DEA[[1]]
for (i in 1:length(LA_DEA_HDO)) {
  LA_DEA_HDO_data <- data.frame(LA_DEA_HDO[[i]])
  
  file_name <- paste0("LumA_", i, "DEA_HDO_result.xlsx")
  write.xlsx(LA_DEA_HDO_data, file = file_name, rowNames = FALSE)
}
LA_DEA_DisGeNET<-miRSM_LumA_DEA[[2]]
for (i in 1:length(LA_DEA_DisGeNET)) {
  LA_DEA_DisGeNET_data <- data.frame(LA_DEA_DisGeNET[[i]])
  
  file_name <- paste0("LumA_", i, "DEA_DisGeNET_result.xlsx")
  write.xlsx(LA_DEA_DisGeNET_data, file = file_name, rowNames = FALSE)
}
LA_DEA_NCG<-miRSM_LumA_DEA[[3]]
for (i in 1:length(LA_DEA_NCG)) {
  LA_DEA_NCG_data <- data.frame(LA_DEA_NCG[[i]])
  
  file_name <- paste0("LumA_", i, "DEA_NCG_result.xlsx")
  write.xlsx(LA_DEA_NCG_data, file = file_name, rowNames = FALSE)
}
# Her2
Her_DEA_HDO<-miRSM_Her_DEA[[1]]
for (i in 1:length(Her_DEA_HDO)) {
  Her_DEA_HDO_data <- data.frame(Her_DEA_HDO[[i]])
  
  file_name <- paste0("Her_", i, "DEA_HDO_result.xlsx")
  write.xlsx(Her_DEA_HDO_data, file = file_name, rowNames = FALSE)
}
Her_DEA_DisGeNET<-miRSM_Her_DEA[[2]]
for (i in 1:length(Her_DEA_DisGeNET)) {
  Her_DEA_DisGeNET_data <- data.frame(Her_DEA_DisGeNET[[i]])
  
  file_name <- paste0("Her_", i, "DEA_DisGeNET_result.xlsx")
  write.xlsx(Her_DEA_DisGeNET_data, file = file_name, rowNames = FALSE)
}
Her_DEA_NCG<-miRSM_Her_DEA[[3]]
for (i in 1:length(Her_DEA_NCG)) {
  Her_DEA_NCG_data <- data.frame(Her_DEA_NCG[[i]])
  
  file_name <- paste0("Her_", i, "DEA_NCG_result.xlsx")
  write.xlsx(Her_DEA_NCG_data, file = file_name, rowNames = FALSE)
}
# LumB
LB_DEA_HDO<-miRSM_LumB_DEA[[1]]
for (i in 1:length(LB_DEA_HDO)) {
  LumB_DEA_HDO_data <- data.frame(LB_DEA_HDO[[i]])
  
  file_name <- paste0("LumB_", i, "DEA_HDO_result.xlsx")
  write.xlsx(LumB_DEA_HDO_data, file = file_name, rowNames = FALSE)
}
LB_DEA_DisGeNET<-miRSM_LumB_DEA[[2]]
for (i in 1:length(LB_DEA_DisGeNET)) {
  LumB_DEA_DisGeNET_data <- data.frame(LB_DEA_DisGeNET[[i]])
  
  file_name <- paste0("LumB_", i, "DEA_DisGeNET_result.xlsx")
  write.xlsx(LumB_DEA_DisGeNET_data, file = file_name, rowNames = FALSE)
}
LB_DEA_NCG<-miRSM_LumB_DEA[[3]]
for (i in 1:length(LB_DEA_NCG)) {
  LumB_DEA_NCG_data <- data.frame(LB_DEA_NCG[[i]])
  
  file_name <- paste0("LumB_", i, "DEA_NCG_result.xlsx")
  write.xlsx(LumB_DEA_NCG_data, file = file_name, rowNames = FALSE)
}
# Basal
Ba_DEA_HDO<-miRSM_Basal_DEA[[1]]
for (i in 1:length(Ba_DEA_HDO)) {
  Basal_DEA_HDO_data <- data.frame(Ba_DEA_HDO[[i]])
  
  file_name <- paste0("Basal_", i, "DEA_HDO_result.xlsx")
  write.xlsx(Basal_DEA_HDO_data, file = file_name, rowNames = FALSE)
}
Ba_DEA_DisGeNET<-miRSM_Basal_DEA[[2]]
for (i in 1:length(Ba_DEA_DisGeNET)) {
  Basal_DEA_DisGeNET_data <- data.frame(Ba_DEA_DisGeNET[[i]])
  
  file_name <- paste0("Basal_", i, "DEA_DisGeNET_result.xlsx")
  write.xlsx(Basal_DEA_DisGeNET_data, file = file_name, rowNames = FALSE)
}
Ba_DEA_NCG<-miRSM_Basal_DEA[[3]]
for (i in 1:length(Ba_DEA_NCG)) {
  Basal_DEA_NCG_data <- data.frame(Ba_DEA_NCG[[i]])
  
  file_name <- paste0("Basal_", i, "DEA_NCG_result.xlsx")
  write.xlsx(Basal_DEA_NCG_data, file = file_name, rowNames = FALSE)
}

#### 分布分析
## 提取每个乳腺癌特异性ceRNA模块的共享miRNA
# Nor_share: Normal亚型各ceRNA模块共享的miRNA
# LumA_share: LumA亚型各ceRNA模块共享的miRNA
# Her_share: Her2亚型各ceRNA模块共享的miRNA
# LumB_share: LumB亚型各ceRNA模块共享的miRNA
# Basal_share: Basal亚型各ceRNA模块共享的miRNA
Nor_share<-share_miRs(miRNA_se,miRTar_se,Normal_ceRNA)
LumA_share<-share_miRs(miRNA_se,miRTar_se,LumA_ceRNA)
Her_share<-share_miRs(miRNA_se,miRTar_se,Her2_ceRNA)
LumB_share<-share_miRs(miRNA_se,miRTar_se,LumB_ceRNA)
Basal_share<-share_miRs(miRNA_se,miRTar_se,Basal_ceRNA)

Normal_share_upset<-upset(fromList(Nor_share),
                          sets=names(Nor_share),
                          point.size=3.5,
                          line.size=1.5,
                          mainbar.y.label="Number of shared miRNAs",
                          main.bar.color = 'black',
                          matrix.color="black",
                          sets.x.label="Set size",
                          sets.bar.color=c('red','blue'),
                          order.by = "freq",
                          text.scale=c(3,3,2,2.5,2.5,3)
)
Normal_share_upset
LumA_share_upset<-upset(fromList(LumA_share),
                        sets=names(LumA_share),
                        point.size=3.5,
                        line.size=1.5,
                        mainbar.y.label="Number of shared miRNAs",
                        main.bar.color = 'black',
                        matrix.color="black",
                        sets.x.label="Set size",
                        sets.bar.color=c('pink','red', 'orange','yellow','green','cyan','blue','purple','magenta'),
                        order.by = "freq",
                        text.scale=c(3.2,3,2,2.5,2.5,3)
)
LumA_share_upset
Her2_share_upset<-upset(fromList(Her_share),
                        sets=names(Her_share),
                        point.size=3.5,
                        line.size=1.5,
                        mainbar.y.label="Number of shared miRNAs",
                        main.bar.color = 'black',
                        matrix.color="black",
                        sets.x.label="Set size",
                        sets.bar.color=c('red','green','blue'),
                        order.by = "freq",
                        text.scale=c(3.2,3,2,2.5,2.5,3)
)
Her2_share_upset
LumB_share_upset<-upset(fromList(LumB_share),
                        sets=names(LumB_share),
                        point.size=3.5,
                        line.size=1.5,
                        mainbar.y.label="Number of shared miRNAs",
                        main.bar.color = 'black',
                        matrix.color="black",
                        sets.x.label="Set size",
                        sets.bar.color=c('red', 'orange','green','blue'),
                        order.by = "freq",
                        text.scale=c(3.2,3,2,2.5,2.5,3)
)
LumB_share_upset
Basal_share_upset<-upset(fromList(Basal_share),
                         sets=names(Basal_share),
                         point.size=3.5,
                         line.size=1.5,
                         mainbar.y.label="Number of shared miRNAs",
                         main.bar.color = 'black',
                         matrix.color="black",
                         sets.x.label="Set size",
                         sets.bar.color=c('red', 'orange','green','blue'),
                         order.by = "freq",
                         text.scale=c(3.2,3,2,2.5,2.5,3)
)
Basal_share_upset

####生存分析
lnc<-assay(lncRNA_se)
mR<-assay(mRNA_se)
miRNA<-assay(miRNA_se)
exp<-cbind(miRNA,mR,lnc)
row<-BRCA_sur$sample
row_index<-match(row,rownames(exp))
expdata<-exp[row_index,]

Normal_list <- list()
for (i in 1:length(Normal_ceRNA)) {
  Normal_list[[i]] <- c(Normal_ceRNA[[i]]$ceRNA, Normal_ceRNA[[i]]$mRNA)
}

LumA_list <- list()
for (i in 1:length(LumA_ceRNA)) {
  LumA_list[[i]] <- c(LumA_ceRNA[[i]]$ceRNA, LumA_ceRNA[[i]]$mRNA)
}

Her2_list <- list()
for (i in 1:length(Her2_ceRNA)) {
  Her2_list[[i]] <- c(Her2_ceRNA[[i]]$ceRNA, Her2_ceRNA[[i]]$mRNA)
}

LumB_list <- list()
for (i in 1:length(LumB_ceRNA)) {
  LumB_list[[i]] <- c(LumB_ceRNA[[i]]$ceRNA, LumB_ceRNA[[i]]$mRNA)
}

Basal_list <- list()
for (i in 1:length(Basal_ceRNA)) {
  Basal_list[[i]] <- c(Basal_ceRNA[[i]]$ceRNA, Basal_ceRNA[[i]]$mRNA)
}
Nor_sur<-moduleSurvival(Normal_list , expdata ,BRCA_sur,devidePercentage=.5, plot = TRUE)
LumA_sur<-moduleSurvival(LumA_list , expdata ,BRCA_sur,devidePercentage=.5, plot = TRUE)
Her_sur<-moduleSurvival(Her2_list , expdata ,BRCA_sur,devidePercentage=.5, plot = TRUE)
LumB_sur<-moduleSurvival(LumB_list , expdata ,BRCA_sur,devidePercentage=.5, plot = TRUE)
Basal_sur<-moduleSurvival(Basal_list , expdata ,BRCA_sur,devidePercentage=.5, plot = TRUE)


# sur_senlin：生存分析结果数据框
tabletext <- cbind(
  c("Module", sur_senlin$Module),
  c("Chi-square",round(sur_senlin$`Chi-square`,2)),
  c("HR (95% CI)", paste0(round(sur_senlin$`HR (95% CI)`, 2), " (", round(sur_senlin$HRlow95, 2), " - ", round(sur_senlin$HRup95, 2), ")")),
  c("p-value", sur_senlin$`p-value`)
)
library(forestplot)
forest_plot <- forestplot(
  labeltext = tabletext,
  #title = paste0(gene_name, " – Overall Survival"),
  #title_gp = grid::gpar(fontsize = 14, fontface = "bold", just = "left"),
  mean = c(NA, sur_senlin$`HR (95% CI)`),
  lower = c(NA, sur_senlin$HRlow95),
  upper = c(NA, sur_senlin$HRup95),
  hrzl_lines = list(`1` = grid::gpar(lwd = 2, col = "black"),
                    `2` = grid::gpar(lwd = 2, col = "black"),
                    `24` = grid::gpar(lwd = 2, col = "black")),
  is.summary = rep(FALSE, nrow(sur_senlin) + 1),
  zero = 1,
  boxsize = 0.2,
  lineheight = unit(0.7, "cm"),
  xlog = FALSE,
  col = fpColors(box = "#4DBBD5", lines = "black", zero = "gray50"),
  lwd.ci = 2.5,
  ci.vertices = TRUE,
  ci.vertices.height = 0.02,
  clip = c(0, 30),
  xticks = seq(0, 30, by = 5),
  graph.pos = 4,
  graphwidth = unit(6, "cm"),
  txt_gp = forestplot::fpTxtGp(ticks = grid::gpar(cex = 0.9)),
)
forest_plot

####多分类分析
library(e1071)
library(plyr)
library(mldr)
library(utiml)

BRCA_Exp<-rbind(lncREXP_final,mREXP_final)
BRCA_Exp<-t(BRCA_Exp)

subtype_type <- cbind(rownames(subtype_col), subtype_col)
colnames(subtype_type)[2]="type"
colnames(subtype_type)[1]="sample"
rownames(subtype_type) <- NULL
##整理亚型特异性ceRNA模块列表
names(Normal_list)[1] <- "Normal-ceRM1"
names(Normal_list)[2] <- "Normal-ceRM2"

names(LumA_list)[1]<-"LumA-ceRM1"
names(LumA_list)[2]<-"LumA-ceRM2"
names(LumA_list)[3]<-"LumA-ceRM3"
names(LumA_list)[4]<-"LumA-ceRM4"
names(LumA_list)[5]<-"LumA-ceRM5"
names(LumA_list)[6]<-"LumA-ceRM6"
names(LumA_list)[7]<-"LumA-ceRM7"
names(LumA_list)[8]<-"LumA-ceRM8"
names(LumA_list)[9]<-"LumA-ceRM9"

names(Her2_list)[1]<-"Her2-ceRM1"
names(Her2_list)[2]<-"Her2-ceRM2"
names(Her2_list)[3]<-"Her2-ceRM3"

names(LumB_list)[1]<-"LumB-ceRM1"
names(LumB_list)[2]<-"LumB-ceRM2"
names(LumB_list)[3]<-"LumB-ceRM3"
names(LumB_list)[4]<-"LumB-ceRM4"

names(Basal_list)[1]<-"Basal-ceRM1"
names(Basal_list)[2]<-"Basal-ceRM2"
names(Basal_list)[3]<-"Basal-ceRM3"
names(Basal_list)[4]<-"Basal-ceRM4"

# ceRExp: ceRNA基因表达矩阵
# mRExp: mRNA基因表达矩阵
# BRCA_subtype: 带有样本乳腺癌亚型注释的数据框架
# Modulelist: 乳腺癌亚型特异性ceRNA模块列表集合
# Output: module_classifier为各个ceRNA模块的多类分类分析结果矩阵
module.classify <- function(ceRExp, mRExp, BRCA_subtype, Modulelist, method = "br", base.algorith = "SVM", cv.folds = 10,
                            cv.sampling = "stratified", cv.seed = 12345) {
  
  module_ceRExp <- lapply(seq_along(Modulelist), function(i) ceRExp[, which(colnames(ceRExp) %in% Modulelist[[i]])])
  module_mRExp <- lapply(seq_along(Modulelist), function(i) mRExp[, which(colnames(mRExp) %in% Modulelist[[i]])])
  Basal <- as.numeric(BRCA_subtype[, 2] == "Basal")
  Her2 <- as.numeric(BRCA_subtype[, 2] == "Her2")
  LumA <- as.numeric(BRCA_subtype[, 2] == "LumA")
  LumB <- as.numeric(BRCA_subtype[, 2] == "LumB")
  Normal <- as.numeric(BRCA_subtype[, 2] == "Normal")
  module_classify <- list()
  
  for (i in seq_along(Modulelist)){
    
    temp <- as.data.frame(cbind(module_ceRExp[[i]], module_mRExp[[i]], Basal, Her2, LumA, LumB, Normal))
    Indices <- ncol(temp)
    temp_mldr <- mldr_from_dataframe(temp, labelIndices = c(Indices-4, Indices-3, Indices-2, Indices-1, Indices), name = "TEMPMLDR")
    temp_res <- cv(temp_mldr, method = method, base.algorith = base.algorith, cv.folds = cv.folds,
                   cv.sampling = cv.sampling, cv.seed = cv.seed)
    module_classify[[i]] <- temp_res
    
  }
  
  return(module_classify)
}

lncR<-t(lncREXP_final)
mR<-t(mREXP_final)
Modulelist<-c(Normal_list,LumA_list,Her2_list,LumB_list,Basal_list)

moduel_Classification_baseline<-do.call(cbind,module.classify(lncR, mR, subtype_type, Modulelist, method = "baseline"))
moduel_Classification<-do.call(cbind,module.classify(lncR, mR, subtype_type, Modulelist, method = "br"))
colnames(moduel_Classification)=names(Modulelist)
colnames(moduel_Classification_baseline)=names(Modulelist)

####免疫浸润分析
library(GSVA)
Normal_exp_list <- list()
for (miRSM_name in names(Normal_ceRNA)) {
  lnc <- Normal_ceRNA[[miRSM_name]]$ceRNA
  lnc_in <- match(lnc, rownames(lncREXP_final))
  lncRExp <- lncREXP_final[lnc_in,]
  que <- any(is.na(lncRExp))
  
  mR <- Normal_ceRNA[[miRSM_name]]$mRNA
  mR_in <- match(mR, rownames(mREXP_final))
  mRExp <- mREXP_final[mR_in,]
  
  exp <- rbind(lncRExp, mRExp)
  
  Normal_exp_list[[miRSM_name]] <- exp
}
Normal_exp <- Normal_exp_list[[1]]
for (i in 2:length(Normal_exp_list)) {
  Normal_exp <- rbind(Normal_exp, Normal_exp_list[[i]])
}
# Normal_es: 正常ceRNA模块富集评分矩阵
Normal_para<-gsvaParam(Normal_exp,Normal_list)
Normal_es<-gsva(Normal_para,verbose=TRUE)
Normal_es_data<-as.data.frame(Normal_es)
rownames(Normal_es_data)=names(Normal_list)

# 包含每个LumA亚型ceRNA模块中基因表达矩阵的列表
LumA_exp_list <- list()
for (miRSM_name in names(LumA_ceRNA)) {
  Lu_lnc1 <- LumA_ceRNA[[miRSM_name]]$ceRNA
  Lu_lnc1_in <- match(Lu_lnc1, rownames(lncREXP_final))
  Lu_lncRExp1 <- lncREXP_final[Lu_lnc1_in,]
  que <- any(is.na(Lu_lncRExp1))
  
  Lu_mR1 <- LumA_ceRNA[[miRSM_name]]$mRNA
  Lu_mR1_in <- match(Lu_mR1, rownames(mREXP_final))
  Lu_mRExp1 <- mREXP_final[Lu_mR1_in,]
  
  Lu_exp1 <- rbind(Lu_lncRExp1, Lu_mRExp1)
  
  LumA_exp_list[[miRSM_name]] <- Lu_exp1
}

LumA_exp <- LumA_exp_list[[1]]
for (i in 2:length(LumA_exp_list)) {
  LumA_exp <- rbind(LumA_exp, LumA_exp_list[[i]])
}
# LumA_es: LumA亚型ceRNA模块富集评分矩阵
LumA_para<-gsvaParam(LumA_exp,LumA_list)
LumA_es<-gsva(LumA_para,verbose=TRUE)
LumA_es_data<-as.data.frame(LumA_es)
rownames(LumA_es_data)=names(LumA_list)

Her_exp_list <- list()
for (miRSM_name in names(Her2_ceRNA)) {
  Her_lnc1 <- Her2_ceRNA[[miRSM_name]]$ceRNA
  Her_lnc1_in <- match(Her_lnc1, rownames(lncREXP_final))
  Her_lncRExp1 <- lncREXP_final[Her_lnc1_in,]
  
  Her_mR1 <- Her2_ceRNA[[miRSM_name]]$mRNA
  Her_mR1_in <- match(Her_mR1, rownames(mREXP_final))
  Her_mRExp1 <- mREXP_final[Her_mR1_in,]
  
  Her_exp <- rbind(Her_lncRExp1, Her_mRExp1)
  
  Her_exp_list[[miRSM_name]] <- Her_exp
}

Her_exp <- Her_exp_list[[1]]
for (i in 2:length(Her_exp_list)) {
  Her_exp <- rbind(Her_exp, Her_exp_list[[i]])
}
# Her2_es:Her2亚型ceRNA模块富集评分矩阵
Her_para<-gsvaParam(Her_exp,Her2_list)
Her2_es<-gsva(Her_para,verbose=TRUE)
Her2_es_data<-as.data.frame(Her2_es)
rownames(Her2_es_data)=names(Her2_list)

LumB_exp_list <- list()
for (miRSM_name in names(LumB_ceRNA)) {
  lnc <- LumB_ceRNA[[miRSM_name]]$ceRNA
  lnc_in <- match(lnc, rownames(lncREXP_final))
  lncRExp <- lncREXP_final[lnc_in,]
  
  mR <- LumB_ceRNA[[miRSM_name]]$mRNA
  mR_in <- match(mR, rownames(mREXP_final))
  mRExp <- mREXP_final[mR_in,]
  
  exp <- rbind(lncRExp, mRExp)
  
  LumB_exp_list[[miRSM_name]] <- exp
}

LumB_exp <- LumB_exp_list[[1]]
for (i in 2:length(LumB_exp_list)) {
  LumB_exp <- rbind(LumB_exp, LumB_exp_list[[i]])
}
# LumB_es: LumB亚型ceRNA模块富集评分矩阵
LumB_para<-gsvaParam(LumB_exp,LumB_list)
LumB_es<-gsva(LumB_para,verbose=TRUE)
LumB_es_data<-as.data.frame(LumB_es)
rownames(LumB_es_data)=names(LumB_list)

Basal_exp_list <- list()
for (miRSM_name in names(Basal_ceRNA)) {
  lnc <- Basal_ceRNA[[miRSM_name]]$ceRNA
  lnc_in <- match(lnc, rownames(lncREXP_final))
  lncRExp <- lncREXP_final[lnc_in,]
  
  mR <- Basal_ceRNA[[miRSM_name]]$mRNA
  mR_in <- match(mR, rownames(mREXP_final))
  mRExp <- mREXP_final[mR_in,]
  
  exp <- rbind(lncRExp, mRExp)
  
  Basal_exp_list[[miRSM_name]] <- exp
}

Basal_exp <- Basal_exp_list[[1]]
for (i in 2:length(Basal_exp_list)) {
  Basal_exp <- rbind(Basal_exp, Basal_exp_list[[i]])
}
# Basal_es: Basal亚型ceRNA模块的富集评分矩阵
Basal_para<-gsvaParam(Basal_exp,Basal_list)
Basal_es<-gsva(Basal_para,verbose=TRUE)
Basal_es_data<-as.data.frame(Basal_es)
rownames(Basal_es_data)=names(Basal_list)

# module_es: ceRNA模块富集评分矩阵
module_es<-rbind(Normal_es_data,LumA_es_data,Her2_es_data,LumB_es_data,Basal_es_data)

library(readxl)
imm_data <- read_excel("Immune cell markers.xlsx", col_names = TRUE)
imm_data<-as.data.frame(imm_data)
immune_geneset <- lapply(1:ncol(imm_data), function(col_index) {
  imm_data[, col_index]
})
for (i in 1:length(immune_geneset)) {
  immune_geneset[[i]] <- na.omit(immune_geneset[[i]])
}
names(immune_geneset)<-colnames(imm_data)

EXP<-rbind(mREXP_final,lncREXP_final)

# immune_geneset: 24种免疫细胞标志物的列表
# exp_geneSet: 24种免疫细胞的富集评分矩阵
exp_Param<-gsvaParam(EXP,immune_geneset,assay = NA_character_,annotation = NA_character_,
                     minSize = 1,maxSize = Inf,kcdf = "Gaussian",tau = 1,maxDiff = TRUE,absRanking = FALSE)
exp_geneSet<- gsva(exp_Param, verbose = TRUE)
exp_geneSet_t<-t(exp_geneSet)
module_es_t<-t(module_es)
## Calculate the correlation of two enrichment score matrixs
module_geneset_cor <- cor(module_es_t, exp_geneSet_t, use = "everything", method = "pearson")
library(pheatmap)
pheatmap(module_geneset_cor,cluster_rows = F,cluster_cols = F,
         color=colorRampPalette(c("#63baf8","white","#ff7c9c"))(100),
         show_colnames = T,
         show_rownames = T,
         fontsize = 13
)

################################################组织特异性cceRNA模块推断###################################################
####乳腺癌组织数据下载及预处理
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(limma)
library(edgeR)
library(org.Hs.eg.db)
library(clusterProfiler)
library(igraph)
library(tidyverse)

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab"
)
GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)
save(clinical.BCRtab.all,file = 'clinical_all.RData')

clinical <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")

# 查询RNA-seq基因表达数据
query_rna <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",  # 或者 "HTSeq - Counts"
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)
# 下载数据
GDCdownload(query_rna, method = "api", directory = "GDCdata")

RNAexp <- GDCprepare(query_rna,save = T,save.filename = "brcaExp.rda")

##加载RNA数据
load("brcaExp.rda")
#####mRNA
mrna <- data[data@rowRanges$gene_type=="protein_coding",]
assay_list <- assayNames(mrna)
mRNAMatrix <- assay(mrna, "unstranded")
gene_info <- rowData(mrna)
print(head(gene_info))
# 将矩阵转换为数据框，保留原始列名
mRNA <- data.frame(mRNAMatrix, check.names = FALSE)
# 添加基因名
mRNA$gene_name <- gene_info$gene_name
# 对相同基因名的表达量求和（解决一个基因多个转录本的问题）
mRNA_agg <- aggregate(mRNA, . ~ gene_name, "sum")
# 设置行名为基因名
rownames(mRNA_agg) <- mRNA_agg$gene_name
# 移除gene_name列（现在已经是行名了）
mRNA_agg <- mRNA_agg[, -which(colnames(mRNA_agg) == "gene_name")]
# 检查是否有缺失值
cat("缺失值数量:", sum(is.na(mRNA_agg)), "\n")
# 检查零表达基因
zero_genes <- rowSums(mRNA_agg == 0) == ncol(mRNA_agg)
cat("在所有样本中零表达的基因数量:", sum(zero_genes), "\n")
mRNA_agg_count<-mRNA_agg
save(mRNA_agg_count,file = 'mRNA_agg_count.RData')

#####lncRNA
lncrna <- data[data@rowRanges$gene_type=="lncRNA",]
lncRNAMatrix <- assay(lncrna, "unstranded")
gene_info <- rowData(lncrna)
print(head(gene_info))
# 将矩阵转换为数据框，保留原始列名
lncRNA <- data.frame(lncRNAMatrix, check.names = FALSE)
# 添加基因名
lncRNA$gene_name <- gene_info$gene_name
# 对相同基因名的表达量求和（解决一个基因多个转录本的问题）
lncRNA_agg <- aggregate(lncRNA, . ~ gene_name, "sum")
# 设置行名为基因名
rownames(lncRNA_agg) <- lncRNA_agg$gene_name
# 移除gene_name列（现在已经是行名了）
lncRNA_agg <- lncRNA_agg[, -which(colnames(lncRNA_agg) == "gene_name")]
# 检查是否有缺失值
cat("缺失值数量:", sum(is.na(lncRNA_agg)), "\n")
# 检查零表达基因
zero_genes <- rowSums(lncRNA_agg == 0) == ncol(lncRNA_agg)
cat("在所有样本中零表达的基因数量:", sum(zero_genes), "\n")
lncRNA_agg_count<-lncRNA_agg
save(lncRNA_agg_count,file = 'lncRNA_agg_count.RData')

###fpkm
mrna_fpkm <- assay(mrna, "fpkm_unstrand")
gene_info_mrna <- rowData(mrna)
# 转换为数据框并去重
mrna_fpkm_df <- data.frame(mrna_fpkm, check.names = FALSE)
mrna_fpkm_df$gene_name <- gene_info_mrna$gene_name

mrna_fpkm_agg <- aggregate(mrna_fpkm_df, . ~ gene_name, "sum")
rownames(mrna_fpkm_agg) <- mrna_fpkm_agg$gene_name
mrna_fpkm_agg <- mrna_fpkm_agg[, -which(colnames(mrna_fpkm_agg) == "gene_name")]

cat("原始mRNA FPKM矩阵:", dim(mrna_fpkm_agg), "\n")

#lncRNA
lncrna_fpkm_raw <- assay(lncrna, "fpkm_unstrand")
gene_info_lncrna <- rowData(lncrna)

# 转换为数据框并去重
lncrna_fpkm_df <- data.frame(lncrna_fpkm_raw, check.names = FALSE)
lncrna_fpkm_df$gene_name <- gene_info_lncrna$gene_name

lncrna_fpkm_agg <- aggregate(lncrna_fpkm_df, . ~ gene_name, "sum")
rownames(lncrna_fpkm_agg) <- lncrna_fpkm_agg$gene_name
lncrna_fpkm_agg <- lncrna_fpkm_agg[, -which(colnames(lncrna_fpkm_agg) == "gene_name")]

cat("原始lncRNA FPKM矩阵:", dim(lncrna_fpkm_agg), "\n")

# 检查重叠基因
overlap_genes <- intersect(rownames(lncRNA_agg_count), rownames(mRNA_agg_count))
genes_to_remove_from_lncRNA <- overlap_genes  # 6个全部删除
genes_to_remove_from_mRNA <- c("ARMCX5-GPRASP2", "SFTA3")  # 只删除这2个
lncRNA_count_cleaned <- lncRNA_agg_count[!rownames(lncRNA_agg_count) %in% genes_to_remove_from_lncRNA, ]
mRNA_count_cleaned <- mRNA_agg_count[!rownames(mRNA_agg_count) %in% genes_to_remove_from_mRNA, ]
lncRNA_fpkm_cleaned <- lncrna_fpkm_agg[!rownames(lncrna_fpkm_agg) %in% genes_to_remove_from_lncRNA, ]
mRNA_fpkm_cleaned <- mrna_fpkm_agg[!rownames(mrna_fpkm_agg) %in% genes_to_remove_from_mRNA, ]
save(lncRNA_count_cleaned,mRNA_count_cleaned,file='lncRNA_mRNA_count_cleaned.RData')
save(lncRNA_fpkm_cleaned,mRNA_fpkm_cleaned,file = 'lncRNA_mRNA_fpkm_cleaned.RData')

##样本去重
colnames(mRNA_count_cleaned) <- substr(colnames(mRNA_count_cleaned), 1, 16)
colnames(lncRNA_count_cleaned) <- substr(colnames(lncRNA_count_cleaned), 1, 16)
colnames(lncRNA_fpkm_cleaned) <- substr(colnames(lncRNA_fpkm_cleaned), 1, 16)
colnames(mRNA_fpkm_cleaned) <- substr(colnames(mRNA_fpkm_cleaned), 1, 16)


#' 处理TCGA表达矩阵，识别样本类型并准备配对分析数据
#'
#' @param expression_matrix 表达矩阵，行是基因，列是样本(TCGA barcode格式)
#' @param handle_duplicates 是否处理重复的肿瘤样本(默认TRUE，取均值合并)
#' @param require_paired 是否只保留既有肿瘤又有正常样本的患者(默认TRUE)
#' @param tumor_code 肿瘤样本标识代码(默认"01A"，TCGA标准)
#' @param normal_code 正常样本标识代码(默认"11A"，TCGA标准)
#' @param verbose 是否显示处理过程信息(默认TRUE)
#'
#' @return 一个列表，包含处理后的各种数据
#' @export
#'
#' @examples
#' # 基本用法

library(dplyr)
library(tibble)

process_tcga_expression_pair <- function(expression_matrix, 
                                         handle_duplicates = TRUE,
                                         tumor_code = "01A",
                                         normal_code = "11A",
                                         verbose = TRUE) {
  
  # 辅助函数：从barcode提取样本类型
  extract_sample_type <- function(barcodes, tumor_code, normal_code) {
    sample_codes <- substr(barcodes, 14, 16)
    ifelse(
      sample_codes == tumor_code, "Tumor",
      ifelse(sample_codes == normal_code, "Normal", "Other")
    )
  }
  
  # 辅助函数：从barcode提取患者ID（前12位）
  extract_patient_id <- function(barcodes) {
    substr(barcodes, 1, 12)
  }
  
  if (verbose) {
    cat("开始处理TCGA表达矩阵...\n")
    cat("原始矩阵维度:", nrow(expression_matrix), "个基因 ×", ncol(expression_matrix), "个样本\n")
  }
  
  # 1. 识别样本类型
  sample_types <- extract_sample_type(colnames(expression_matrix), tumor_code, normal_code)
  
  if (verbose) {
    cat("样本类型分布:\n")
    print(table(sample_types))
  }
  
  # 2. 分离肿瘤和正常样本
  tumor_samples <- expression_matrix[, sample_types == "Tumor", drop = FALSE]
  normal_samples <- expression_matrix[, sample_types == "Normal", drop = FALSE]
  
  tumor_patients <- extract_patient_id(colnames(tumor_samples))
  normal_patients <- extract_patient_id(colnames(normal_samples))
  
  if (verbose) {
    cat("肿瘤样本数:", ncol(tumor_samples), "\n")
    cat("正常样本数:", ncol(normal_samples), "\n")
  }
  
  # 3. 处理重复的肿瘤样本（如果需要）
  if (handle_duplicates && ncol(tumor_samples) > 0) {
    
    # 检查是否有重复患者
    duplicated_patients <- tumor_patients[duplicated(tumor_patients)]
    unique_duplicated <- unique(duplicated_patients)
    
    if (length(unique_duplicated) > 0) {
      if (verbose) {
        cat("发现", length(unique_duplicated), "名患者有多个肿瘤样本，正在合并...\n")
      }
      
      # 创建包含完整样本名的数据框
      tumor_with_names <- as.data.frame(t(tumor_samples))
      tumor_with_names$patient_id <- tumor_patients
      tumor_with_names$full_sample_name <- colnames(tumor_samples)
      
      # 计算每个患者的最大值
      tumor_max <- tumor_with_names %>%
        group_by(patient_id) %>%
        summarise(
          across(-full_sample_name, ~ max(.x, na.rm = TRUE)),  
          representative_sample = full_sample_name[1],
          .groups = 'drop'
        ) %>%
        as.data.frame()
      
      # 准备最终矩阵
      representative_samples <- tumor_max$representative_sample
      tumor_max_final <- tumor_max %>%
        select(-representative_sample) %>%
        column_to_rownames("patient_id") %>%
        t() %>%
        as.matrix()
      
      colnames(tumor_max_final) <- representative_samples
      tumor_samples <- tumor_max_final
      tumor_patients <- extract_patient_id(colnames(tumor_samples))
      
      if (verbose) {
        cat("已合并重复样本，现肿瘤样本数:", ncol(tumor_samples), "\n")
      }
    }
  }
  
  # 4. 找出既有肿瘤又有正常样本的患者
  paired_patients <- intersect(tumor_patients, normal_patients)
  
  if (length(paired_patients) == 0) {
    warning("未找到配对患者，返回空矩阵")
    return(matrix(nrow = nrow(expression_matrix), ncol = 0))
  }
  
  if (verbose) {
    cat("配对患者数量:", length(paired_patients), "\n")
  }
  
  # 提取配对样本
  tumor_paired <- tumor_samples[, tumor_patients %in% paired_patients, drop = FALSE]
  normal_paired <- normal_samples[, normal_patients %in% paired_patients, drop = FALSE]
  
  # 确保样本顺序一致（按患者ID排序）
  tumor_patients_paired <- extract_patient_id(colnames(tumor_paired))
  normal_patients_paired <- extract_patient_id(colnames(normal_paired))
  
  tumor_paired <- tumor_paired[, order(tumor_patients_paired), drop = FALSE]
  normal_paired <- normal_paired[, order(normal_patients_paired), drop = FALSE]
  
  # 5. 合并为配对表达矩阵
  paired_expression <- cbind(tumor_paired, normal_paired)
  
  if (verbose) {
    cat("最终配对表达矩阵维度:", nrow(paired_expression), "个基因 ×", ncol(paired_expression), "个样本\n")
    cat("处理完成！\n")
  }
  
  return(paired_expression)
}

mRNA_countexp<-process_tcga_expression_pair(mRNA_count_cleaned)
lncRNA_countexp <- process_tcga_expression_pair(lncRNA_count_cleaned)
save(mRNA_countexp,lncRNA_countexp,file = 'mRNA_lncRNA_countexp.RData')
mRNA_fpkmexp<-process_tcga_expression_pair(mRNA_fpkm_cleaned)
lncRNA_fpkmtexp <- process_tcga_expression_pair(lncRNA_fpkm_cleaned)
save(mRNA_fpkmexp,lncRNA_fpkmtexp,file = 'mRNA_lncRNA_fpkmexp.RData')

#######下载miRNA数据
# 查询miRNA表达数据
query_mirna <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "miRNA Expression Quantification",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

# 下载miRNA数据
GDCdownload(query_mirna, method = "api", directory = "GDCdata_miRNA", files.per.chunk = 50)
# 处理miRNA数据
mirna_data <- GDCprepare(query_mirna, directory = "GDCdata_miRNA", save = TRUE, save.filename = "mirna_data.rda")

##加载miRNA数据
load("mirna_data.rda")
#miRNA
count_cols <- grep("read_count", colnames(data), value = TRUE)
mirna <- data.frame(data[, count_cols])
# 设置行名为miRNA_ID
rownames(mirna) <- data$miRNA_ID
# 简化列名（去掉前缀）
colnames(mirna) <- gsub("read_count_", "", count_cols)
mirna1 <- mirna
colnames(mirna1) <- substr(colnames(mirna1), 1, 16)

rpm_cols <- grep("reads_per_million", colnames(data), value = TRUE)
mirna_rpm <- data.frame(data[, rpm_cols])
# 设置行名为miRNA_ID
rownames(mirna_rpm) <- data$miRNA_ID

# 简化列名（去掉前缀）
colnames(mirna_rpm) <- gsub("reads_per_million_miRNA_mapped_", "", rpm_cols)
colnames(mirna_rpm) <- substr(colnames(mirna_rpm), 1, 16)

miRNA_countexp<-process_tcga_expression_pair(mirna1)
miRNA_rpmexp<-process_tcga_expression_pair(mirna_rpm)
save(miRNA_countexp,miRNA_rpmexp,file = 'miRNA_rpm_countexp.RData')

mR_index<-match(colnames(miRNA_countexp),colnames(mRNA_countexp))
mRNAexp_count<-mRNA_countexp[,mR_index]
lncRNAexp_count<-lncRNA_countexp[,mR_index]
save(mRNAexp_count,lncRNAexp_count,file = 'mR_lncRNAexp_count.RData')
mRNAexp_fpkm<-mRNA_fpkmexp[,mR_index]
lncRNAexp_fpkm<-lncRNA_fpkmtexp[,mR_index]
save(mRNAexp_fpkm,lncRNAexp_fpkm,file = 'mR_lncRNAexp_fpkm.RData')

#####################差异表达分析
####limma
library(limma)
library(ggplot2) #用于绘制火山图
library(DESeq2)
library(edgeR)

gp <- substring(colnames(mRNAexp_count), 14, 15)
table(gp)
brca_tumor <- mRNAexp_count[, as.numeric(gp) < 10]
brca_normal <- mRNAexp_count[, as.numeric(gp) >= 10]
group <- c(rep('tumor', ncol(brca_tumor)), rep('normal', ncol(brca_normal)))
group <- factor(group, levels = c("normal", "tumor"))

####DESeq2（Differential Expression using DESeq2）
colData <- data.frame(row.names = colnames(mRNAexp_count),
                      group = group)
colData$group <- factor(colData$group, levels = c("normal", "tumor"))
head(colData)
dds <- DESeqDataSetFromMatrix(countData = mRNAexp_count, # 表达矩阵
                              colData = colData,        # 表达矩阵列名和分组信息的对应关系
                              design = ~ group)         # group为colData中的group，也就是分组信息
# 查看一下构建好的矩阵
head(dds)
# 进行差异表达分析
dds <- DESeq(dds)
# 提取差异表达结果，进行对比，这里contrast参数指定了对比的组别
# contrast参数必须写成下面三个元素的向量格式，且顺序不能反
res <- results(dds, contrast = c("group", rev(levels(group))))
# 按照padj（调整后的p值）的大小对差异结果进行排序（只有DESeq2需要，limma和edgeR会自动排好）
resOrdered <- res[order(res$padj), ]

# 将差异表达结果转换为数据框
DEG <- as.data.frame(resOrdered)

# 去除缺失值，如果没有这一步，一些表达量很低的基因计算后会出现NA，给后续分析和绘图带来麻烦，远离麻烦！
mRNA_DEG_deseq <- na.omit(DEG)

# 将处理后的差异表达结果保存为R数据文件
save(mRNA_DEG_deseq, file = 'mRNA_DEG_deseq.Rdata')

#' 使用DESeq2进行差异表达分析
#'
#' @param count_matrix 表达矩阵，行是基因，列是样本（必须是整数counts）
#' @param group 分组向量，长度等于样本数，通常为"normal"和"tumor"
#' @param reference_level 参考组（分母），默认为"normal"
#'
#' @return 差异表达结果数据框（已按padj排序并去除NA）
#' @export
#'
#' @examples
#' # 基本用法
#' deg_results <- run_deseq2_simple(mRNAexp_count, group)
#'
#' # 指定参考组为"tumor"
#' deg_results <- run_deseq2_simple(mRNAexp_count, group, reference_level = "tumor")

run_deseq2_simple <- function(count_matrix, 
                              group, 
                              reference_level = "normal") {
  
  # 加载DESeq2包
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("请先安装DESeq2包: BiocManager::install('DESeq2')")
  }
  library(DESeq2)
  
  # 检查分组信息
  if (!is.factor(group)) {
    group <- factor(group)
  }
  
  # 设置参考水平
  if (!reference_level %in% levels(group)) {
    stop("参考组 '", reference_level, "' 不在分组水平中")
  }
  group <- relevel(group, ref = reference_level)
  
  # 创建colData
  colData <- data.frame(row.names = colnames(count_matrix),
                        group = group)
  
  # 创建DESeq2对象并运行分析
  dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                colData = colData,
                                design = ~ group)
  
  dds <- DESeq(dds)
  
  # 提取结果
  contrast_levels <- rev(levels(group))
  res <- results(dds, contrast = c("group", contrast_levels))
  
  # 按padj排序并去除NA
  resOrdered <- res[order(res$padj), ]
  DEG <- as.data.frame(resOrdered)
  DEG_clean <- na.omit(DEG)
  
  return(DEG_clean)
}

##lncRNA
lncRNA_DEG_deseq <- run_deseq2_simple(lncRNAexp_count, group)
miRNA_DEG_deseq <- run_deseq2_simple(miRNA_countexp, group)
save(mRNA_DEG_deseq,lncRNA_DEG_deseq,miRNA_DEG_deseq,file = 'mR_lncR_miRNA_DEG_deseq.RData')

logFC = 1
P.Value = 0.05
k1 <- (mRNA_DEG_deseq$pvalue < P.Value) & (mRNA_DEG_deseq$log2FoldChange < -logFC)
k2 <- (mRNA_DEG_deseq$pvalue < P.Value) & (mRNA_DEG_deseq$log2FoldChange > logFC)
mRNA_DEG_deseq <- mutate(mRNA_DEG_deseq, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(mRNA_DEG_deseq$change)

k1 <- (lncRNA_DEG_deseq$pvalue < P.Value) & (lncRNA_DEG_deseq$log2FoldChange < -logFC)
k2 <- (lncRNA_DEG_deseq$pvalue < P.Value) & (lncRNA_DEG_deseq$log2FoldChange > logFC)
lncRNA_DEG_deseq <- mutate(lncRNA_DEG_deseq, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(lncRNA_DEG_deseq$change)

k1 <- (miRNA_DEG_deseq$pvalue < P.Value) & (miRNA_DEG_deseq$log2FoldChange < -logFC)
k2 <- (miRNA_DEG_deseq$pvalue < P.Value) & (miRNA_DEG_deseq$log2FoldChange > logFC)
miRNA_DEG_deseq <- mutate(miRNA_DEG_deseq, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(miRNA_DEG_deseq$change)

#######edgeR（Empirical Analysis of Digital Gene Expression Data in R）
# 创建 DGEList 对象，用于存储基因表达数据和组信息，还是使用原始计数矩阵
d <- DGEList(counts = mRNAexp_count, group = group)

# 根据每个基因在每个样本中的 CPM（Counts Per Million）值去除低表达基因
keep <- rowSums(cpm(d) > 1) >= 2

# 或者自动过滤，去除低表达基因
# keep <- filterByExpr(d)

# 从 DGEList 对象中筛选出符合条件的基因
d <- d[keep, , keep.lib.sizes = FALSE]

# 更新样本的库大小信息
d$samples$lib.size <- colSums(d$counts)

# 归一化，TMM 方法
d <- calcNormFactors(d)
dge = d
# 创建设计矩阵，用于指定差异分析模型
design <- model.matrix(~0 + factor(group))
rownames(design) <- colnames(dge)
colnames(design) <- levels(factor(group))

# 估计数据的离散度 —— common离散度、trended离散度、tagwise离散度
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

# 在估计的模型基础上进行 广义线性模型 (GLM) 拟合
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, contrast = c(-1, 1))

# 从 LRT 计算结果中获取前 nrow(dge) 个顶部差异表达基因
nrDEG <- topTags(lrt, n = nrow(dge))

# 将差异表达基因结果转换为数据框形式
DEG_edgeR <- as.data.frame(nrDEG)

logFC = 1
P.Value = 0.05
k1 <- (DEG_edgeR$PValue < P.Value) & (DEG_edgeR$logFC < -logFC)
k2 <- (DEG_edgeR$PValue < P.Value) & (DEG_edgeR$logFC > logFC)
DEG_edgeR <- mutate(DEG_edgeR, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(DEG_edgeR$change)
mRNA_DEG_edgeR<-DEG_edgeR

run_edger_analysis <- function(count_matrix, 
                               group, 
                               logFC_threshold = 1, 
                               pvalue_threshold = 0.05,
                               use_filterByExpr = TRUE) {
  
  # 检查并加载edgeR包
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("请先安装edgeR包: BiocManager::install('edgeR')")
  }
  library(edgeR)
  library(dplyr)
  
  # 创建DGEList对象
  d <- DGEList(counts = count_matrix, group = group)
  
  # 过滤低表达基因
  if (use_filterByExpr) {
    keep <- filterByExpr(d)
  } else {
    keep <- rowSums(cpm(d) > 1) >= 2
  }
  
  d <- d[keep, , keep.lib.sizes = FALSE]
  d$samples$lib.size <- colSums(d$counts)
  
  # 归一化
  d <- calcNormFactors(d)
  
  # 创建设计矩阵
  design <- model.matrix(~0 + factor(group))
  rownames(design) <- colnames(d)
  colnames(design) <- levels(factor(group))
  
  # 估计离散度
  d <- estimateGLMCommonDisp(d, design)
  d <- estimateGLMTrendedDisp(d, design)
  d <- estimateGLMTagwiseDisp(d, design)
  
  # 拟合GLM模型并进行LRT检验
  fit <- glmFit(d, design)
  lrt <- glmLRT(fit, contrast = c(-1, 1))  # 假设只有两组，比较第二组vs第一组
  
  # 提取结果
  nrDEG <- topTags(lrt, n = nrow(d))
  DEG_edgeR <- as.data.frame(nrDEG)
  
  # 标记上下调基因
  k_down <- (DEG_edgeR$PValue < pvalue_threshold) & (DEG_edgeR$logFC < -logFC_threshold)
  k_up <- (DEG_edgeR$PValue < pvalue_threshold) & (DEG_edgeR$logFC > logFC_threshold)
  
  DEG_edgeR <- DEG_edgeR %>%
    mutate(change = case_when(
      k_up ~ "up",
      k_down ~ "down",
      TRUE ~ "stable"
    ))
  
  return(DEG_edgeR)
}

lncRNA_DEG_edgeR <- run_edger_analysis(lncRNAexp_count, group)
table(lncRNA_DEG_edgeR$change)
miRNA_DEG_edgeR <- run_edger_analysis(miRNA_countexp, group)
table(miRNA_DEG_edgeR$change)
save(mRNA_DEG_edgeR,lncRNA_DEG_edgeR,miRNA_DEG_edgeR,file = 'mR_lnR_miRNA_DEG_edgeR.RData')


######limma
design <- model.matrix(~0 + factor(group))
colnames(design) = levels(factor(group))
rownames(design) = colnames(mRNAexp_count)

# 创建 DGEList 对象
dge <- DGEList(counts = mRNAexp_count)

# 归一化，得到的归一化系数被用作文库大小的缩放系数
dge <- calcNormFactors(dge)

# 使用 voom 方法进行标准化
v <- voom(dge, design, plot = TRUE, normalize = "quantile")

# 如果是芯片数据、TPM数据或已标准化的数据，不需要再进行标准化，可直接从这里开始进行差异分析

# 使用线性模型进行拟合
fit <- lmFit(v, design)

# 和上面两个包一样，需要说明是谁比谁
con <- paste(rev(levels(group)), collapse = "-")
con
cont.matrix <- makeContrasts(contrasts = c(con), levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# 获取差异表达基因结果
tempOutput <- topTable(fit2, coef = con, n = Inf)
mRNA_DEG_limma_voom <- na.omit(tempOutput)
k1 <- (mRNA_DEG_limma_voom$P.Value < P.Value) & (mRNA_DEG_limma_voom$logFC < -logFC)
k2 <- (mRNA_DEG_limma_voom$P.Value < P.Value) & (mRNA_DEG_limma_voom$logFC > logFC)
mRNA_DEG_limma_voom <- mutate(mRNA_DEG_limma_voom, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(mRNA_DEG_limma_voom$change)


#' 使用limma-voom进行差异表达分析
#'
#' @param count_matrix 表达矩阵，行是基因，列是样本
#' @param group 分组向量，通常为"normal"和"tumor"
#' @param normalize_method 归一化方法，默认为"TMM"
#' @param plot_voom 是否绘制voom图，默认为FALSE
#' @param logFC_threshold log2倍数变化阈值，默认为1
#' @param pvalue_threshold P值阈值，默认为0.05
#'
#' @return 差异表达结果数据框，包含change列标记上下调
#' @export
#'
#' @examples
#' # 基本用法
#' DEG_limma <- run_limma_voom(mRNAexp_count, group)

run_limma_voom <- function(count_matrix, 
                           group, 
                           normalize_method = "TMM",
                           plot_voom = FALSE,
                           logFC_threshold = 1,
                           pvalue_threshold = 0.05) {
  
  # 检查并加载必要的包
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("请先安装limma包: BiocManager::install('limma')")
  }
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("请先安装edgeR包: BiocManager::install('edgeR')")
  }
  library(limma)
  library(edgeR)
  library(dplyr)
  
  # 创建DGEList对象
  dge <- DGEList(counts = count_matrix)
  
  # 归一化
  dge <- calcNormFactors(dge, method = normalize_method)
  
  # 创建设计矩阵
  design <- model.matrix(~0 + factor(group))
  rownames(design) <- colnames(dge)
  colnames(design) <- levels(factor(group))
  
  # 使用voom进行标准化
  v <- voom(dge, design, plot = plot_voom, normalize = "quantile")
  
  # 拟合线性模型
  fit <- lmFit(v, design)
  
  # 创建对比矩阵
  con <- paste(rev(levels(factor(group))), collapse = "-")
  cont.matrix <- makeContrasts(contrasts = c(con), levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  
  # 获取差异表达结果
  tempOutput <- topTable(fit2, coef = con, n = Inf)
  DEG_limma <- na.omit(tempOutput)
  
  # 标记上下调基因
  k_down <- (DEG_limma$P.Value < pvalue_threshold) & (DEG_limma$logFC < -logFC_threshold)
  k_up <- (DEG_limma$P.Value < pvalue_threshold) & (DEG_limma$logFC > logFC_threshold)
  
  DEG_limma <- DEG_limma %>%
    mutate(change = case_when(
      k_up ~ "up",
      k_down ~ "down",
      TRUE ~ "stable"
    ))
  
  return(DEG_limma)
}

mRNA_DEG_limma <- run_limma_voom(mRNAexp_count, group)
table(mRNA_DEG_limma$change)
lncRNA_DEG_limma <- run_limma_voom(lncRNAexp_count, group)
table(lncRNA_DEG_limma$change)
miRNA_DEG_limma <- run_limma_voom(miRNA_countexp, group)
table(miRNA_DEG_limma$change)
save(mRNA_DEG_limma,lncRNA_DEG_limma,miRNA_DEG_limma,file = 'mR_lncR_miRNA_DEG_limma.RData')

##########################################################两两交集的总和###################################################
# 函数：获取三种方法两两交集的并集（不考虑上下调方向）
get_shared_diff_genes_no_direction <- function(deseq_res, edger_res, limma_res, rna_type = "RNA") {
  
  cat("\n", strrep("=", 60), "\n")
  cat(rna_type, "差异表达基因 - 至少两种方法共有（不考虑方向）\n")
  cat(strrep("=", 60), "\n")
  
  # 1. 提取各方法的差异基因（不考虑方向，只要不是"stable"就是差异基因）
  deseq_diff <- rownames(deseq_res[deseq_res$change != "stable", ])
  edger_diff <- rownames(edger_res[edger_res$change != "stable", ])
  limma_diff <- rownames(limma_res[limma_res$change != "stable", ])
  
  # 2. 计算两两交集
  cat("\n--- 计算两两交集 ---\n")
  
  deseq_edger <- intersect(deseq_diff, edger_diff)
  cat("DESeq2 和 edgeR 的交集:", length(deseq_edger), "个基因\n")
  
  deseq_limma <- intersect(deseq_diff, limma_diff)
  cat("DESeq2 和 limma 的交集:", length(deseq_limma), "个基因\n")
  
  edger_limma <- intersect(edger_diff, limma_diff)
  cat("edgeR 和 limma 的交集:", length(edger_limma), "个基因\n")
  
  # 3. 计算所有两两交集的并集（至少两种方法共有的基因）
  cat("\n--- 至少两种方法共有的基因 ---\n")
  
  # 取所有两两交集的并集
  shared_genes <- unique(c(deseq_edger, deseq_limma, edger_limma))
  cat("至少两种方法共有的基因总数:", length(shared_genes), "\n")
  
  # 4. 计算每种方法支持的次数
  cat("\n--- 方法支持统计 ---\n")
  
  # 创建支持矩阵
  support_matrix <- data.frame(
    gene = shared_genes,
    stringsAsFactors = FALSE
  )
  
  # 添加方法支持信息
  support_matrix$in_DESeq2 <- ifelse(support_matrix$gene %in% deseq_diff, "yes", "no")
  support_matrix$in_edgeR <- ifelse(support_matrix$gene %in% edger_diff, "yes", "no")
  support_matrix$in_limma <- ifelse(support_matrix$gene %in% limma_diff, "yes", "no")
  
  # 计算支持方法数量
  support_matrix$support_count <- apply(
    support_matrix[, c("in_DESeq2", "in_edgeR", "in_limma")], 
    1, 
    function(x) sum(x == "yes")
  )
  
  # 统计支持方法数量的分布
  support_dist <- table(support_matrix$support_count)
  cat("\n基因支持方法数量分布:\n")
  for (i in 2:3) {
    if (i == 2) {
      cat("支持2种方法的基因:", ifelse(is.na(support_dist[as.character(i)]), 0, support_dist[as.character(i)]), "\n")
    } else if (i == 3) {
      cat("支持3种方法的基因:", ifelse(is.na(support_dist[as.character(i)]), 0, support_dist[as.character(i)]), "\n")
    }
  }
  
  # 5. 计算各方法独有的差异基因
  cat("\n--- 各方法独有的差异基因 ---\n")
  
  deseq_only <- setdiff(deseq_diff, union(edger_diff, limma_diff))
  edger_only <- setdiff(edger_diff, union(deseq_diff, limma_diff))
  limma_only <- setdiff(limma_diff, union(deseq_diff, edger_diff))
  
  cat("DESeq2 独有:", length(deseq_only), "个基因\n")
  cat("edgeR 独有:", length(edger_only), "个基因\n")
  cat("limma 独有:", length(limma_only), "个基因\n")
  
  # 6. 创建交集组合信息
  cat("\n--- 交集组合详情 ---\n")
  
  # 计算各种交集组合的基因
  only_deseq_edger <- setdiff(deseq_edger, limma_diff)
  only_deseq_limma <- setdiff(deseq_limma, edger_diff)
  only_edger_limma <- setdiff(edger_limma, deseq_diff)
  all_three <- intersect(intersect(deseq_diff, edger_diff), limma_diff)
  
  cat("仅 DESeq2 和 edgeR 共有（不包括limma）:", length(only_deseq_edger), "\n")
  cat("仅 DESeq2 和 limma 共有（不包括edgeR）:", length(only_deseq_limma), "\n")
  cat("仅 edgeR 和 limma 共有（不包括DESeq2）:", length(only_edger_limma), "\n")
  cat("三种方法共有:", length(all_three), "\n")
  
  # 7. 统计信息汇总
  cat("\n--- 统计汇总 ---\n")
  
  total_diff_genes <- unique(c(deseq_diff, edger_diff, limma_diff))
  cat("DESeq2 差异基因数:", length(deseq_diff), "\n")
  cat("edgeR 差异基因数:", length(edger_diff), "\n")
  cat("limma 差异基因数:", length(limma_diff), "\n")
  cat("所有方法差异基因总数（去重）:", length(total_diff_genes), "\n")
  cat("\n至少两种方法共有的基因占比:", 
      round(length(shared_genes) / length(total_diff_genes) * 100, 2), "%\n")
  
  # 8. 返回结果
  result <- list(
    # 主要结果：至少两种方法共有的基因
    shared_genes = shared_genes,
    
    # 支持矩阵
    support_matrix = support_matrix,
    
    # 各方法独有的基因
    unique_genes = list(
      DESeq2_only = deseq_only,
      edgeR_only = edger_only,
      limma_only = limma_only
    ),
    
    # 详细交集组合
    intersection_details = list(
      only_deseq_edger = only_deseq_edger,
      only_deseq_limma = only_deseq_limma,
      only_edger_limma = only_edger_limma,
      all_three = all_three
    ),
    
    # 原始差异基因列表
    all_diff_genes = list(
      DESeq2 = deseq_diff,
      edgeR = edger_diff,
      limma = limma_diff
    ),
    
    # 统计信息
    statistics = list(
      DESeq2_total = length(deseq_diff),
      edgeR_total = length(edger_diff),
      limma_total = length(limma_diff),
      total_unique = length(total_diff_genes),
      shared_genes_count = length(shared_genes),
      shared_percentage = round(length(shared_genes) / length(total_diff_genes) * 100, 2),
      support_2_methods = sum(support_matrix$support_count == 2),
      support_3_methods = sum(support_matrix$support_count == 3),
      DESeq2_only = length(deseq_only),
      edgeR_only = length(edger_only),
      limma_only = length(limma_only),
      only_deseq_edger = length(only_deseq_edger),
      only_deseq_limma = length(only_deseq_limma),
      only_edger_limma = length(only_edger_limma),
      all_three = length(all_three)
    ),
    
    # 参数信息
    parameters = list(
      rna_type = rna_type,
      analysis_type = "shared_diff_genes_no_direction",
      criteria = "至少两种方法共有的差异基因，不考虑上下调方向"
    )
  )
  
  # 9. 保存到文件
  output_file <- paste0(rna_type, "_shared_diff_genes_no_direction.RData")
  save(result, file = output_file)
  cat("\n结果已保存到:", output_file, "\n")
  
  # 10. 输出文本文件
  # 保存共享基因列表
  if (length(shared_genes) > 0) {
    write.table(shared_genes, 
                file = paste0(rna_type, "_shared_genes.txt"),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    write.csv(support_matrix, 
              file = paste0(rna_type, "_shared_genes_details.csv"),
              row.names = FALSE)
  }
  
  # 保存独有基因列表
  if (length(deseq_only) > 0) {
    write.table(deseq_only, 
                file = paste0(rna_type, "_deseq_only_genes.txt"),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  if (length(edger_only) > 0) {
    write.table(edger_only, 
                file = paste0(rna_type, "_edger_only_genes.txt"),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  if (length(limma_only) > 0) {
    write.table(limma_only, 
                file = paste0(rna_type, "_limma_only_genes.txt"),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  # 11. 输出建议
  cat("\n", strrep("-", 60), "\n")
  cat("分析建议:\n")
  cat("1. 使用 'shared_genes' 获取至少两种方法共有的差异基因\n")
  cat("2. 使用 'support_matrix' 查看每个基因被哪些方法支持\n")
  cat("3. 共享基因已保存到文本文件，可用于后续富集分析\n")
  cat("4. 共有", result$statistics$shared_genes_count, "个基因被至少两种方法共同检测到\n")
  cat(strrep("-", 60), "\n")
  
  return(result)
}

mRNA_shared <- get_shared_diff_genes_no_direction(
  deseq_res = mRNA_DEG_deseq,
  edger_res = mRNA_DEG_edgeR,
  limma_res = mRNA_DEG_limma,
  rna_type = "mRNA"
)

lncRNA_shared <- get_shared_diff_genes_no_direction(
  deseq_res = lncRNA_DEG_deseq,
  edger_res = lncRNA_DEG_edgeR,
  limma_res = lncRNA_DEG_limma,
  rna_type = "lncRNA"
)

miRNA_shared <- get_shared_diff_genes_no_direction(
  deseq_res = miRNA_DEG_deseq,
  edger_res = miRNA_DEG_edgeR,
  limma_res = miRNA_DEG_limma,
  rna_type = "miRNA"
)

mR_index<-match(mRNA_shared[["shared_genes"]],rownames(mRNAexp_fpkm))
mR_shared_fpkm<-mRNAexp_fpkm[mR_index,]

lncR_index<-match(lncRNA_shared[["shared_genes"]],rownames(lncRNAexp_fpkm))
lncR_shared_fpkm<-lncRNAexp_fpkm[lncR_index,]

miR_index<-match(miRNA_shared[["shared_genes"]],rownames(miRNA_rpmexp))
miR_shared_rpm<-miRNA_rpmexp[miR_index,]
#miRNA ID转化
library(miRBaseConverter)
library(dplyr)
miR_row_shared <- rownames(miR_shared_rpm)
miR_mature_trans <- miRNA_PrecursorToMature(miR_row_shared, version = "v21")

# 2. 将数据框转换为长格式，同时考虑Mature1和Mature2
miR_mature_long <- miR_mature_trans %>%
  pivot_longer(
    cols = c(Mature1, Mature2),
    names_to = "MatureType",
    values_to = "MatureID"
  ) %>%
  filter(!is.na(MatureID)) %>%  # 移除NA值
  distinct()  # 去重

# 3. 合并表达矩阵
merged_data <- miR_shared_rpm %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Precursor") %>%
  inner_join(miR_mature_long, by = c("Precursor" = "OriginalName"))

# 4. 统计每个成熟miRNA的来源前体数量
mature_summary <- merged_data %>%
  group_by(MatureID) %>%
  summarise(
    n_precursors = n_distinct(Precursor),
    .groups = 'drop'
  )

#5. 如果有多个前体对应同一个成熟miRNA，可以选择：求和
miR_final_max <- merged_data %>%
  group_by(MatureID) %>%
  summarise(across(
    where(is.numeric),
    ~ max(.x, na.rm = TRUE)  # 取最大值而不是平均值
  ))

miR_rpm_shared <- as.matrix(miR_final_max[, -1])
rownames(miR_rpm_shared)<-miR_final_max$MatureID

miR_shared_se<-SummarizedExperiment(t(miR_rpm_shared))
mR_shared_se<-SummarizedExperiment(t(mR_shared_fpkm))
lncR_shared_se<-SummarizedExperiment(t(lncR_shared_fpkm))
miRtar_se<-SummarizedExperiment(miR_Target)
save(miR_shared_se,mR_shared_se,lncR_shared_se,miRtar_se,file = 'share_sedata.RData')

########################################乳腺癌组织特异性ceRNA模块推断#######################################
library(miRSM)
set.seed(1234)
Tumor_modulegenes_WGCNA_shared <- module_WGCNA(lncR_shared_se[index[[1]], ],mR_shared_se[index[[1]], ])
Tumor_miRSM_WGCNA_SDC_shared <- miRSM(miR_shared_se[index[[1]], ], lncR_shared_se[index[[1]], ], mR_shared_se[index[[1]], ], miRtar_se,
                                      Tumor_modulegenes_WGCNA_shared,
                                      num_shared_miRNAs = 3, pvalue.cutoff = 0.05,
                                      method = "SDC")
save(Tumor_miRSM_WGCNA_SDC_shared,file='Tumor_miRSM_WGCNA_SDC_shared.RData')

Normal_modulegenes_WGCNA_shared  <- module_WGCNA(lncR_shared_se[index[[2]], ],mR_shared_se[index[[2]], ])
Normal_miRSM_WGCNA_SDC_shared  <- miRSM(miR_shared_se[index[[2]], ], lncR_shared_se[index[[2]], ], mR_shared_se[index[[2]], ], miRtar_se,
                                        Normal_modulegenes_WGCNA_shared,
                                        num_shared_miRNAs = 3, pvalue.cutoff = 0.05,
                                        method = "SDC")
save(Normal_miRSM_WGCNA_SDC_shared,file='Normal_miRSM_WGCNA_SDC_shared.RData')

##################################乳腺癌组织特异性ceRNA模块功能解析##############################################
sample_names <- colnames(mRNAexp_count)
Tumor <- grep('0[1-9]', substr(sample_names, 14, 16))
Normal <- grep('1[0-9]', substr(sample_names, 14, 16))
group_filtered <- factor(c(rep("Tumor", length(Tumor)), rep("Normal", length(Normal))))
print(table(group_filtered))

deg_res <- diff_analysis(
  exprset = mRNAexp_count,
  is_count = TRUE,
  logFC_cut = 1.5,
  pvalue_cut = 0.01,
  adjpvalue_cut = 0.05,
  group = group_filtered, # 正确设置分组信息
) 

##将肿瘤和正常样本的结果合并
ceRNA_result<-list("Normal" = Normal_miRSM_WGCNA_SDC_shared$`miRNA sponge modules`,
                   "Tumor" = Tumor_miRSM_WGCNA_SDC_shared$`miRNA sponge modules`)

sample_mapping <- list(
  "Normal" = "Normal",
  "Tumor" = "Tumor"
)

#####查看重叠模块
count_overlap_modules <- function(mod_list1, mod_list2) {
  # 创建模块标识符的函数
  create_id <- function(mod) {
    paste0(
      paste(sort(mod$ceRNA), collapse = "|"),
      "##",
      paste(sort(mod$mRNA), collapse = "|")
    )
  }
  
  # 生成标识符
  ids1 <- sapply(mod_list1, create_id)
  ids2 <- sapply(mod_list2, create_id)
  
  # 计算重叠
  overlap <- intersect(ids1, ids2)
  
  return(list(
    list1_count = length(mod_list1),
    list2_count = length(mod_list2),
    overlap_count = length(overlap),
    overlap_ids = overlap
  ))
}

overlap_stats <- count_overlap_modules(
  Normal_miRSM_WGCNA_SDC_shared$`miRNA sponge modules`,
  Tumor_miRSM_WGCNA_SDC_shared$`miRNA sponge modules`
)


##相似性分析
library(miRSM)
sim_Normal_Tumor<-module_group_sim(Normal_miRSM_WGCNA_SDC_shared$`miRNA sponge modules`, Tumor_miRSM_WGCNA_SDC_shared$`miRNA sponge modules`, sim.method = "Simpson")

##富集分析
library(openxlsx)
library(miRSM)
#功能富集分析
miRSM_Normal_FEA <- module_FA(Normal_miRSM_WGCNA_SDC_shared$`miRNA sponge modules`, Analysis.type = "FEA")
Nor_FEA_BP<-miRSM_Normal_FEA[[1]]
for (i in 1:length(Nor_FEA_BP)) {
  Nor_FEA_BP_data <- data.frame(Nor_FEA_BP[[i]])
  
  file_name <- paste0("Normal_", i, "FEA_BP_result.xlsx")
  write.xlsx(Nor_FEA_BP_data, file = file_name, rowNames = FALSE)
}
Nor_FEA_KEGG<-miRSM_Normal_FEA[[2]]
for (i in 1:length(Nor_FEA_KEGG)) {
  Nor_FEA_KEGG_data <- data.frame(Nor_FEA_KEGG[[i]])
  
  file_name <- paste0("Normal_", i, "FEA_KEGG_result.xlsx")
  write.xlsx(Nor_FEA_KEGG_data, file = file_name, rowNames = FALSE)
}
Nor_FEA_Reactome<-miRSM_Normal_FEA[[3]]
for (i in 1:length(Nor_FEA_Reactome)) {
  Nor_FEA_Reactome_data <- data.frame(Nor_FEA_Reactome[[i]])
  
  file_name <- paste0("Normal_", i, "FEA_Reactome_result.xlsx")
  write.xlsx(Nor_FEA_Reactome_data, file = file_name, rowNames = FALSE)
}

miRSM_Tumor_FEA <- module_FA(Tumor_miRSM_WGCNA_SDC_shared$`miRNA sponge modules`, Analysis.type = "FEA")
Tum_FEA_BP<-miRSM_Tumor_FEA[[1]]
for (i in 1:length(Tum_FEA_BP)) {
  Tum_FEA_BP_data <- data.frame(Tum_FEA_BP[[i]])
  
  file_name <- paste0("Tumor_", i, "FEA_BP_result.xlsx")
  write.xlsx(Tum_FEA_BP_data, file = file_name, rowNames = FALSE)
}
Tum_FEA_KEGG<-miRSM_Tumor_FEA[[2]]
for (i in 1:length(Tum_FEA_KEGG)) {
  Tum_FEA_KEGG_data <- data.frame(Tum_FEA_KEGG[[i]])
  
  file_name <- paste0("Tumor_", i, "FEA_KEGG_result.xlsx")
  write.xlsx(Tum_FEA_KEGG_data, file = file_name, rowNames = FALSE)
}
Tum_FEA_Reactome<-miRSM_Tumor_FEA[[3]]
for (i in 1:length(Tum_FEA_Reactome)) {
  Tum_FEA_Reactome_data <- data.frame(Tum_FEA_Reactome[[i]])
  
  file_name <- paste0("Tumor_", i, "FEA_Reactome_result.xlsx")
  write.xlsx(Tum_FEA_Reactome_data, file = file_name, rowNames = FALSE)
}
#疾病富集分析
miRSM_Normal_DEA <- module_FA(Normal_miRSM_WGCNA_SDC_shared$`miRNA sponge modules`, Analysis.type = "DEA")
Nor_DEA_HDO<-miRSM_Normal_DEA[[1]]
for (i in 1:length(Nor_DEA_HDO)) {
  Nor_DEA_HDO_data <- data.frame(Nor_DEA_HDO[[i]])
  
  file_name <- paste0("Normal_", i, "DEA_HDO_result.xlsx")
  write.xlsx(Nor_DEA_HDO_data, file = file_name, rowNames = FALSE)
}
Nor_DEA_DisGeNET<-miRSM_Normal_DEA[[2]]
for (i in 1:length(Nor_DEA_DisGeNET)) {
  Nor_DEA_DisGeNET_data <- data.frame(Nor_DEA_DisGeNET[[i]])
  
  file_name <- paste0("Normal_", i, "DEA_DisGeNET_result.xlsx")
  write.xlsx(Nor_DEA_DisGeNET_data, file = file_name, rowNames = FALSE)
}
Nor_DEA_NCG<-miRSM_Normal_DEA[[3]]
for (i in 1:length(Nor_DEA_NCG)) {
  Nor_DEA_NCG_data <- data.frame(Nor_DEA_NCG[[i]])
  
  file_name <- paste0("Normal_", i, "DEA_NCG_result.xlsx")
  write.xlsx(Nor_DEA_NCG_data, file = file_name, rowNames = FALSE)
}

miRSM_Tumor_DEA <- module_FA(Tumor_miRSM_WGCNA_SDC_shared$`miRNA sponge modules`, Analysis.type = "DEA")
Tum_DEA_HDO<-miRSM_Tumor_DEA[[1]]
for (i in 1:length(Tum_DEA_HDO)) {
  Tum_DEA_HDO_data <- data.frame(Tum_DEA_HDO[[i]])
  
  file_name <- paste0("Tumor_", i, "DEA_HDO_result.xlsx")
  write.xlsx(Tum_DEA_HDO_data, file = file_name, rowNames = FALSE)
}
Tum_DEA_DisGeNET<-miRSM_Tumor_DEA[[2]]
for (i in 1:length(Tum_DEA_DisGeNET)) {
  Tum_DEA_DisGeNET_data <- data.frame(Tum_DEA_DisGeNET[[i]])
  
  file_name <- paste0("Tumor_", i, "DEA_DisGeNET_result.xlsx")
  write.xlsx(Tum_DEA_DisGeNET_data, file = file_name, rowNames = FALSE)
}
Tum_DEA_NCG<-miRSM_Tumor_DEA[[3]]
for (i in 1:length(Tum_DEA_NCG)) {
  Tum_DEA_NCG_data <- data.frame(Tum_DEA_NCG[[i]])
  
  file_name <- paste0("Tumor_", i, "DEA_NCG_result.xlsx")
  write.xlsx(Tum_DEA_NCG_data, file = file_name, rowNames = FALSE)
}
#富集通路数汇总
library(readxl)
library(ggplot2)
gongneng_term<-read_excel("富集分析/功能富集分析通路数量.xlsx",col_names = TRUE)
ggplot(data = gongneng_term, aes(x = Modules, y = `Number of significantly enriched terms`, fill = `sample type`)) +
  geom_bar(stat = "identity", width = 0.8, position = position_dodge(width = 0.9)) +
  geom_text(aes(label = `Number of significantly enriched terms`), vjust = -0.5, color = "black", size = 4.5) +  # 设置字体为新罗马字体
  facet_grid(rows = vars(ontology)) +
  scale_y_continuous(limits = c(0, max(gongneng_term$`Number of significantly enriched terms`) * 1.1)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14))

jibing_term<-read_excel("富集分析/疾病富集分析通路数量 .xlsx",col_names = TRUE)

ggplot(data = jibing_term, aes(x = Modules, y = `Number of significantly enriched terms`, fill = `sample type`)) +
  geom_bar(stat = "identity", width = 0.8, position = position_dodge(width = 0.9)) +
  geom_text(aes(label = `Number of significantly enriched terms`), vjust = -0.5, color = "black", size = 4.5) +  # 设置字体为新罗马字体
  facet_grid(rows = vars(ontology)) +
  scale_y_continuous(limits = c(0, max(jibing_term$`Number of significantly enriched terms`) * 1.1)) +
  theme_bw() +
  theme( # 去掉所有背景网格线
    axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),  # 设置字体为新罗马字体
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14))

#分布分析
library(UpSetR)
Nor_share<-share_miRs(miR_shared_se,miRtar_se,ceRNA_result[["Normal"]])
Tumor_share<-share_miRs(miR_shared_se,miRtar_se,ceRNA_result[["Tumor"]])

Nor_tu_share<-c(Nor_share,Tumor_share)
names(Nor_tu_share)[1] <- "Normal-ceRM1"
names(Nor_tu_share)[2] <- "Normal-ceRM2"
names(Nor_tu_share)[3] <- "Normal-ceRM3"
names(Nor_tu_share)[4] <- "Normal-ceRM4"
names(Nor_tu_share)[5] <- "Tumor-ceRM1"
Normal_tur_share_upset<-upset(fromList(Nor_tu_share),
                              sets=names(Nor_tu_share),
                              point.size=3.5,
                              line.size=1.5,
                              mainbar.y.label="Number of shared miRNAs",
                              main.bar.color = 'black',
                              matrix.color="black",
                              sets.x.label="Set size",
                              sets.bar.color=c('red', 'orange','yellow','green','blue'),
                              order.by = "freq",
                              text.scale=c(3,3,2,2.5,2.5,3)
)
Normal_tur_share_upset




#RNA分布
RNA_data <- read_excel("正常和肿瘤样本模块的RNA分布情况.xlsx" ,sheet="Normal_Tur")
RNA_plot<-ggplot(RNA_data, aes(fill = RNA类型,
                               y = 基因数量, x = ceRNA模块))+
  geom_bar(position = "stack", stat = "identity", width = 0.4)+
  geom_text(aes(label = 基因数量), position = position_stack(vjust = 0.5), color = "black", size = 4) +
  ggtitle("正常样本")+
  #scale_fill_manual(values = c("pink", "skyblue")) +  # 修改颜色
  theme(plot.title = element_text(hjust = 0.5))
RNA_plot


##生存分析
survival_matched$sample <- id_mapping_df$full_sample_id[matched_indices]
lnc<-assay(lncR_shared_se)
mR<-assay(mR_shared_se)
miRNA<-assay(miR_shared_se)
exp<-cbind(mR,lnc)

# 注意：只保留肿瘤样本（样本类型代码为‘01’）
all_sample_ids <- rownames(lnc)

# 创建映射数据框
id_map <- data.frame(
  full_sample_id = all_sample_ids,
  patient_id = substr(all_sample_ids, 1, 12),
  stringsAsFactors = FALSE
)

# 筛选肿瘤样本：第14-15位字符为‘01’
id_map <- id_map[substr(id_map$full_sample_id, 14, 15) == "01", ]

# 检查映射表：一个患者ID是否对应了多个肿瘤样本？
dup_patients <- id_map$patient_id[duplicated(id_map$patient_id)]
if(length(dup_patients) > 0) {
  cat("注意：发现", length(dup_patients), "个患者有多个肿瘤样本。\n")
  cat("示例：\n")
  print(head(id_map[id_map$patient_id %in% dup_patients, ]))
  # 默认处理：保留每个患者的第一个肿瘤样本
  id_map <- id_map[!duplicated(id_map$patient_id), ]
}

# 2. 将生存数据中的患者ID匹配替换为完整的肿瘤样本ID
# 使用match函数进行精确匹配
match_idx <- match(final_survival_df$sample, id_map$patient_id)

# 3. 创建新的生存数据框，替换ID并移除匹配失败的行
survival_matched <- final_survival_df
survival_matched$sample <- id_map$full_sample_id[match_idx]  # 替换为完整样本ID
survival_matched <- survival_matched[!is.na(survival_matched$sample), ]  # 移除未匹配到的
##将生存时间单位改为月份
survival_matched$time <- round(survival_matched$time / 30.44, 2)

row<-survival_matched$sample
row_index<-match(row,rownames(exp))
expdata<-exp[row_index,]

Tumor_list <- list()
for (i in 1:length(Tumor_miRSM_WGCNA_SDC_shared$`miRNA sponge modules`)) {
  Tumor_list[[i]] <- c(Tumor_miRSM_WGCNA_SDC_shared$`miRNA sponge modules`[[i]]$ceRNA, Tumor_miRSM_WGCNA_SDC_shared$`miRNA sponge modules`[[i]]$mRNA)
}
library(miRspongeR)

# 加载必要的包
library(survival)
library(survminer)

# 假设已得到以下对象：
# Tumor_list, expdata, survival_matched

# 提取模块基因
module_genes <- Tumor_list[[1]]

# 提取模块基因表达数据
exp_module <- expdata[, colnames(expdata) %in% module_genes, drop = FALSE]

# 构建 Cox 模型数据框
cox_data <- data.frame(
  time = survival_matched$time,
  status = survival_matched$status,
  exp_module
)

# 拟合多因素 Cox 模型
mm <- coxph(Surv(time, status) ~ ., data = cox_data)

# 计算风险评分
risk_score <- predict(mm, type = "risk")

# 按中位数分组（确保因子顺序为 "High" 在前，"Low" 在后，与 moduleSurvival 一致）
group <- ifelse(risk_score > median(risk_score), "High", "Low")
group <- factor(group, levels = c("High", "Low"))  # 重要：High 为第一组

# 创建生存数据框
surv_df <- data.frame(
  time = survival_matched$time,
  status = survival_matched$status,
  group = group
)

# 拟合 Kaplan-Meier 曲线
fit <- survfit(Surv(time, status) ~ group, data = surv_df)

# 计算 Log-rank 检验
sdf <- survdiff(Surv(time, status) ~ group, data = surv_df)

# 提取 Log-rank p 值
logrank_p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)

# 按照 moduleSurvival 的方法计算 HR 和 95% CI
obs <- sdf$obs
exp <- sdf$exp
HR <- (obs[1] / exp[1]) / (obs[2] / exp[2])  # High 组为第1组，Low 组为第2组
hr_low <- exp(log(HR) - qnorm(0.975) * sqrt(1/exp[1] + 1/exp[2]))
hr_up  <- exp(log(HR) + qnorm(0.975) * sqrt(1/exp[1] + 1/exp[2]))

# 获取各组样本数
n_high <- sdf$n[1]
n_low <- sdf$n[2]

# 使用 ggsurvplot 绘图
gg <- ggsurvplot(
  fit,
  data = surv_df,
  pval = FALSE,                # 不自动添加 p 值（我们手动添加）
  conf.int = TRUE,              # 显示置信区间
  risk.table = FALSE,             # 显示风险表
  risk.table.col = "strata",
  palette = c("#E7B800", "#2E9FDF"),  # 高风险黄色，低风险蓝色
  xlab = "Time (months)",
  ylab = "Survival probability",
  legend.title = "Risk group",
  legend.labs = c(paste0("High (n=", n_high, ")"), 
                  paste0("Low (n=", n_low, ")")),
  ggtheme = theme_classic()
)

# 在图上添加 HR 和 p 值标注
gg$plot <- gg$plot + 
  annotate("text", 
           x = max(surv_df$time) * 0.7, 
           y = 0.2,
           label = paste0("HR = ", round(HR, 2), 
                          " (", round(hr_low, 2), "-", round(hr_up, 2), ")\n",
                          "Log-rank p = ", format.pval(logrank_p, digits = 3)),
           hjust = 0, size = 4, color = "black")

# 显示图形
print(gg)

colnames(Tu_sur)[3]<-"HR (95% CI)"
Tu_sur<-as.data.frame(Tu_sur)
Tu_sur$Module <- names(Tumor_list)
all_cols <- colnames(Tu_sur)
# 将Module列移到最前面，其他列保持原有顺序
new_order <- c("Module", setdiff(all_cols, "Module"))
Tu_sur <- Tu_sur[, new_order]

tabletext <- cbind(
  c("Module", Tu_sur$Module),
  c("Chi-square",round(Tu_sur$`Chi-square`,2)),
  c("HR (95% CI)", paste0(round(Tu_sur$`HR (95% CI)`, 2), " (", round(Tu_sur$HRlow95, 2), " - ", round(Tu_sur$HRup95, 2), ")")),
  c("p-value", Tu_sur$`p-value`)
)
library(forestplot)
forest_plot <- forestplot(
  labeltext = tabletext,
  mean = c(NA, Tu_sur$`HR (95% CI)`),
  lower = c(NA, Tu_sur$HRlow95),
  upper = c(NA, Tu_sur$HRup95),
  hrzl_lines = list(`1` = grid::gpar(lwd = 2, col = "black"),
                    `2` = grid::gpar(lwd = 2, col = "black"),
                    `3` = grid::gpar(lwd = 2, col = "black")),
  is.summary = rep(FALSE, nrow(Tu_sur) + 1),
  zero = 1,
  boxsize = 0.2,
  lineheight = unit(0.7, "cm"),
  xlog = FALSE,
  col = fpColors(box = "#4DBBD5", lines = "black", zero = "gray50"),
  lwd.ci = 2.5,
  ci.vertices = TRUE,
  ci.vertices.height = 0.02,
  clip = c(0, 30),
  xticks = seq(0, 30, by = 5),
  graph.pos = 4,
  graphwidth = unit(6, "cm"),
  txt_gp = forestplot::fpTxtGp(ticks = grid::gpar(cex = 0.9)),
)
forest_plot



######二分类分析
library(e1071)
library(plyr)
library(mldr)
library(utiml)

Normal_list <- list()
for (i in 1:length(Normal_miRSM_WGCNA_SDC_shared$`miRNA sponge modules`)) {
  Normal_list[[i]] <- c(Normal_miRSM_WGCNA_SDC_shared$`miRNA sponge modules`[[i]]$ceRNA, Normal_miRSM_WGCNA_SDC_shared$`miRNA sponge modules`[[i]]$mRNA)
}
Tumor_list <- list()
for (i in 1:length(Tumor_miRSM_WGCNA_SDC_shared$`miRNA sponge modules`)) {
  Tumor_list[[i]] <- c(Tumor_miRSM_WGCNA_SDC_shared$`miRNA sponge modules`[[i]]$ceRNA, Tumor_miRSM_WGCNA_SDC_shared$`miRNA sponge modules`[[i]]$mRNA)
}
Modulelist<-c(Tumor_list,Normal_list)

sample_names <- rownames(lnc)  # 或 rownames(mRExp)

# 根据TCGA命名规则提取样本类型
# 第14-15位表示样本类型：01-09表示肿瘤，10-19表示正常
sample_labels <- ifelse(substr(sample_names, 14, 15) %in% sprintf("%02d", 1:9), 
                        "Tumor", "Normal")

# 转换为因子
sample_labels <- factor(sample_labels, levels = c("Normal", "Tumor"))


# 加载必要的库
library(e1071)      # SVM
library(caret)      # 交叉验证和模型评估
library(pROC)       # ROC曲线
library(ggplot2)    # 可视化
library(dplyr)      # 数据处理
library(tidyr)      # 数据整理

# 假设你已经有了以下数据：
# Modulelist: ceRNA模块基因列表
# ceRExp: lncRNA表达矩阵（行：样本，列：lncRNA）
# mRExp: mRNA表达矩阵（行：样本，列：mRNA）
# 样本标签：tumor vs normal

# 1. 准备样本标签（这里假设你已经有标签数据）
# 请根据你的实际情况调整
# sample_labels <- factor(c(rep("Tumor", n_tumor), rep("Normal", n_normal)))

# 2. 整合所有ceRNA模块的基因
all_ceRNA_genes <- unique(unlist(Modulelist))

# 3. 从表达矩阵中提取ceRNA基因的表达数据
# 假设ceRExp和mRExp的行名是样本ID，列名是基因名
# 首先合并lncRNA和mRNA表达矩阵
# 注意：确保两个矩阵的样本顺序一致
combined_exp <- cbind(lnc, mR)

# 提取ceRNA模块中存在的基因
ceRNA_exp <- combined_exp[, colnames(combined_exp) %in% all_ceRNA_genes]


# 4. 准备完整的数据集
full_data <- data.frame(
  Class = sample_labels,  # 样本标签
  ceRNA_exp              # 表达数据
)

# 5. 对每个模块分别进行分析的函数
analyze_module <- function(module_name, module_genes, exp_matrix, labels) {
  cat("\n分析模块:", module_name, "\n")
  cat("基因数量:", length(module_genes), "\n")
  
  # 提取该模块的基因
  module_exp <- exp_matrix[, colnames(exp_matrix) %in% module_genes]
  
  # 如果基因数量太少，跳过
  if(ncol(module_exp) < 5) {
    cat("基因数量太少，跳过此模块\n")
    return(NULL)
  }
  
  # 准备数据
  data_df <- data.frame(
    Class = labels,
    module_exp
  )
  
  # 移除有缺失值的样本
  data_df <- na.omit(data_df)
  
  # 设置十折交叉验证
  set.seed(123)
  folds <- createFolds(data_df$Class, k = 10, list = TRUE, returnTrain = FALSE)
  
  # 存储每折的结果
  results <- list()
  predictions <- data.frame()
  
  for(fold in 1:10) {
    cat("正在进行第", fold, "折交叉验证...\n")
    
    # 划分训练集和测试集
    test_idx <- folds[[fold]]
    train_idx <- setdiff(1:nrow(data_df), test_idx)
    
    train_data <- data_df[train_idx, ]
    test_data <- data_df[test_idx, ]
    
    # 标准化数据（基于训练集）
    preproc <- preProcess(train_data[, -1], method = c("center", "scale"))
    train_scaled <- predict(preproc, train_data)
    test_scaled <- predict(preproc, test_data)
    
    # 训练SVM模型
    # 可以调整参数，这里使用默认的径向基核函数
    svm_model <- svm(Class ~ ., 
                     data = train_scaled,
                     kernel = "radial",
                     probability = TRUE)
    
    # 预测
    pred <- predict(svm_model, test_scaled[, -1], probability = TRUE)
    pred_prob <- attr(pred, "probabilities")[, "Tumor"]  # 假设Tumor是正类
    
    # 保存结果
    fold_results <- data.frame(
      Fold = fold,
      True = test_scaled$Class,
      Predicted = pred,
      Probability = pred_prob,
      Module = module_name
    )
    
    predictions <- rbind(predictions, fold_results)
    
    # 计算本折的性能指标
    cm <- confusionMatrix(pred, test_scaled$Class)
    
    results[[fold]] <- list(
      Accuracy = cm$overall["Accuracy"],
      Sensitivity = cm$byClass["Sensitivity"],
      Specificity = cm$byClass["Specificity"],
      F1 = cm$byClass["F1"],
      AUC = NA  # 将在后面计算
    )
  }
  
  # 计算AUC
  roc_obj <- roc(predictions$True, predictions$Probability, levels = c("Normal", "Tumor"))
  auc_value <- auc(roc_obj)
  
  # 汇总结果
  summary_results <- data.frame(
    Module = module_name,
    Accuracy = mean(sapply(results, function(x) x$Accuracy), na.rm = TRUE),
    Sensitivity = mean(sapply(results, function(x) x$Sensitivity), na.rm = TRUE),
    Specificity = mean(sapply(results, function(x) x$Specificity), na.rm = TRUE),
    F1_Score = mean(sapply(results, function(x) x$F1), na.rm = TRUE),
    AUC = auc_value,
    Num_Genes = ncol(module_exp),
    Num_Samples = nrow(data_df)
  )
  
  return(list(
    summary = summary_results,
    predictions = predictions,
    roc = roc_obj,
    module_genes = colnames(module_exp)
  ))
}

# 6. 对所有模块进行分析
all_results <- list()

for(module_name in names(Modulelist)) {
  module_result <- analyze_module(
    module_name = module_name,
    module_genes = Modulelist[[module_name]],
    exp_matrix = ceRNA_exp,
    labels = sample_labels
  )
  
  if(!is.null(module_result)) {
    all_results[[module_name]] <- module_result
  }
}

# 7. 汇总所有模块的结果
summary_table <- do.call(rbind, lapply(all_results, function(x) x$summary))
save(summary_table,file = '分类分析/classify summary_table.RData')
library(writexl)
library(tibble)
# 保存为Excel文件
# 将矩阵转换为数据框，并将行名添加为一列
classdf <- as.data.frame(summary_table)
classdf <- rownames_to_column(classdf, var = "RowNames")

# 写入Excel
write_xlsx(classdf, path = "分类分析/classify.xlsx")

# 绘制ROC曲线
roc_plot_data <- data.frame()
for(module_name in names(all_results)) {
  roc_df <- data.frame(
    Module = module_name,
    FPR = 1 - all_results[[module_name]]$roc$specificities,
    TPR = all_results[[module_name]]$roc$sensitivities
  )
  roc_plot_data <- rbind(roc_plot_data, roc_df)
}
# 首先加载RColorBrewer包
library(RColorBrewer)

# 查看所有可用的调色板
display.brewer.all()

####绘制ROC曲线
roc_curve_plot <- ggplot(roc_plot_data, aes(x = FPR, y = TPR, color = Module)) +
  geom_line(linewidth = 1.2) +
  geom_abline(linetype = "dashed", color = "gray") +
  labs(title = "各ceRNA模块ROC曲线",
       x = "假阳性率 (1 - Specificity)",
       y = "真阳性率 (Sensitivity)") +
  theme_minimal()
#scale_color_brewer(palette = "Set2")

print(roc_curve_plot)

performance_plot <- summary_table %>%
  select(Module, Accuracy, Sensitivity, Specificity, AUC) %>%
  gather(key = "Metric", value = "Value", -Module) %>%
  ggplot(aes(x = Module, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "各ceRNA模块分类性能比较",
       x = "模块",
       y = "性能指标值") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(performance_plot)

######免疫浸润分析                                                                                                                                                             ##免疫浸润分析
library(GSVA)
library(SummarizedExperiment)
lncRa<-assay(lncR_shared_se)
lnct<-t(lncRa)
mRa<-assay(mR_shared_se)
mRt<-t(mRa)

names(Normal_list)[1] <- "Normal-ceRM1"
names(Normal_list)[2] <- "Normal-ceRM2"
names(Normal_list)[3] <- "Normal-ceRM3"
names(Normal_list)[4] <- "Normal-ceRM4"

names(Tumor_list)[1] <- "Tumor-ceRM1"
Modulelist<-c(Tumor_list,Normal_list)

Normal_exp_list <- list()
for (miRSM_name in names(Normal_miRSM_WGCNA_SDC_shared[["miRNA sponge modules"]])) {
  lnc <- Normal_miRSM_WGCNA_SDC_shared[["miRNA sponge modules"]][[miRSM_name]]$ceRNA
  lnc_in <- match(lnc, rownames(lnct))
  lncRExp <- lnct[lnc_in,]
  que <- any(is.na(lncRExp))
  
  mR <- Normal_miRSM_WGCNA_SDC_shared[["miRNA sponge modules"]][[miRSM_name]]$mRNA
  mR_in <- match(mR, rownames(mRt))
  mRExp <- mRt[mR_in,]
  
  exp <- rbind(lncRExp, mRExp)
  
  Normal_exp_list[[miRSM_name]] <- exp
}
Normal_exp <- Normal_exp_list[[1]]
for (i in 2:length(Normal_exp_list)) {
  Normal_exp <- rbind(Normal_exp, Normal_exp_list[[i]])
}
Normal_para<-gsvaParam(Normal_exp,Normal_list)
Normal_es<-gsva(Normal_para,verbose=TRUE)
Normal_es_data<-as.data.frame(Normal_es)

mdule_exp_list<-c(Normal_exp_list,Tumor_exp_list)
mdule_exp<-exp_list[[1]]
for (i in 2:length(mdule_exp_list)) {
  mdule_exp <- rbind(mdule_exp, mdule_exp_list[[i]])
}
mdule_exp_para<-gsvaParam(mdule_exp,Modulelist)
mdule_exp_es<-gsva(exp_para,verbose=TRUE)
module_es_t<-t(mdule_exp_es)


library(readxl)
imm_data <- read_excel("E:/experiment/data download/BRCA/Scene/data/Immune cell markers.xlsx", col_names = TRUE)
imm_data<-as.data.frame(imm_data)

immune_geneset <- lapply(1:ncol(imm_data), function(col_index) {
  imm_data[, col_index]
})
for (i in 1:length(immune_geneset)) {
  immune_geneset[[i]] <- na.omit(immune_geneset[[i]])
}
names(immune_geneset)<-colnames(imm_data)

EXP<-rbind(mRt,lnct)

exp_Param<-gsvaParam(EXP,immune_geneset,assay = NA_character_,annotation = NA_character_,
                     minSize = 1,maxSize = Inf,kcd = "Gaussian",tau = 1,maxDiff = TRUE,absRanking = FALSE)
exp_geneSet<- gsva(exp_Param, verbose = TRUE)
exp_geneSet_t<-t(exp_geneSet)
module_es_t<-t(module_es)
## Calculate the correlation of two enrichment score matrixs
module_geneset_cor <- cor(module_es_t, exp_geneSet_t, use = "everything", method = "pearson")
library(pheatmap)
pheatmap(module_geneset_cor,cluster_rows = F,cluster_cols = F,
         color=colorRampPalette(c("#63baf8","white","#ff7c9c"))(100),
         show_colnames = T,
         show_rownames = T,
         fontsize = 13
)
library(writexl)
library(tibble)
# 保存为Excel文件
# 将矩阵转换为数据框，并将行名添加为一列
df <- as.data.frame(module_geneset_cor)
df <- rownames_to_column(df, var = "RowNames")

# 写入Excel
write_xlsx(df, path = "免疫浸润/module_geneset_correlation.xlsx")