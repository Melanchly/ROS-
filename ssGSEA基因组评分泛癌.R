setwd("D:/Dim/TCGAplot/UQC")
rm(list = ls())
source("../ssGSEA-Pan-Cancer.R")

cox_genes1 <- c(
  "COX6A1", "COX4I2", "COX7A2", "COX6B1", "COX7C", "COX4I1", "COX7B",
  "COX5B", "COX6A2", "COX6B2", "COX7A1", "COX6C","COX7B2", "COX8A",
  "COX5A", "COX8C", "MT-CO2", "MT-CO1", "MT-CO3"
)

cox_genes2 = c(
  "UQCRB","UQCRQ","UQCRC1","UQCRC2","MT-CYB",
  "CYC1","UQCRFS1","UQCRH","UQCR10","UQCR11"
)
cox_genes3 = c(
  "NDUFAB1","NDUFA1","NDUFA2","NDUFA3","NDUFA5","NDUFA6","NDUFA7",
  "NDUFA8","NDUFA9","NDUFA10","NDUFA11","NDUFA12","NDUFA13",
  "NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7",
  "NDUFB8","NDUFB9","NDUFB10","NDUFB11",
  "NDUFC1","NDUFC2",
  "NDUFS4","NDUFS5","NDUFS6","NDUFV3"
)


# cox_genes2 <- c(
#   "GSTO2", "GSTP1", "GSTZ1", "GSTT2B", "GSTM1", "GSTM5", "GSTM3",
#   "GSTO1", "HPGDS", "GSTM4", "GSTA4", "GSTA3", "GSTA5", "GSTK1",
#   "GSTM2", "GSTA1", "GSTA2", "GSTT4"
# )
# 
# 
# cox_genes3 <- c(
#   "ATP5F1D", "ATP5F1B", "ATP5PB", "ATP5F1E", "ATP5IF1", "ATP5MC2",
#   "ATP5F1A", "ATP5MC3", "ATP5PF", "ATP5MC1", "ATP5F1C", "ATP5MG",
#   "ATP5PD", "ATP5ME", "MT-ATP6", "MT-ATP8", "ATP5MF", "ATP5PO"
# )
# 
# cox_genes4 <- c(
#   "TRPM4", "TRPM2", "TRPM8", "TRPC4", "TRPC3", "TRPV3",
#   "MLKL", "TRPC4AP", "TRPC6", "TRPC1", "TRPV5", "TRPM1",
#   "TRPV2", "TRPM5", "MCOLN1", "TRPM7", "TRPV4", "RIPK3",
#   "RIPK1", "TRPV6", "MCOLN3", "TRPA1", "MCOLN2", "TRPM3",
#   "TRPM6", "TRPV1"
# )
# 
# cox_genes5 <- c("MAPK1", "MAPK3", "MAPK4", "MAPK6", "MAPK7", "MAPK8",
#   "MAPK9", "MAPK10", "MAPK11", "MAPK12", "MAPK13", "MAPK14",
#   "MAPK15", "NLK"
# )


cox_list <- list(cox_genes2, cox_genes3)

# cox_list <- list(cox_genes)
for (cox_genes in cox_list) {
  ssGsea_Score(tpm, cox_genes)
}
# 执行计算基因组ssGsea得分##############
ssGsea_Score(tpm, cox_genes)


#基因家族差异表达LogFC########
plot_diff_heatmap(rt)
#TCGA+GTEx
cox_genes <- c(
  "COX6A1", "COX4I2", "COX7A2", "COX6B1", "COX7C", "COX4I1", "COX7B",
  "COX5B", "COX6A2", "COX6B2", "COX7A1", "COX6C","COX7B2", "COX8A",
  "COX5A", "COX8C", "MT-CO2", "MT-CO1", "MT-CO3"
)
cox_genes = c(
  "UQCRB","UQCRQ","UQCRC1","UQCRC2","MT-CYB",
  "CYC1","UQCRFS1","UQCRH","UQCR10","UQCR11"
)
cox_genes = c(
  "NDUFAB1","NDUFA1","NDUFA2","NDUFA3","NDUFA5","NDUFA6","NDUFA7",
  "NDUFA8","NDUFA9","NDUFA10","NDUFA11","NDUFA12","NDUFA13",
  "NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7",
  "NDUFB8","NDUFB9","NDUFB10","NDUFB11",
  "NDUFC1","NDUFC2",
  "NDUFS4","NDUFS5","NDUFS6","NDUFV3"
)

#TCGA+GTEx
#基因家族差异表达LogFC########
plot_diff_heatmap1(dir_name)




plot_OS_KM_expression(rt)
inputFile <- paste0(dir_name, "_expTime.txt")
rt_cox = read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
plot_cox_analysis_forest(rt_cox)# 调用函数获取每个基因的COX.txt文件


#COX生存分析热图########
folder_path <- "生存分析(COX回归)"
plot_cox_analysis_forest_picture(folder_path)#循环上一步的cox.txt绘制每个基因的森林图
plot_cox_heatmap1(folder_path)#无统计
plot_cox_heatmap2(folder_path)#有统计


#突变汇总########
folder_path <- "../突变汇总"
Alterations(folder_path)



#CNV和甲基化与mRNA的关系###############
analysis_target <- "CNV"
analysis_target <- "Methylation"
CNV_Methylation_mRNA(analysis_target)
#CNV_Methylation_mRNA1(analysis_target)


#GSVA得分正常-肿瘤对照箱线图######
ssGsea_Score_diff(df)
ssGsea_Score_diff1(df)

#GSVA得分纯肿瘤趋势箱线图######
ssGsea_Score_diff_Tumor(df)


#GSVA泛癌TMB######  
ssGsea_Score_pan_cancer_TMB(df)


#GSVA泛癌MSI######  
ssGsea_Score_pan_cancer_MSI(df)


#基因集的泛癌种与免疫相关基因的相关性分析#########
#基因集的泛癌种与ICG的相关性#############
ssGsea_Score_pan_cancer_ICG(df,tpm)


#基因集的泛癌种与趋势化因子的相关性#############
ssGsea_Score_pan_cancer_chemokine(df,tpm)


#基因集的泛癌种与趋化因子受体的相关性#############
ssGsea_Score_pan_cancer_chemokine_receptor(df,tpm)


#基因集的泛癌种与免疫刺激剂的相关性############
ssGsea_Score_pan_cancer_immune_stimulator(df,tpm)


#基因集的泛癌种与免疫抑制剂的相关性############
ssGsea_Score_pan_cancer_immune_inhibitor(df,tpm)


#基因集的泛癌种与免疫细胞比率相关性########
ssGsea_Score_pan_cancer_immunecellratio(df)


#基因集的泛癌种与免疫评分的相关性############
ssGsea_Score_pan_cancer_immunescore(df)


#基因集的泛癌种与TME_Pathaway通路评分的相关性############
ssGsea_Score_pan_cancer_TME_GSVA(df)


#基因集的泛癌种与HALLMARK通路评分的相关性############
ssGsea_Score_pan_cancer_HALLMARK_GSVA(df)


#泛癌森林图######
Targets <- c("PFI", "DFI", "DSS", "OS")
for (Target in Targets) {
  ssGsea_Score_pan_cancer_fores(df,Target)
}


#基于基因集的癌症类型特异性表达分析##################
#按临床信息分组的基于基因组的样本配对肿瘤-正常箱线图（一张图）######
ssGsea_Score_paired_boxplot(df)



#免疫应答######
file_path1 <- paste0("D:/Dim/TCGAplot/免疫治疗队列/", dir_name, "-GSE91061_OS_ssGSEA_scores.Rdata")
file_path2 <- paste0("D:/Dim/TCGAplot/免疫治疗队列/", dir_name, "-GSE135222_PFS_ssGSEA_scores.Rdata")
file_path3 <- paste0("D:/Dim/TCGAplot/免疫治疗队列/", dir_name, "-GSE78220_OS_ssGSEA_scores.Rdata")
file_path4 <- paste0("D:/Dim/TCGAplot/免疫治疗队列/", dir_name, "-PMID32472114_OS_ssGSEA_scores.Rdata")
vector <- list(file_path1, file_path2, file_path3, file_path4)
for (file_path in vector) {
  plot_km(file_path) 
}

#免疫应答1######
vector <- list(file_path1)
vector <- list(file_path2)
vector <- list(file_path3)
vector <- list(file_path4)
for (file_path in vector) {
  plot_km1(file_path) 
}



#按临床信息分组的基于基因组的表达分析######
cancer="LUAD"
ssGsea_Score_Tumor_Normal(df,cancer)


#按临床信息分组的基于基因组的样本配对肿瘤-正常箱线图（单肿瘤）######
cancer <- "LUAD"
ssGsea_Score_one_paired_boxplot(df,cancer)


#基因集泛癌ROC曲线#########
cancer <- "LUAD"
ssGsea_Score_ROC(df,cancer)


#基因集泛癌生存KM###############
cancer <- "COAD"
ssGsea_Score_KM(df,cancer)
