
library(reshape2)
library(dplyr)
library(ggplot2)
library(dplyr)
rm(list = ls())
load("D:/Dim/TCGAplot/TCGA_GTEx_pancancer_mrna_pheno.rdata")

gene_families <- list(
  
  COX = c(
    "COXFA4","COXFA4L2","COXFA4L3","COX4I1","COX4I2",
    "COX5A","COX5B","COX6A1","COX6A2","COX6B1","COX6B2","COX6C",
    "COX7A1","COX7A2","COX7B","COX7B2","COX7C","COX8A","COX8C",
    "MT-CO1","MT-CO2","MT-CO3"
  ),
  GST = c(
    "GSTA1","GSTA2","GSTA3","GSTA4","GSTA5","GSTA6P","GSTA7P","GSTA8P","GSTA9P",
    "GSTA10P","GSTA11P","GSTA12P",
    "GSTK1",
    "GSTM1","GSTM2","GSTM3","GSTM4","GSTM5",
    "GSTO1","GSTO2","GSTO3P",
    "GSTP1",
    "HPGDS",
    "GSTT1","GSTT2","GSTT2B","GSTT3P","GSTT4",
    "GSTZ1"
  ),
  ATP5 = c(
    "ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D","ATP5F1E",
    "ATP5MC1","ATP5MC2","ATP5MC3","ATP5ME","ATP5MF","ATP5MG","ATP5MJ","ATP5MK",
    "MT-ATP6","MT-ATP8",
    "ATP5PB","ATP5PD","ATP5PF","ATP5PO","ATP5IF1"
  ),
  MAPK = c(
    "MAPK1","MAPK3","MAPK4","MAPK6","MAPK7",
    "MAPK8","MAPK9","MAPK10","MAPK11","MAPK12",
    "MAPK13","MAPK14","MAPK15","NLK"
  ),
  NDU = c("MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6",
          "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS7", "NDUFS8", "NDUFV1", "NDUFV2"
  ),
  UQC = c(
    "UQCRB","UQCRQ","UQCRC1","UQCRC2","MT-CYB",
    "CYC1","UQCRFS1","UQCRH","UQCR10","UQCR11"
  ),
  SLC26 = c(
    "SLC26A1","SLC26A2","SLC26A3","SLC26A4","SLC26A5","SLC26A6",
    "SLC26A7","SLC26A8","SLC26A9","SLC26A10P","SLC26A11"
  ),
  ASM = c(    "ASMT", "ASMTL", "AS3MT", "CARNMT1", "COMT", "COMTD1", "COQ3", 
              "COQ5", "GAMT", "GNMT", "HNMT", "INMT", "TMT1A", "TMT1B", 
              "NNMT", "PNMT", "TOMT", "TPMT"),
  NFE = c(
    "ATF1","ATF2","ATF3","ATF4","ATF5","ATF6","ATF6B","ATF7",
    "BACH1","BACH2",
    "BATF","BATF2","BATF3",
    "CEBPA","CEBPB","CEBPD","CEBPE","CEBPG",
    "CREBL2","CREBRF","CREBZF","CREB1","CREB5","CREM",
    "DDIT3",
    "FOS","FOSB","FOSL1","FOSL2",
    "JDP2","JUN","JUNB","JUND",
    "NFE2","NFE2L1","NFE2L2","NFE2L3",
    "NFILZ","NFIL3",
    "XBP1"
  ),
  SLC25 = c(
    "SLC25A1","SLC25A2","SLC25A3","SLC25A4","SLC25A5","SLC25A6",
    "UCP1","UCP2","UCP3",
    "SLC25A10","SLC25A11","SLC25A12","SLC25A13","SLC25A14","SLC25A15",
    "SLC25A16","SLC25A17","SLC25A18","SLC25A19","SLC25A20","SLC25A21",
    "SLC25A22","SLC25A23","SLC25A24","SLC25A25","SLC25A26","SLC25A27",
    "SLC25A28","SLC25A29","SLC25A30","SLC25A31","SLC25A32","SLC25A33",
    "SLC25A34","SLC25A35","SLC25A36","SLC25A37","SLC25A38","SLC25A39",
    "SLC25A40","SLC25A41","SLC25A42","SLC25A43","SLC25A44","SLC25A45",
    "SLC25A46","SLC25A47","SLC25A48",
    "MTCH1","MTCH2","SLC25A51","SLC25A52","SLC25A53"
  ),
  
  NDUF = c(
    "NDUFAB1","NDUFA1","NDUFA2","NDUFA3","NDUFA5","NDUFA6","NDUFA7",
    "NDUFA8","NDUFA9","NDUFA10","NDUFA11","NDUFA12","NDUFA13",
    "NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7",
    "NDUFB8","NDUFB9","NDUFB10","NDUFB11",
    "NDUFC1","NDUFC2",
    "NDUFS4","NDUFS5","NDUFS6","NDUFV3"
  ),
  CYP = c(
    "CYP2AB1P","CYP2AC1P","CYP2A6","CYP2A7","CYP2A7P1","CYP2A13",
    "CYP2B6","CYP2B7P",
    "CYP2C8","CYP2C9","CYP2C18","CYP2C19","CYP2C23P","CYP2C56P",
    "CYP2C58P","CYP2C59P","CYP2C60P","CYP2C61P","CYP2C63P","CYP2C64P",
    "CYP2C115P",
    "CYP2D6","CYP2D7","CYP2D8P",
    "CYP2E1",
    "CYP2F1","CYP2F2P",
    "CYP2G1P","CYP2G2P",
    "CYP2J2",
    "CYP2R1",
    "CYP2S1",
    "CYP2T1P","CYP2T3P",
    "CYP2U1",
    "CYP2W1"
  ),
  AKR = c(
    "AKR1A1","AKR1B1","AKR1B10","AKR1B15",
    "AKR1C1","AKR1C2","AKR1C3","AKR1C4","AKR1C5P","AKR1C6P","AKR1C7P","AKR1C8",
    "AKR1D1","AKR1E2",
    "KCNAB1","KCNAB2","KCNAB3",
    "AKR7A2","AKR7A3"
  ),
  EPH = c(
    "ENOPH1","EPHX2","LHPP","NANP","PHOSPHO1","PHOSPHO2",
    "PMM1","PMM2","PNKP","PSPH","PUDP"
  ),
  ARN = c(
    "AHR","AHRR","ARNT","BMAL1","BMAL2","ARNT2","CLOCK","EPAS1","HIF1A","HIF3A",
    "KCNH1","KCNH2","KCNH3","KCNH4","KCNH5","KCNH6","KCNH7","KCNH8",
    "NCOA1","NCOA2","NCOA3",
    "NPAS1","NPAS2","NPAS3","NPAS4",
    "PASD1","PASK",
    "PDE8A","PDE8B",
    "PER1","PER2","PER3",
    "SIM1","SIM2"
  )
)

rt=tcga_gtex_mrna_pheno
str(rt)
class(rt)
calculate_family_logFC_percent <- function(rt, gene_families) {
  # 改列名统一
  colnames(rt)[2] <- "Cancer"
  colnames(rt)[4] <- "Group"
  
  # 提取基因列
  expr <- rt[, 5:ncol(rt)]
  gene_names <- colnames(expr)
  
  # 结果矩阵
  cancers <- unique(rt$Cancer)
  logFC_matrix <- matrix(NA, nrow = length(gene_names), ncol = length(cancers))
  pval_matrix <- matrix(NA, nrow = length(gene_names), ncol = length(cancers))
  rownames(logFC_matrix) <- gene_names
  rownames(pval_matrix) <- gene_names
  colnames(logFC_matrix) <- cancers
  colnames(pval_matrix) <- cancers
  
  # 计算logFC和p值
  for (gene in gene_names) {
    for (cancer in cancers) {
      tmp <- rt[rt$Cancer == cancer, c("Group", gene)]
      if (length(unique(tmp$Group)) < 2) next
      tumor <- tmp[tmp$Group == "TCGA_tumor", gene]
      normal <- tmp[tmp$Group != "TCGA_tumor", gene]
      if (length(tumor) < 3 || length(normal) < 3) next
      
      logFC_matrix[gene, cancer] <- mean(tumor, na.rm = TRUE) - mean(normal, na.rm = TRUE)
      pval_matrix[gene, cancer] <- wilcox.test(tumor, normal)$p.value
    }
  }
  
  # 计算家族结果
  family_logFC_percent <- data.frame(
    Family = names(gene_families),
    Total = NA,
    SigPosCount = NA,
    Proportion = NA,
    GeneCount = NA
  )
  
  for (i in seq_along(gene_families)) {
    fam_genes <- intersect(rownames(logFC_matrix), gene_families[[i]])
    Total <- sum(!is.na(logFC_matrix[fam_genes, ]))
    SigPosCount <- sum(logFC_matrix[fam_genes, ] > 0 & pval_matrix[fam_genes, ] < 0.05, na.rm = TRUE)
    Proportion <- SigPosCount / Total
    GeneCount <- length(fam_genes)
    
    family_logFC_percent[i, c("Total", "SigPosCount", "Proportion", "GeneCount")] <- 
      c(Total, SigPosCount, Proportion, GeneCount)
  }
  
  return(family_logFC_percent)
}

# 使用示例
family_logFC_percent <- calculate_family_logFC_percent(rt, gene_families)
print(family_logFC_percent)
family_logFC_percent1=family_logFC_percent
family_logFC_percent=family_logFC_percent1

library(dplyr)

# 假设 family_logFC_percent 已经计算好
family_logFC_percent <- family_logFC_percent %>%
  mutate(
    AdjProp = (SigPosCount + 1) / (Total + 2),
    Quality = AdjProp / sqrt(GeneCount)
  ) %>%
  arrange(desc(Quality))

# 查看排序后的结果
print(family_logFC_percent)

# 如果要选择前三个优质家族
top_families <- head(family_logFC_percent, 3)
print(top_families)





library(dplyr)

family_logFC_percent <- family_logFC_percent %>%
  mutate(
    Score = (Proportion * Total) / sqrt(GeneCount)
  ) %>%
  arrange(desc(Score))

# 查看排序结果
print(family_logFC_percent)

# 选择前三个综合得分最高的家族
top_families <- head(family_logFC_percent, 3)
print(top_families)

