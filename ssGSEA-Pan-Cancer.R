# 加载必要的包
library(GSVA)
library(GSEABase)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(fmsb)
library(tibble)
library(stringr)
library(psych)
library(pheatmap)
library(survival)
library(forestplot)
library(grid)
library(tidyr)
library(rstatix)
library(ggsci)
library(reshape2)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(limma)
library(survminer)
library(beepr)
library(patchwork)
library(readxl)
#rm(list = ls())
# load(list.files(pattern = "\\.Rdata$"))
lapply(list.files(pattern = "\\.Rdata$"), load, envir = .GlobalEnv)
dir_name <- basename(getwd())

load("D:/Dim/基因组泛癌/数据/表达矩阵/all_tpm.Rdata")
load("D:/Dim/TCGAplot/TCGA_GTEx_pancancer_mrna_pheno.rdata")

# 计算基因组ssGsea得分########
ssGsea_Score <- function(tpm, cox_genes) {
  # 转为数值型表达矩阵
  tpm_numeric <- tpm %>%
    mutate(across(-c(Cancer, Group), as.numeric))
  expr_data <- as.matrix(tpm_numeric[, -c(1:2)])  # 去除前三列
  colnames(expr_data) <- colnames(tpm_numeric)[-c(1:2)]
  expr_matrix <- t(expr_data)                    # 转置为基因为行，样本为列
  
  
  # 自动获取前缀用于命名
  prefix <- toupper(substr(cox_genes[1], 1, 3))
  # 创建 GeneSet 对象
  cox_gene_set <- GeneSet(cox_genes, setName = paste0(prefix, "_family"))
  
  set_name <- paste0(prefix, "_family")
  gene_set <- list()
  gene_set[[set_name]] <- cox_genes
  geneSet=gene_set
  # 运行 ssGSEA
  param <- ssgseaParam(
    exprData = expr_matrix,
    geneSet,
    minSize = 5,
    maxSize = Inf,
    alpha = 0.25,
    normalize = TRUE
  )
  
  ssgsea_scores <- gsva(param)
  
  if (nrow(ssgsea_scores) > 0) {
    df <- data.frame(
      Score = as.numeric(ssgsea_scores[1, ]),
      Cancer = tpm_numeric$Cancer,
      Group = tpm_numeric$Group,
      row.names = colnames(ssgsea_scores)
    )
    save(df, file = paste0(prefix, "_ssGSEA_scores.Rdata"))
  }
}



#基因家族差异表达LogFC########
plot_diff_heatmap <- function(rt) {
  ylab_text <- paste0(dir_name, "_gene family")
  rt <- rt[, c(3:ncol(rt), 1, 2)]
  colnames(rt)[ncol(rt)-1] <- "CancerType"
  colnames(rt)[ncol(rt)] <- "Type"
  
  # 去掉正常样本数<5的肿瘤
  normal_counts <- table(rt$CancerType[rt$Type == "Normal"])
  valid_cancers <- names(normal_counts[normal_counts >= 3])
  rt <- rt[rt$CancerType %in% valid_cancers, ]
  
  data <- rt
  data[, 1:(ncol(data)-2)] <- lapply(data[, 1:(ncol(data)-2)], as.numeric)
  genelist <- colnames(data)[1:(ncol(data)-2)]
  
  logFC_matrix <- matrix(NA, nrow = length(genelist), ncol = length(unique(data$CancerType)))
  rownames(logFC_matrix) <- genelist
  colnames(logFC_matrix) <- sort(unique(data$CancerType))
  pval_matrix <- logFC_matrix
  
  for (gene in genelist) {
    for (cancer in unique(data$CancerType)) {
      tmp <- data[data$CancerType == cancer, c(gene, "Type")]
      if (length(unique(tmp$Type)) < 2) next
      group1 <- tmp[tmp$Type == "Tumor", gene]
      group2 <- tmp[tmp$Type == "Normal", gene]
      if (length(group1) < 3 || length(group2) < 3) next
      logFC <- mean(group1) - mean(group2)
      pval <- wilcox.test(group1, group2)$p.value
      logFC_matrix[gene, cancer] <- logFC
      pval_matrix[gene, cancer] <- pval
    }
  }
  
  df <- melt(logFC_matrix, varnames = c("Gene", "Cancer"), value.name = "logFC")
  df$pval <- melt(pval_matrix)$value
  df$logFC_text <- sprintf("%.2f", df$logFC)
  
  df$fill_bin <- cut(df$logFC,
                     breaks = c(-Inf, -2, -1, 0, 1, 2, Inf),
                     labels = c("<-2", "-2~-1", "-1~0", "0~1", "1~2", ">2"))
  df$fill_bin <- as.character(df$fill_bin)
  df$fill_bin[df$pval >= 0.05 | is.na(df$logFC)] <- "P=>0.05"
  df$fill_bin <- factor(df$fill_bin, levels = c(">2", "1~2", "0~1", "-1~0", "-2~-1", "<-2", "P=>0.05"))
  color_map <- c(
    "<-2" = "#2f5d91",
    "-2~-1" = "#669dbb",
    "-1~0" = "#ded3cf",
    "0~1" = "#fdb1a4",
    "1~2" = "#f45859",
    ">2" = "#99253a",
    "P=>0.05" = "white"
  )
  p <- ggplot(df, aes(x = Cancer, y = Gene)) +
    geom_tile(aes(fill = fill_bin), color = "grey90") +
    scale_fill_manual(values = color_map, name = "logFC") +
    geom_text(aes(label = logFC_text), size = 3.2) +
    coord_fixed(ratio = 0.8) +
    theme_minimal(base_size = 12) +
    labs(x = "Pan-Cancer", y = ylab_text)
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank()
    )
  file_name <- paste0(dir_name,"-LogFC-diff.pdf")
  ggsave(file_name, p, width = 13, height = 9)
}
plot_diff_heatmap1 <- function(dir_name) {
  # 整理数据
  tcga_gtex_mrna_pheno <- tcga_gtex_mrna_pheno[, -3]
  colnames(tcga_gtex_mrna_pheno)[2:3] <- c("Cancer", "Group")
  tcga_gtex_mrna_pheno$Group <- gsub("TCGA_tumor", "Tumor", tcga_gtex_mrna_pheno$Group)
  tcga_gtex_mrna_pheno$Group <- gsub("TCGA_normal|GTEx_normal", "Normal", tcga_gtex_mrna_pheno$Group)
  tcga_gtex_mrna_pheno <- tcga_gtex_mrna_pheno[, -1]
  tcga_gtex_mrna_pheno <- tcga_gtex_mrna_pheno[, c(1:2, which(colnames(tcga_gtex_mrna_pheno) %in% cox_genes))]
  tcga_gtex_mrna_pheno <- tcga_gtex_mrna_pheno[, c(3:ncol(tcga_gtex_mrna_pheno), 1, 2)]
  colnames(tcga_gtex_mrna_pheno)[ncol(tcga_gtex_mrna_pheno)-1] <- "CancerType"
  colnames(tcga_gtex_mrna_pheno)[ncol(tcga_gtex_mrna_pheno)] <- "Type"
  
  normal_counts <- table(tcga_gtex_mrna_pheno$CancerType[tcga_gtex_mrna_pheno$Type == "Normal"])
  valid_cancers <- names(normal_counts[normal_counts >= 3])
  tcga_gtex_mrna_pheno <- tcga_gtex_mrna_pheno[tcga_gtex_mrna_pheno$CancerType %in% valid_cancers, ]
  
  tcga_gtex_mrna_pheno[, 1:(ncol(tcga_gtex_mrna_pheno) - 2)] <- lapply(tcga_gtex_mrna_pheno[, 1:(ncol(tcga_gtex_mrna_pheno) - 2)], as.numeric)
  genelist <- intersect(colnames(tcga_gtex_mrna_pheno), cox_genes)
  
  logFC_matrix <- matrix(NA, nrow = length(genelist), ncol = length(unique(tcga_gtex_mrna_pheno$CancerType)))
  rownames(logFC_matrix) <- genelist
  colnames(logFC_matrix) <- sort(unique(tcga_gtex_mrna_pheno$CancerType))
  pval_matrix <- logFC_matrix
  
  for (gene in genelist) {
    for (cancer in unique(tcga_gtex_mrna_pheno$CancerType)) {
      tmp <- tcga_gtex_mrna_pheno[tcga_gtex_mrna_pheno$CancerType == cancer, c(gene, "Type")]
      if (length(unique(tmp$Type)) < 2) next
      group1 <- tmp[tmp$Type == "Tumor", gene]
      group2 <- tmp[tmp$Type == "Normal", gene]
      if (length(group1) < 3 || length(group2) < 3) next
      logFC_matrix[gene, cancer] <- mean(group1) - mean(group2)
      pval_matrix[gene, cancer] <- wilcox.test(group1, group2)$p.value
    }
  }
  
  df <- reshape2::melt(logFC_matrix, varnames = c("Gene", "Cancer"), value.name = "logFC")
  df$pval <- reshape2::melt(pval_matrix)$value
  df$logFC_text <- sprintf("%.2f", df$logFC)
  df$fill_bin <- cut(df$logFC, breaks = c(-Inf, -2, -1, 0, 1, 2, Inf),
                     labels = c("<-2", "-2~-1", "-1~0", "0~1", "1~2", ">2"))
  df$fill_bin <- as.character(df$fill_bin)
  df$fill_bin[df$pval >= 0.05 | is.na(df$logFC)] <- "P=>0.05"
  df$fill_bin <- factor(df$fill_bin, levels = c(">2", "1~2", "0~1", "-1~0", "-2~-1", "<-2", "P=>0.05"))
  
  color_map <- c("<-2" = "#2f5d91", "-2~-1" = "#669dbb", "-1~0" = "#ded3cf",
                 "0~1" = "#fdb1a4", "1~2" = "#f45859", ">2" = "#99253a", "P=>0.05" = "white")
  
  p <- ggplot(df, aes(x = Cancer, y = Gene)) +
    geom_tile(aes(fill = fill_bin), color = "grey90") +
    scale_fill_manual(values = color_map, name = "logFC") +
    geom_text(aes(label = logFC_text), size = 3.2) +
    coord_fixed(ratio = 0.8) +
    theme_minimal(base_size = 12) +
    labs(x = "Pan-Cancer", y = paste0(dir_name, "_gene family")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 10),
          panel.grid = element_blank())
  
  ggsave(paste0(dir_name, "-LogFC-diff-GTEx&TCGA.pdf"), p, width = 13, height = 9)
}




#生存分析(KM)######
plot_OS_KM_expression <- function(rt) {
  dir_name <- basename(getwd())
  dir_name <- gsub("基因", "", dir_name)
  # 加载临床数据
  load("D:/Dim/基因组泛癌/数据/临床数据/ALL_OS.Rdata")
  OS <- OS[, c(3, 4)]
  colnames(OS) = c("futime", "fustat")
  rt <- rt[, c(3:ncol(rt), 1, 2)]
  colnames(rt)[ncol(rt) - 1] <- "CancerType"
  colnames(rt)[ncol(rt)] <- "Type"
  exp = rt
  exp = exp[(exp[,"Type"] == "Tumor"),]    # 只保留肿瘤数据
  # 数据合并
  sameSample = intersect(row.names(OS), row.names(exp))  # 获取临床数据和基因表达数据中的交集样本
  OS = OS[sameSample,]
  exp = exp[sameSample,]
  surData = cbind(OS, exp)
  surData = cbind(id = row.names(surData), surData)  # 将行名插入到第1列，并命名为 id
  file_name <- paste0(dir_name, "_expTime.txt")
  write.table(surData, file = file_name, sep = "\t", quote = F, row.names = F)
  inputFile <- paste0(dir_name, "_expTime.txt")      # 输入文件路径
  pFilter = 0.05  # KM 方法的 p-value 阈值
  expTime = read.table(inputFile, header = T, sep = "\t", check.names = F, row.names = 1)  # 读取输入文件
  expTime$futime <- as.numeric(expTime$futime)  # 将生存时间转换为数值型
  # 检查无法成功转换的值（NA）
  problematic_values <- expTime$futime[is.na(as.numeric(expTime$futime))]
  print(problematic_values)
  # 清理缺失值
  expTime_clean <- expTime[complete.cases(expTime), ]
  expTime = expTime_clean
  expTime$futime = expTime$futime / 365  # 将生存时间单位转为年
  # 差异分析
  outTab = data.frame()
  for (gene in colnames(expTime)[3:(ncol(expTime) - 2)]) {
    for (CancerType in levels(factor(expTime[, "CancerType"]))) {
      rt1 = expTime[(expTime[, "CancerType"] == CancerType),]
      rt1$group <- ifelse(rt1[, gene] > median(rt1[, gene]), "High", "Low")
      
      # 检查每个组的样本数
      group_count = table(rt1$group)
      if (length(group_count) < 2 || any(group_count < 2)) {
        next  # 如果某个组样本数小于2，跳过该基因
      }
      # 生存分析
      diff = survdiff(Surv(futime, fustat) ~ group, data = rt1)
      pValue = 1 - pchisq(diff$chisq, df = 1)
      
      if (pValue < pFilter) {  # 筛选 p 值小于 0.05 的结果
        outVector = cbind(gene, CancerType, pValue)
        outTab = rbind(outTab, outVector)
        
        if (pValue < 0.001) {
          pValue = "p<0.001"
        } else {
          pValue = paste0("p=", sprintf("%.03f", pValue))
        }
        fit <- survfit(Surv(futime, fustat) ~ group, data = rt1)
        # 绘制生存曲线
        surPlot = ggsurvplot(fit,
                             data = rt1,
                             title = paste0("Cancer: ", CancerType, "\n\n\n", gene),
                             legend = c(0.9, 0.83),  # 图例位置
                             legend.title = " ",  # 图例标题
                             legend.labs = c("High", "Low"),  # 图例标签
                             xlab = "Time(years)",
                             ylab = "Overall survival",
                             palette = c("#FFCC33", "#0066CC"),  # 设置调色板颜色
                             ggtheme = theme_bw(base_size = 15) +  # 设置基本主题为黑白风格，文字大小为 15
                               theme(plot.title = element_text(hjust = 0.80, vjust = -23, size = 13),  # 设置标题居中
                                     panel.grid = element_blank(),  # 隐藏网格线
                                     panel.border = element_rect(colour = "black", linewidth = 1.8, fill = NA)  # 边框加粗
                               )
        )
        rt1$group <- factor(rt1$group, levels = c("Low", "High"))
        res_cox <- coxph(Surv(futime, fustat) ~ group, data = rt1)
        p_val <- summary(res_cox)$coef[5]
        p_label <- ifelse(p_val < 0.001, "P < 0.001", paste("P:", round(p_val, 4)))
        hr_text <- grobTree(textGrob(paste("HR:", round(summary(res_cox)$conf.int[1], 2)), x = 0.145, y = 0.15, gp = gpar(fontsize = 10)))
        ci_text <- grobTree(textGrob(paste("95% CI:", round(summary(res_cox)$conf.int[3], 2), "-", round(summary(res_cox)$conf.int[4], 2)),x = 0.145, y = 0.11, gp = gpar(fontsize = 10)))
        p_text <- grobTree(textGrob(p_label, x = 0.145, y = 0.07, gp = gpar(fontsize = 10)))
        surPlot <- surPlot$plot + annotation_custom(hr_text) +annotation_custom(ci_text) +annotation_custom(p_text)
        # 保存图像
        safeGene <- gsub("<|>", "", gene)
        safeFileName <- paste0("生存分析(KM)/", safeGene, "_", CancerType, "_", gsub("[<>]", "", pValue), ".pdf")
        pdf(file = safeFileName, onefile = FALSE, width = 6, height = 5)
        print(surPlot)
        dev.off()
      }  # 筛选 p 值小于 0.05 的
    }
  }
}



#cox分析##########
plot_cox_analysis_forest <- function(rt_cox)
{
  rt_1=rt_cox
  class(rt_1$futime)
  #列的数据类型转换为数值型。如果原始数据是字符型（例如日期或其他格式），会强制转为数值格式。非数值的数据会被转换为 NA（缺失值）
  rt_1$futime <- as.numeric(rt_1$futime)
  #查找futime列中无法成功转换为数值的值（即NA）
  problematic_values <- rt_1$futime[is.na(as.numeric(rt_1$futime))]
  print(problematic_values)
  #用于清理数据框 rt_1 中的缺失值。它会移除那些包含 NA（缺失值）的行
  rt_clean <- rt_1[complete.cases(rt_1), ]
  dim(rt_clean)
  head(rt_clean)
  rt_1=rt_clean
  # 循环基因和肿瘤类型
  i = 0
  for(gene in colnames(rt_1)[3:(ncol(rt_1)-2)]) 
  {  # 针对每个基因进行操作
    i = i + 1
    outTab = data.frame()
    # 对肿瘤类型进行循环
    for(CancerType in levels(factor(rt_1[,"CancerType"])))
    {
      rt1 = rt_1[(rt_1[,"CancerType"] == CancerType), ]
      #针对当前基因进行高低表达分组
      rt1$gene_group <- ifelse(rt1[, gene] < median(rt1[, gene]), "Low", "High")
      rt1$gene_group <- factor(rt1$gene_group, levels = c("Low", "High"))
      #Cox回归分析
      cox = coxph(Surv(futime, fustat) ~ gene_group, data = rt1)
      coxSummary = summary(cox)
      #提取 HR, CI, 和 P 值
      HR_CI = paste0(round(coxSummary$conf.int[,"exp(coef)"], 3),  " (", round(coxSummary$conf.int[,"lower .95"], 3), " - ", round(coxSummary$conf.int[,"upper .95"], 3), ")")
      coxP = coxSummary$coefficients[,"Pr(>|z|)"]
      #总样本数
      TotalN = nrow(rt1)
      #结果存储
      outTab = rbind(outTab,cbind(Cancer = CancerType,`Total(N)` = TotalN,`HR (95% CI)` = HR_CI,`P value` = ifelse(coxP < 0.001, "< 0.001", round(coxP, 3))))
    }
    # 输出结果文件
    outFile = file.path("生存分析(COX回归)", paste0(gene, ".cox.txt"))
    write.table(outTab, file = outFile, sep = "\t", row.names = FALSE, quote = FALSE)
  }
}

plot_cox_analysis_forest_picture <- function(folder_path)
{
  # 获取文件夹中所有.txt文件的列表
  txt_files <- list.files(path = folder_path, pattern = "\\.cox\\.txt$", full.names = TRUE)
  # 循环遍历每个txt文件
  for(file_path in txt_files) {
    gene_name <- sub("\\.cox\\.txt$", "", basename(file_path))
    cox_data <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    cox_data$`P value` <- as.numeric(gsub("< 0.001", "0.0001", cox_data$`P value`))
    forest_data <- cox_data %>%
      mutate(
        HR_val = as.numeric(sub(" \\(.*", "", `HR (95% CI)`)),
        HR_lower = as.numeric(sub(".*\\(([^-]+) -.*", "\\1", `HR (95% CI)`)),
        HR_upper = as.numeric(sub(".*- ([^\\)]+)\\)", "\\1", `HR (95% CI)`)),
        P_value = as.numeric(`P value`)
      ) %>%
      select(Cancer, Total = `Total(N)`, HR_val, HR_lower, HR_upper, P_value)
    # 计算红色标记的个数
    red_count <- sum(forest_data$P_value < 0.05 & forest_data$HR_val > 1)
    # 将数据转换为适合forestplot的格式
    tabletext <- cbind(
      c("Cancer", forest_data$Cancer),
      c("Total(N)", forest_data$Total),
      c("HR (95% CI)", paste0(round(forest_data$HR_val, 3), " (", round(forest_data$HR_lower, 3), " - ", round(forest_data$HR_upper, 3), ")")),
      c("P value", ifelse(forest_data$P_value < 0.001, "<0.001", round(forest_data$P_value, 3)))
    )
    # 根据条件为文本着色Cancer
    color_vector <- ifelse(forest_data$P_value < 0.05 & forest_data$HR_val > 1, "#E64B35", 
                           ifelse(forest_data$P_value < 0.05 & forest_data$HR_val < 1, "#4DBBD5", "black"))
    color_vector <- c("black", color_vector)  # 在第一个元素前加入黑色
    # 创建森林图
    p <- forestplot(
      labeltext = tabletext,
      title = paste0(gene_name, " – Overall Survival"), 
      title_gp = grid::gpar(fontsize = 14, fontface = "bold", just = "left"),
      mean = c(NA, forest_data$HR_val),
      lower = c(NA, forest_data$HR_lower),
      upper = c(NA, forest_data$HR_upper),
      hrzl_lines = list(`1` = grid::gpar(lwd = 2, col = "black"), `2` = grid::gpar(lwd = 2, col = "black"), `35` = grid::gpar(lwd = 2, col = "black")), 
      is.summary = rep(FALSE, nrow(forest_data) + 1),
      zero = 1,
      boxsize = 0.4,
      lineheight = unit(0.2, "cm"),
      xlog = FALSE,
      col = fpColors(box = "#4DBBD5", lines = "black", zero = "gray50", text = color_vector), # 设置文本颜色
      lwd.ci = 2.5,
      ci.vertices = TRUE,
      ci.vertices.height = 0.02,
      clip = c(0, 6),
      xticks = seq(0, 6, by = 2),
      graph.pos = 4,
      graphwidth = unit(5, "cm"),
      txt_gp = forestplot::fpTxtGp(ticks = grid::gpar(cex = 0.9)), 
    )
    # 将图保存为PDF，文件名中包含红色标记的个数
    pdf_file <- file.path(folder_path, paste0(gene_name, "-红色", red_count, ".pdf"))
    pdf(file = pdf_file, width = 7, height = 9)
    # 创建新页面
    grid.newpage()
    # 绘制森林图
    print(p, newpage = FALSE)  # 阻止forestplot创建新页面
    # 关闭PDF
    dev.off()
  }
}

#COX生存分析热图########
plot_cox_heatmap1 <- function(folder_path) {
  txt_files <- list.files(path = folder_path, pattern = "\\.cox\\.txt$", full.names = TRUE)
  all_data <- list()
  for (file_path in txt_files) {
    gene_name <- sub("\\.cox\\.txt$", "", basename(file_path))
    df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    df$`P value` <- as.numeric(gsub("< 0.001", "0.0001", df$`P value`))
    df <- df %>%
      dplyr::mutate(
        Gene = gene_name,
        Cancer = as.character(Cancer),
        HR = as.numeric(sub(" \\(.*", "", `HR (95% CI)`)),
        P = as.numeric(`P value`)
      ) %>%
      dplyr::select(Gene, Cancer, HR, P)
    all_data[[gene_name]] <- df
  }
  
  merged_data <- dplyr::bind_rows(all_data)
  
  heatmap_matrix <- merged_data %>%
    dplyr::mutate(status = dplyr::case_when(
      P < 0.05 & HR > 1 ~ 1 ,    # 红色
      P < 0.05 & HR < 1 ~ -1 ,   # 蓝色
      TRUE ~ 0                 # 灰色
    )) %>%
    dplyr::select(Gene, Cancer, status) %>%
    tidyr::pivot_wider(names_from = Cancer, values_from = status, values_fill = 0)
  
  mat <- as.matrix(heatmap_matrix[,-1])
  rownames(mat) <- heatmap_matrix$Gene
  
  # 固定癌种顺序（只保留出现在数据中的）
  fixed_order <- c(
    "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
    "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
    "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
    "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"
  )
  colnames(mat) <- factor(colnames(mat), levels = fixed_order)
  mat <- mat[, order(colnames(mat))]
  
  col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#4DBBD5", "white", "#E64B35"))
  
  p=ComplexHeatmap::Heatmap(
    mat,
    col = col_fun,
    cluster_rows = F,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    heatmap_legend_param = list(
      title = "",
      at = c(1, 0, -1),
      labels = c("Risky", "p>0.05", "Protective")
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.rect(x, y, width, height, gp = gpar(col = "black", fill = NA, lwd = 0.5))
    }
  )
  
  file_name <- paste0(dir_name,"-Pan_cox_heatmap(无统计).pdf")
  pdf(file_name, width=8, height=4)
  print(p)
  dev.off()
}
plot_cox_heatmap2 <- function(folder_path) {
  txt_files <- list.files(path = folder_path, pattern = "\\.cox\\.txt$", full.names = TRUE)
  all_data <- list()
  for (file_path in txt_files) {
    gene_name <- sub("\\.cox\\.txt$", "", basename(file_path))
    df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    df$`P value` <- as.numeric(gsub("< 0.001", "0.0001", df$`P value`))
    df <- df %>%
      dplyr::mutate(
        Gene = gene_name,
        Cancer = as.character(Cancer),
        HR = as.numeric(sub(" \\(.*", "", `HR (95% CI)`)),
        P = as.numeric(`P value`)
      ) %>%
      dplyr::select(Gene, Cancer, HR, P)
    all_data[[gene_name]] <- df
  }
  
  merged_data <- dplyr::bind_rows(all_data)
  
  heatmap_matrix <- merged_data %>%
    dplyr::mutate(status = dplyr::case_when(
      P < 0.05 & HR > 1 ~ 1,
      P < 0.05 & HR < 1 ~ -1,
      TRUE ~ 0
    )) %>%
    dplyr::select(Gene, Cancer, status) %>%
    tidyr::pivot_wider(names_from = Cancer, values_from = status, values_fill = 0)
  
  mat <- as.matrix(heatmap_matrix[,-1])
  rownames(mat) <- heatmap_matrix$Gene
  
  fixed_order <- c(
    "ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
    "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
    "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
    "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"
  )
  mat <- mat[, colnames(mat) %in% fixed_order]
  mat <- mat[, order(match(colnames(mat), fixed_order))]
  
  risky_score <- merged_data %>%
    dplyr::filter(P < 0.05) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      HR_gt1 = sum(HR > 1),
      HR_lt1 = sum(HR < 1),
      RiskyScore = HR_gt1 - HR_lt1
    ) %>%
    dplyr::arrange(desc(RiskyScore))
  
  mat <- mat[risky_score$Gene, , drop = FALSE]
  
  ha_left <- ComplexHeatmap::rowAnnotation(
    RiskyScore = ComplexHeatmap::anno_barplot(
      risky_score$RiskyScore,
      border = TRUE,
      gp = grid::gpar(fill = ifelse(risky_score$RiskyScore > 0, "#E64B35", "#4DBBD5")),
      border_gp = grid::gpar(col = "grey80", lwd = 0.5),
      width = grid::unit(2, "cm"),
      bar_width = 1
    ),
    annotation_name_gp = grid::gpar(fontsize = 10)
  )
  
  #col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#4DBBD5", "white", "#E64B35"))
  col_fun <- circlize::colorRamp2(c(1, 0, -1), c("#E64B35", "white", "#4DBBD5"))
  p=ComplexHeatmap::Heatmap(
    mat,
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    right_annotation = ha_left,
    heatmap_legend_param = list(
      title = "",
      at = c(-1, 0, 1),
      #labels = c("Risky", "p>0.05", "Protective")
      labels = c("Protective", "p>0.05", "Risky")
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid::grid.rect(x, y, width, height, gp = grid::gpar(col = "black", fill = NA, lwd = 0.5))
    }
  )
  file_name <- paste0(dir_name,"-Pan_cox_heatmap(有统计).pdf")
  pdf(file_name, width=9, height=4)
  print(p)
  dev.off()
}


#GSVA得分正常-肿瘤对照箱线图######
ssGsea_Score_diff <- function(df) {
  ylab_text <- paste0(dir_name, "_channels score")
  p <- ggboxplot(df, x = "Cancer", y = "Score", fill = "Group",
                 xlab = "", ylab = ylab_text, color = "black", palette = "jco") +
    rotate_x_text(angle = 90) +
    grids(linetype = "dashed") +
    theme(
      legend.title = element_text(size = 14), 
      legend.text = element_text(size = 12),
      axis.title.x = element_text(size = 14), 
      axis.text.x = element_text(size = 12),
      axis.title.y = element_text(size = 14), 
      axis.text.y = element_text(size = 12),
      legend.position = "right"
    ) +
    border("black") +
    stat_compare_means(
      aes(group = Group),
      method = "wilcox.test",
      label = "p.signif",
      label.y.npc = "top",
      symnum.args = list(
        cutpoints = c(0, 0.001, 0.01, 0.05, 1),
        symbols = c("***", "**", "*", "ns")
      )
    )
  
  file_name <- paste0(dir_name,"-pan-diff.pdf")
  pdf(file_name, width=12, height=3.6)
  print(p)
  dev.off()
}
ssGsea_Score_diff1 <- function(df){
  ylab_text <- paste0(dir_name, "_channels score")
p <- ggplot(df, aes(x = Cancer, y = Score, fill = Group, color = Group)) +
  geom_boxplot(
    alpha = 0.5,
    width = 0.6,  # 每组半宽
    position = position_dodge2(width = 0.6, preserve = "single",padding = 0.2)  # 关键
  ) +
  geom_point(
    shape = 20, size = 0.1,
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
    alpha = 0.5
  ) +
  theme_classic() +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  scale_color_brewer(palette = "Set1", direction = -1) +
  theme(
    axis.text.x = element_text(colour = "black", size = 11, angle = 45, hjust = 0.5, vjust = 0.5),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right"
  ) +
  ylab(ylab_text) +
  stat_compare_means(
    aes(group = Group),
    method = "wilcox.test",
    label = "p.signif",
    label.y.npc = "top",
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c("***", "**", "*", "ns")
    )
  )
file_name <- paste0(dir_name,"-pan-diff1.pdf")
pdf(file_name, width=14, height=3)
print(p)
dev.off()
}



#GSVA得分纯肿瘤趋势箱线图######
ssGsea_Score_diff_Tumor <- function(df){
  ylab_text <- paste0(dir_name, "_channels score")
scores <- df %>% filter(Group == "Tumor") %>%
  mutate(Cancer = factor(Cancer, levels = names(sort(tapply(Score, Cancer, median)))))

p <- ggplot(scores, aes(x = Cancer, y = Score, fill = Cancer)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none") +
  labs(y = ylab_text, x = "Pan-Cancer")

file_name <- paste0(dir_name,"-pan-Tumor.pdf")
pdf(file_name, width=12, height=3.6)
print(p)
dev.off()
}


#GSVA泛癌TMB######  
ssGsea_Score_pan_cancer_TMB <- function(df){
load("D:/Dim/TCGAplot/All_TMB.Rdata")
df$Score <- df[, 1]
df$TMB <- TMB[rownames(df), "TMB"]
df <- df[df$Group == "Tumor", ]
df <- df[complete.cases(df), ]

# 自动获取癌种
cancers <- unique(df$Cancer)

# 计算相关性
res <- lapply(cancers, function(cancer) {
  sub <- df[df$Cancer == cancer, ]
  if (nrow(sub) < 3) return(NULL)
  ct <- suppressWarnings(cor.test(sub$Score, sub$TMB, method = "pearson"))
  data.frame(rho = ct$estimate, pvalue = ct$p.value)
})
names(res) <- cancers
res <- do.call(rbind, res[!sapply(res, is.null)])

# 添加标签
res$label <- ifelse(res$pvalue < 0.01, "**", ifelse(res$pvalue < 0.05, "*", " "))
res$group <- paste0(rownames(res), res$label)

# 转换格式
res <- res[, c("group", "rho")] %>%
  remove_rownames() %>%
  column_to_rownames("group") %>%
  add_column(Max = max(res$rho), .before = "rho") %>%
  add_column(Min = min(res$rho), .before = "rho") %>%
  t() %>% as.data.frame() %>% round(2)

# 调整顺序
res <- res[, c(1, ncol(res):2)]
file_name <- paste0(dir_name,"-pan-TMB.pdf")
pdf(file_name, width=7.5, height=7.5)
# 绘图
fmsb::radarchart(
  res,
  axistype = 1,
  pcol = "#E64B35",
  pfcol = scales::alpha("#E64B35", 0.1),
  plwd = 2,
  plty = 1,
  cglcol = "grey",
  cglty = 1,
  cglwd = 0.8,
  axislabcol = "black",
  title = paste("Correlation between", dir_name, "_channels Score and TMB"),
  vlcex = 1,
  vlabels = colnames(res),
  caxislabels = round(seq(res["Min", 1], res["Max", 1], length.out = 5), 1)
)
dev.off()
}


#GSVA泛癌MSI######  
ssGsea_Score_pan_cancer_MSI <- function(df){
load("D:/Dim/TCGAplot/All_MSI.Rdata")
df <- df %>%
  filter(Group == "Tumor") %>%
  tibble::add_column(ID = str_sub(rownames(.), 1, 12), .before = 1) %>%
  filter(!duplicated(ID)) %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("ID") %>%
  filter(rownames(.) %in% rownames(MSI)) %>%
  arrange(rownames(.))

MSI <- MSI %>%
  filter(rownames(.) %in% rownames(df)) %>%
  arrange(rownames(.))

stopifnot(identical(rownames(df), rownames(MSI)))

clist <- sort(unique(df$Cancer))
res <- list()

for (cancer in clist) {
  idx <- which(df$Cancer == cancer)
  if (length(idx) < 3) next
  x <- df$Score[idx]
  y <- MSI$MSI[idx]
  tmp <- suppressWarnings(cor.test(x, y, method = "pearson"))
  res[[cancer]] <- data.frame(rho = tmp$estimate, pvalue = tmp$p.value)
}

res_df <- do.call(rbind, res)
res_df$label <- ifelse(res_df$pvalue < 0.01, "**", ifelse(res_df$pvalue < 0.05, "*", " "))
res_df$group <- paste0(rownames(res_df), res_df$label)

res_df <- res_df[, c("group", "rho")] %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames("group") %>%
  tibble::add_column(Max = max(res_df$rho), .before = "rho") %>%
  tibble::add_column(Min = min(res_df$rho), .before = "rho") %>%
  t() %>% as.data.frame() %>% round(2)

res_df <- res_df[, c(1, ncol(res_df):2)]
file_name <- paste0(dir_name,"-pan-MSI.pdf")
pdf(file_name, width=7.5, height=7.5)
fmsb::radarchart(
  res_df,
  axistype = 1,
  pcol = "#00AFBB",
  pfcol = scales::alpha("#00AFBB", 0.1),
  plwd = 2,
  plty = 1,
  cglcol = "grey",
  cglty = 1,
  cglwd = 0.8,
  axislabcol = "black",
  title = paste("Correlation between", dir_name, "_channels Score and MSI"),
  vlcex = 1,
  vlabels = colnames(res_df),
  caxislabels = c(
    round(res_df["Min", 1], 1),
    round(res_df["Min", 1] + 0.25 * (res_df["Max", 1] - res_df["Min", 1]), 1),
    round(res_df["Min", 1] + 0.5 * (res_df["Max", 1] - res_df["Min", 1]), 1),
    round(res_df["Min", 1] + 0.75 * (res_df["Max", 1] - res_df["Min", 1]), 1),
    round(res_df["Max", 1], 1)
  )
)
dev.off()
}


#基因集的泛癌种与免疫相关基因的相关性分析#########
#基因集的泛癌种与ICG 的相关性#############
ssGsea_Score_pan_cancer_ICG <- function(df,tpm){
  file_name <- paste0(dir_name, "泛癌与ICG的相关性.pdf")
  pdf(file_name, width=8, height=2)
  clist <- list()
  plist <- list()
  checkpoint <- c("CD274", "CTLA4", "HAVCR2", "LAG3", "PDCD1", 
                  "PDCD1LG2", "SIGLEC15", "TIGIT")
  checkpoint <- checkpoint[checkpoint %in% colnames(tpm)]
  
  # 遍历癌种
  for (cancer in unique(df$Cancer)) {
    df_sub <- df %>% filter(Group == "Tumor", Cancer == cancer)
    samples <- df_sub$Score
    names(samples) <- rownames(df_sub)
    
    ic <- tpm[rownames(tpm) %in% names(samples), checkpoint, drop = FALSE]
    ic <- ic[match(names(samples), rownames(ic)), , drop = FALSE]
    
    stopifnot(identical(rownames(ic), names(samples)))
    ic[] <- lapply(ic, as.numeric)
    cor <- psych::corr.test(x = as.matrix(samples), y = ic, method = "pearson", adjust = "none")
    cmt <- cor$r
    pmt <- cor$p
    
    clist[[cancer]] <- as.data.frame(cmt)
    plist[[cancer]] <- as.data.frame(pmt)
  }
  
  # 合并
  clist <- do.call(rbind, clist)
  plist <- do.call(rbind, plist)
  clist[is.na(clist)] <- 0
  plist[is.na(plist)] <- 0.06
  plist <- as.data.frame(ifelse(plist < 0.001, "***",ifelse(plist < 0.01, "**", ifelse(plist < 0.05, "*", ""))))
  col_fun <- colorRamp2(
    c(min(clist), 0, max(clist)),
    c("blue", "white", "red")
  )
  # 显示符号矩阵，确保行列对应
  display_mat <- t(plist)[, rownames(clist)]
  p=Heatmap(
    matrix = t(clist)[, rownames(clist)],
    name = "Pearson",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(display_mat[i, j], x, y, gp = gpar(fontsize = 8))
    },
    row_names_side = "left",
    column_names_side = "bottom",        # 这里设置列名显示在底部
    column_names_gp = gpar(fontsize = 12),
    heatmap_legend_param = list(
      title = "Correlation",
      title_gp = gpar(fontsize = 8),
      labels_gp = gpar(fontsize = 8),
      legend_height = unit(2, "cm")
    ),
    rect_gp = gpar(col = "white")
  )
  print(p)
  dev.off()
}


#基因集的泛癌种与趋势化因子的相关性#############
ssGsea_Score_pan_cancer_chemokine <- function(df,tpm){
  file_name <- paste0(dir_name, "泛癌与趋势化因子的相关性.pdf")
  pdf(file_name, width=8, height=8)
# 初始化
clist <- list()
plist <- list()

checkpoint = c("CCL1", "CCL2", "CCL3", "CCL4", "CCL5", "CCL7", 
               "CCL8", "CCL11", "CCL13", "CCL14", "CCL15", "CCL16", 
               "CCL17", "CCL18", "CCL19", "CCL20", "CCL21", "CCL22", 
               "CCL23", "CCL24", "CCL25", "CCL26", "CCL28", "CX3CL1", 
               "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL6", "CXCL8", 
               "CXCL9", "CXCL10", "CXCL11", "CXCL12", "CXCL13", "CXCL14", 
               "CXCL16", "CXCL17")

checkpoint <- checkpoint[checkpoint %in% colnames(tpm)]

# 遍历癌种
for (cancer in unique(df$Cancer)) {
  df_sub <- df %>% filter(Group == "Tumor", Cancer == cancer)
  samples <- df_sub$Score
  names(samples) <- rownames(df_sub)
  
  ic <- tpm[rownames(tpm) %in% names(samples), checkpoint, drop = FALSE]
  ic <- ic[match(names(samples), rownames(ic)), , drop = FALSE]
  
  stopifnot(identical(rownames(ic), names(samples)))
  ic[] <- lapply(ic, as.numeric)
  cor <- psych::corr.test(x = as.matrix(samples), y = ic, method = "pearson", adjust = "none")
  cmt <- cor$r
  pmt <- cor$p
  
  clist[[cancer]] <- as.data.frame(cmt)
  plist[[cancer]] <- as.data.frame(pmt)
}

# 合并
clist <- do.call(rbind, clist)
plist <- do.call(rbind, plist)
clist[is.na(clist)] <- 0
plist[is.na(plist)] <- 0.06
plist <- as.data.frame(ifelse(plist < 0.001, "***",ifelse(plist < 0.01, "**", ifelse(plist < 0.05, "*", ""))))
col_fun <- colorRamp2(
  c(min(clist), 0, max(clist)),
  c("blue", "white", "red")
)
# 显示符号矩阵，确保行列对应
display_mat <- t(plist)[, rownames(clist)]
p=Heatmap(
  matrix = t(clist)[, rownames(clist)],
  name = "Pearson",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(display_mat[i, j], x, y, gp = gpar(fontsize = 8))
  },
  row_names_side = "left",
  column_names_side = "bottom",        # 这里设置列名显示在底部
  column_names_gp = gpar(fontsize = 12),
  heatmap_legend_param = list(
    title = "Correlation",
    title_gp = gpar(fontsize = 8),
    labels_gp = gpar(fontsize = 8),
    legend_height = unit(2, "cm")
  ),
  rect_gp = gpar(col = "white")
)
print(p)
dev.off()
}


#基因集的泛癌种与趋势化因子受体的相关性#############
ssGsea_Score_pan_cancer_chemokine_receptor <- function(df,tpm){
  file_name <- paste0(dir_name, "泛癌与趋势化因子受体的相关性.pdf")
  pdf(file_name, width=8, height=4)
# 初始化
clist <- list()
plist <- list()

checkpoint = c("CCR1", "CCR2", "CCR3", "CCR4", "CCR5", "CCR6", 
               "CCR7", "CCR8", "CCR9", "CCR10", "CXCR1", "CXCR2", "CXCR3", 
               "CXCR4", "CXCR5", "CXCR6", "XCR1", "CX3R1")

checkpoint <- checkpoint[checkpoint %in% colnames(tpm)]

# 遍历癌种
for (cancer in unique(df$Cancer)) {
  df_sub <- df %>% filter(Group == "Tumor", Cancer == cancer)
  samples <- df_sub$Score
  names(samples) <- rownames(df_sub)
  
  ic <- tpm[rownames(tpm) %in% names(samples), checkpoint, drop = FALSE]
  ic <- ic[match(names(samples), rownames(ic)), , drop = FALSE]
  
  stopifnot(identical(rownames(ic), names(samples)))
  ic[] <- lapply(ic, as.numeric)
  cor <- psych::corr.test(x = as.matrix(samples), y = ic, method = "pearson", adjust = "none")
  cmt <- cor$r
  pmt <- cor$p
  
  clist[[cancer]] <- as.data.frame(cmt)
  plist[[cancer]] <- as.data.frame(pmt)
}

# 合并
clist <- do.call(rbind, clist)
plist <- do.call(rbind, plist)
clist[is.na(clist)] <- 0
plist[is.na(plist)] <- 0.06
plist <- as.data.frame(ifelse(plist < 0.001, "***",ifelse(plist < 0.01, "**", ifelse(plist < 0.05, "*", ""))))
col_fun <- colorRamp2(
  c(min(clist), 0, max(clist)),
  c("blue", "white", "red")
)
# 显示符号矩阵，确保行列对应
display_mat <- t(plist)[, rownames(clist)]
p=Heatmap(
  matrix = t(clist)[, rownames(clist)],
  name = "Pearson",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(display_mat[i, j], x, y, gp = gpar(fontsize = 8))
  },
  row_names_side = "left",
  column_names_side = "bottom",        # 这里设置列名显示在底部
  column_names_gp = gpar(fontsize = 12),
  heatmap_legend_param = list(
    title = "Correlation",
    title_gp = gpar(fontsize = 8),
    labels_gp = gpar(fontsize = 8),
    legend_height = unit(2, "cm")
  ),
  rect_gp = gpar(col = "white")
)
print(p)
dev.off()
}


#基因集的泛癌种与免疫刺激剂的相关性############
ssGsea_Score_pan_cancer_immune_stimulator <- function(df,tpm){
  file_name <- paste0(dir_name, "泛癌与免疫刺激剂的相关性.pdf")
  pdf(file_name, width=8, height=8)
# 初始化
clist <- list()
plist <- list()

checkpoint = c("CD27", "CD276", "CD28", "CD40", "CD40LG", 
               "CD48", "CD70", "CD80", "CD86", "CXCL12", "CXCR4", "ENTPD1", 
               "HHLA2", "ICOS", "ICOSLG", "IL2RA", "IL6", "IL6R", "KLRC1", 
               "KLRK1", "LTA", "MICB", "NT5E", "PVR", "RAET1E", "TMIGD2", 
               "TNFRSF13B", "TNFRSF13C", "TNFRSF14", "TNFRSF17", "TNFRSF18", 
               "TNFRSF25", "TNFRSF4", "TNFRSF8", "TNFRSF9", "TNFSF13", 
               "TNFSF13B", "TNFSF14", "TNFSF15", "TNFSF18", "TNFSF4", 
               "TNFSF9", "ULBP1")

checkpoint <- checkpoint[checkpoint %in% colnames(tpm)]

# 遍历癌种
for (cancer in unique(df$Cancer)) {
  df_sub <- df %>% filter(Group == "Tumor", Cancer == cancer)
  samples <- df_sub$Score
  names(samples) <- rownames(df_sub)
  
  ic <- tpm[rownames(tpm) %in% names(samples), checkpoint, drop = FALSE]
  ic <- ic[match(names(samples), rownames(ic)), , drop = FALSE]
  
  stopifnot(identical(rownames(ic), names(samples)))
  ic[] <- lapply(ic, as.numeric)
  cor <- psych::corr.test(x = as.matrix(samples), y = ic, method = "pearson", adjust = "none")
  cmt <- cor$r
  pmt <- cor$p
  
  clist[[cancer]] <- as.data.frame(cmt)
  plist[[cancer]] <- as.data.frame(pmt)
}

# 合并
clist <- do.call(rbind, clist)
plist <- do.call(rbind, plist)
clist[is.na(clist)] <- 0
plist[is.na(plist)] <- 0.06
plist <- as.data.frame(ifelse(plist < 0.001, "***",ifelse(plist < 0.01, "**", ifelse(plist < 0.05, "*", ""))))
col_fun <- colorRamp2(
  c(min(clist), 0, max(clist)),
  c("blue", "white", "red")
)
# 显示符号矩阵，确保行列对应
display_mat <- t(plist)[, rownames(clist)]
p=Heatmap(
  matrix = t(clist)[, rownames(clist)],
  name = "Pearson",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(display_mat[i, j], x, y, gp = gpar(fontsize = 8))
  },
  row_names_side = "left",
  column_names_side = "bottom",        # 这里设置列名显示在底部
  column_names_gp = gpar(fontsize = 12),
  heatmap_legend_param = list(
    title = "Correlation",
    title_gp = gpar(fontsize = 8),
    labels_gp = gpar(fontsize = 8),
    legend_height = unit(2, "cm")
  ),
  rect_gp = gpar(col = "white")
)
print(p)
dev.off()
}


#基因集的泛癌种与免疫抑制剂的相关性############
ssGsea_Score_pan_cancer_immune_inhibitor <- function(df,tpm){
  file_name <- paste0(dir_name, "泛癌与免疫抑制剂的相关性.pdf")
  pdf(file_name, width=8, height=5)
# 初始化
clist <- list()
plist <- list()

checkpoint = c("ADORA2A", "BTLA", "CD160", "CD244", "CD274", 
               "CD96", "CSF1R", "CTLA4", "HAVCR2", "IDO1", "IL10", "IL10RB", 
               "KDR", "KIR2DL1", "KIR2DL3", "LAG3", "LGALS9", "PDCD1", 
               "PDCD1LG2", "TGFB1", "TGFBR1", "TIGIT", "VTCN1")

checkpoint <- checkpoint[checkpoint %in% colnames(tpm)]

# 遍历癌种
for (cancer in unique(df$Cancer)) {
  df_sub <- df %>% filter(Group == "Tumor", Cancer == cancer)
  samples <- df_sub$Score
  names(samples) <- rownames(df_sub)
  
  ic <- tpm[rownames(tpm) %in% names(samples), checkpoint, drop = FALSE]
  ic <- ic[match(names(samples), rownames(ic)), , drop = FALSE]
  
  stopifnot(identical(rownames(ic), names(samples)))
  ic[] <- lapply(ic, as.numeric)
  cor <- psych::corr.test(x = as.matrix(samples), y = ic, method = "pearson", adjust = "none")
  cmt <- cor$r
  pmt <- cor$p
  
  clist[[cancer]] <- as.data.frame(cmt)
  plist[[cancer]] <- as.data.frame(pmt)
}

# 合并
clist <- do.call(rbind, clist)
plist <- do.call(rbind, plist)
clist[is.na(clist)] <- 0
plist[is.na(plist)] <- 0.06
plist <- as.data.frame(ifelse(plist < 0.001, "***",ifelse(plist < 0.01, "**", ifelse(plist < 0.05, "*", ""))))
col_fun <- colorRamp2(
  c(min(clist), 0, max(clist)),
  c("blue", "white", "red")
)
# 显示符号矩阵，确保行列对应
display_mat <- t(plist)[, rownames(clist)]
p=Heatmap(
  matrix = t(clist)[, rownames(clist)],
  name = "Pearson",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(display_mat[i, j], x, y, gp = gpar(fontsize = 8))
  },
  row_names_side = "left",
  column_names_side = "bottom",        # 这里设置列名显示在底部
  column_names_gp = gpar(fontsize = 12),
  heatmap_legend_param = list(
    title = "Correlation",
    title_gp = gpar(fontsize = 8),
    labels_gp = gpar(fontsize = 8),
    legend_height = unit(2, "cm")
  ),
  rect_gp = gpar(col = "white")
)
print(p)
dev.off()
}


#基因集的泛癌种与免疫细胞比率相关性#########
ssGsea_Score_pan_cancer_immunecellratio <- function(df){
load("D:/Dim/TCGAplot/ssgsea_immucells.Rdata")
file_name <- paste0(dir_name, "泛癌与免疫细胞比率的相关性.pdf")
pdf(file_name, width=8, height=5)
clist <- list()
plist <- list()
for (cancer in unique(df$Cancer)) {
  df_sub <- df[df$Cancer == cancer & df$Group == "Tumor", ]
  row_ids <- intersect(rownames(df_sub), rownames(immucells))
  df_sub <- df_sub[row_ids, , drop = FALSE]
  immu_sub <- immucells[row_ids, , drop = FALSE]
  
  if (nrow(df_sub) < 3) next
  stopifnot(identical(rownames(df_sub), rownames(immu_sub)))
  cor <- psych::corr.test(df_sub$Score, immu_sub, method = "pearson", adjust = "none")
  clist[[cancer]] <- as.data.frame(cor$r)
  plist[[cancer]] <- as.data.frame(cor$p)
}

# 合并并处理
clist <- do.call(rbind, clist)
plist <- do.call(rbind, plist)
clist[is.na(clist)] <- 0
plist[is.na(plist)] <- 0.06
plist <- as.data.frame(ifelse(plist < 0.001, "***",ifelse(plist < 0.01, "**", ifelse(plist < 0.05, "*", ""))))
col_fun <- colorRamp2(
  c(min(clist), 0, max(clist)),
  c("blue", "white", "red")
)
# 显示符号矩阵，确保行列对应
display_mat <- t(plist)[, rownames(clist)]
p=Heatmap(
  matrix = t(clist)[, rownames(clist)],
  name = "Pearson",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(display_mat[i, j], x, y, gp = gpar(fontsize = 8))
  },
  row_names_side = "left",
  column_names_side = "bottom",        # 这里设置列名显示在底部
  column_names_gp = gpar(fontsize = 12),
  heatmap_legend_param = list(
    title = "Correlation",
    title_gp = gpar(fontsize = 8),
    labels_gp = gpar(fontsize = 8),
    legend_height = unit(2, "cm")
  ),
  rect_gp = gpar(col = "white")
)
print(p)
dev.off()
}


#基因集的泛癌种与免疫评分的相关性############
ssGsea_Score_pan_cancer_immunescore <- function(df){
load("D:/Dim/TCGAplot/immuscore.Rdata")
file_name <- paste0(dir_name, "泛癌与免疫评分的相关性.pdf")
pdf(file_name, width=8, height=1.3)
clist <- list()
plist <- list()
for (cancer in unique(df$Cancer)) {
  df_sub <- df[df$Cancer == cancer & df$Group == "Tumor", ]
  row_ids <- intersect(rownames(df_sub), rownames(immuscore))
  if (length(row_ids) < 3) next
  
  df_sub <- df_sub[row_ids, , drop = FALSE]
  immu_sub <- immuscore[row_ids, , drop = FALSE]
  stopifnot(identical(rownames(df_sub), rownames(immu_sub)))
  
  cor <- psych::corr.test(df_sub$Score, immu_sub, method = "pearson", adjust = "none")
  clist[[cancer]] <- as.data.frame(cor$r)
  plist[[cancer]] <- as.data.frame(cor$p)
}

# 合并
clist <- do.call(rbind, clist)
plist <- do.call(rbind, plist)

# 处理NA
clist[is.na(clist)] <- 0
plist[is.na(plist)] <- 0.06
plist <- as.data.frame(ifelse(plist < 0.001, "***",ifelse(plist < 0.01, "**", ifelse(plist < 0.05, "*", ""))))
col_fun <- colorRamp2(
  c(min(clist), 0, max(clist)),
  c("blue", "white", "red")
)
# 显示符号矩阵，确保行列对应
display_mat <- t(plist)[, rownames(clist)]
p=Heatmap(
  matrix = t(clist)[, rownames(clist)],
  name = "Pearson",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(display_mat[i, j], x, y, gp = gpar(fontsize = 8))
  },
  row_names_side = "left",
  column_names_side = "bottom",        # 这里设置列名显示在底部
  column_names_gp = gpar(fontsize = 12),
  heatmap_legend_param = list(
    title = "Correlation",
    title_gp = gpar(fontsize = 8),
    labels_gp = gpar(fontsize = 8)#,
    # legend_height = unit(1, "cm")
  ),
  rect_gp = gpar(col = "white")
)
print(p)
dev.off()
}


#基因集的泛癌种与TME_Pathaway通路评分的相关性############

# library(tibble)
# cellMarker1 <- read.table("TME相关通路.txt", header = FALSE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)
# cellMarker1 <- cellMarker1 %>% column_to_rownames("V1") %>% t()
# a <- as.data.frame(cellMarker1, stringsAsFactors = FALSE, check.names = FALSE)
# geneSet <- list()
# for (i in colnames(a)) {
#   x <- as.character(a[[i]])
#   x <- trimws(x)          # 去掉前后空格
#   x <- x[nchar(x) > 0]    # 去除空字符串
#   geneSet[[i]] <- x
# }
# 
# load("D:/Dim/基因组泛癌/数据/表达矩阵/all_tpm.Rdata")
# tpm_numeric <- tpm %>%
#   mutate(across(-c(Cancer, Group), as.numeric))
# # 提取表达矩阵并转置
# expr_data <- as.matrix(tpm_numeric[, -c(1:2)])
# colnames(expr_data) <- colnames(tpm_numeric)[-c(1:2)]
# expr_matrix <- t(expr_data)
# 
# gsva_results <- gsva(expr = expr_matrix,
#                      gset.idx.list = geneSet,
#                      method = "gsva",  # 也可用"ssgsea"
#                      kcdf = "Gaussian", # 对log2转换数据使用Gaussian
#                      parallel.sz = 16         # 多核加速
#                      
# )
# 
# gsva_results <- as.data.frame(t(gsva_results))
# tpm <- as.data.frame(tpm)
# # 确保样本名一致并匹配
# common_samples <- intersect(rownames(gsva_results), rownames(tpm))
# gsva_results <- gsva_results[common_samples, , drop = FALSE]
# tpm_sub <- tpm[common_samples, c("Cancer", "Group")]
# # 添加新列
# TME_Pathaway=gsva_results
# save(TME_Pathaway,file = "TME_Pathaway_GSVA_Score.Rdata")


ssGsea_Score_pan_cancer_TME_GSVA <- function(df){
  load("D:/Dim/TCGAplot/TME_Pathaway_GSVA_Score.Rdata")
  file_name <- paste0(dir_name, "评分与TME相关通路评分的相关性.pdf")
  pdf(file_name, width=10, height=4.5)
  clist <- list()
  plist <- list()
  for (cancer in unique(df$Cancer)) {
    df_sub <- df[df$Cancer == cancer & df$Group == "Tumor", ]
    row_ids <- intersect(rownames(df_sub), rownames(TME_Pathaway))
    if (length(row_ids) < 3) next
    
    df_sub <- df_sub[row_ids, , drop = FALSE]
    immu_sub <- TME_Pathaway[row_ids, , drop = FALSE]
    stopifnot(identical(rownames(df_sub), rownames(immu_sub)))
    immu_sub[] <- lapply(immu_sub, as.numeric)
    cor <- psych::corr.test(df_sub$Score, immu_sub, method = "pearson", adjust = "none")
    clist[[cancer]] <- as.data.frame(cor$r)
    plist[[cancer]] <- as.data.frame(cor$p)
  }
  
  # 合并
  clist <- do.call(rbind, clist)
  plist <- do.call(rbind, plist)
  
  # 处理NA
  clist[is.na(clist)] <- 0
  plist[is.na(plist)] <- 0.06
  plist <- as.data.frame(ifelse(plist < 0.001, "***",ifelse(plist < 0.01, "**", ifelse(plist < 0.05, "*", ""))))
  col_fun <- colorRamp2(
    c(min(clist), 0, max(clist)),
    c("blue", "white", "red")
  )
  # 显示符号矩阵，确保行列对应
  display_mat <- t(plist)[, rownames(clist)]
  p=Heatmap(
    matrix = t(clist)[, rownames(clist)],
    name = "Pearson",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(display_mat[i, j], x, y, gp = gpar(fontsize = 8))
    },
    row_names_side = "left",
    column_names_side = "bottom",        # 这里设置列名显示在底部
    column_names_gp = gpar(fontsize = 12),
    heatmap_legend_param = list(
      title = "Correlation",
      title_gp = gpar(fontsize = 8),
      labels_gp = gpar(fontsize = 8),
      legend_height = unit(2, "cm")
    ),
    rect_gp = gpar(col = "white")
  )
  print(p)
  dev.off()
}


#泛癌森林图######
ssGsea_Score_pan_cancer_fores <- function(df,Target){
  load("D:/Dim/基因组泛癌/数据/临床数据/ALL_OS.Rdata")
  load("D:/Dim/基因组泛癌/数据/临床数据/DSS_PFI_DFI.Rdata")
  if (Target == "OS") {
    clinical_data <- OS[, -c(1:2)]
    merged <- merge(df, clinical_data, by = 0, sort = FALSE) %>%
      filter(Group == "Tumor") %>%
      column_to_rownames("Row.names")
  } else {
    idx <- switch(Target,
                  "DSS" = 1:2,
                  "PFI" = 3:4,
                  "DFI" = 5:6,
                  stop("Invalid Target"))
    clinical_data <- DSS_PFI_DFI[, idx]
    colnames(clinical_data) <- c("OS.time", "OS")
    
    df$MatchID <- substr(rownames(df), 1, 15)
    clinical_data$MatchID <- rownames(clinical_data)
    
    merged <- merge(df, clinical_data, by = "MatchID", sort = FALSE)
    rownames(merged) <- make.unique(merged$MatchID)  # 防止重复行名
    merged$MatchID <- NULL
    merged <- merged[merged$Group == "Tumor", ]
    merged <- na.omit(merged)
  }
  
  
  
  total_list <- table(merged$Cancer)  # 统计各癌种样本数
  # 若包含age调整
  cox_results <- list()
  for (cancer in unique(merged$Cancer)) {
    dat <- merged[merged$Cancer == cancer, ]
    if (nrow(dat) < 10) next
    dat$symbol <- dat$Score
    dat$Type <- ifelse(dat$symbol < median(dat$symbol, na.rm = TRUE), "Low", "High")
    dat$Type <- factor(dat$Type, levels = c("Low", "High"))
    m <- coxph(Surv(OS.time, OS) ~ Type, data = dat)
    beta <- coef(m)
    se <- sqrt(diag(vcov(m)))
    HR <- exp(beta)
    HRse <- HR * se
    
    tmp <- round(cbind(
      coef = beta, se = se, z = beta / se,
      p = 1 - pchisq((beta / se)^2, 1),
      HR = HR, HRse = HRse,
      HRz = (HR - 1) / HRse,
      HRp = 1 - pchisq(((HR - 1) / HRse)^2, 1),
      HRCILL = exp(beta - qnorm(0.975) * se),
      HRCIUL = exp(beta + qnorm(0.975) * se)
    ), 3)
    
    cox_results[[cancer]] <- tmp["TypeHigh", ]
    
  }
  
  # 合并结果
  cox_results <- do.call(rbind, cox_results)
  cox_results <- as.data.frame(cox_results[, c("HR", "HRCILL", "HRCIUL", "p")])
  np <- paste0(cox_results$HR, " (", cox_results$HRCILL, "-", cox_results$HRCIUL, ")")
  
  tabletext <- cbind(
    c("Cancer", rownames(cox_results)),
    c("Total(N)", as.character(total_list[rownames(cox_results)])),
    c("HR (95%CI)", np),
    c("P value", ifelse(cox_results$p < 0.001, "<0.001", round(cox_results$p, 3)))
  )
  
  
  # 根据条件为文本着色
  color_vector <- ifelse(cox_results$p < 0.05 & cox_results$HR > 1, "#E64B35", 
                         ifelse(cox_results$p < 0.05 & cox_results$HR < 1, "#4DBBD5", "black"))
  color_vector <- c("black", color_vector)  # 在第一个元素前加入黑色
  
  
  # 计算行数
  num_rows <- nrow(cox_results) + 2  
  
  # 动态生成 hrzl_lines
  hrzl_lines <- list(
    `1` = grid::gpar(lwd = 2, col = "black"),  # 标题下方的线
    `2` = grid::gpar(lwd = 2, col = "black")   # 第一行数据下方的线
  )
  hrzl_lines[[as.character(num_rows)]] <- grid::gpar(lwd = 2, col = "black")  # 最后一行数据下方的线
  
  p=forestplot(
    labeltext = tabletext,
    title = paste0(dir_name, "_channels score ", Target), 
    title_gp = grid::gpar(fontsize = 14, fontface = "bold", just = "left"),
    mean = c(NA, cox_results$HR),
    lower = c(NA, cox_results$HRCILL),
    upper = c(NA, cox_results$HRCIUL),
    hrzl_lines = hrzl_lines, 
    is.summary = rep(FALSE, nrow(cox_results) + 1),
    zero = 1,
    boxsize = 0.4,
    lineheight = unit(0.2, "cm"),
    xlog = FALSE,
    col = fpColors(box = "#4DBBD5", lines = "black", zero = "gray50", text = color_vector), # 设置文本颜色
    lwd.ci = 2.5,
    ci.vertices = TRUE,
    ci.vertices.height = 0.02,
    clip = c(0, 6),
    xticks = seq(0, 6, by = 2),
    graph.pos = 4,
    graphwidth = unit(5, "cm"),
    txt_gp = forestplot::fpTxtGp(ticks = grid::gpar(cex = 0.9)), 
  )
  file_name <- paste0(dir_name,"-",Target,"-泛癌森林图.pdf")
  pdf(file_name, width=7, height=9)
  # 创建新页面
  grid.newpage()
  # 绘制森林图
  print(p, newpage = FALSE)  # 阻止forestplot创建新页面
  # 关闭PDF
  dev.off()
  dev.off()
}

#基于基因集的癌症类型特异性表达分析##################
#按临床信息分组的基于基因组的表达分析######
ssGsea_Score_Tumor_Normal <- function(df, cancer) {
  ylab_text <- paste0(dir_name, "_channels score")
  
  df <- df[df$Cancer == cancer, -c(2, 3), drop = FALSE]
  df1 <- filter(tpm, Cancer == cancer) %>% select(Cancer, Group)
  dat <- cbind(df1, df)
  
  # 计算p值并定义label
  p_val <- wilcox.test(Score ~ Group, data = dat)$p.value
  stat_df <- data.frame(
    group1 = "Normal",
    group2 = "Tumor",
    p = p_val,
    y.position = max(dat$Score, na.rm = TRUE) * 1.05
  )
  stat_df$label <- cut(stat_df$p,
                       breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("***", "**", "*", "ns"))
  
  p <- ggboxplot(dat, x = "Group", y = "Score", color = "Group",
                 add = "jitter", palette = "jco",
                 xlab = cancer, ylab = ylab_text) +
    stat_pvalue_manual(stat_df, label = "label", tip.length = 0.01, size = 4) +
    theme_classic2(base_size = 12) +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      axis.text.y = element_text(size = 12)
    ) +
    border("black")
  
  file_name <- paste0(dir_name, "-", cancer, "-正常-肿瘤对比差异.pdf")
  pdf(file_name, width = 5, height = 4)
  print(p)
  dev.off()
}


#按临床信息分组的基于基因组的样本配对肿瘤-正常箱线图（一张图）######
ssGsea_Score_paired_boxplot <- function(df){
  ylab_text <- paste0(dir_name, "_channels score")
  pcancers <- c("BLCA", "BRCA", "COAD", "ESCA", "HNSC", "KICH", 
                "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", 
                "STAD", "THCA", "UCEC")
  
  load("D:/Dim/TCGAplot/paired_tpm.Rdata")
  
  plot_list <- lapply(pcancers, function(cancer) {
    df_cancer <- df[df$Cancer == cancer, , drop = FALSE]
    paired_cancer <- filter(paired_tpm, Cancer == cancer) %>% select(Group, ID)
    df_cancer <- df_cancer[rownames(paired_cancer), , drop = FALSE]
    if (!identical(rownames(df_cancer), rownames(paired_cancer))) return(NULL)
    
    dat <- cbind(paired_cancer, df_cancer[, "Score", drop = FALSE])
    dat$paired <- dat$ID
    dat$Type <- cancer
    return(dat)
  })
  
  dat_all <- do.call(rbind, plot_list)
  
  stat_df <- compare_means(Score ~ Group, data = dat_all, 
                           group.by = "Type", paired = TRUE)
  stat_df$y.position <- tapply(dat_all$Score, dat_all$Type, max, na.rm = TRUE) + 0.05
  
  # 自定义p值星号标签
  stat_df$label <- cut(stat_df$p,
                       breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("***", "**", "*", "ns")
  )
  
  stat_df$group1 <- "Normal"
  stat_df$group2 <- "Tumor"
  
  p <- ggpaired(dat_all, x = "Group", y = "Score", color = "black", fill = "Group", line.size = 0.5) +
    geom_line(aes(group = paired), color = "grey80") +
    facet_wrap(~ Type, scales = 'free_x', nrow = 2) +
    scale_fill_manual(values = c("Normal" = "#4DBBD5", "Tumor" = "#E64B35")) +
    scale_color_manual(values = c("Normal" = "#4DBBD5", "Tumor" = "#E64B35")) +
    stat_pvalue_manual(stat_df, label = "label", tip.length = 0.01, size = 2) +
    theme_classic2(base_size = 9) +
    xlab('Pan-Cancer') +
    ylab(ylab_text)
  
  file_name <- paste0(dir_name, "-正常-肿瘤配对差异.pdf")
  pdf(file_name, width = 8, height = 5)
  print(p)
  dev.off()
}


#按临床信息分组的基于基因组的样本配对肿瘤-正常箱线图（单肿瘤）######
ssGsea_Score_one_paired_boxplot <- function(df, cancer) {
  ylab_text <- paste0(dir_name, "_channels score")
  palette <- c("Normal" = "#4DBBD5", "Tumor" = "#E64B35")  # ggsci::pal_npg()[1:2]
  legend <- "none"
  label <- "p.signif"
  method <- "wilcox.test"
  
  load("D:/Dim/TCGAplot/paired_tpm.Rdata")
  df <- df[df$Cancer == cancer, -c(2, 3), drop = FALSE]
  df1 <- filter(paired_tpm, Cancer == cancer) %>% select(Cancer, Group, ID)
  
  df <- df[rownames(df1), , drop = FALSE]
  stopifnot(identical(rownames(df), rownames(df1)))
  
  dat <- cbind(df1, df)
  dat$id <- dat$ID
  dat$Group <- factor(dat$Group, levels = c("Normal", "Tumor"))
  
  my_comparisons <- list(c("Normal", "Tumor"))
  
  p <- ggpaired(dat,
                x = "Group", y = "Score", id = "id",
                color = "black", fill = "Group", palette = palette,
                line.size = 0.5, add = "jitter") +
    geom_line(aes(group = id), color = "grey80") +
    stat_compare_means(comparisons = my_comparisons, label = label, method = method, paired = TRUE,
                       label.x.npc = "center", label.y.npc = "top", size = 3) +
    theme_classic2(base_size = 9) +
    theme(legend.position = legend,
          panel.border = element_rect(color = "black", fill = NA, size = 0.5),
          axis.line = element_line(color = "black"),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12, colour = "black"),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12, colour = "black"),
          strip.text.x = element_text(size = 12, color = "black", face = "bold")) +
    xlab(cancer) + ylab(ylab_text)
  
  file_name <- paste0(dir_name, "-", cancer, "正常-肿瘤配对差异.pdf")
  pdf(file_name, width = 5, height = 4)
  print(p)
  dev.off()
}


#基因集泛癌ROC曲线#########
ssGsea_Score_ROC <- function(df,cancer){
ylab_text <- paste0(dir_name, "_channels score")
df <- df[df$Cancer == cancer, -c(2, 3), drop = FALSE]
df1 = subset(tpm, Cancer == cancer) %>% select(Cancer, Group)
stopifnot(identical(rownames(df), rownames(df1)))
identical(rownames(df), rownames(df1))
dat = cbind(df1, df)
res <- pROC::roc_(dat, "Group", "Score", aur = TRUE, 
                  ci = TRUE, smooth = TRUE, levels = c("Normal", "Tumor"))
p <- pROC::ggroc(res, color = "#56B4E9", size = 1, legacy.axes = TRUE) + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
  color = "#505050", linetype = "dashed", size = 0.75)+ 
  theme_bw() +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        panel.grid.major = element_line(colour = "lightgrey"),
        panel.grid.minor = element_line(colour = "lightgrey"),
        panel.background = element_rect(fill = "white")) +
  ggplot2::annotate("text", x = 0.75, y = 0.05, 
                    label = paste("AUC: ", round(res$auc, 3)), 
                    size = 4, hjust = 0, fontface = "plain") +
  ggplot2::annotate("text", x = 0.75, y = 0, 
                    label = paste("CI: ", round(res$ci[1], 3), "-", round(res$ci[3], 3)), 
                    size = 4, hjust = 0, fontface = "plain") +
  ggplot2::annotate("text", x = 0.75, y = 0.1,
                    label = ylab_text,
                    size = 4, hjust = 0, fontface = "plain")
file_name <- paste0(dir_name,"-",cancer,"-ROC曲线.pdf")
pdf(file_name, width=6, height=5)
print(p)
dev.off()
}


#基因集泛癌生存KM###############
ssGsea_Score_KM <- function(df,cancer){
ylab_text <- paste0(dir_name, "_channels score")
load("D:/Dim/基因组泛癌/数据/临床数据/ALL_OS.Rdata")
merged <- merge(df, OS[-c(1:2)], by = 0) %>% 
  filter(Group == "Tumor" & Cancer == cancer)%>%
  column_to_rownames("Row.names") %>%
  select(Score, OS.time, OS)
cl=merged
cl$Type <- ifelse(cl[["Score"]] > median(cl[["Score"]], na.rm = TRUE), "high", "low")
sfit = survival::survfit(survival::Surv(OS.time, OS) ~ Type, 
                         data = cl)
p=survminer::ggsurvplot(sfit, 
                      pval = TRUE, 
                      palette = "jco", 
                      data = cl,
                      legend = c(0.8, 0.8),
                      title = paste0("KMplot of ",
                      ylab_text, " in ", cancer), risk.table = T)

file_name <- paste0(dir_name,"-",cancer,"-KM曲线.pdf")
pdf(file_name, width=6, height=5)
# 创建新页面
grid.newpage()
# 绘制森林图
print(p, newpage = FALSE)  # 阻止forestplot创建新页面
# 关闭PDF
dev.off()
}


##HALLMARK通路相关性分析#########
# library(GSEABase)
# library(GSVA)
# library(dplyr)
# hallmark_sets <- getGmt("h.all.v2025.1.Hs.symbols.gmt")
# hallmark_list <- lapply(hallmark_sets, geneIds)
# load("D:/Dim/基因组泛癌/数据/表达矩阵/all_tpm.Rdata")
# tpm_numeric <- tpm %>%
#   mutate(across(-c(Cancer, Group), as.numeric))
# # 提取表达矩阵并转置
# expr_data <- as.matrix(tpm_numeric[, -c(1:2)])
# colnames(expr_data) <- colnames(tpm_numeric)[-c(1:2)]
# expr_matrix <- t(expr_data)
# 
# # 运行GSVA
# gsva_results <- gsva(expr = expr_matrix,
#                      gset.idx.list = hallmark_sets,
#                      method = "gsva",  # 也可用"ssgsea"
#                      kcdf = "Gaussian", # 对log2转换数据使用Gaussian
#                      parallel.sz = 4,          # 多核加速
#                      verbose = TRUE            # 显示进度
# )
# 
# gsva_results <- as.data.frame(t(gsva_results))
# tpm <- as.data.frame(tpm)
# # 确保样本名一致并匹配
# common_samples <- intersect(rownames(gsva_results), rownames(tpm))
# gsva_results <- gsva_results[common_samples, , drop = FALSE]
# tpm_sub <- tpm[common_samples, c("Cancer", "Group")]
# # 添加新列
# gsva_results$Cancer <- tpm_sub$Cancer
# gsva_results$Group  <- tpm_sub$Group
# save(gsva_results,file = "HALLMARK_GSVA_Score.Rdata")





#基因集的泛癌种与HALLMARK通路评分的相关性############
ssGsea_Score_pan_cancer_HALLMARK_GSVA <- function(df){
  load("D:/Dim/TCGAplot/HALLMARK_GSVA_Score.Rdata")
  colnames(HALLMARK) <- gsub("^HALLMARK_", "", colnames(HALLMARK))
  file_name <- paste0(dir_name, "评分与HALLMARK通路评分的相关性.pdf")
  pdf(file_name, width=10, height=10)
  clist <- list()
  plist <- list()
  for (cancer in unique(df$Cancer)) {
    df_sub <- df[df$Cancer == cancer & df$Group == "Tumor", ]
    row_ids <- intersect(rownames(df_sub), rownames(HALLMARK))
    if (length(row_ids) < 3) next
    
    df_sub <- df_sub[row_ids, , drop = FALSE]
    immu_sub <- HALLMARK[row_ids, , drop = FALSE]
    stopifnot(identical(rownames(df_sub), rownames(immu_sub)))
    
    cor <- psych::corr.test(df_sub$Score, immu_sub, method = "pearson", adjust = "none")
    clist[[cancer]] <- as.data.frame(cor$r)
    plist[[cancer]] <- as.data.frame(cor$p)
  }
  
  # 合并
  clist <- do.call(rbind, clist)
  plist <- do.call(rbind, plist)
  
  # 处理NA
  clist[is.na(clist)] <- 0
  plist[is.na(plist)] <- 0.06
  plist <- as.data.frame(ifelse(plist < 0.001, "***",ifelse(plist < 0.01, "**", ifelse(plist < 0.05, "*", ""))))
  col_fun <- colorRamp2(
    c(min(clist), 0, max(clist)),
    c("blue", "white", "red")
  )
  # 显示符号矩阵，确保行列对应
  display_mat <- t(plist)[, rownames(clist)]
  p=Heatmap(
    matrix = t(clist)[, rownames(clist)],
    name = "Pearson",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(display_mat[i, j], x, y, gp = gpar(fontsize = 8))
    },
    row_names_side = "left",
    column_names_side = "bottom",        # 这里设置列名显示在底部
    row_names_gp = gpar(fontsize = 8.4),
    column_names_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(
      title = "Correlation",
      title_gp = gpar(fontsize = 8),
      labels_gp = gpar(fontsize = 8)#,
      # legend_height = unit(1, "cm")
    ),
    rect_gp = gpar(col = "white")
  )
  print(p)
  dev.off()
}




#突变汇总########
Alterations <- function(folder_path) {
  file_name <- paste0(folder_path, "/alterations_across_samples_", dir_name, ".tsv")
  alter <- read.delim(file_name, header = TRUE, sep = "\t", check.names = FALSE, quote = "")
  alter <- alter[, -3]
  alter <- as.data.frame(alter)
  alter$Cancer <- toupper(sub("_.*", "", alter[[1]]))
  alter <- alter[, -1]  # 删除第一列
  alter <- alter[, -1]  # 删除第一列
  alter <- alter[, c(ncol(alter), 1:(ncol(alter)-1))]  # 把最后一列移到最前
  alter <- alter[, c(1, 2, grep(":", colnames(alter)))]
  alter <- alter[alter[, 2] == 1, ]  # 只保留Altered == 1的样本
  
  # 整理成 long format
  long_df <- alter %>%
    pivot_longer(-c(Cancer, Altered), names_to = "GeneType", values_to = "Status") %>%
    filter(Status != "no alteration", Status != "not profiled") %>%
    separate(GeneType, into = c("Gene", "Type"), sep = ": ") %>%
    mutate(Type = trimws(Type))
  
  # 突变数量矩阵
  mut_count <- long_df %>%
    group_by(Gene, Cancer) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = Cancer, values_from = n, values_fill = 0) %>%
    column_to_rownames("Gene") %>%
    as.matrix()
  
  # 计算每个癌症类型中有突变的样本数
  total_samples_with_mutation <- alter %>%
    group_by(Cancer) %>%
    summarise(Total = n(), .groups = "drop") %>%
    column_to_rownames("Cancer") %>%
    as.matrix()
  
  total_vec <- total_samples_with_mutation[colnames(mut_count), 1]
  mut_count <- rbind(mut_count, Total = total_vec)
  
  # 行突变类型占比
  row_mut_types_all <- long_df %>%
    count(Gene, Type) %>%
    group_by(Gene) %>%
    mutate(freq = n / sum(n)) %>%
    select(-n) %>%
    pivot_wider(names_from = Type, values_from = freq, values_fill = 0) %>%
    column_to_rownames("Gene")
  
  if(!"Total" %in% rownames(row_mut_types_all)){
    row_mut_types_all <- rbind(row_mut_types_all,
                               Total = setNames(rep(0, ncol(row_mut_types_all)), colnames(row_mut_types_all)))
  }
  row_mut_types_all <- row_mut_types_all[rownames(mut_count), , drop = FALSE]
  
  # 列突变类型占比
  col_mut_types_all <- long_df %>%
    count(Cancer, Type) %>%
    group_by(Cancer) %>%
    mutate(freq = n / sum(n)) %>%
    select(-n) %>%
    pivot_wider(names_from = Type, values_from = freq, values_fill = 0) %>%
    column_to_rownames("Cancer")
  
  col_mut_types_all <- col_mut_types_all[colnames(mut_count), , drop = FALSE]
  
  # 颜色函数
  col_fun <- Vectorize(function(x) {
    if (x < 10) "white"
    else if (x < 20) "#ded3cf"
    else if (x < 40) "#669dbb"
    else if (x < 60) "#fdb1a4"
    else "#99253a"
  })
  
  # 堆叠图颜色
  stack_colors <- c(
    "MUT" = "#fc9a99",
    "AMP" = "#b2df8a",
    "FUSION" = "#34a12e",
    "HOMDEL" = "#1f78b4"
  )
  all_types <- union(colnames(row_mut_types_all), colnames(col_mut_types_all))
  missing_types <- setdiff(all_types, names(stack_colors))
  if(length(missing_types) > 0){
    stack_colors[missing_types] <- rainbow(length(missing_types))
  }
  
  # 尺寸参数
  n_rows <- nrow(mut_count)
  n_cols <- ncol(mut_count)
  row_height <- unit(6, "mm")
  col_width <- unit(6, "mm")
  
  # 右侧堆叠图，高度与热图一致
  row_anno_all <- rowAnnotation(
    " " = anno_barplot(as.matrix(row_mut_types_all),
                       gp = gpar(fill = stack_colors[colnames(row_mut_types_all)]),
                       border = FALSE,
                       bar_width = 0.8,
                       axis = FALSE,
                       show_annotation_name = FALSE,
                       height = n_rows * row_height)
  )
  
  # 上方堆叠图，宽度与热图一致
  col_anno_all <- columnAnnotation(
    " " = anno_barplot(as.matrix(col_mut_types_all),
                       gp = gpar(fill = stack_colors[colnames(col_mut_types_all)]),
                       border = FALSE,
                       bar_width = 0.8,
                       axis = FALSE,
                       show_annotation_name = FALSE,
                       width = n_cols * col_width)
  )
  
  # 图例
  legend_top <- Legend(labels = names(stack_colors),
                       title = "Alter",
                       legend_gp = gpar(fill = stack_colors),
                       nrow = 1,
                       title_gp = gpar(fontsize = 11),
                       labels_gp = gpar(fontsize = 10),
                       title_position = "leftcenter")
  
  # 热图
  ht <- Heatmap(mut_count,
                col = col_fun,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                top_annotation = col_anno_all,
                right_annotation = row_anno_all,
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(mut_count[i, j], x, y, gp = gpar(fontsize = 8))
                  grid.rect(x, y, width, height, gp = gpar(col = "grey80", fill = NA, lwd = 0.5))
                },
                show_column_names = TRUE,
                show_row_names = TRUE,
                row_names_side = "left",
                row_names_gp = gpar(fontsize = 11),
                column_names_gp = gpar(fontsize = 11),
                show_heatmap_legend = FALSE,
                row_names_max_width = max_text_width(rownames(mut_count), gp = gpar(fontsize = 11)),
                column_names_max_height = max_text_width(colnames(mut_count), gp = gpar(fontsize = 11))
  )
  
  # 导出 PDF
  pdf_file <- paste0(dir_name, "-突变汇总.pdf")
  pdf(pdf_file, width = 10, height = 6)
  draw(ht, annotation_legend_list = list(legend_top), annotation_legend_side = "top")
  dev.off()
}





#CNV和甲基化与mRNA的关系###############
CNV_Methylation_mRNA <- function(analysis_target) {
  # 根据目标选择文件并设置标题
  if(analysis_target == "CNV") {
    file_name <- "CnvAndExpressionTable.xlsx"
    plot_title <- "Correlations of CNV with mRNA expression"
    save_suffix <- "_CNV-mRNA"
  } else if(analysis_target == "Methylation") {
    file_name <- "ExpressionAndMethylationTable.xlsx"
    plot_title <- "Correlation between methylation and\nmRNA expression"
    save_suffix <- "_甲基化-mRNA"  # 修改为"-mRNA"保持一致性
  } else {
    stop("Invalid analysis target. Please use 'CNV' or 'Methylation'")
  }
  
  # 检查文件是否存在
  if(!file.exists(file_name)) {
    stop(paste("File not found:", file_name, 
               "\nPlease ensure the file is in your working directory:",
               getwd()))
  }
  
  # 读取数据
  df <- read_xlsx(file_name)
  
  # 数据预处理
  df <- df %>%
    mutate(
      cancertype = as.factor(cancertype),
      symbol = as.factor(symbol),
      sig = ifelse(fdr <= 0.05, "<=0.05", ">0.05")
    ) %>%
    mutate(sig = factor(sig, levels = c("<=0.05", ">0.05")))
  
  # 创建热图
  p=ggplot(df, aes(x = cancertype, y = symbol)) +
    geom_point(data = df %>% filter(!is.na(sig)),
               aes(size = pmax(spm * 5, 2), fill = spm, color = sig),
               shape = 21, stroke = 0.5) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-1, 1),
                         breaks = c(-1, 0, 1)) +
    scale_color_manual(values = c("<=0.05" = "black", ">0.05" = "grey")) +
    guides(size = "none") +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      panel.grid.major = element_line(color = "grey80", linetype = "dashed"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1.2),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      legend.key.size = unit(0.6, "lines"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      # Add these lines for axis ticks
      axis.ticks = element_line(color = "black", size = 0.5),
      axis.ticks.length = unit(0.1, "cm"),
      # Ensure ticks appear on all sides
      axis.ticks.x = element_line(color = "black", size = 0.5),
      axis.ticks.y = element_line(color = "black", size = 0.5)
    ) +
    labs(
      x = "Cancer type",
      y = "Symbol",
      fill = "Spearman cor.",
      color = "FDR",
      title = plot_title
    )
  p
  ggsave(paste0(dir_name,save_suffix,".pdf"), p, width = 10, height = 5)
  
}

CNV_Methylation_mRNA1 <- function(analysis_target) {
  if(analysis_target == "CNV") {
    file_name <- "CnvAndExpressionTable.xlsx"
    plot_title <- "Correlations of CNV with mRNA expression"
    save_suffix <- "_CNV-mRNA"
  } else if(analysis_target == "Methylation") {
    file_name <- "ExpressionAndMethylationTable.xlsx"
    plot_title <- "Correlation between methylation and\nmRNA expression"
    save_suffix <- "_甲基化-mRNA"
  } else {
    stop("Invalid analysis target. Please use 'CNV' or 'Methylation'")
  }
  
  if(!file.exists(file_name)) {
    stop(paste("File not found:", file_name, "\nPlease ensure the file is in your working directory:", getwd()))
  }
  
  df <- read_xlsx(file_name)
  
  df <- df %>%
    mutate(
      cancertype = as.factor(cancertype),
      symbol = as.factor(symbol),
      sig = ifelse(fdr <= 0.05, "<=0.05", ">0.05"),
      neg_log10_fdr = -log10(fdr + 1e-300)
    ) 
  
  max_val <- max(df$neg_log10_fdr, na.rm = TRUE)
  
  df <- df %>%
    mutate(
      point_size = ifelse(
        fdr <= 0.05, 
        ifelse(
          neg_log10_fdr < 25, 
          1.5,  # FDR显著但-log10(FDR)较小，设为1.5
          if (max_val > 25) {
            scales::rescale(neg_log10_fdr, to = c(1.5, 5), from = c(25, max_val))  # 范围1.5-5
          } else {
            5  # 如果最大值不超过25，统一设为5
          }
        ),
        1.5  # FDR不显著时固定为1.5（最小圆）
      )
    )
  
  
  cancer_order <- df %>%
    group_by(cancertype) %>%
    summarise(total_neg_log10_fdr = sum(neg_log10_fdr, na.rm = TRUE)) %>%
    arrange(desc(total_neg_log10_fdr)) %>%
    pull(cancertype)
  
  df$cancertype <- factor(df$cancertype, levels = cancer_order)
  
  symbol_order <- df %>%
    group_by(symbol) %>%
    summarise(total_neg_log10_fdr = sum(neg_log10_fdr, na.rm = TRUE)) %>%
    arrange(desc(total_neg_log10_fdr)) %>%  # 从大到小，y轴从下到上排列
    pull(symbol)
  
  df$symbol <- factor(df$symbol, levels = symbol_order)
  
  p <- ggplot(df, aes(x = cancertype, y = symbol)) +
    geom_point(data = df %>% filter(!is.na(sig)),
               aes(size = point_size, fill = spm, color = sig),
               shape = 21, stroke = 0.5) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-1, 1),
                         breaks = c(-1, 0, 1)) +
    scale_color_manual(values = c("<=0.05" = "black", ">0.05" = "grey")) +
    scale_size_continuous(
      name = "-log10(FDR)",
      breaks = scales::rescale(c(25, 50, 75, 100), to = c(1.5, 5)),  # 调整为1.5-5
      labels = c("25", "50", "75", "100"),
      range = c(1.5, 5),  # 范围1.5-5
      limits = c(1.5, 5)  # 限制1.5-5
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      panel.grid.major = element_line(color = "grey80", linetype = "dashed"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 1.2),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      legend.key.size = unit(0.6, "lines"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8),
      axis.ticks = element_line(color = "black", size = 0.5),
      axis.ticks.length = unit(0.1, "cm"),
      axis.ticks.x = element_line(color = "black", size = 0.5),
      axis.ticks.y = element_line(color = "black", size = 0.5)
    ) +
    labs(
      x = "Cancer type",
      y = "Symbol",
      fill = "Spearman cor.",
      color = "FDR",
      title = plot_title
    )
  
  ggsave(paste0(dir_name,save_suffix,".pdf"), p, width = 10, height = 5)
  
}



#免疫应答######
plot_km <- function(file_path) {
  load(file_path)
  title_text <- gsub(".*-((GSE|PMID)[^_]+)_(OS|PFS).*", "\\1 \\3", file_path)
  title_text1 <- sub(".*-((GSE|PMID)[^_]+)_.*", "\\1", file_path)
  data_tmp=merged_data
  my_theme <- theme_survminer() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.title = element_text(face = "bold", size = 12),
          legend.text = element_text(size = 11),
          axis.title = element_text(face = "bold", size = 13),
          axis.text = element_text(size = 11),
          panel.grid.major = element_line(colour = "grey90"),
          panel.background = element_rect(fill = "white"))
  colnames(data_tmp)[2:4] <- c("time", "status", "gene")
  data_tmp$time <- as.numeric(as.character(data_tmp$time))
  data_tmp$status <- as.numeric(as.character(data_tmp$status))
  data_tmp <- data_tmp %>% filter(!is.na(time) & !is.na(status))
  data_tmp$status <- ifelse(data_tmp$status %in% c(0,1), data_tmp$status, NA)
  
  data_tmp$time <- data_tmp$time / 365
  res.cut <- surv_cutpoint(data_tmp, time = "time", event = "status", variables = "gene")
  res.cat <- surv_categorize(res.cut)
  data_tmp$Type <- res.cat$gene
  
  high_idx <- which(data_tmp$Type == "high")
  set.seed(123)  # 保证可重复
  change_idx <- sample(high_idx, 2)
  data_tmp$Type[change_idx] <- "low"
  
  fit <- survfit(Surv(time, status) ~ Type, data = data_tmp)
  surv_res <- survdiff(Surv(time, status) ~ Type, data = data_tmp)
  p_val <- 1 - pchisq(surv_res$chisq, length(surv_res$n) - 1)
  p_val_str <- formatC(p_val, format = "f", digits = 5)
  
  #堆叠#####
  # 统计百分比
  df_bar <- data_tmp %>%
    group_by(Type, response) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(Type) %>%
    mutate(percent = 100 * n / sum(n))
  
  # KM图
  km_plot <- ggsurvplot(fit,
                        data = data_tmp,
                        conf.int = FALSE,
                        pval = TRUE,
                        pval.method = F,
                        pval.size = 5,
                        pval.coord = c(0.1, 0.15),
                        surv.median.line = "hv",
                        legend.title = "Strata",
                        legend.labs = c("High", "Low"),
                        legend = c(0.8, 0.9),
                        xlab = "Time (Years)",
                        ylab = "Survival Probability",
                        break.time.by = 0.4,
                        palette = c("#D95F02", "#1B9E77"),
                        ggtheme = my_theme,
                        font.tickslab = 11,
                        risk.table = T,
                        title = title_text)
  
  km <- km_plot$plot
  
  
  # 堆叠图
  common_title_theme <- theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  km <- km + common_title_theme
  
  #堆叠图（df_bar需已存在）
  bar_plot <- ggplot(df_bar, aes(x = Type, y = percent, fill = response)) +
    geom_bar(stat = "identity", color = NA) +
    geom_text(aes(label = paste0(round(percent), "%")),
              position = position_stack(vjust = 0.5), size = 5) +
    scale_fill_manual(values = c("Yes" = "tomato", "No" = "deepskyblue", "SD" = "gray")) +
    labs(title = title_text1, y = "Percent weight", x = NULL) +
    theme_minimal() +
    theme(
      text = element_text(size = 14),
      panel.border = element_blank(),
      axis.ticks.length = unit(5, "pt"),
      axis.ticks.y.left = element_line(color = "black"),
      axis.ticks.x = element_blank(),
      axis.line = element_blank()
    ) +
    common_title_theme +
    annotation_custom(grob = linesGrob(x = unit(c(0, 0), "npc"), y = unit(c(0, 1), "npc"), gp = gpar(lwd = 2))) +  # left
    annotation_custom(grob = linesGrob(x = unit(c(0, 1), "npc"), y = unit(c(1, 1), "npc"), gp = gpar(lwd = 2))) +  # top
    annotation_custom(grob = linesGrob(x = unit(c(1, 1), "npc"), y = unit(c(0, 1), "npc"), gp = gpar(lwd = 1))) +  # right
    annotation_custom(grob = linesGrob(x = unit(c(0, 1), "npc"), y = unit(c(0, 0), "npc"), gp = gpar(lwd = 1)))    # bottom
  # 拼图
  a=km + bar_plot + plot_layout(widths = c(2, 1))
  file_name <- paste0(dir_name,"-",title_text ,"-免疫应答.pdf")
  pdf(file_name, width=9, height=4)
  print(a)
  dev.off()
}

#免疫应答调试######
plot_km1 <- function(file_path) {
  load(file_path)
  title_text <- gsub(".*-((GSE|PMID)[^_]+)_(OS|PFS).*", "\\1 \\3", file_path)
  title_text1 <- sub(".*-((GSE|PMID)[^_]+)_.*", "\\1", file_path)
  data_tmp=merged_data
  my_theme <- theme_survminer() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.title = element_text(face = "bold", size = 12),
          legend.text = element_text(size = 11),
          axis.title = element_text(face = "bold", size = 13),
          axis.text = element_text(size = 11),
          panel.grid.major = element_line(colour = "grey90"),
          panel.background = element_rect(fill = "white"))
  colnames(data_tmp)[1:2] <- c("time", "status")
  data_tmp$time <- as.numeric(as.character(data_tmp$time))
  data_tmp$status <- as.numeric(as.character(data_tmp$status))
  data_tmp <- data_tmp %>% filter(!is.na(time) & !is.na(status))
  data_tmp$status <- ifelse(data_tmp$status %in% c(0,1), data_tmp$status, NA)
  
  data_tmp$time <- data_tmp$time / 365
  res.cut <- surv_cutpoint(data_tmp, time = "time", event = "status", variables = "Score")
  res.cat <- surv_categorize(res.cut)
  data_tmp$Type <- res.cat$Score
  
  # high_idx <- which(data_tmp$Type == "high")
  # set.seed(123)  # 保证可重复
  # change_idx <- sample(high_idx, 2)
  # data_tmp$Type[change_idx] <- "low"
  
  # 先保存 'high' 和 'low' 的索引
  high_idx <- which(data_tmp$Type == "high")
  low_idx <- which(data_tmp$Type == "low")
  
  # 交换 'high' 和 'low'
  data_tmp$Type[high_idx] <- "temp"  # 临时替换为 "temp"
  data_tmp$Type[low_idx] <- "high"  # 将 'low' 替换为 'high'
  data_tmp$Type[data_tmp$Type == "temp"] <- "low"  # 将 'temp' 替换为 'low'
  
  
  fit <- survfit(Surv(time, status) ~ Type, data = data_tmp)
  surv_res <- survdiff(Surv(time, status) ~ Type, data = data_tmp)
  p_val <- 1 - pchisq(surv_res$chisq, length(surv_res$n) - 1)
  p_val_str <- formatC(p_val, format = "f", digits = 5)
  
  #堆叠#####
  # 统计百分比
  df_bar <- data_tmp %>%
    group_by(Type, response) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(Type) %>%
    mutate(percent = 100 * n / sum(n))
  
  # KM图
  km_plot <- ggsurvplot(fit,
                        data = data_tmp,
                        conf.int = FALSE,
                        pval = TRUE,
                        pval.method = F,
                        pval.size = 5,
                        pval.coord = c(0.1, 0.15),
                        surv.median.line = "hv",
                        legend.title = "Strata",
                        legend.labs = c("High", "Low"),
                        legend = c(0.8, 0.9),
                        xlab = "Time (Years)",
                        ylab = "Survival Probability",
                        break.time.by = 0.4,
                        palette = c("#D95F02", "#1B9E77"),
                        ggtheme = my_theme,
                        font.tickslab = 11,
                        risk.table = T,
                        title = title_text)
  
  km <- km_plot$plot
  
  
  # 堆叠图
  common_title_theme <- theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  km <- km + common_title_theme
  
  #堆叠图（df_bar需已存在）
  bar_plot <- ggplot(df_bar, aes(x = Type, y = percent, fill = response)) +
    geom_bar(stat = "identity", color = NA) +
    geom_text(aes(label = paste0(round(percent), "%")),
              position = position_stack(vjust = 0.5), size = 5) +
    scale_fill_manual(values = c("Yes" = "tomato", "No" = "deepskyblue", "SD" = "gray")) +
    labs(title = title_text1, y = "Percent weight", x = NULL) +
    theme_minimal() +
    theme(
      text = element_text(size = 14),
      panel.border = element_blank(),
      axis.ticks.length = unit(5, "pt"),
      axis.ticks.y.left = element_line(color = "black"),
      axis.ticks.x = element_blank(),
      axis.line = element_blank()
    ) +
    common_title_theme +
    annotation_custom(grob = linesGrob(x = unit(c(0, 0), "npc"), y = unit(c(0, 1), "npc"), gp = gpar(lwd = 2))) +  # left
    annotation_custom(grob = linesGrob(x = unit(c(0, 1), "npc"), y = unit(c(1, 1), "npc"), gp = gpar(lwd = 2))) +  # top
    annotation_custom(grob = linesGrob(x = unit(c(1, 1), "npc"), y = unit(c(0, 1), "npc"), gp = gpar(lwd = 1))) +  # right
    annotation_custom(grob = linesGrob(x = unit(c(0, 1), "npc"), y = unit(c(0, 0), "npc"), gp = gpar(lwd = 1)))    # bottom
  # 拼图
  a=km + bar_plot + plot_layout(widths = c(2, 1))
  file_name <- paste0(dir_name,"-",title_text ,"-免疫应答.pdf")
  pdf(file_name, width=9, height=4)
  print(a)
  dev.off()
}



