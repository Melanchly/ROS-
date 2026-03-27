
setwd("D:/Dim/TCGAplot/免疫治疗队列")
library(GSEABase)
library(GSVA)
library(matrixStats)

rm(list = ls())
data <- readRDS("Melanoma-GSE78220.Response.Rds")
clin <- readRDS("Melanoma-GSE78220.Response  clin.Rds")

data <- readRDS("RCC-Braun_2020.Response.Rds")
clin <- readRDS("RCC-Braun_2020.Response clin.Rds")

data <- readRDS("Melanoma-GSE91061.Response.Rds")
clin <- readRDS("Melanoma-GSE91061.Response clin.Rds")

data <- readRDS("NSCLC_GSE135222.Response.Rds")
clin <- readRDS("NSCLC_GSE135222.Response clin.Rds")



cox_genes1 <- c(
  "COX6A1", "COX4I2", "COX7A2", "COX6B1", "COX7C", "COX4I1", "COX7B",
  "COX5B", "COX6A2", "COX6B2", "COX7A1", "COX6C","COX7B2", "COX8A",
  "COX5A", "COX8C", "MT-CO2", "MT-CO1", "MT-CO3"
)

cox_genes2 <- c(
  "GSTO2", "GSTP1", "GSTZ1", "GSTT2B", "GSTM1", "GSTM5", "GSTM3",
  "GSTO1", "HPGDS", "GSTM4", "GSTA4", "GSTA3", "GSTA5", "GSTK1",
  "GSTM2", "GSTA1", "GSTA2", "GSTT4"
)


cox_genes31 <- c(
  "ATP5F1D", "ATP5F1B", "ATP5PB", "ATP5F1E", "ATP5IF1", "ATP5MC2",
  "ATP5F1A", "ATP5MC3", "ATP5PF", "ATP5MC1", "ATP5F1C", "ATP5MG",
  "ATP5PD", "ATP5ME", "MT-ATP6", "MT-ATP8", "ATP5MF", "ATP5PO"
)



cox_genes4 <- c(
  "TRPM4", "TRPM2", "TRPM8", "TRPC4", "TRPC3", "TRPV3",
  "MLKL", "TRPC4AP", "TRPC6", "TRPC1", "TRPV5", "TRPM1",
  "TRPV2", "TRPM5", "MCOLN1", "TRPM7", "TRPV4", "RIPK3",
  "RIPK1", "TRPV6", "MCOLN3", "TRPA1", "MCOLN2", "TRPM3",
  "TRPM6", "TRPV1"
)
cox_genes1 <- c("MAPK1", "MAPK3", "MAPK4", "MAPK6", "MAPK7", "MAPK8",
               "MAPK9", "MAPK10", "MAPK11", "MAPK12", "MAPK13", "MAPK14",
               "MAPK15", "NLK"
)



cox_genes1 = c(
  "UQCRB","UQCRQ","UQCRC1","UQCRC2","MT-CYB",
  "CYC1","UQCRFS1","UQCRH","UQCR10","UQCR11"
)
cox_genes2 = c(
  "NDUFAB1","NDUFA1","NDUFA2","NDUFA3","NDUFA5","NDUFA6","NDUFA7",
  "NDUFA8","NDUFA9","NDUFA10","NDUFA11","NDUFA12","NDUFA13",
  "NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7",
  "NDUFB8","NDUFB9","NDUFB10","NDUFB11",
  "NDUFC1","NDUFC2",
  "NDUFS4","NDUFS5","NDUFS6","NDUFV3"
)



cox_list <- list(cox_genes1, cox_genes21, cox_genes41, cox_genes31)
cox_list <- list(cox_genes2)

abcssGsea_Score <- function(data, cox_genes, clin) {
  
  # 数据预处理
  data <- as.data.frame(data)
  data <- data[!is.na(data[, 1]), , drop = FALSE]
  rownames(data) <- data[, 1]        # 将第一列设为行名
  data <- data[, -1, drop = FALSE]   # 删除第一列
  
  # FPKM转TPM函数
  fpkmToTpm <- function(data) {
    exp(log(data) - log(sum(data, na.rm = TRUE)) + log(1e6))
  }
  
  # 应用转换并处理异常值
  eset_tpm <- apply(data, 2, function(x) {
    x[is.na(x)] <- 0
    x[is.infinite(x)] <- 0
    if(sum(x, na.rm = TRUE) == 0) return(x)
    fpkmToTpm(x)
  })
  
  # 检查数据有效性
  if (all(eset_tpm == 0) | all(is.na(eset_tpm))) {
    warning("转换后的TPM数据全为0或NA，跳过该基因集")
    return(NULL)
  }
  
  # 对数转换（如果需要）
  tpm_max <- max(eset_tpm, na.rm = TRUE)
  tpm_min <- min(eset_tpm, na.rm = TRUE)
  
  cat(paste("TPM range:", tpm_min, "-", tpm_max, "\n"))
  
  if (tpm_max > 100 && tpm_min >= 0) {
    eset_tpm <- log2(eset_tpm + 1)
    cat("Applied log2(TPM + 1) transformation\n")
  }
  
  data <- round(eset_tpm, 2)
  data <- as.data.frame(data)
  
  # 准备表达矩阵
  expr_data <- as.matrix(data)
  
  # 自动获取前缀用于命名
  prefix <- toupper(substr(cox_genes[1], 1, 3))
  
  # 创建基因集列表
  gene_sets <- list(cox_genes)
  names(gene_sets) <- paste0(prefix, "_family")
  
  # 运行ssGSEA
  tryCatch({
    param <- ssgseaParam(
      exprData = expr_data,
      geneSets = gene_sets,
      minSize = 5,        # 最小基因集大小
      maxSize = Inf,      # 最大基因集大小
      alpha = 0.25,       # ssGSEA权重参数
      normalize = TRUE    # 标准化分数
    )
    
    ssgsea_scores <- gsva(param)
    ssgsea_scores <- as.data.frame(t(ssgsea_scores))
    colnames(ssgsea_scores)[1] <- "Score"
    
    # 临床数据处理
    clin_df <- as.data.frame(clin)
    rownames(clin_df) <- clin_df[, 1]
    clin_df <- clin_df[, -1, drop = FALSE]
    
    # 检查并匹配临床数据列名
    pfs_time_col <- "overall survival (days)"
    pfs_event_col <- "vital status"
    response_col <- "response"  # 新增：response列名
    
    required_cols <- c(pfs_time_col, pfs_event_col, response_col)
    
    if (all(required_cols %in% colnames(clin_df))) {
      # 提取需要的列：生存时间、生存状态和response
      pfs_data <- clin_df[, required_cols, drop = FALSE]
      
      # 重命名列
      colnames(pfs_data) <- c("OS.time", "OS", "response")
      
      # 转换生存状态
      pfs_data$OS <- ifelse(toupper(pfs_data$OS) %in% c("ALIVE", "LIVING", "0", "FALSE"), 0, 1)
      
      # 转换response列：CR/PR替换为Yes，其他替换为No
      pfs_data$response <- ifelse(toupper(trimws(pfs_data$response)) %in% c("CR", "PR", "R"), "Yes", "No")
      cat(paste("Response distribution:", 
                sum(pfs_data$response == "Yes"), "Yes,", 
                sum(pfs_data$response == "No"), "No\n"))
      
      # 匹配样本
      common_rows <- intersect(rownames(pfs_data), rownames(ssgsea_scores))
      
      if (length(common_rows) > 0) {
        merged_data <- cbind(pfs_data[common_rows, , drop = FALSE], 
                             ssgsea_scores[common_rows, , drop = FALSE])
        
        # 保存结果
        output_file <- paste0(prefix, "-PMID32472114_OS_ssGSEA_scores.Rdata")
        save(merged_data, file = output_file)
        cat(paste("Successfully saved results to", output_file, "\n"))
        
        return(merged_data)
      } else {
        warning("No common samples between expression and clinical data")
        return(NULL)
      }
    } else {
      warning(paste("Required clinical columns not found. Available columns:", 
                    paste(colnames(clin_df), collapse = ", ")))
      return(NULL)
    }
    
  }, error = function(e) {
    message(paste("Error in ssGSEA for", prefix, ":", e$message))
    return(NULL)
  })
}



for (cox_genes in cox_list) {
  abcssGsea_Score(data, cox_genes,clin)
}






library(limma)
library(survival)
library(survminer)
library(beepr)
library(dplyr)
library(tidyr)
rm(list = ls())
data_tmp=OS
my_theme <- theme_survminer() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 11),
        axis.title = element_text(face = "bold", size = 13),
        axis.text = element_text(size = 11),
        panel.grid.major = element_line(colour = "grey90"),
        panel.background = element_rect(fill = "white"))

colnames(data_tmp)[1:3] <- c("time", "status", "gene")
data_tmp$time <- as.numeric(as.character(data_tmp$time))
data_tmp$status <- as.numeric(as.character(data_tmp$status))
data_tmp <- data_tmp %>% filter(!is.na(time) & !is.na(status))
data_tmp$status <- ifelse(data_tmp$status %in% c(0,1), data_tmp$status, NA)

data_tmp$time <- data_tmp$time / 365

res.cut <- surv_cutpoint(data_tmp, time = "time", event = "status", variables = "gene")
res.cat <- surv_categorize(res.cut)
data_tmp$Type <- res.cat$gene

fit <- survfit(Surv(time, status) ~ Type, data = data_tmp)
surv_res <- survdiff(Surv(time, status) ~ Type, data = data_tmp)
p_val <- 1 - pchisq(surv_res$chisq, length(surv_res$n) - 1)
p_val_str <- formatC(p_val, format = "f", digits = 5)

surPlot <- ggsurvplot(fit,
                      data = data_tmp,
                      conf.int = FALSE,
                      pval = TRUE,
                      pval.method = TRUE,
                      pval.size = 5,
                      pval.coord = c(0.1, 0.15),
                      surv.median.line = "hv",
                      legend.title = "Gene Expression",
                      legend.labs = c("High", "Low"),
                      legend = c(0.8, 0.9),
                      xlab = "Time (Years)",
                      ylab = "Overall Survival Probability",
                      break.time.by = 1,
                      palette = c("#D95F02", "#1B9E77"),
                      ggtheme = my_theme,
                      font.tickslab = 11,
                      risk.table = T)
print(surPlot)








#GSE135222######
rm(list = ls())
library(GEOquery)
library(GSVA)
library(GSEABase)
gset <- getGEO('GSE135222', destdir = '.', AnnotGPL = FALSE, getGPL = FALSE)
pheno <- pData(gset[[1]])
library(data.table)
# 解压并读取
expr <- fread("GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz", data.table = FALSE)
library(tidyverse)
library(stringr)
expr=read_tsv('GSE135222_GEO_RNA-seq_omicslab_exp.tsv')
gtf_v22 <- read_tsv(file = "gencode.gene.info.v22.tsv") 
gtf_v22$gene_id=unlist(lapply(gtf_v22$gene_id, function(x){  strsplit(x,'[.]')[[1]][1]}))
table(duplicated(gtf_v22$gene_id))
expr$gene_id=unlist(lapply(expr$gene_id,function(x){  strsplit(x,'[.]')[[1]][1]}))
expr <- inner_join(gtf_v22,expr, by = "gene_id") %>% dplyr::select(-1)
expr=expr[,-c(2:11)]
table(duplicated(expr$gene_name))
expr <- aggregate(.~ gene_name, expr, mean)###对基因名去重取均值
expr <- expr %>% column_to_rownames("gene_name")
fpkmToTpm <- function(fpkm)##FPKM转为TPM
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
expr <- apply(expr, 2, fpkmToTpm)
expr=log2(expr+1)
expr <- round(expr, 2)
cox_genes <- c(
  "TRPM4", "TRPM2", "TRPM8", "TRPC4", "TRPC3", "TRPV3",
  "MLKL", "TRPC4AP", "TRPC6", "TRPC1", "TRPV5", "TRPM1",
  "TRPV2", "TRPM5", "MCOLN1", "TRPM7", "TRPV4", "RIPK3",
  "RIPK1", "TRPV6", "MCOLN3", "TRPA1", "MCOLN2", "TRPM3",
  "TRPM6", "TRPV1"
)
cox_genes <- c(
  "COX6A1", "COX4I2", "COX7A2", "COX6B1", "COX7C", "COX4I1", "COX7B",
  "COX5B", "COX6A2", "COX6B2", "COX7A1", "COX6C","COX7B2", "COX8A",
  "COX5A", "COX8C", "MT-CO2", "MT-CO1", "MT-CO3"
)

cox_genes <- c(
  "GSTO2", "GSTP1", "GSTZ1", "GSTT2B", "GSTM1", "GSTM5", "GSTM3",
  "GSTO1", "HPGDS", "GSTM4", "GSTA4", "GSTA3", "GSTA5", "GSTK1",
  "GSTM2", "GSTA1", "GSTA2", "GSTT4"
)


cox_genes <- c(
  "ATP5F1D", "ATP5F1B", "ATP5PB", "ATP5F1E", "ATP5IF1", "ATP5MC2",
  "ATP5F1A", "ATP5MC3", "ATP5PF", "ATP5MC1", "ATP5F1C", "ATP5MG",
  "ATP5PD", "ATP5ME", "MT-ATP6", "MT-ATP8", "ATP5MF", "ATP5PO"
)

expr_data <- as.matrix(expr)

# 自动获取前缀用于命名
prefix <- toupper(substr(cox_genes[1], 1, 3))
# 创建 GeneSet 对象
cox_gene_set <- GeneSet(cox_genes, setName = paste0(prefix, "_family"))

# 运行 ssGSEA
param <- ssgseaParam(
  exprData = expr_data,
  geneSets = setNames(list(cox_genes), paste0(prefix, "_family")),
  minSize = 5,
  maxSize = Inf,
  alpha = 0.25,
  normalize = TRUE
)
ssgsea_scores <- gsva(param)
ssgsea_scores <- as.data.frame(t(ssgsea_scores))
colnames(ssgsea_scores)[1] <- "Score"
# Subset the pheno data (keeping columns 1, 43, 44)
pheno <- pheno[, c(1, 43, 44)]

# 删除 pheno 第一列中的所有空格
pheno[,1] <- gsub("\\s+", "", pheno[,1])

rownames(ssgsea_scores) <- rownames(pheno)[match(rownames(ssgsea_scores), pheno[[1]])]
pheno <- pheno[, -1, drop = FALSE]
# Rename the columns
colnames(pheno) <- c("OS.time", "OS")
# 检查样本名是否一致
# 去除 pheno 行名中的所有空格
common_samples <- intersect(rownames(pheno), rownames(ssgsea_scores))
Score <- ssgsea_scores[common_samples, ]       # 筛选共有的样本
pheno_matched <- pheno[common_samples, ]   # 同样筛选
# 合并（cbind 按行名自动对齐）
merged_data <- cbind(pheno_matched, Score)
clin=as.data.frame(clin)
rownames(clin) <- clin[,1]
clin <- clin[,-1, drop = FALSE]
clin <- clin[, 5, drop = FALSE]  # 保持为单列数据框
merged_data <- cbind(clin, merged_data)
save(merged_data,file = "GST-GSE135222.Rdata")









library(limma)
library(survival)
library(survminer)
library(beepr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
load("D:/Dim/TCGAplot/免疫治疗队列/MAP-GSE91061_OS_ssGSEA_scores.Rdata")
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

# 画图
ggplot(df_bar, aes(x = Type, y = percent, fill = response)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(percent), "%")), 
            position = position_stack(vjust = 0.5), size = 5) +
  scale_fill_manual(values = c("PD" = "tomato", "PR" = "deepskyblue", "SD" = "gray")) +
  labs(title = "GSE135222", y = "Percent weight", x = NULL) +
  theme_minimal() +
  theme(text = element_text(size = 14))


# KM图
km_plot <- ggsurvplot(fit,
                      data = data_tmp,
                      conf.int = FALSE,
                      pval = TRUE,
                      pval.method = F,
                      pval.size = 5,
                      pval.coord = c(0.1, 0.15),
                      surv.median.line = "hv",
                      legend.title = "Strta",
                      legend.labs = c("High", "Low"),
                      legend = c(0.8, 0.9),
                      xlab = "Time (Years)",
                      ylab = "Overall Survival Probability",
                      break.time.by = 1,
                      palette = c("#D95F02", "#1B9E77"),
                      ggtheme = my_theme,
                      font.tickslab = 11,
                      risk.table = T,
                      title = "Kaplan-Meier Survival Curve")

km <- km_plot$plot


# 堆叠图
bar_plot <- ggplot(df_bar, aes(x = Type, y = percent, fill = response)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(percent), "%")), 
            position = position_stack(vjust = 0.5), size = 5) +
  scale_fill_manual(values = c("PD" = "tomato", "PR" = "deepskyblue", "SD" = "gray")) +
  labs(title = "GSE135222", y = "Percent weight", x = NULL) +
  theme_minimal() +
  theme(text = element_text(size = 14))

# 拼图
km + bar_plot + plot_layout(widths = c(2, 1))






