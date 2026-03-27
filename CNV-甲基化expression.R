
library(ggplot2)
library(readxl)
library(dplyr)
dir_name <- basename(getwd())
# 设置分析目标
analysis_target <- "CNV"
analysis_target <- "Methylation"
CNV-Methylation-mRNA <- function(analysis_target) {
# 根据目标选择文件并设置标题
if(analysis_target == "CNV") {
  file_name <- "CnvAndExpressionTable.xlsx"
  plot_title <- "Correlations of CNV with mRNA expression"
  save_suffix <- "_CNV-mRNA.pdf"
} else if(analysis_target == "Methylation") {
  file_name <- "ExpressionAndMethylationTable.xlsx"
  plot_title <- "Correlation between methylation and\nmRNA expression"
  save_suffix <- "_甲基化-mRNA.pdf"  # 修改为"-mRNA"保持一致性
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

ssGsea_Score(analysis_target)






library(ggplot2)
library(readxl)
library(dplyr)
library(scales)
analysis_target <- "CNV"
analysis_target <- "Methylation"

if(analysis_target == "CNV") {
  file_name <- "CnvAndExpressionTable.xlsx"
  plot_title <- "Correlations of CNV with mRNA expression"
  save_suffix <- "_CNV-mRNA.pdf"
} else if(analysis_target == "Methylation") {
  file_name <- "ExpressionAndMethylationTable.xlsx"
  plot_title <- "Correlation between methylation and\nmRNA expression"
  save_suffix <- "_甲基化-mRNA.pdf"
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
    sig = factor(sig, levels = c("<=0.05", ">0.05")),
    point_size = ifelse(
      neg_log10_fdr < 25, 3,                 # 小点改大一点
      if (max_val > 25) {
        scales::rescale(neg_log10_fdr, to = c(3, 6), from = c(25, max_val))  # 大点范围扩大
      } else {
        6
      }
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
    breaks = scales::rescale(c(25, 50, 75, 100), to = c(3, 6)),
    labels = c("25", "50", "75", "100"),
    range = c(3, 6),
    limits = c(3, 6)
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

print(p)









