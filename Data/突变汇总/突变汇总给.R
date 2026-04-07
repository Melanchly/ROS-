
setwd("D:/Dim/TCGAplot/TRP")
rm(list = ls())
alter <- read.delim("alterations_across_samples_TRP.tsv", header = TRUE, sep = "\t", check.names = FALSE, quote = "")
alter <- alter[, -3]
alter=as.data.frame(alter)
alter$Cancer <- toupper(sub("_.*", "", alter[[1]]))
alter <- alter[, -1]  # 删除第一列
alter <- alter[, -1]  # 删除第一列
alter <- alter[, c(ncol(alter), 1:(ncol(alter)-1))]  # 把最后一列移到最前
alter <- alter[, -c(3:28)]
alter <- alter[alter[, 2] == 1, ]

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(patchwork)
library(tibble)
# 假设你的原始数据是 alter，先整理成 long format
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

# 每列总突变数
total_cancer <- colSums(mut_count)
mut_count <- rbind(mut_count, Total = total_cancer)

# 设置颜色
col_fun <- colorRamp2(c(0, max(mut_count)), c("white", "red"))


col_fun <- Vectorize(function(x) {
  
  if (x < 10) "white"
  else if (x < 20) "#ded3cf"
  else if (x < 40) "#669dbb"
  else if (x < 60) "#fdb1a4"
  else "#99253a"
})



p = main_ht <- Heatmap(mut_count,
                       col = col_fun,
                       cluster_rows = FALSE,
                       cluster_columns = FALSE,
                       cell_fun = function(j, i, x, y, width, height, fill) {
                         grid.text(mut_count[i, j], x, y, gp = gpar(fontsize = 8))
                         grid.rect(x, y, width, height, gp = gpar(col = "grey80", fill = NA, lwd = 0.5))
                       },
                       show_column_names = TRUE,
                       show_row_names = TRUE,
                       row_names_side = "left",
                       row_names_gp = gpar(fontsize = 8),
                       column_names_gp = gpar(fontsize = 8),
                       show_heatmap_legend = FALSE
)

print(p)












