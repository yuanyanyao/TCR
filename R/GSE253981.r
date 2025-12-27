# %% [模块：提取黑质TCR]
library(tidyverse)

tcr_files <- list.files("/home/yuanyanyao/Downloads/PD_TCR_Integration/data/raw/GSE253981/", pattern = "\\.filtJS\\.tsv\\.gz$", full.names = TRUE)

# 定义一个读取函数
read_sn_tcr <- function(file) {
  sample_id <- basename(file) %>% str_remove(".filtJS.tsv.gz")
  df <- read_tsv(file, show_col_types = FALSE) %>%
    select(cdr3, v_gene, j_gene, count) %>%
    mutate(sample = sample_id, tissue = "SN")
  return(df)
}

# 合并所有黑质克隆
sn_tcr_master <- map_df(tcr_files, read_sn_tcr)
# 去重，保留唯一序列清单用于比对
sn_unique_clones <- sn_tcr_master %>% distinct(cdr3)