library(Seurat)
library(tidyverse)
options(future.globals.maxSize = 8000 * 1024^2) # 允许 R 使用更多全局变量空间
# 1. 设置路径
base_path <- "/home/yuanyanyao/Downloads/PD_TCR_Integration/data/raw/Zenodo3993994/"
counts_file <- paste0(base_path, "counts.rds")
norm_file <- paste0(base_path, "normalized.rds")
meta_file <- paste0(base_path, "metadata.txt")

# 2. 读取数据
message("Reading files... This may take a minute.")
counts <- readRDS(counts_file)
norm_data <- readRDS(norm_file) # 保持其为稀疏矩阵格式
meta <- read.table(meta_file, header = TRUE, sep = "\t", row.names = 1)

# 3. 创建 Seurat 对象 (适配 v5)
# 直接用 counts 创建对象
pd_integrated <- CreateSeuratObject(counts = counts, meta.data = meta, project = "PD_TCR")

# 4. 关键修复：将归一化数据存入 Seurat v5 的 layer
# 在 v5 中，我们直接赋值给 layer，不要用 as.matrix()!!
# 如果 norm_data 已经是矩阵格式，Seurat 会自动处理
pd_integrated[["RNA"]]$data <- norm_data 

# 5. 验证结构
message("Object Structure:")
print(pd_integrated)
# %%
unique(pd_integrated@meta.data)
table(pd_integrated$cellType2)
# 6. 提取 T 细胞子集并保存
# 检查一下 metadata 里的 cell_type 列名是否正确


saveRDS(pd_integrated, "/home/yuanyanyao/Downloads/PD_TCR_Integration/data/processed/Zenodo_3993994_Base.rds")
message("Done! Saved to processed folder.")