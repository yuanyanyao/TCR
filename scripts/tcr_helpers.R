if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# 安装核心框架
BiocManager::install(c("Seurat", "scRepertoire", "immunarch", "Harmony"))

# 安装可视化增强包
install.packages(c("ggalluvial", "treemap", "pheatmap"))