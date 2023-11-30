print("该脚本用于将subset中的新的注释细胞类型替换到原rds中")

library(Seurat)
library('getopt')
para<-matrix(c(
    'help','h',0,"logical",
    'newrds','n',1,"character",
    'rawrds','r',1,"character",
    'outfile','o',1,"character"
),byrow=TRUE,ncol=4)
opt <- getopt(para,debug=F)


raw_rds<-readRDS(opt$rawrds)
new_rds<-readRDS(opt$newrds)

#将原rds文件中的labels复制一份到new_labels中
raw_rds@meta.data$new_labels <- raw_rds@meta.data$labels

#
sub_barcode <- rownames(new_rds@meta.data)
sub_celltype <- new_rds@meta.data$predicted.celltype.l2

#将新的注释细胞类型替换到原rds中
raw_rds@meta.data[sub_barcode, "new_labels"]  <- sub_celltype

pdf(opt$outfile,w=14,h=6)
DimPlot(raw_rds, reduction = "umap", group.by = "new_labels",label = TRUE, repel = TRUE, label.size = 2.5)
dev.off()