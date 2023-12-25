print("该脚本用于将rds中的细胞类型替换成二级细胞类型")


library(Seurat)

all_rds = "/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/02_Spatial/4_RNA_ST_2/input/PanSCs/P22082903.diff_PRO.rds"
sub_rds = "/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/02_Spatial/4_RNA_ST_2/input/FibroBlasts/P22082903.diff_PRO.rds"
out_rds = "/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/02_Spatial/4_RNA_ST_2/result/PanSCs.rds"

a <- readRDS(all_rds)
b <- readRDS(sub_rds)

aa = data.frame(name=as.character(a$cluster_standard))
rownames(aa) <- rownames(a@meta.data)
bb = data.frame(name=as.character(b$cluster_standard))
rownames(bb) <- rownames(b@meta.data)

for (ia in rownames(aa)) {
  if (ia %in% rownames(bb)) {
    aa$name[rownames(aa) == ia] <- bb$name[rownames(bb) == ia]
  }
}

a$cluster_standard <- as.factor(aa$name)
saveRDS(a, out_rds)

