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
#将细胞类型重新命名
#方法1
rr@meta.data[rr@meta.data$cluster_standard=="Conventional dendritic cells",]["cluster_standard"] = "Conventional_dendritic"
#方法2：
rds2 <- subset(sce, subset= SingleR_Prediction=="Skeletal muscle")
a1 <- rownames(rds2@meta.data["SingleR_Prediction"])
sce@meta.data[a1,"SingleR_Prediction"]
