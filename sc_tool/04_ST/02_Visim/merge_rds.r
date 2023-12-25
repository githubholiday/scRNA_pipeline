print("将各个细胞类型的rds文件合并，并将其进行标准化处理，用于后续空间数据的注释")
library(Seurat)
B <- readRDS("/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/02_Spatial/4_RNA_ST_2/result/rna//B.rds")
T_rds <- readRDS("/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/02_Spatial/4_RNA_ST_2/result/rna//T.rds")
ECs.rds <- readRDS("/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/02_Spatial/4_RNA_ST_2/result/rna//ECs.rds")
EpithelialCells.rds <- readRDS("/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/02_Spatial/4_RNA_ST_2/result/rna/EpithelialCells.rds")
MPs.rds <- readRDS("/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/02_Spatial/4_RNA_ST_2/result/rna/MPs.rds")
PanSCs.rds <- readRDS("/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/02_Spatial/4_RNA_ST_2/result/rna/PanSCs.rds")
PlasmaCells.rds <- readRDS("/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/02_Spatial/4_RNA_ST_2/result/rna/PlasmaCells.rds")
combine <- merge(B, y=c(T_rds,ECs.rds,EpithelialCells.rds,MPs.rds,PanSCs.rds,PlasmaCells.rds))

outfile <- "/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/02_Spatial/4_RNA_ST_2/result/rna/merge.rds"
saveRDS( combine, outfile )

rds <- readRDS( outfile )
rds <- NormalizeData(rds, normalization.method = "LogNormalize", scale.factor = 10000)
rds <- FindVariableFeatures(rds, selection.method = "vst", nfeatures = 20)
rds <- ScaleData(rds)
rds <- RunPCA(rds,features= VariableFeatures(object = rds))
rds <- RunUMAP(rds,, reduction = "pca",dims = 1:10)
rds <- FindNeighbors(rds, dims = 1:10)
rds <- FindClusters(rds, resolution =0.6)

out_pdf = "/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/02_Spatial/4_RNA_ST_2/result/rna/rna.dimplpt.pdf"
pdf(out_pdf)
DimPlot(rds, group.by="cluster_standard",label=F)
dev.off()


saveRDS(rds, "/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/02_Spatial/4_RNA_ST_2/result/rna/rna.rds")
