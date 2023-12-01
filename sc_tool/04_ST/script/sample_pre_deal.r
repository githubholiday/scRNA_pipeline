print("按照Normalize->FindVariable->ScaleData->RunPCA->FindNeighbors->FindClusters->DimPlot等，主要是前面的步骤，将rds文件进行初步处理，再进行后面的绘图等")

library('getopt')
para<- matrix(c(
        'help',         'h',    0,      "logical",
        'rds',          'r',    1,      "character",
        'prefix',       'p',    1,      "character",
        'outdir',       'o',    1,      "character"
),byrow=TRUE,ncol=4)
#===========================================================
opt <- getopt(para,debug=FALSE)

library(Seurat)
rds <- readRDS( opt$rds)
rds <- NormalizeData(rds, normalization.method = "LogNormalize", scale.factor = 10000)
rds <- FindVariableFeatures(rds, selection.method = "vst", nfeatures = 20)
rds <- ScaleData(rds)
rds <- RunPCA(rds,features= VariableFeatures(object = rds))
rds <- RunUMAP(rds,, reduction = "pca",dims = 1:10)
rds <- FindNeighbors(rds, dims = 1:10)
rds <- FindClusters(rds, resolution =0.6)

out_pre = paste(opt$outdir, opt$prefix, sep="/")
out_pdf = paste(out_pre,".dimplpt.pdf",sep="")
pdf(out_pdf)
DimPlot(rds, group.by="cluster_standard",label=F)
dev.off()

out_pre = paste(opt$outdir, opt$prefix, sep="/")
saveRDS(rds, paste(out_pre, ".rds", sep=""))