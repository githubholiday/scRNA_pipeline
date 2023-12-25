rm(list = ls()) library(Seurat) 
# devtools::install_github('satijalab/seurat-data') 
library(SeuratData) AvailableData() 
# InstallData("pbmc3k") # (89.4 MB) data("pbmc3k") 
exprMat <- as.matrix(pbmc3k@assays$RNA@data) 
dim(exprMat) 
exprMat[1:4,1:4] 
cellInfo <- pbmc3k@meta.data[,c(4,2,3)] 
colnames(cellInfo)=c('CellType', 'nGene' ,'nUMI') 
#head(cellInfo) 
table(cellInfo$CellType) 
### Initialize settings library(SCENIC) 
# 保证cisTarget_databases 文件夹下面有下载好2个1G的文件 
scenicOptions <- initializeScenic(org="hgnc", dbDir="cisTarget_databases", nCores=1) 
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
### Co-expression network 
genesKept <- geneFiltering(exprMat, scenicOptions) 
exprMat_filtered <- exprMat[genesKept, ] 
exprMat_filtered[1:4,1:4] 
dim(exprMat_filtered) 
runCorrelation(exprMat_filtered, scenicOptions) 
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions) 
### Build and score the GRN 
exprMat_log <- log2(exprMat+1) 
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] 
# Toy run settings 
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions) 
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) 
# Toy run settings library(doParallel) 
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log ) 
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions) 
tsneAUC(scenicOptions, aucType="AUC") 
# choose settings #运算结果保存方式 
export2loom(scenicOptions, exprMat) 
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
#查看转录因子富集结果 
rm(list = ls()) 
library(Seurat) 
library(SCENIC) 
library(doParallel) 
scenicOptions=readRDS(file="int/scenicOptions.Rds") 
#每个motif的数量 
as.data.frame(sort(table(motifEnrichment_selfMotifs_wGenes$highlightedTFs),decreasing = T)) 
#可视化 
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="IRF7"] viewMotifs(tableSubset)