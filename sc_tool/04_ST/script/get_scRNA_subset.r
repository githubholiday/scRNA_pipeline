print("第一步按照样本信息将不同细胞类型的rds中的样本提取出来，并进行简单质控")
print("第一步按照样本信息将不同细胞类型的rds中的样本提取出来，并进行简单质控")
print("该脚本用于将单细胞转录组不同细胞类型中的特定样本数据提取出来，保存到rds文件中,并将细胞等信息输出到表达中")

library(Seurat)
library('getopt')
para<-matrix(c(
    'help','h',0,"logical",
    'rds','n',1,"character",
    'outdir','o',1,"character"
),byrow=TRUE,ncol=4)
opt <- getopt(para,debug=F)


rds<-readRDS(opt$rds)
#提取以下三个样本的数据保存到新的rds文件中
mysubset <-  subset(x=rds, subset=(sample=="ASCP_1"|sample=="PDAC_P2"|sample=="PDAC_P3"))
meta<-data.frame(cell_id=rownames(mysubset@meta.data),cell_cluster=mysubset@meta.data$raw_cluster,cell_type=mysubset@meta.data$cluster_standard,sample=mysubset@meta.data$sample)
colnames(meta) <- c('cell_id','cell_cluster','cell_type','sample')
outfile <- paste(opt$outdir,"cell_type.tsv",sep='/')
write.table(meta,outfile,sep='\t',quote=F,row.names=F)

outrds <- paste(opt$outdir,"subset.rds",sep='/')
unique(mysubset@meta.data$sample)
saveRDS( mysubset, outrds)