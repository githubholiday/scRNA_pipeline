library(Seurat)
library('getopt')
para<-matrix(c(
    'help','h',0,"logical",
    'rdsfile','r',1,"character",
    'outfile','o',1,"character"
),byrow=TRUE,ncol=4)
opt <- getopt(para,debug=F)


rds<-readRDS(opt$rdsfile)

#sub<-subset(rds,downsample = opt$num)
meta<-data.frame(cell_id=rownames(rds@meta.data),cell_cluster=rds@meta.data$seurat_clusters,cell_type=rds@meta.data$labels,sample=rds@meta.data$orig.ident)

colnames(meta) <- c('cell_id','cell_cluster','cell_type','sample')
write.table(meta,opt$outfile,sep='\t',quote=F,row.names=F)