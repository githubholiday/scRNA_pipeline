print("第二步将不同细胞类型的subset rds文件合并到一起")

library(Seurat)
library('getopt')
para<-matrix(c(
    'help','h',0,"logical",
    'rds1','r1',1,"character",
    'rds2','r2',1,"character",
    'outfile','o',1,"character"
),byrow=TRUE,ncol=4)
opt <- getopt(para,debug=F)

rds1 <- readRDS( opt$rds1)
rds2 <- readRDS( opt$rds2)
combine <- merge( rds1, y=rds2 )
saveRDS( combine, opt$outfile )