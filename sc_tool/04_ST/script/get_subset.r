
library(Seurat)
library('getopt')
para<-matrix(c(
    'help','h',0,"logical",
    'rds','r',1,"character",
    'sample','s',1,"character"
    'outrds','o',1,"character"
),byrow=TRUE,ncol=4)
opt <- getopt(para,debug=F)


rds<-readRDS(opt$rds)
epi_subset <-subset(x=rds, subset=(term=="mEpi"|term=="vEpi"))
epi_subset@meta.data$orig.ident <- opt$sample
saveRDS( epi_subset, opt$outrds)
