library('getopt')
para<- matrix(c(
    'help', 'h',0,"logical",
    'infile','i',1,"character",
    'object', 's',1,"character",
    'outfile', 'o',1,"character"
),byrow=TRUE,ncol=4)

opt <- getopt(para,debug=FALSE)


print_usage <-  function(para=NULL){
    cat(getopt(para,usage=TRUE))
    cat("
    该脚本用于将rda文件中的基因表达量进行输出"
    )
    q(status=1)
}

library(Seurat)
sce = load(opt$infile)
sce
mydata <- opt$object
#write.csv(t(as.matrix(mormal2@assays$RNA@counts)),file = opt$outfile)#样本3_RNA
write.csv(t(as.matrix(plasma@assays$RNA@counts)),file = opt$outfile) #样本5——RNA