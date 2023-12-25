library('getopt')
rgs=commandArgs(T)
para<-matrix(c(
    'help','h',0,"logical",
    'rdsfile','rds',1,"character",
    'pathway_gene','pg',1,"character",
    'type','t',1,"character",
    'outd','od',1,"character",
    'outpdf','op',1,"character",
    'outtable','ot',1,"character"
),byrow=TRUE,ncol=4)
opt <- getopt(para,debug=FALSE)

if ( is.null(opt$rdsfile)){ print_usage }

print_usage <- function(para=NULL){
    cat(getopt(para,usage=TRUE))
    cat("
    Usage example:
    Rscript GSVA.r -rds file.rds -pg pathway_gene -c 1,2,3 -g TB,TE -t T -sample T1E,T1B,T2E,T2B -o outdir/ -name 'GSVA'
    Options:
    --helphNULLget this help
    --rdsfilerdscharacterthe rds file for all cells[forced]
    --pathway_genepgcharacterpathway and their corresponding gene list separated by '|'[forced]
    --outdodcharacter the outdir for result [forced]
    --outpdfopcharacter the output pdf for result [forced]
    --outtableotcharacter the output table for result [forced]
    --type cluster or group? [forced]
    \n")
    q(status=1)
}

#if ( is.null(opt$rdsfile)){ print_usage }
if ( is.null(opt$cluster)){ opt$cluster <- 'all' }
if ( is.null(opt$group)){ opt$group <- 'defult' }

require(Seurat)
library(GSVA)
library(pheatmap)

immune.combined <-readRDS(opt$rdsfile)
path_gene <-read.table(opt$pathway_gene,header=F,sep='\t')

pathway <-path_gene$V1
gene_all <-rownames(immune.combined)
geneset<-list()
for (i in path_gene$V1){
tmp_list <-path_gene[path_gene$V1==i,]$V2
tmp_list<-unlist(strsplit(as.character(tmp_list[1]),split = "|",fixed=T))
tmp_list_a <-list(intersect(gene_all,tmp_list))
geneset<-c(geneset,tmp_list_a)
}

names(geneset) <-path_gene$V1
Idents(immune.combined)<-immune.combined@meta.data[opt$type]
averages <- AverageExpression(object = immune.combined, assays ='RNA',return.seurat = F)

data_file <-paste(opt$outd,opt$outtable,sep='/')
pdf_file<-paste(opt$outd,opt$outpdf,sep='/')
exp<-averages$RNA
es.max_tmp <- gsva(as.matrix(exp), geneset,kcdf="Poisson", verbose=FALSE, parallel.sz=32)
write.table(es.max_tmp,data_file,quote=F,sep='\t')
pheatmap(es.max_tmp,show_rownames=T,show_colnames=T,cluster_cols = T,cluster_rows= T,treeheight_row=0,cellheight=12,fontsize=8,border=FALSE,color = colorRampPalette(c('blue', "white", 'red'))(50),main ='',filename = pdf_file,width = 10,height=6.18)


cat('恭喜您，完成了GSVA的所有分析\n')

