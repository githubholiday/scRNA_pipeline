library('getopt')
rgs=commandArgs(T)
para<-matrix(c(
    'help','h',0,"logical",
    'rdafile','r',1,"character",
    'database','d',1,"character",
    'outfile','o',1,"character"
),byrow=TRUE,ncol=4)
opt <- getopt(para,debug=FALSE)

#if ( is.null(opt$rdafile)){ print_usage }

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

immune.combined <- load(opt$rdafile)
sub <- normal2
#获取表达量信息
expMat<-GetAssayData(sub,slot="data")
mydata<-as.matrix(expMat)
annotation_col<-data.frame(Type=Idents(sub))

#读取数据库
geneSets<-GSEABase::getGmt(opt$database)

#计算GSVA
res_es <- gsva(mydata, geneSets,min.sz=10, max.sz=500, kcdf="Gaussian", verbose=FALSE, method = "gsva",parallel.sz=32)

write.table(t(res_es),opt$outfile,quote=F,sep='\t')


#绘图
#meta <- as.data.frame(sub@meta.data[,c('orig.ident',"labels")])
#meta <- meta %>%arrange(meta$labels) 
#table(meta$labels) #展示每个细胞类型的细胞数量

#data <- res_es
#data <- data[,rownames(meta)]
#identical(colnames(data),rownames(meta))




#pheatmap(es.max_tmp,show_rownames=T,show_colnames=T,cluster_cols = T,cluster_rows= T,treeheight_row=0,cellheight=12,fontsize=8,border=FALSE,color = colorRampPalette(c('blue', "white", 'red'))(50),main ='',filename = pdf_file,width = 10,height=6.18)


cat('恭喜您，完成了GSVA的所有分析\n')

