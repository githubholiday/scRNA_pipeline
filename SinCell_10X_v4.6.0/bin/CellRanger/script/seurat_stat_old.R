library('getopt')
para<-matrix(c(
	'help',	'h',	0,	"logical",
	'datadir',	'd',	1,	"character",
	'outdir',	'o',	1,	"character",
	'sample',	's',	1,	"character",
	'mincell',	'mc',	2,	"character",
	'mingene',	'mg',	2,	"character",
	'mt',	'mt',	2,	"character",
	'testmethod',	'tm',	2,	"character"
),byrow=TRUE,ncol=4)
opt <- getopt(para,debug=FALSE)
#=============
print_usage <- function(para=NULL){
	cat(getopt(para,usage=TRUE))
	cat("
	Usage example:
	Rscript securat_parameter.R -d hg19/ -p ANproject -o outdir/
	Options:
	--help		h	NULL		get this help
	--datadir	d	character	the cellRanger output datadir [forced]
	--sample	s	character	SampleName[forced]
	--outdir	o	character	output file dir [forced]
	--mincell	mc	character	Filter minium cell numbers [default: 3]
	--mingene	mg	character	Filter minium gene numbers [default: 200]
	--mt	mt	character	mitochondria Gene symbol [default : MT]
	--testmethod	tm	character	test method for findMakers [default: negbinom]
	\n")
	q(status=1)
}
if ( is.null(opt$datadir) || is.null(opt$sample) || is.null(opt$outdir)){ print_usage } 
if ( is.null(opt$mincell))	{ opt$mincell <- c("0") }
if ( is.null(opt$mingene))	{ opt$mingene <- c("0") }
if ( is.null(opt$mt))	{ opt$mt <- c("^MT-") }

#======================
library(Seurat)
library(dplyr)
library(Matrix)
library(magrittr)

prefix <- paste(opt$outdir,opt$sample,sep='/')
## creat object
expression_matrix <- Read10X(data.dir = opt$datadir) ## 数据读入路径
object <- CreateSeuratObject(raw.data = expression_matrix, min.cells = opt$mincell, min.genes = opt$mingene, project = opt$sample) ### 设置最小细胞中表达的基因和最小基因表达数的细胞

#head(object@raw.data[,1:3])
#head(object@data[,1:3])
#head(object@meta.data[,1:3])
#head(object@ident)
nGene <- nrow(object@meta.data)
print(nGene)
file <- paste(prefix,'all_cell.csv',sep='_')
print(file)
write.csv(object@meta.data,file,quote=F)
file <- paste(prefix,'all_UMI.csv',sep='_')
print(file)
write.csv(as.matrix(object@data),file,quote=F)
##QC
#mito.genes
#mito.genes <- grep(pattern = opt$mt, ignore.case = TRUE, x = rownames(x = object@data), value = TRUE)
#print(mito.genes)
#percent.mito <- Matrix::colSums(object@raw.data[mito.genes, ]) / Matrix::colSums(object@raw.data)
#object <- AddMetaData(object = object, metadata = percent.mito, col.name = "percent.mito")
# Expression Violin plot
graph <- paste(prefix,'Gene_UMI_Vln.pdf',sep='_')
unlink('Rplots.pdf')
pdf(graph,width=12,height=8) ## 注意画图输出
VlnPlot(object = object, features.plot = c("nGene", "nUMI"), nCol = 2)
print(graph)
dev.off()
# UMI correlation scatter plot
graph<-paste(prefix,'Gene_UMI_Scatter.pdf',sep='_')
pdf(graph,width=12,height=8)
par(mfrow = c(1, 1)) 
#GenePlot(object = object, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = object, gene1 = "nUMI", gene2 = "nGene")
print(graph)
dev.off()
#nGene = FetchData(object=object, vars.all="nGene")
#nUMI = FetchData(object=object, vars.all="nUMI")

