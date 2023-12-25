library('getopt')
para<-matrix(c(
	'help',	'h',	0,	"logical",
	'datadir',	'd',	1,	"character",
	'outdir',	'o',	1,	"character",
	'sample',	's',	1,	"character",
	'datatype', '',	2,	"character",
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
	--datatype		t	character	10X or 10XV2 [default: matrix]
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
if ( is.null(opt$testmethod))	{ opt$testmethod <- c("negbinom") }
if ( is.null(opt$datatype))	{ opt$datatype <- c("matrix") }
if ( is.null(opt$mt))	{ opt$mt <- c("^MT-") }


#======================
library(Seurat)
library(dplyr)
library(Matrix)
library(magrittr)

prefix <- paste(opt$outdir,opt$sample,sep='/')
## creat object

print(opt$datatype)
if (opt$datatype == 'matrix'){
	expression_matrix <- Read10X(data.dir = opt$datadir, gene.column = 2, unique.features = TRUE) ## 数据读入路径
} else if ( opt$datatype == 'h5') {
   input_matrix_dir <- Read10X_h5(filename = opt$datadir,use.names = TRUE, unique.features = TRUE) ## 数据读入路径
}else{
	print_usage
}
object <- CreateSeuratObject(counts = expression_matrix, min.cells = opt$mincell, min.features = opt$mingene, project = opt$sample) ### 设置最小细胞中表达的基因和最小基因表达数的细胞默认不进行过滤

nGene <- nrow(object@meta.data)
print(nGene)
file <- paste(prefix,'all_cell.csv',sep='_')
print(file)
### 解决生成特大矩阵问题 ###
### 用稀疏矩阵替换as.matrix ###
as_matrix <- function(mat){

  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
    
  for (i in seq_along(val)){
      tmp[row_pos[i],col_pos[i]] <- val[i]
  }
    
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

cell_data  <-  data.frame(data.frame(Gene=rownames(object@meta.data)),object@meta.data)
data.table::fwrite(data.table::data.table(cell_data),file,quote=F)
file <- paste(prefix,'all_UMI.csv',sep='_')
print(file)
raw_data<-as_matrix(object@assays$RNA@data)
raw_data  <-  data.frame(data.frame(Gene=rownames(raw_data)),raw_data)
data.table::fwrite(data.table::data.table(raw_data),file,quote=F,sep = ",",col.names = TRUE,row.names = FALSE)
#write.csv(as.matrix(object@assays$RNA@data),file,quote=F)
##QC
#mito.genes
#mito.genes <- grep(pattern = opt$mt, ignore.case = TRUE, x = rownames(x = object@data), value = TRUE)
#print(mito.genes)
#percent.mito <- Matrix::colSums(object@raw.data[mito.genes, ]) / Matrix::colSums(object@raw.data)
#object <- AddMetaData(object = object, metadata = percent.mito, col.name = "percent.mito")
# Expression Violin plot

#计算线粒体基因比例
#object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
graph <- paste(prefix,'Gene_UMI_Vln.pdf',sep='_')
unlink('Rplots.pdf')
pdf(graph,width=12,height=8) ## 注意画图输出
par(mfrow=c(1,2))
#VlnPlot(object = object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), nCol = 3, combine = "False")
plot1 <- VlnPlot(object = object, features = "nFeature_RNA", ncol = 1, combine = "False")
plot2 <- VlnPlot(object = object, features = "nCount_RNA", ncol = 1, combine = "False")
#plot3 <- VlnPlot(object = object, features = "percent.mt", nCol = 1, combine = "False")
CombinePlots(plots=c(plot1,plot2), ncol=2)
print(graph)
dev.off()
# UMI correlation scatter plot
graph<-paste(prefix,'Gene_UMI_Scatter.pdf',sep='_')
pdf(graph,width=12,height=8)
par(mfrow = c(1, 1)) 
#GenePlot(object = object, gene1 = "nUMI", gene2 = "percent.mito")
FeatureScatter(object = object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(graph)
dev.off()
#nGene = FetchData(object=object, vars.all="nGene")
#nUMI = FetchData(object=object, vars.all="nUMI")

