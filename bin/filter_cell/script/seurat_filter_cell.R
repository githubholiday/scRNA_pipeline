library('getopt')
para<-matrix(c(
	'help',	'h',	0,	"logical",
	'datadir',	'd',	1,	"character",
	'outdir',	'o',	1,	"character",
	'sample',	's',	1,	"character",
	'minumi',	'mu',	2,	"integer",
	'mincell',	'mc',	2,	"integer",
	'min',	'min',	2,	"integer",
	'max',	'max',	2,	"integer",
	'mtp',	'mtp',	2,	"integer",
	'hbp',	'hbp',	2,	"integer",
	'species',	'sp',	1,	"character",
	'mtt',	'mtt',	1,	"character",
	'doublet',	'doub',	1,	"character"
),byrow=TRUE,ncol=4)
opt <- getopt(para,debug=FALSE)
#=============
print_usage <- function(para=NULL){
	cat(getopt(para,usage=TRUE))
	cat("
	Usage example:
	Rscript securat_filter_cell.R -d datadir -s sample -o outdir --species
	Options:
	--help		h	NULL		get this help
	--datadir	d	character	the cellRanger output datadir [forced]
	--sample	s	character	SampleName[forced]
	--outdir	o	character	output file dir [forced]
	--minumi	mu	character	Filter genes have minium umi count [default: 200]
	--mincell	mc	character	Filter genes expession based on minium cell numbers [default: 3]
	--min	min	character	Filter cells based on minium gene numbers [default: 200]
	--max	max	character	 Filter cells based on max gene numbers [default: 10000]
	--mtp	mtp	character	Filter cells based on mito percent [default:20]
	--hbp	hbp	character	Filter cells based on HB percent [default:5]
	--mtt	mtt	character	mitochondria Gene symbol []
	--species	sp	character	species[human:9606 mouse:10090]
	--doublet	doub	character	doublet_cell file
	\n")
	q(status=1)
}
if ( is.null(opt$datadir) || is.null(opt$sample) || is.null(opt$outdir) || is.null(opt$species)){ print_usage(para) } 
if ( is.null(opt$mincell))	{ opt$mincell <- 3 }
if ( is.null(opt$minumi))	{ opt$minumi <- 200 }
if ( is.null(opt$min))	{ opt$min <- 200 }
if ( is.null(opt$max))	{ opt$max <- 10000 }
if ( is.null(opt$mtp))	{ opt$mtp <- 20 }
if ( is.null(opt$hbp))	{ opt$hbp <- 5}

#======================
library(Seurat)
library(dplyr)
library(Matrix)
library(magrittr)
library(configr)


prefix <- paste(opt$outdir,opt$sample,sep='/')
## creat object
expression_matrix <- Read10X(data.dir = opt$datadir) ## 数据读入路径
object <- CreateSeuratObject(counts = expression_matrix, project = opt$sample, min.cells = opt$mincell)#, min.features = opt$minumi) ### 设置最小细胞中表达的基因和最小基因表达数的细胞
nGene <- nrow(object@meta.data)
print(nGene)

doublet_cell <- read.csv(opt$doublet, header = T)
doublet_score <- doublet_cell[match(rownames(object@meta.data), doublet_cell$Barcode), 3]
object@meta.data$DoubletCell <- doublet_score

#threshold <- unname(quantile(object@meta.data$nFeature_RNA,0.25)) #细胞中表达基因数目的下四分位数，probs参数可调
#print(threshold)
#threshold <-max(threshold,200)
file <- paste(prefix,'filter_UMI.csv',sep='_')
print(file)
species<<-opt$species
s.features<-cc.genes$s.genes
g2m.features<-cc.genes$g2m.genes
if (species == '9606'){ 
  HB.genes_total <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") # 人类血液常见红细胞基因
  HB_m <- match(HB.genes_total,rownames(object@assays$RNA))
  HB.genes <- rownames(object@assays$RNA)[HB_m]
  HB.genes <- HB.genes[!is.na(HB.genes)]
  object[["percent.HB"]]<-PercentageFeatureSet(object,features=HB.genes)
 }else if (species == '10090'){
   HB.genes_total <- c("Hbb-bt","Hbb-bs","Hbb-bh2","Hbb-bh1","Hbb-y","Hba-x","Hba-a1","Hba-a2") # 小鼠血液常见红细胞基因
   HB_m <- match(HB.genes_total,rownames(object@assays$RNA))
   HB.genes <- rownames(object@assays$RNA)[HB_m]
   HB.genes <- HB.genes[!is.na(HB.genes)]
   object[["percent.HB"]]<-PercentageFeatureSet(object,features=HB.genes)
   library(stringr)
   s.features <- str_to_title(cc.genes$s.genes)
   g2m.features <- str_to_title(cc.genes$g2m.genes)
 }else{
   object@meta.data$percent.HB<-rep(0,length(rownames(object@meta.data)))
}
A<<-opt$min
B<<-opt$max
mt_percent<<-opt$mtp
HB_percent<<-opt$hbp
print(A)
print(B)
print(mt_percent)
print(HB_percent)
mt_gene_index <- grep('^MT-',rownames(object@assays$RNA),ignore.case = TRUE) ##不区分大小写,检索基因中前3个字符串为mt-的基因作为线粒体基因
MT_genes <-rownames(object@assays$RNA)[mt_gene_index]
print(MT_genes)
#object[["percent.mt"]] <- PercentageFeatureSet(object, features=MT_genes)
if ( is.null(MT_genes) && is.null(mitoName)){ ##如果没有线粒体基因,且没给前缀名，就退出
	q(status=1)
}else if (is.null(MT_genes)){
	mitoName<<-opt$mtt
	object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = paste("^", mitoName,"-",sep=""))
}else{
	object[["percent.mt"]] <- PercentageFeatureSet(object, features=MT_genes)
}


#object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = paste("^", mitoName,"-",sep=""))
print('start filter')
file_stat <- paste(prefix,'filter_stat_cell.csv',sep='_')
removed_cells<-data.frame(sample=0,total_cell=0,high_quality_cell=0,low_nFeature=0,high_nFeature=0,high_MT=0,high_HB=0,all_filtered_cell=0)
del_cells<-data.frame(sample=opt$sample,total_cell=as.character(nrow(object@meta.data)),high_quality_cell=as.character(nrow(subset(object@meta.data,nFeature_RNA > A & nFeature_RNA < B & percent.mt < mt_percent & percent.HB < HB_percent))),low_nFeature=as.character(nrow(subset(object@meta.data,nFeature_RNA <= A))),high_nFeature=as.character(nrow(subset(object@meta.data,nFeature_RNA >= B))),high_MT=as.character(nrow(subset(object@meta.data,percent.mt >= mt_percent))),high_HB=as.character(nrow(subset(object@meta.data,percent.HB >= HB_percent))),all_filtered_cell=as.character(nrow(subset(object@meta.data,nFeature_RNA <= A | percent.mt >= mt_percent|nFeature_RNA >= B|percent.HB >= HB_percent))))
removed_cells<<-rbind(removed_cells,del_cells)
write.csv(removed_cells[-1,],file_stat,quote=F,row.names=F)
object <- subset(object,  subset = nFeature_RNA > A & nFeature_RNA < B & percent.mt < mt_percent & percent.HB < HB_percent)
object
## 改为fwrite写入
#raw_data<-as.matrix(object@assays$RNA@counts)
#raw_data  <-  data.frame(data.frame(Gene=rownames(raw_data)),raw_data)
#data.table::fwrite(data.table::data.table(raw_data),file = file,sep = ",",col.names = TRUE,row.names = FALSE,quote=F)
file <- paste(prefix,'filter_cell.csv',sep='_')
print(file)
b <- object@meta.data
b2 <- data.frame(barcode=rownames(b), b)
#write.csv(object@meta.data,file,quote=F)
write.table(b2,file,quote=F,sep="\t",row.names=F)
#cell_for_heatmap = round(length(object@meta.data[,1])*0.05)
#print(cell_for_heatmap)

##QC
# Expression Violin plot
gene_umi_mito <- paste(prefix,'Gene_UMI_mito_percent.pdf',sep='_')
#unlink('Rplots.pdf')
pdf(gene_umi_mito,width=16,height=8) ## 注意画图输出
#VlnPlot(object = object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), nCol = 3)
plot1 <- VlnPlot(object = object, features = "nFeature_RNA", ncol = 1, combine = "False")
plot2 <- VlnPlot(object = object, features = "nCount_RNA", ncol = 1, combine = "False")
plot3 <- VlnPlot(object = object, features = "percent.mt", ncol = 1, combine = "False")
CombinePlots(plots=c(plot1,plot2,plot3), ncol=3)
print(gene_umi_mito)
dev.off()
# UMI correlation scatter plot
#graph<-paste(prefix,'GenePlot.pdf',sep='_')
#pdf(graph,width=12,height=8)
#par(mfrow = c(1, 2)) 
#plot1 <- FeatureScatter(object = object, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(object = object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#CombinePlots(plots = list(plot1, plot2))
#print(graph)
#dev.off()
#object <- NormalizeData(object = object, normalization.method = ini.list$Para$normalization.method, scale.factor = as.numeric(ini.list$Para$scale.factor), verbose = FALSE)
#object<- FindVariableFeatures(object = object, selection.method = opt$findvariablemethod, nfeatures = as.numeric(opt$findvariablefeatures))

#细胞周期

object <- NormalizeData(object = object, verbose = FALSE)
print(s.features)
print(g2m.features)
object<-CellCycleScoring(object, s.features=s.features, g2m.features=g2m.features)

file <- paste(prefix,'filter_cell.rds',sep='_')
saveRDS(object, file)
