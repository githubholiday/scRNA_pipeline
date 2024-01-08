#名称：scatter_ggplot2.r
#作者：刘慧玲
#邮箱：huilingliu@genome.cn
#时间：20180604
#版本：v0.0.1
#用途：利用ggplot2绘制散点图，按照不同需求进行着色
###说明：
#程序开发环境 R-3.2.2
#主要功能为绘制有颜色区分的散点图
#===========================================================
library('getopt')
para<- matrix(c(
	'help',	'h',	0,	"logical",
	'infile',	'i',	1,	"character",
	'output',	'o',	1,	"character"
),byrow=TRUE,ncol=4)
opt <- getopt(para,debug=FALSE)
#===========================================================
print_usage <- function(para=NULL){
	cat(getopt(para,usage=TRUE))
	cat("
	Input format:
	Barcode	TSNE-1	TSNE-2	Cluster group
	AAACCTGAGAACAATC-1	-17.084590815492913	-14.682643460038516	Cluster1	KO
	AAACCTGAGGTAGCTG-1	6.520995680485729	-27.409997612070427	Cluster3	KO
	AAACCTGAGTGTCCAT-1	-44.83264455613646	5.7341402920770825	Cluster6	KO
	Usage example:
	Rscript this.r -i tsne_cluster.xls -o sample
	Options:
	--help		h	NULL		get this help
	--infile	i	character	the data file [forced]
	--output	o	character	PDF file name [forced]
	\n")
	q(status=1)	
}
#===========================================================
if ( !is.null(opt$help) )	{ print_usage(para) }
if ( is.null(opt$infile) )	{ cat("Please input the data file to draw plot...\n\n") ; print_usage(para)}
if ( is.null(opt$output) )	{ opt$output <- c("") }
#===========================================================
#library("AnnoColor")
#[1] "#4DBBD5FF" "#E64B35FF" "#00A087FF" "#3C5488FF" "#F39B7FFF" "#8491B4FF" "#91D1C2FF" "#DC0000FF"
#library(RColorBrewer)
library(ggplot2)
library(Cairo)
library(grid)
library(gridExtra)
data <- read.table(file=opt$infile,header=T,row.names=1,sep='\t',stringsAsFactors = FALSE,  check.names = FALSE,quote = "")
#print(head(data$'group'))
# sample 
outfile <- paste(opt$output,'group.pdf',sep='_')
Cairo(file = outfile, type='pdf',unit='cm', width=16, height=16)
p1 <- ggplot(data,aes(x=data$'TSNE-1', y=data$'TSNE-2', colour=factor(data$group)))+ geom_point(size=0.8)+
	coord_fixed(ratio = 1) + theme_bw(base_size = 12, base_family = "Arial") +  
	theme(panel.grid=element_blank(), legend.background = element_rect(colour = NA),
	legend.position='right', legend.margin=margin(rep(0,4)),
	legend.title = element_blank(),legend.text = element_text(face="bold", color="black", size=12),
	axis.text.x = element_text(color="black", size=12),
	axis.text.y = element_text(color="black", size=12),
	axis.title.x = element_text(face="bold", color="black", size=12),
	axis.title.y = element_text(face="bold", color="black", size=12))+
	labs(x='tSNE 1', y='tSNE 2')
grid.arrange(p1)
dev.off()
# cluster
outfile <- paste(opt$output,'cluster.pdf',sep='_')
Cairo(file = outfile, type='pdf',unit='cm', width=16, height=16)
p2 <- ggplot(data, aes(x=data$'TSNE-1', y=data$'TSNE-2', colour=factor(data$Cluster))) + geom_point(size=0.1)+
	coord_fixed(ratio = 1) + theme_bw(base_size = 12, base_family = "Arial") +  
	theme(panel.grid=element_blank(), legend.background = element_rect(colour = NA),
	legend.position='right', legend.margin=margin(rep(0,4)),
	legend.title = element_blank(),legend.text = element_text(face="bold", color="black", size=12),
	axis.text.x = element_text(color="black", size=12),
	axis.text.y = element_text(color="black", size=12),
	axis.title.x = element_text(face="bold", color="black", size=12),
	axis.title.y = element_text(face="bold", color="black", size=12))+
	labs(x='tSNE 1', y='tSNE 2')
grid.arrange(p2)
dev.off()
outfile <- paste(opt$output,'merge.pdf',sep='_')
Cairo(file = outfile, type='pdf',unit='cm', width=33, height=16)
grid.arrange(p1, p2, ncol=2, respect=FALSE)
dev.off()
