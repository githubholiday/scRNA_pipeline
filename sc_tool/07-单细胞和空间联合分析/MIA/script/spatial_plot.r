

library('getopt')
para<- matrix(c(
	'help',	'h',	0,	"logical",
	'infile',	'i',	1,	"character",
	'rds',	'r',	1,	"character",
    'outfile',	'o',	1,	"character"

),byrow=TRUE,ncol=4)
opt <- getopt(para,debug=FALSE)
print_usage <- function(para=NULL){
	cat(getopt(para,usage=TRUE))
	cat("
	==============================================
	Usage example:
	Rscript this.r -p prefix -r RNArds -o outdir -c config.ini -s spacerds
	Options:
	--help	h	NULL		get this help
	--infile	i	character	the file of MIA enrichment and group column [forced]
	--outpre	o	character	the prefix of outfile[forced]
	\n")
	q(status=1)
}
library(Seurat)
st <- readRDS(opt$rds)
lable_names <- read.csv(opt$infile,header=T,sep="\t")
labels <- lable_names[match(as.numeric(as.character(st@active.ident)), lable_names[,1]),2]

st@meta.data$group <-  labels
p1 <- SpatialPlot(st, group.by="group")
pdf(opt$outfile)
print(p1)
dev.off()
