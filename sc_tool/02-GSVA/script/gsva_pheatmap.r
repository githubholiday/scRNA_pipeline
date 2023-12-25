library('getopt')
library('pheatmap')
para<- matrix(c(
	'help',	'h',	0,	"logical",
	'infile',	'i',	1,	"character",
	'outfile',	'o',	1,	"character"
),byrow=TRUE,ncol=4)

opt <- getopt(para,debug=FALSE)
#===========================================================
print_usage <- function(para=NULL){
	cat(getopt(para,usage=TRUE))
	q(status=1)	
}

data <- read.table(file=opt$infile,header=T,row.names=1,sep='\t',stringsAsFactors = FALSE,  check.names = FALSE,quote = "")

data = as.matrix(data)
order <- order(data[,1])

data_t <- data[order,]
pdf( opt$outfile, width=20, height=10)
pheatmap(data_t,scale="none",border=FALSE,color=colorRampPalette(c("steelblue","white","firebrick3"))(100),cellheight=20,cellwidth=20,fontsize=8)
dev.off()