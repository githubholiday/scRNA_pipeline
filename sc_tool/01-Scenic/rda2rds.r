library('getopt')
library(Seurat)
para<- matrix(c(
	'help',	'h',	0,	"logical",
	'infile',	'i',	1,	"character",
	'outfile',	'o',	2,	"character" 
),byrow=TRUE,ncol=4)
opt <- getopt(para,debug=FALSE)
sce<-load(opt$infile)
sce
saveRDS(plasma,file=opt$outfile) 
#saveRDS(normal2,file=opt$outfile) 