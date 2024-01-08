library(getopt)

para <- matrix(c(
    "help",   "h",	0,	"logical",
    "prefix",	"p",	1,	"character",
    "geneinfo",	"g",	1,	"character",
    "species",   "s",   1, "character",
    "maplist,",   "m",   1, "character",
    "kegg_dir,",   "k",   1, "character", 
    "column,",   "n",   1, "integer", 
    "limit,",   "l",   1, "double", 
    "config,",   "c",   1, "character"), byrow = TRUE, ncol = 4)

opt <- getopt(para, debug = FALSE)
print_usage <- function(para = NULL) {
    cat(getopt(para, usage = TRUE))
    cat("
    Options:
    --help		h	NULL		get this help
    --prefix	p	dir		    dir to cantain all the new map[forced]
    --geneinfo 	g	file	    input file that the first column is genelist, other info is mapped info maps[forced]
    --species   s   character   species, homo_sapiens or mus_musculus [forced]
    --maplist	m	file		the first line is map ID [forced]
    --column	n	int	        the column that used to draw map[default is 2]
    --config	c	file		OrgDB list 
    --kegg_dir	k	dir	    	kegg_map dir  [optional]
    \n")
    q(status = 1)
}
if (!is.null(opt$help))	{
    print_usage(para)
}
if (is.null(opt$geneinfo)) {
    cat("Please input the gene info data  ...\n\n")
     print_usage(para)
}
if (is.null(opt$prefix)) {
    cat("Please speicfify the output prefix...\n\n")
    print_usage(para)
}
if (is.null(opt$maplist)) {
    cat("Please speicfify the maplist prefix ...\n\n")
    print_usage(para)
}
if (is.null(opt$kegg_dir)) {
    cat("Please speicfify the kegg_dir ...\n\n")
    print_usage(para)
}

path=commandArgs(trailingOnly=F)
if (is.null(opt$config)) {
    script_dir=dirname(sub("--file=","",path[grep("--file",path)]))
    config = paste(c(script_dir,"/OrgDb.list"),collapse="")
}else{
	config =opt$config
}
column=1
if (is.null(opt$column)) {
	column=1
}else{
	column =opt$column
}

species="no"
kegg_species="no"
if (is.null(opt$species)){
    cat("Please speicfify the speceise ...\n\n")
    print_usage(para)
}else{
	orgdb = read.table(config,T)
	species =orgdb$OrgDb[which(orgdb$species==opt$species)]
	kegg_species =orgdb$kegg[which(orgdb$species==opt$species)]
	print(species)
    print(kegg_species)
	if(length(species) == 0){ 
        print("[Error] your species is wrong, please choose right species based on config " )
        q(status=1)
    }
}
library(dplyr)
library(pathview)
library(clusterProfiler)


## 读取物种信息

#/annoroad/data1/bioinfo/PMO/suyalei/software/Anaconda/minconda3/envs/r_4/bin/R
#在外部传入前就先提取好，不在该脚本里处理文件
id_output=paste(opt$prefix,"gene.id.xls",sep=".")

input_table <- read.table(opt$geneinfo,T,sep="\t")
map_table <- read.table(opt$maplist,T,sep="\t")
if (length(map_table) <1){
    print("there is no map info to get pathview map")
    q(status=0)
}
if(column == 1){
	input_table$Info=1
}else{
	names(input_table)[column] <- "Info"
}
names(input_table)[1] <- "Gene"
names(map_table)[1] <- "ID"
input_table_new <- distinct(input_table, Gene, .keep_all=T)
ids = bitr(input_table_new$Gene, fromType="ENSEMBL", toType=c("ENTREZID","SYMBOL"), OrgDb = species)
write.table(ids,id_output,quote = FALSE , sep="\t" , row.names=FALSE)


ID_geneinfo = merge(ids,input_table_new,by.x="ENSEMBL",by.y="Gene")
gene_data=ID_geneinfo$Info

gene_max = max(abs(gene_data))
if (is.null(opt$limit)) {
	gene_max = max(abs(gene_data))
}else{
	gene_max = opt$limit
}
attr(gene_data,"names") = ID_geneinfo$ENTREZID
newID=gsub("map","",map_table$ID)
for (i in newID){
    tryCatch(
        {
        pathview(gene.data = gene_data,
            pathway.id = i,
            species    = kegg_species, 
            out.suffix =  opt$prefix,
            kegg.native = T, same.layer = T, kegg.dir=opt$kegg_dir,
            limit      = list(gene=gene_max))
        },
            warning = function(w) {message("There is something warnings, maybe Negative values!")},
            error = function(e) {message("There is something wrong, please Enter numeric value!")}
    )
}
