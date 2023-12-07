library(getopt)

para <- matrix(c(
    "help",   "h",	0,	"logical",
    "output",	"o",	1,	"character",
    "input",	"i",	1,	"character",
    "species",   "s",   2, "character",
    "term2gene,",   "g",   2, "character",
    "term2name,",   "n",   2, "character", 
    "head,",   "t",   1, "logical", 
    "maxgene,",   "m",   1, "integer", 
    "pvalue,",   "p",   1, "double", 
    "config,",   "c",   1, "character"), byrow = TRUE, ncol = 4)

opt <- getopt(para, debug = FALSE)
print_usage <- function(para = NULL) {
    cat(getopt(para, usage = TRUE))
    cat("
    Options:
    --help		h	NULL		get this help
    --output	o	file		outfile for enrichment report[forced]
    --input 	i	file	input file that the first column is genelist [forced]
    --species   s   character   species, homo_sapiens or mus_musculus [optional]
    --term2gene	g	file	go term of genelist  [optional]
    --term2name	n	file	go term of genename  [optional]
    --config	c	file	go term of genename  [optional]
    --padj	p	float	pvalue threadshold,default padj <0.05 is significant  [optional]
    --head		t	NULL	whehter the tabel has title T or F [optional,default is T]
    --maxgene	m	int	maximal size of genes annotated for testing [optional,default is 1000]
    \n")
    q(status = 1)
}
if (!is.null(opt$help))	{
    print_usage(para)
}
if (is.null(opt$input)) {
    cat("Please input the data file1 ...\n\n")
     print_usage(para)
}
if (is.null(opt$output)) {
    cat("Please speicfify the output prefix ...\n\n")
    print_usage(para)
}
if (is.null(opt$head)) {
	head = TRUE
}else{
	head = opt$head
}
if (is.null(opt$maxgene)) {
	maxGene =1000
}else{
	maxGene = opt$maxgene
}
if (is.null(opt$head)) {
	head = TRUE
}else{
	head = opt$head
}
if (is.null(opt$padj)) {
	padj = 0.05
}else{
	padj = opt$padj
}
path=commandArgs(trailingOnly=F)
if (is.null(opt$config)) {
    script_dir=dirname(sub("--file=","",path[grep("--file",path)]))
    config = paste(c(script_dir,"/OrgDb.list"),collapse="")
}else{
	config =opt$config
}
std="yes"
species="no"
if (is.null(opt$species)){
	std = "no"
}else{
	orgdb = read.table(config,T)
	species =orgdb$OrgDb[which(orgdb$species==opt$species)]
	print(species)
	if(length(species) == 0){ std ="no"}
}
if (std =="no"){
	if (is.null(opt$term2gene) || is.null(opt$term2name))  {
    	cat("Please input species or provide term2gene and term2name to condact local go enrichment...\n\n")
    	print_usage(para)
	}
}

library(readr)
library(dplyr)
library(clusterProfiler)
#library(enrichplot)
get_go <- function(x){
	result=data.frame("go"=c(),"gene"=c())
	for (i in 1:nrow(x)){
		for (j in 2:length(x[i,])){
			new=data.frame("go"=c(x[i,j]),"gene"=c(x[i,1]))
			result = rbind(result,new)
		}
	}
	result
}

## 读取物种信息
#/annoroad/data1/bioinfo/PMO/suyalei/software/Anaconda/minconda3/envs/r_4/bin/R
#在外部传入前就先提取好，不在该脚本里处理文件
input_table <- read_tsv(file = opt$input,head)
#de_table <- filter(input_table, Significant == "yes")
names(input_table)[1] <- "Gene"
##基因列去重
input_table <- distinct(input_table, Gene)

if(std == "yes"){
	print("go enrichment starts with standard database")
	ego <- enrichGO(gene         = input_table$Gene,
                OrgDb         = species,
                keyType       = "ENSEMBL",
                ont           = "ALL",
                pAdjustMethod = "BH",
                maxGSSize = maxGene,
                pvalueCutoff  = 1,
                qvalueCutoff  = 1)
	ego@result %>% arrange(ONTOLOGY,p.adjust) %>%mutate(Significant=ifelse(p.adjust<padj,"yes","no"))->result
	write.table(result,opt$output,
            quote = FALSE , sep="\t" , row.names=FALSE)
}else{
	print("go enrichment starts with local database: go term2gene and term2name")
	term2gene <- read.table(opt$term2gene,header=F,sep="\t")
	term2name <- read.table(opt$term2name, header=T ,sep="\t",quote="")
	a= intersect(input_table$Gene,term2gene[,2])
	if (length(a)==0){
		print("Gene in input has no Go info in golist")
		q()
	}
	ego <- enricher(gene=input_table$Gene,  ## 基因集合
           pvalueCutoff = 1, 
           pAdjustMethod = "BH", 
           qvalueCutoff = 1, 
           minGSSize = 2,
           maxGSSize = maxGene,
           TERM2GENE = term2gene, ## GO  gene
           TERM2NAME = term2name  ## GO  name
       )
	result=filter(ego@result,Description!="NA")
	result_new = merge(result,term2name,by.x="ID",by.y="GO",all.x=T)
	namelist=c("ONTOLOGY",names(result))
	result_new[,namelist] %>% arrange(ONTOLOGY,p.adjust) %>%mutate(Significant=ifelse(p.adjust<padj,"yes","no"))->result_new
	write.table(result_new,opt$output,
            quote = FALSE , sep="\t" , row.names=FALSE)
}

print("go enrichment finished")

