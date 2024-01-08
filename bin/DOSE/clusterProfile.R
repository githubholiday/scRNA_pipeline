#!Rscript
.libPaths(new="/annoroad/data1/bioinfo/PROJECT/big_Commercial/Cooperation/B_TET/B_TET-024/std/result/changjiyang/yaomengcheng/anaconda3/Anaconda3/lib/R/library")
#library(clusterProfiler)
#library(DOSE)
#library(org.Hs.eg.db)
#library(pathview)
#library(ReactomePA)
#library(dplyr)
library('getopt')
para<- matrix(c(
	'help',	'h',	0,	"logical",
	'diffgene',	'g',	1,	"character",
	'header',	't',	1,	"logical",
	'species',	's',	1,	"character",
	'prefix',	'f',	1,	"character",
	'foldup',	'u',	1,	"numeric",
	'folddown',	'd',	1,	"numeric",
	'pvalueCutoff',	'p',	1,	"numeric",
	'qvalueCutoff',	'q',	1,	"numeric",
	'outdir',	'o',	1,	"character",
	'bindir',	'b',	1,	"character"	
),byrow=TRUE,ncol=4)
opts <- getopt(para)
print_usage <-function(para=NULL){
	cat(getopt(para,usage=TRUE))
	cat("
使用方法：
Rscript analysi.R -g gene.txt -t F/T -s hsa -f prefix -u 1.5 -d 0.3 -p 1 -q 1 -b /
参数说明：
--help	-h	NULL	get this help
--diffgene	-g	character	the differential gene file	#差异基因
--header	-t	character	T or F	#文件是否有表头
--species	-s	character	the name of species	#物种缩写hsa/mmu/cbr
--prefix	-f	character	the prefix of output file	#文件名前缀
--foldup	-u	character	Foldchanger Up	#用以过滤差异倍数UP
--folddown	-d	character	Foldchanger Down	#用以过滤差异倍数Down 
--pvalueCutoff	-p	character	P value cutoff	#富集P阈值
--qvalueCutoff	-q	character	q value cutoff  #富集q阈值
--outdir	-o	character	Dir of result	#输出结果路径	
--bindir	-b	character	Bindir	#Bin路径
\n")
	q(status=1)
}
if (!is.null(opts$help)){print_usage(para)};if (is.null(opts$diffgene)){print_usage(para)}
if (is.null(opts$species)){print_usage(para)};if (is.null(opts$prefix)){print_usage(para)};if (is.null(opts$outdir)){print_usage(para)}
mkdirs <- function(outdir,fp) {
        if(!file.exists(file.path(outdir,fp))) {
                dir.create(file.path(outdir,fp))
        }else{
                        print(paste(fp,"Dir already exists!",sep="     "))
                        #unlink(file.path(outdir,fp), recursive=TRUE)
                        #dir.create(file.path(outdir,fp))
                }
}
diffgene <- opts$diffgene
header <- opts$header
species <- opts$species
prefix <- opts$prefix
foldup <- opts$foldup
folddown <- opts$folddown
pvalueCutoff <- opts$pvalueCutoff
qvalueCutoff <- opts$qvalueCutoff
filter <- opts$filter
outdir <- opts$outdir 
bindir <- opts$bindir
print(diffgene)
print(species)
print(prefix)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Bt.eg.db)
library(org.Rn.eg.db)
library(pathview)
library(ReactomePA)
#library(dplyr)
library(enrichplot)

diffgene<-read.csv(file=diffgene,header=header,sep="\t")##读取差异基因文件第一列必须为ENSEMBL，第二列为change
if  ( species =='hsa' ) {
list=select(org.Hs.eg.db,keys=keys(org.Hs.eg.db,keytype = "ENSEMBL"),columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
}
if ( species =='mmu' ) {
list=select(org.Mm.eg.db,keys=keys(org.Mm.eg.db,keytype = "ENSEMBL"),columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
}
if ( species =='cbr' ) {
list=select(org.Bt.eg.db,keys=keys(org.Bt.eg.db,keytype = "ENSEMBL"),columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
}
if ( species =='rat' ) {
list=select(org.Rn.eg.db,keys=keys(org.Rn.eg.db,keytype = "ENSEMBL"),columns = c("ENTREZID","SYMBOL"), keytype="ENSEMBL")
}
colnames(diffgene)<-c("ENSEMBL","FoldChange")
new_diff<-merge(diffgene,list,by="ENSEMBL",all.x=F)[,c("ENTREZID","FoldChange")]
head(new_diff)
diff<-new_diff[,2]
names(diff) <- as.character(new_diff[,1])
diff <- sort(diff, decreasing = TRUE) ## 基因降序排序
if ( foldup == 0 ){
	gene <- names(diff) 
}else{
	gene <- names(diff)[diff >= foldup | diff <= folddown ]
	}
####进行 wikipathway analysis
#mkdirs(outdir,prefix)
#setwd(paste(outdir,prefix,sep='/'))
if (species=='hsa') {
wpgmtfile<-read.gmt(paste(bindir,"database/wikipathways-gmt-Homo_sapiens.gmt",sep='/'))
org.db <- org.Hs.eg.db }
if (species=='mmu') { 
wpgmtfile<-read.gmt(paste(bindir,"database/wikipathways-gmt-Mus_musculus.gmt",sep='/'))
org.db <- org.Mm.eg.db }
if (species=='cbr') {
wpgmtfile<-read.gmt(paste(bindir,"database/wikipathways-gmt-Bos_taurus.gmt",sep='/'))
org.db <- org.Bt.eg.db }
if (species=='rat'){
wpgmtfile<-read.gmt(paste(bindir,"database/wikipathways-gmt-Rattus_norvegicus.gmt",sep='/'))
org.db <- org.Rn.eg.db }
wp2gene<-wpgmtfile
library(dplyr)
wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
##富集wikipathway
ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name,pvalueCutoff=pvalueCutoff,qvalueCutoff=qvalueCutoff) ##为了获取更多的信息建议设置较大的阈值
ewp2 <- GSEA(diff, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE,pvalueCutoff=pvalueCutoff)
ewp <- setReadable(ewp, org.db, keyType = "ENTREZID")
ewp2 <- setReadable(ewp2, org.db, keyType = "ENTREZID")
enrich_number<-nrow(ewp)
enrichmin<-min(10,enrich_number)
gsea_number<-nrow(ewp2)
gseamin<-min(10,gsea_number)
print(enrichmin)
print(gseamin)
#保存富集分析结果并画图
mkdirs(paste(outdir,prefix,sep='/'),'WikiPathway')
setwd(paste(outdir,prefix,'WikiPathway',sep='/'))
file <- paste(prefix,'enrichPathway.xls',sep='_')
write.table(ewp,file=file,sep="\t",row.names=F,quote=F)
print(file)
file <- paste(prefix,'GSEAPathway.xls',sep='_')
write.table(ewp2,file=file,sep="\t",row.names=F,quote=F)
print(file)
## plot bar+dot+cnetplot+heatmap+gseaplot
graph <- paste(prefix,'enrichPathway_bar.pdf',sep='_')
pdf(graph,w=12, h=8)
print(graph)
barplot(ewp,showCategory=enrichmin)
dev.off()
graph <- paste(prefix,'enrichPathway_dot.pdf',sep='_')
pdf(graph,w=12, h=8)
print(graph)
dotplot(ewp,showCategory=enrichmin)
dev.off()
graph <- paste(prefix,'enrichPathway_net.pdf',sep='_')
pdf(graph,w=12, h=8)
print(graph)
cnetplot(ewp, node_label="category")
dev.off()
graph <- paste(prefix,'enrichPathway_heatmap.pdf',sep='_')
pdf(graph,w=12, h=8)
print(graph)
heatplot(ewp,foldChange=diff,showCategory=enrichmin)
dev.off()
##绘制前10个pathway
if(gseamin == 0){
	print('No GSEA result!')
}else{
graph <- paste(prefix,'GSEAPathway.pdf',sep='_')
pdf(graph,w=12, h=8)
for(i in 1:gseamin){
p<-gseaplot2(ewp2, ewp2[i,1], title=ewp2[i,2])
print(p)
}
dev.off()
print("GSEA analysis finish!")
}
setwd(outdir)
