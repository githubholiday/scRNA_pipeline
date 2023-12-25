#名称：heatmap_pheatmap.r
#作者：解飞
#邮箱：feixie@genome.cn
#时间：20190515
#版本：v0.0.1
#用途：利用pheatmap包绘制热图
###说明：
#程序开发环境 R-3.2.2
#主要功能为绘制热图
#通过输入参数可以自定义是否对数据进行标准化，是否按行或列进行聚类等
#输入文件要求有行名，列名
#===========================================================
library('getopt')
para<- matrix(c(
	'help',	'h',	0,	"logical",
	'infile',	'i',	1,	"character",
	'outfile',	'o',	1,	"character",
	'scale',	's',	2,	"character",
	'color',	'l',	2,	"character",
	'merge',	'm',	2,	"character",
	'poltfile',	'f',	2,	"character",
	'cclst',	'C',	2,	"logical",
	'rclst',	'R',	2,	"logical",
	'cname',	'c',	2,	"logical",
	'rname',	'r',	2,	"logical",
	'Canno',	'a',	2,	"character",
	'Ranno',	'n',	2,	"character",
	'order',	'd',	2,	"character",
	'rgaps',	'p',	2,	"character",
	'cgaps',	'G',	2,	"character",
	'wcell',	'W',	2,	"numeric",
	'hcell',	'H',	2,	"numeric",
	'ctree',	't',	2,	"numeric",
	'rtree',	'T',	2,	"numeric",
	'fsize',	'F',	2,	"numeric",
	'border',	'B',	2,	"logical"
),byrow=TRUE,ncol=4)

opt <- getopt(para,debug=FALSE)
#===========================================================
print_usage <- function(para=NULL){
	cat(getopt(para,usage=TRUE))
	q(status=1)	
}

#===========================================================
if ( !is.null(opt$help) )	{ print_usage(para) }
if ( is.null(opt$infile) )	{ cat("Please input the data file to draw plot...\n\n") ; print_usage(para)} #作图文件
if ( is.null(opt$outfile) )	{ opt$outfile <- c("heatmap.pdf") }# 输出文件
if ( is.null(opt$scale))	{ opt$scale <- c("row") }	#归一化方式
if ( is.null(opt$color))	{ opt$color <- c("blue,white,red") } #作图颜色
#cname <- ifelse(is.null(opt$cname),F,T) #显示列名
#rname <- ifelse(is.null(opt$rname),F,T) # 显示行名
#cclst <- ifelse(is.null(opt$cclst),FALSE,TRUE) # 列聚类
#rclst <- ifelse(is.null(opt$rclst),FALSE,TRUE) # 行聚类
#wcell <- ifelse(is.null(opt$wcell),15,opt$wcell) #cell 宽
#hcell <- ifelse(is.null(opt$hcell),10,opt$hcell) #cell 高
#ctree <- ifelse(is.null(opt$ctree),40,opt$ctree) #列树高
#rtree <- ifelse(is.null(opt$rtree),60,opt$rtree) #列树高
#fsize <- ifelse(is.null(opt$fsize),8,opt$fsize) #字体大小
#border <-ifelse(is.null(opt$border),FALSE,TRUE) #是否显示每个热图格子的边框
#===========================================================
library(pheatmap)
library(RColorBrewer)
library(reshape2)

data <- read.table(file=opt$infile,header=T,row.names=1,sep='\t',stringsAsFactors = FALSE,  check.names = FALSE,quote = "")

data = as.matrix(data)

# 如果有分组文件
if (!is.null(opt$merge)){	
	merge <- read.table(file=opt$merge,header=F,sep='\t',stringsAsFactors = FALSE,  check.names = FALSE,quote = "")
	merge <- merge[,c(ncol(merge)-1,ncol(merge))]
	colnames(merge) <- c("sample","cluster")
	ex <- data.frame(rownames(data))
	colnames(ex) <- "name"
	a <- as.vector(unique(merge$cluster))
	for (l in a){
		if(length(merge[merge$cluster==l,]$sample)>1){
			b <- data.frame(apply(data[,as.vector(merge[merge$cluster==l,]$sample)],1,mean))
		}else{
			b <-data.frame(data[,as.vector(merge[merge$cluster==l,]$sample)])
		}
		colnames(b) <- l
		b$name <- rownames(b)
		ex <- merge(ex,b,by="name",sort=F)
	}
	rownames(ex) <- ex$name
	ex$name <- NULL
	data=ex
}	


# 判断文件中是否有NA
if ( any(is.na(data)) ){
	cat("There is NA in infile. Please check your data. \n");
	q(status=1)
}



# 判断一行值或一列值是否完全相同，不能归一化
if (!is.null(opt$scale)){
	if(opt$scale !="none"){
		if(opt$scale=="row"){
			C <-apply(data,1,function(x) sd(x))
			if(!all(C>0)){
				num <- which(C==0,arr.ind = TRUE)
				data <- data[-num,]
				cat("第",num,"行数据完全相同，已删除并进行后续绘图。\n")
#				cat("某一行数据完全相同，无法按行归一化，请选择其他归一化方式！ \n");
#				q(status=1)
			}
		}
		if(opt$scale=="column"){
			C <-apply(data,2,function(x) sd(x))
			if(!all(C>0)){
				num <- which(C==0,arr.ind = TRUE)
				data <- data[,-num]
				cat("第",num,"列数据完全相同，已删除并进行后续绘图。\n")
#				cat("某一列数据完全相同，无法按列归一化，请选择其他归一化方式！ \n");
#				q(status=1)
			}
		}
	}
}

cnum <- ncol(data)# 列数
rnum <- nrow(data)#行数
color = strsplit(opt$color,split=',')
color = unlist(color)
colorRamp = colorRampPalette(color)(256)

#处理gaps 参数
rgaps <-c()
cgaps <-c()
if(!is.null(opt$rgaps)){
	rgaps = unlist(strsplit(opt$rgaps,split=',|，'))
	rgaps = as.numeric(rgaps)	
}
if(!is.null(opt$cgaps)){
	cgaps =unlist(strsplit(opt$cgaps,split=',|，')) # 列gaps
	cgaps =as.numeric(cgaps)
}
if(cnum==1){cgaps=c()}
if(rnum==1){rgaps=c()}


#处理样本顺序，返回排序后的data
getOrder<-function() { 
	list <- read.table(file=opt$order,header=F,stringsAsFactors = FALSE,  check.names = FALSE,quote = "")
	list <- list$V1
	data <- data[,list]
	return (data)
}

#处理列注释文件
getColAnno <-function(){
	#两列或多列，第一列为样本，后面每列是注释信息
	#	##      	 anno1	anno2
	#	## Test1      CT1    1
	#	## Test2      CT2    2
	#	## Test3      CT1    3
	#	## Test4      CT2    4
	#	## Test5      CT1    5
	#	## Test6      CT2    1
	anno <- read.table(file=opt$Canno,header=T,sep='\t',stringsAsFactors = FALSE,row.names=1,check.names = FALSE,quote = "")
	return(anno)
}

#处理行注释文件
getRowAnno <-function(){
	#两列，第一列为行名，第二列是注释信息
	#   ##           pathway
	#   ## Gene1      A
	#   ## Gene2      A
    #   ## Gene3      A
    #   ## Gene4      A
    #   ## Gene5      B
    #   ## Gene6      C
	Ranno <- read.table(file=opt$Ranno,header=T,sep='\t',stringsAsFactors = FALSE,row.names=1,check.names = FALSE,quote = "")
	return (Ranno)
}

# 保存画图的数据顺序
getPoltReorder <-function(cluster_col=FALSE,cluster_rows=FALSE){
	# 不对列做聚类
	if(cluster_col==FALSE){
		order_col = seq(from = 1,to=cnum,by=1)
	}else{
		order_col = p$tree_col$order    #记录热图的列排序
	}
	# 不对行做聚类
	if(cluster_rows==FALSE){
		order_row = seq(from = 1,to=rnum,by=1)
	}else{
		order_row = p$tree_row$order  #记录热图的行排序
	}
	datat = data.frame(data[order_row,order_col])   # 按照热图的顺序，重新排原始数据
	datat = data.frame(rownames(datat),datat,check.names =F)  # 将行名加到表格数据中
	colnames(datat)[1] = "ID"
	write.table(datat,file=opt$poltfile ,row.names=FALSE,quote = FALSE,sep='\t')  #输出结果，按照热图中的顺序
}

if( !is.null(opt$order)){
	#有控制样本顺序文件传入
	data<-getOrder() #按顺序重排
	cclst=F #列上不做聚类，按排序好的顺序输出
}else{
	cclst=opt$cclst	
}

#unlink('Rplots.pdf')
if( !is.null(opt$Canno) && !is.null(opt$Ranno) ){
	#行和列注释文件都存在
	canno <-getColAnno()
	ranno <-getRowAnno()
	p = pheatmap(data,color=colorRamp, scale=opt$scale,cluster_cols=cclst, cluster_rows=opt$rclst, show_rownames=opt$rname, show_colnames=opt$cname,treeheight_col=opt$ctree, treeheight_row=opt$rtree, fontsize=opt$fsize, filename=opt$outfile,gaps_row =rgaps, gaps_col =cgaps, annotation_col = canno,annotation_row = ranno,border=opt$border)
}else if(!is.null(opt$Canno) && is.null(opt$Ranno)){
	#只有列注释文件
	canno <-getColAnno()
	p = pheatmap(data,color=colorRamp, scale=opt$scale, cluster_cols=cclst, cluster_rows=opt$rclst, show_rownames=opt$rname, show_colnames=opt$cname,  treeheight_col=opt$ctree, treeheight_row=opt$rtree, fontsize=opt$fsize, filename=opt$outfile,annotation_col = canno,gaps_row =rgaps, gaps_col =cgaps,border=opt$border)
}else if(is.null(opt$Canno) && !is.null(opt$Ranno)){
	#只有行注释文件
	ranno <-getRowAnno()
	p = pheatmap(data, scale=opt$scale,color=colorRamp, cluster_cols=cclst, cluster_rows=opt$rclst, show_rownames=opt$rname, show_colnames=opt$cname,  treeheight_col=opt$ctree, treeheight_row=opt$rtree, fontsize=opt$fsize, filename=opt$outfile,annotation_row = ranno,gaps_row =rgaps, gaps_col =cgaps,border=opt$border)
}else{
	#不做注释
	p = pheatmap(data, scale=opt$scale, color=colorRamp,cluster_cols=cclst, cluster_rows=opt$rclst, show_rownames=opt$rname, show_colnames=opt$cname, treeheight_col=opt$ctree, treeheight_row=opt$rtree, fontsize=opt$fsize, filename=opt$outfile,gaps_row =rgaps, gaps_col =cgaps,border=opt$border)	
}

print(opt$outfile)
#输出作图的数据顺序
getPoltReorder(cluster_col=cclst, cluster_rows=opt$rclst)
dev.off(opt$outfile)
