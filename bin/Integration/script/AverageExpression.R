###Seurat云工具
library(Seurat)
library(ggplot2)
library(getopt)

args = commandArgs(trailingOnly=T)
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.path <- normalizePath(dirname(script))

options = matrix(c(
  "help","h","0","logical","help",
  "file_rds","f","1","character","input_file_rds",
  "outdir","o","1","character","outdir",
  "gene_list","l","2","character","gene to show",                            #参数类型：0：不需要参数，1：必须提供参数，2：可选参数
  "scale","s","2","integer","scale or not"
),ncol = 5, byrow = T)                                                        #matrix中5列分别为（长参数，短参数，参数类型，参数数据类型，参数说明）

opts = getopt(options)

##-----------检查参数是否给定
if (!is.null(opts$help) || is.null(opts$file_rds) || is.null(opts$outdir)){     #在输入help或没有给运行参数时，不执行脚本内容
  cat(paste(getopt(options, usage = T), "\n"))
  q()
}
##-----------检查参数是否正确
if (file.exists(opts$file_rds)){
  time=Sys.time()
  mess=paste(format(time)," - AverageExpression - INFO - the file_rds has been detected")
  cat(mess,"\n")
}else{
  time=Sys.time()
  mess=paste(format(time)," - AverageExpression - INFO - the file_rds has not been detected, please check your input")
  cat(mess,"\n")
}

if (dir.exists(opts$outdir)){
  time=Sys.time()
  mess=paste(format(time)," - AverageExpression - INFO - the dir has been detected")
  cat(mess,"\n")
}else{
  time=Sys.time()
  mess=paste(format(time)," - AverageExpression - INFO - the dir has not been detected, please check your input")
  cat(mess,"\n")
}

##############主体
seurat_obj = readRDS(opts$file_rds)       #读取rds文件，该文件来自于10x项目的中间文件

cluster.averages = AverageExpression(seurat_obj,return.seurat = TRUE)           #计算结果存于cluster.averages[['RNA']]@data中，是matrix

#用标准化后的scale.data(solt默认)降维后绘制
if (!is.null(opts$gene_list)){
  if (file.exists(opts$gene_list)){
    temp_fea = read.table(opts$gene_list,sep = ',')
    fea = array(unlist(temp_fea))}else{
	 cat('list not found, please check\n')
	 q()}
}else{
  fea = unlist(TopFeatures(seurat_obj[['pca']],balanced = TRUE))
}
prefix = substr(basename(opts$file_rds),start=1,nchar(basename(opts$file_rds))-4)
out_csv = paste(opts$outdir,"/AverageExpression.xls",sep = "")
out_pdf = paste(opts$outdir,"/AverageExpression_heatmap.pdf",sep = "")
pdf(file = out_pdf)
if (opts$scale == 'yes'){
  table_out = GetAssayData(object = cluster.averages[['RNA']],slot = 'scale.data')
  #temp_name = paste(prefix,'_cluster_',colnames(table_out),sep = "")
  #colnames(table_out) = temp_name

  write.table(table_out,out_csv,quote = FALSE,sep="\t")
  DoHeatmap(cluster.averages, features = fea,size = 3,draw.lines = FALSE) + scale_fill_gradientn(colors = c("blue", "white", "red"))
}else if(opts$scale == 'no'){
  table_out = GetAssayData(object = cluster.averages[['RNA']],slot = 'data')
  #temp_name = paste(prefix,'_cluster_',colnames(table_out),sep = "")
  #colnames(table_out) = temp_name

  write.table(table_out,out_csv,quote = FALSE,sep="\t")
  DoHeatmap(cluster.averages, slot='data', features = fea,size = 3,draw.lines = FALSE) + scale_fill_gradientn(colors = c("blue", "white", "red"))
}
dev.off()






#传入参数要给一个feature

# cluster.averages                                                                #查看基本情况，是一个seurat类，包含一个assay，有19800个特征和18个样本也就是簇
# 
# cluster.averages@assays                                                         #查看assay数据结构里面有RNA的一条信息
# 
# cluster.averages@assays$RNA                                                     #RNA里有一个data，包含有19800个特征也就是基因，然后分成18类细胞簇
# 
# head(cluster.averages@assays$RNA@counts)                                        #每cluster的细胞数，.表示0
# head(cluster.averages@assays$RNA@data)                                          #每cluster的表达量
# head(cluster.averages@assays$RNA@scale.data)                                    #标准化后的表达量
# head(cluster.averages@assays$RNA@key)
# head(cluster.averages@assays$RNA@assay.orig)
# head(cluster.averages@assays$RNA@var.features)
# head(cluster.averages@assays$RNA@meta.features)
# head(as.matrix(cluster.averages@assays$RNA@meta.features))
# 
# GetAssayData(object = cluster.averages[['RNA']],slot = 'data')[1:10,]
# 
# CellScatter(cluster.averages, cell1 = 1, cell2 = 10)
# ?CellScatter
