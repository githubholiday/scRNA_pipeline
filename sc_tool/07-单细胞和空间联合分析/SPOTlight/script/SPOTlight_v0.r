#!/annoroad/data1/software/install/Miniconda-2023/envs/R4.0/bin/Rscript(203使用)
#作者：张阳
#邮箱：yangzhang@genome.cn
#时间：2023年05月17日 星期三 11时10分59秒
#版本：v0.0.1
#用途：
library('getopt')
para= matrix(c(
    'help',    'h',    0,  "logical",
    'scdata',   'c',    1,  "character",
    'stdata',   't',    1,  "character",
    'prefix',    'p',    1,  "character",
    'configini',   'f',    1,  "character",
    'scriptdir',   's',    1,  "character",
    'outdir',   'o',    1,  "character"
),byrow=TRUE,ncol=4)
opt = getopt(para,debug=FALSE)
print_usage = function(para=NULL){
    cat(getopt(para,usage=TRUE))
    cat("
    Options:
    help    h   NULL        get this help
    scdata  c   character   单细胞10X的rds，需要带有注释
    stdata  t character 空间转录组的rds
    outdir  o character 输出路径
    prefix  p   character   输出结果的前缀
    configini  f character 参数配置文件
    scriptdir  s character source的脚本所在路径
    \n")
    q(status=1)
}

if ( !is.null(opt$help) )	{ print_usage(para) }
if ( is.null(opt$scdata) )	{ cat("Please input the RNA rds data file ...\n\n") ; print_usage(para)}
if ( is.null(opt$stdata) )	{ cat("Please input the space rds data file ...\n\n") ; print_usage(para)}
if ( is.null(opt$prefix) )	{ cat("Please input the prefix for analysis ...\n\n") ; print_usage(para)}
if ( is.null(opt$outdir) )	{ cat("Please give the outdir for analysis ...\n\n") ; print_usage(para)}
if ( is.null(opt$configini) )	{ cat("Please give the configini file for analysis ...\n\n") ; print_usage(para)}
if ( is.null(opt$scriptdir) )	{ cat("Please give the configini file for analysis ...\n\n") ; print_usage(para)}

library(ggplot2)
library(SPOTlight)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)
library(Seurat)
library(dplyr)
library(NMF)
library(ggsci)
library(configr)

###创建目录
mkdirs = function(outdir,fp) {
	if(!file.exists(file.path(outdir,fp))) {
		dir.create(file.path(outdir,fp))
	}else{
		print(paste(fp,"Dir already exists!",sep="     "))
		}
}

scdata = opt$scdata
stdata = opt$stdata
outdir = opt$outdir
prefix = opt$prefix
mkdirs(outdir,prefix)
scriptdir = opt$scriptdir
ini = opt$configini
ini.list = read.config(file = ini)
setwd(paste(outdir,prefix,sep='/'))
cell_number = ini.list$Para$cellnumber
gene_number = ini.list$Para$genenumber
celltype=ini.list$Para$RNAcelltype
seurat_clusters=ini.list$Para$spacecluster
repeat_gene_num=as.integer(ini.list$Para$repeat_gene_num)
mgs_col=ini.list$Para$mgs_col
mgs_col_threshold=as.integer(ini.list$Para$mgs_col_threshold)
# 一、检查celltype是否存在rds中
stdata = readRDS( stdata )
scdata = readRDS( scdata )
if(is.null(scdata@meta.data[,celltype])){
	print(paste("celltype变量在",opt$RNArds,"文件中不存在，流程退出",sep="_"))
	q();
}
#检查seurat_clusters是否存在spacerds中
if(is.null(stdata@meta.data[,seurat_clusters])){
	print(paste("seurat_clusters变量在",opt$spacerds,"文件中不存在，流程退出",sep="_"))
	q();
}
#检查两者的基因名字是否有高于3000的重复
gene_space=row.names(stdata@assays$Spatial@counts)
gene_RNA=row.names(scdata@assays$RNA@counts)
common_gene=intersect(gene_space,gene_RNA)
if(length(common_gene) > repeat_gene_num){
	print(paste("RNA和空间数据集中的gene重复超过 ", repeat_gene_num , " 个，继续分析",sep=""))
}else{
	print(paste("RNA和空间数据集中的gene重复数为" , length(common_gene) , "流程退出",sep=""))
	q();
}

# 二、空转数据情况描述
Idents(stdata)=stdata@meta.data[,seurat_clusters]
pdf( paste( prefix , "_spatial_cluster_umap.pdf" , sep = '') , width = 12 , height = 8)
p1 = DimPlot( stdata, reduction = "umap", label = TRUE )
p2 = SpatialDimPlot( stdata, label = TRUE, label.size = 3 )
p1 + p2
dev.off()

# 三、单细胞数据情况描述
Idents(scdata)=scdata@meta.data[,celltype]
sc = Seurat::SCTransform(scdata,verbose = FALSE) %>% Seurat::RunPCA(., verbose = FALSE) %>% Seurat::RunUMAP(., dims = 1:30,verbose = FALSE)
pdf(paste( prefix , "_singlecell_cluster_umap.pdf" , sep = '' ))
Seurat::DimPlot(sc , group.by = celltype , label = TRUE ) + Seurat::NoLegend()
dev.off()
### 统计细胞数量
cell_num = as.data.frame(table(scdata@meta.data[,celltype]))
names(cell_num) = c("celltype","count")
cell_num = cell_num[order(cell_num$count),]
write.table(cell_num,paste(prefix , "_celltype_count.xls" , sep='' ),sep="\t",quote=F,row.names=F)

# 四、联合分析
sce = as.SingleCellExperiment(scdata , assay='RNA')
sce = logNormCounts(sce)
dec = modelGeneVar(sce)  
# 挑选前3000个基因
hvg = getTopHVGs(dec, n = as.integer(gene_number))
genes = rownames(sce)
colLabels(sce) = colData(sce)[,celltype]
mgs = scoreMarkers(sce, subset.row = genes)
mgs_fil = lapply(names(mgs), function(i) {
    x = mgs[[i]] 
    x = x[x[,mgs_col] > mgs_col_threshold, ] 
    x = x[order(x[,mgs_col], decreasing = TRUE), ]
    x$gene = rownames(x)
    x$cluster = i
    data.frame(x)
})
mgs_df = do.call(rbind, mgs_fil)
idx = split(seq(ncol(sce)), colLabels(sce))
n_cells = as.integer(cell_number)
cs_keep = lapply(idx, function(i) {
    n = length(i)
    if (n < n_cells)
        n_cells = n
    sample(i, n_cells)
})
sce = sce[, unlist(cs_keep)]
# 解卷积
sc_dat = sce@assays@data$counts
sp_dat = stdata@assays$Spatial@counts
res = SPOTlight(
    x = sc_dat,
    y = sp_dat,
    groups = as.character(colLabels(sce)),
    mgs = mgs_df,
    hvg = hvg,
    weight_id = mgs_col,
    group_id = "cluster",
    gene_id = "gene")
saveRDS(res,paste( prefix , "_res.rds" ,sep = '' ))
mat = res$mat
head(mat)[, seq_len(3)]

# 结果1、spot中的细胞类型占比统计，整体结果的饼图展示
source(paste( scriptdir , "/plotSpatialScatterpie.R" , sep = ''))
source(paste( scriptdir , "/utils.R" , sep = ''))
spot_data = data.frame(spot=rownames(mat),mat)
write.table(spot_data , file=paste(prefix,"_all_spot_celltype.xls",sep='') ,sep="\t",quote=F,row.names=F)
spot_data_tmp = t(spot_data[1,])
write.table(spot_data_tmp ,file=paste(prefix,"_all_spot_celltype_example.xls",sep='') ,sep="\t",quote=F,row.names=T,col.names=F)
ct = colnames(mat)
mat[mat < 0.1] = 0
paletteMartin = c("#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d")
pal = colorRampPalette(paletteMartin)(length(ct))
names(pal) = ct
pdf(paste(prefix,"_SpatialScatterpie.pdf",sep=''))
plotSpatialScatterpie(
    x = stdata,
    y = mat,
    img = TRUE,
    cell_types = colnames(mat), # 默认是y的列名，至少两个，类型为字符串
    scatterpie_alpha = 1,  # 1-10没有变化
    pie_scale = 0.5) +    # 控制每一个饼的大小，数值变大容易叠在一起
    scale_fill_manual(
        values = pal,
        breaks = names(pal))
dev.off()

# 结果3、单cluter细胞类型饼图
cluster_name = unique(stdata@meta.data[,seurat_clusters])

for (i in cluster_name){
    stdata_n = subset(stdata,idents=i)
    mat_n = mat[rownames(stdata_n@meta.data),]
    plotSpatialScatterpie(  
        x = stdata_n, 
        y = mat_n,  
        img = TRUE, 
        cell_types = colnames(mat_n), 
        scatterpie_alpha = 1,  
        pie_scale = 0.9) +    
        scale_fill_manual(
            values = pal,
            breaks = names(pal)) +
    ggtitle(paste( seurat_clusters , i , sep = '_')) +
    theme(plot.title = element_text(hjust=0.5))
    ggsave(paste(prefix , "_SpatialScatterpie_",i,".pdf" , sep=''))
}

# 结果4、各细胞类型箱线图
cell_type = c()
mat = res[[1]]
for (i in 1:nrow(mat)){
    tmp = mat[i,]
    max_num = max(tmp)
    type = names( tmp[which(tmp == max_num)])
    cell_type[i] = type
}
cell_type = as.data.frame(cell_type)
result = cbind(mat,cell_type)
cluster = as.data.frame(stdata@meta.data[,seurat_clusters])
colnames(cluster)="cluster"
rownames(cluster)=rownames(stdata@meta.data)
cluster_cell = merge(result, cluster ,by="row.names" , all=T)
colors = c(pal_aaas("default", alpha = 0.8)(10),pal_lancet("lanonc", alpha = 0.8)(9),pal_nejm("default", alpha = 0.8)(8),pal_jama("default", alpha = 0.8)(7),pal_jco("default", alpha = 0.8)(10))
group_color = colors[1:length(cluster_name)]
for (i in ct){
    cluster_cell_i = cluster_cell[,c("cluster",i)]
    colnames(cluster_cell_i) = c("cluster","i")
    p = ggplot(cluster_cell_i, aes(x=cluster, y=i ,color=cluster)) + geom_boxplot() + 
        geom_jitter(shape=16, position = position_jitter(0.2)) +
        scale_color_manual(values = group_color) +  ylab(i) +
        theme_bw()
    ggsave(paste( prefix , "_boxplot_" , i , ".pdf" , sep = '') , p ) 
}

# 结果5、细胞类型相关性矩阵
pdf(paste(prefix,"_CorrelationMatrix.pdf",sep=''))
plotCorrelationMatrix(mat,tl.cex=10)
dev.off()
mat = res[[1]]
mat = mat[, colSums(mat)>0]
mat_cor = cor(mat)
write.table(data.frame(sample=rownames(mat_cor),mat_cor),file=paste(prefix,"_CorrelationMatrix.xls",sep='') ,sep="\t",quote=F,row.names=F)

# 结果6、各细胞类型共位情况-加权网络图
pdf(paste(prefix,"_Interactions.pdf",sep=''))
source(paste(scriptdir, "/plotInteractionsold.R" , sep=''))
graph_ntw = get_spatial_interaction_graph(res[[1]])
edge_importance = E(graph_ntw)$importance
E(graph_ntw)[[]]
plot(graph_ntw,
     edge.width = edge_importance,  # 边的粗细。
     edge.color = "#cde394",  # 边的颜色 
     vertex.color = "#fcb857",  # 圈的颜色
     vertex.frame.color = "white",  # 圈的框的颜色
     vertex.label.color = "black",  # 圈上细胞类型文字的颜色
     layout = layout.circle)
dev.off()

# 结果2、单细胞类型占比组合图
stdata@meta.data = data.frame(res[[1]])
cell_types_all = colnames(stdata@meta.data)
for (i in cell_types_all){
    # pdf( paste( prefix , "_single_celltype_" , i ,".pdf" , sep = '') , width = 12, height = 8 )
    p1 = Seurat::SpatialDimPlot(stdata, label = TRUE, label.size = 3)
    p2 = Seurat::SpatialFeaturePlot(stdata , features = i , alpha = c(0.1, 1) )
    p1 + p2
    # dev.off()
    ggsave(paste( prefix , "_single_celltype_" , i ,".pdf" , sep = '') , p1 + p2 , width = 12, height = 8 )
}
