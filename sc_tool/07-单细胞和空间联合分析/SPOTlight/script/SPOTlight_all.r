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
    'cellnumber',    'n',    1,  "logical",
    'genenumber',   'g',    1,  "character",
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
    cellnumber  n   character   解卷积时每个细胞类型用到的代表细胞数量
    genenumber  g character 选取的高可变基因的数量
    scriptdir  s character source的脚本所在路径
    \n")
    q(status=1)
}
if (is.null(opt$scdata)) {print_usage(para)}

scdata = opt$scdata
stdata = opt$stdata
outdir = opt$outdir
cell_number = opt$cellnumber
gene_number = opt$genenumber
scriptdir = opt$scriptdir

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

# 一、加载本地空转数据
stdata = readRDS( stdata )
pdf( paste( outdir , "/spatial_cluster_umap.pdf" , sep = '') )
p1 <- DimPlot( stdata, reduction = "umap", label = TRUE )
p2 <- SpatialDimPlot( stdata, label = TRUE, label.size = 3 )
p1 + p2
dev.off()
# 二、加载本地单细胞数据
scdata = readRDS( scdata )
sce = Seurat::SCTransform(scdata,verbose = FALSE) %>% Seurat::RunPCA(., verbose = FALSE) %>% Seurat::RunUMAP(., dims = 1:30,verbose = FALSE)
pdf(paste( outdir , "/singlecell_cluster_umap.pdf" , sep = '' ))
Seurat::DimPlot(sce , group.by = "celltype" , label = TRUE ) + Seurat::NoLegend()
dev.off()
# 统计细胞数量
table(sce$celltype) %>% write.table(.,file=paste(outdir,"/cell_number_stat.xls",sep=''),quote=T,sep="\t",row.names=F,col.names=F)

# 三、联合分析
# 1、10x单细胞数据从seurat转成SingleCellExperiment类型数据
sce=as.SingleCellExperiment(scdata,assay='RNA')
# 标准化，方法来自scater包
sce = logNormCounts(sce)
# 模型拟合获得高可变基因，方法来自scran包
dec = modelGeneVar(sce)  
# 挑选前3000个基因
hvg = getTopHVGs(dec, n = gene_number)
# 对marker基因打分
genes = rownames(sce)
colLabels(sce) = colData(sce)$celltype
mgs = scoreMarkers(sce, subset.row = genes)
mgs_fil = lapply(names(mgs), function(i) {
    x = mgs[[i]] # x是该细胞类型下的打分data.frame
    # Filter and keep relevant marker genes, those with AUC > 0.8
    x = x[x$mean.AUC > 0.8, ] # 选其中的一种打分进行筛选
    # Sort the genes from highest to lowest weight
    x = x[order(x$mean.AUC, decreasing = TRUE), ]
    # Add gene and cluster id to the dataframe
    x$gene = rownames(x)
    x$cluster = i
    data.frame(x)
})
# 将mgs_fil这个list合并成为数据框
mgs_df = do.call(rbind, mgs_fil)
# 获取细胞类型及其分组的列表
idx = split(seq(ncol(sce)), sce$celltype)
n_cells = cell_number
cs_keep = lapply(idx, function(i) {
    n = length(i)
    if (n < n_cells)
        n_cells = n
    sample(i, n_cells)
})
sce = sce[, unlist(cs_keep)]
# 筛选完细胞之后，进行联合分析。其中对于各自数据的要求如下，即取出表达矩阵
sc_dat = sce@assays@data$counts
sp_dat = stdata@assays$Spatial@counts
res = SPOTlight(
    x = sc_dat,
    y = sp_dat,
    groups = as.character(sce$celltype),
    mgs = mgs_df,
    hvg = hvg,
    weight_id = "mean.AUC",
    group_id = "cluster",
    gene_id = "gene")
saveRDS(res,paste( outdir , "/res.rds" ,sep = '' ))
# res是一个list，有3个： "mat"    "res_ss" "NMF"
head(mat = res$mat)[, seq_len(3)]
write.table(data.frame(spot=rownames(mat),mat),file=paste(outdir,"/all_spot_celltype.xls",sep='') ,sep="\t",quote=T,row.names=F)
mod = res$NMF
# 联合分析结束，后续为图表展示
# 1、topic profile
pdf(paste(outdir,"/topic_profile.pdf",sep=''))
plotTopicProfiles(
    x = mod,
    y = sce$celltype,
    facet = FALSE,
    min_prop = 0.01,
    ncol = 1) +
    theme(aspect.ratio = 1)
dev.off()

pdf(paste(outdir,"/cell_identity_topic_profile.pdf",sep=''))
plotTopicProfiles(
    x = mod,
    y = sce$celltype,
    facet = TRUE,
    min_prop = 0.01,
    ncol = 6)
dev.off()

sign = basis(mod)
colnames(sign) = paste0("Topic", seq_len(ncol(sign)))
head(sign)
write.table(data.frame(gene=rownames(sign),sign),file=paste(outdir,"/all_topic_profile.xls",sep='') ,sep="\t",quote=T,row.names=F)
# 2、细胞类型相关性矩阵
pdf(paste(outdir,"/CorrelationMatrix.pdf",sep=''))
plotCorrelationMatrix(mat,tl.cex=10)
dev.off()
mat = res[[1]]
mat = mat[, colSums(mat)>0]
mat_cor = cor(mat)
write.table(data.frame(sample=rownames(mat_cor),mat_cor),file=paste(outdir,"/CorrelationMatrix.xls",sep='') ,sep="\t",quote=T,row.names=F)
# 3、各细胞类型共位情况-加权网络图
pdf(paste(outdir,"/Interactions.pdf",sep=''))
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
# 4、空间细胞类型饼图
ct = colnames(mat)
mat[mat < 0.1] = 0
paletteMartin = c("#000000", "#004949", "#009292", "#ff6db6", "#ffb6db", "#490092", "#006ddb", "#b66dff", "#6db6ff", "#b6dbff", "#920000", "#924900", "#db6d00", "#24ff24", "#ffff6d")
pal = colorRampPalette(paletteMartin)(length(ct))
names(pal) = ct
pdf(paste(outdir,"/SpatialScatterpie.pdf",sep=''))
plotSpatialScatterpie(
    x = stdata,
    y = mat,
    cell_types = colnames(mat), # 默认是y的列名，至少两个，类型为字符串
    scatterpie_alpha = 1,  # 1-10没有变化
    pie_scale = 0.2) +    # 控制每一个饼的大小，数值变大容易叠在一起
    scale_fill_manual(
        values = pal,
        breaks = names(pal))
dev.off()
# 5、验证结果
#stdata$res_ss = res[[2]][colnames(stdata)]
#xy = SpatialExperiment::spatialCoords(stdata)
#stdata$x = xy[, 1]
#stdata$y = xy[, 2]
#pdf("5-RSS.pdf")
#ggcells(stdata, aes(x, y, color = res_ss)) +
#    geom_point() +
#    scale_color_viridis_c() +
#    coord_fixed() +
#    theme_bw()
#dev.off()

# 6、单细胞类型百分比展示
stdata@meta.data = data.frame(res[[1]])
cell_types_all = colnames(stdata@meta.data)
for (i in 1:cell_types_all){
    pdf( paste( outdir , "/single_celltype_" , i ,".pdf" , sep = '') , width = 12, height = 8 )
    # pdf(paste( outdir , "/single-celltype.pdf" , sep = '') , height=8*length(cell_types_all) )
    p1 = Seurat::SpatialDimPlot(stdata, label = TRUE, label.size = 3)
    p2 = Seurat::SpatialFeaturePlot(stdata , features = i , alpha = c(0.1, 1) )
    p1 + p2
    dev.off()
}

# 7、细胞在空间上的分布饼图（单个展示）
source(paste( outdir , "/plotSpatialScatterpie.R" , sep = ''))
source(paste( outdir , "/utils.R" , sep = ''))
cluster_num = length(stdata$seurat_clusters[!duplicated(stdata$seurat_clusters)])
Idents(stdata)=stdata$seurat_clusters
for (i in 1:cluster_num){
    stdata_n = subset(stdata,idents=c(i))
    mat_n = mat[rownames(stdata_n@meta.data),]
    plotSpatialScatterpie(  # 只能依赖于source源代码。目前不清楚原因。
        x = stdata_n, 
        y = mat_n,  
        img = TRUE, # 用于加底图
        cell_types = colnames(mat_n), 
        scatterpie_alpha = 1,  
        pie_scale = 0.5) +    
        scale_fill_manual(
            values = pal,
            breaks = names(pal))
    ggsave(paste(outdir , "/SpatialScatterpie_",i,".pdf" , sep=''))
}

# 8、不同细胞类型中cluster的统计及箱线图
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
cluster = as.data.frame(stdata$seurat_clusters)
colnames(cluster)="cluster"
cluster_cell = merge(result, cluster ,by="row.names" , all=T)
write.table( cluster_cell,file = paste(outdir , "/cluster_cell.xls" , sep = ''),sep="\t" , quote=T , row.names=F )
stat = as.data.frame(table(cluster_cell[,c("cell_type","cluster")]))
stat_wide = dcast(stat,cell_type~cluster)
write.table(stat_wide , file = paste( outdir , "/cluster_cell_stat.xls" , sep = '') , sep="\t" , quote=T , row.names=F)
pdf(paste(outdir , "/cluster_cell_stat_dotplot.pdf" , sep = '') , width=ncol(stat_wide) , height=nrow(stat_wide))
ggplot(stat, aes(cluster, cell_type), showCategory=8) +
	geom_point(aes(size=Freq)) +
    theme_bw() 
dev.off()
# 箱式图
colors = c(pal_aaas("default", alpha = 0.8)(10),pal_lancet("lanonc", alpha = 0.8)(9),pal_nejm("default", alpha = 0.8)(8),pal_jama("default", alpha = 0.8)(7),pal_jco("default", alpha = 0.8)(10))
group_color = colors[1:cluster_num]
for (i in colnames(mat)){
    cluster_cell_i = cluster_cell[,c("cluster",i)]
    colnames(cluster_cell_i) = c("cluster","i")
    p = ggplot(cluster_cell_i, aes(x=cluster, y=i ,color=cluster)) + geom_boxplot() + 
        geom_jitter(shape=16, position = position_jitter(0.2)) +
        scale_color_manual(values = group_color) +  ylab(i) +
        theme_bw()
    ggsave(paste( outdir , "/boxplot_" , i , ".pdf" , sep = '') , p ) 
}

# 9、展示最大比例细胞类型（总分）
# 将最大的细胞类型作为该spot的细胞类型，然后映射回stdata
stdata = readRDS( stdata )
stdata@meta.data$orig.spot = rownames(stdata@meta.data)
labers=cluster_cell[match(as.character(stdata@meta.data$orig.spot),as.character(cluster_cell$Row.names)),]$cell_type
labers=c(labers)
stdata$celltype = labers
#将所有细胞展示在切片上
Idents(stdata) = stdata$celltype
SpatialDimPlot(stdata) + theme(legend.position = "right")
ggsave(file=paste( outdir , "/SpatialDimPlot.pdf" , sep = ''),width=10, height=8)

#将每种细胞单独展示
SpatialDimPlot(stdata,cells.highlight = CellsByIdentities(object = stdata, idents = unique(Idents(stdata))),facet.highlight = TRUE)
ggsave(file=paste( outdir , "/SpatialDimPlot_Pit1_eachcell.pdf" , sep = ''),width=10, height=8)
