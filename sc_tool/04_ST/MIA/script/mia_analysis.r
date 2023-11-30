#/annoroad/data1/bioinfo/PMO/yaomengcheng/bk_Anaconda3/envs/monocle3/bin/R
#MIA-单细胞转录组和空间转录组联合分析
library('getopt')
para<- matrix(c(
	'help',	'h',	0,	"logical",
	'prefix',	'p',	1,	"character",
	'RNArds',	'r',	1,	"character",
	'outdir',	'o',	1,	"character",
	'configini',	'c',	1,	"character",
	'spacerds',	's',	1,	"character"
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
	--prefix	p	character	the prefix for outputfiles [forced]
	--RNArds	r	character	RNArdsfile[forced]
	--outdir	o	character	The	resurt of out dir for analysis [forced]
	--configini	c	character	config.ini  [forced]
	--spacerds	s	character	spacerdsfile  [forced]
	\n")
	q(status=1)
}
#===========================================================
if ( !is.null(opt$help) ) { print_usage(para) }
if ( is.null(opt$RNArds) ) { cat("Please input the RNArds data file ...\n\n") ; print_usage(para)}
if ( is.null(opt$prefix) ) { cat("Please input the prefix for analysis ...\n\n") ; print_usage(para)}
if ( is.null(opt$outdir) ) { cat("Please give the outdir for analysis ...\n\n") ; print_usage(para)}
if ( is.null(opt$configini) ) { cat("Please give the configini file for analysis ...\n\n") ; print_usage(para)}
if ( is.null(opt$spacerds) ) { cat("Please input the spacerds data file ...\n\n") ; print_usage(para)}

library(Seurat)
library(tidyverse)
library(configr)
library(cowplot)

###创建目录
mkdirs <- function(outdir,fp) {
	if(!file.exists(file.path(outdir,fp))) {
		dir.create(file.path(outdir,fp))
	}else{
		print(paste(fp,"Dir already exists!",sep="     "))
		}
}

syMIA=function(
    region_specific,
    celltype_specific,
    N,                             #所有基因的数量
    pvalue_log10_neg=20,           #p值相关的阈值
    axis.text.y.left.size=14,      #主panel纵轴的文本大小
    axis.ticks.y.left.length=0.2,  #主panel纵轴的刻度线长度
    legend.title.size=14,          #主panel的legend的title大小
    legend.text.size_uppanel=14,   #上层panel的legend的文本大小
    axis.text.y.left_uppanel=14,   #上层panel的纵轴的文本大小
    plotwidth=18,
    plotheight=25
){
  miares=data.frame()
  
  rnames = unique(as.character(celltype_specific$celltype))
  cnames = unique(as.character(region_specific$region))
  num = length(rnames) * length(cnames)
  seq_n = rep(0, num)
  my = matrix(seq_n, nrow = length(rnames), ncol = length(cnames), dimnames = list(rnames, cnames))
  same_genes = matrix(seq_n, nrow = length(rnames), ncol = length(cnames), dimnames = list(rnames, cnames))

  for (ri in unique(as.character(region_specific$region))) {
    smalldf1=region_specific[region_specific$region == ri,]
    n=length(smalldf1$gene)
    term=c()
    pvalue=c()
    
    for (ci in unique(as.character(celltype_specific$celltype))) {
      one.celltype=celltype_specific[ celltype_specific$celltype %in% ci ,]
      M=length(one.celltype$gene)

      same_gene = smalldf1$gene %in% one.celltype$gene
      k=sum(same_gene)
      same_genes[ci, ri] = paste(smalldf1$gene[same_gene],collapse=',')
      my[ci, ri] = k

      one.pvalue=phyper(k-1,M, N-M, n, lower.tail=FALSE)
      term=append(term,ci)
      pvalue=append(pvalue,one.pvalue)
    }
    
    one.miares=data.frame(region=ri,term=term,pvalue=pvalue)
    one.miares=one.miares%>%arrange(pvalue)
    if(which(unique(as.character(region_specific$region)) %in% ri) == 1) {
      miares = one.miares
    } else {
      miares = miares%>%rbind(one.miares)
    }
  }
  write.table(same_genes, paste(prefix, 'same_genes.xls', sep='_'), quote = F, sep="\t")
  for (i in rownames(same_genes)){
    sam = c()
    for (k in colnames(same_genes)){
      same_gene=same_genes[i,k]
      sam = c(sam, strsplit(same_gene,',')[[1]])
    }
    print(sam)
    stat = as.data.frame(table(sam))
    if ( dim(stat)[1] > 0){
      stat = stat[order(-as.numeric(stat$Freq)),]
      print(i)
      print(stat)
      top_plot = SpatialFeaturePlot(space_rds, features = stat$sam[0:3])
      ggsave(top_plot, file=paste(i, 'top3_gene_FeaturePlot.pdf', sep='_'), width=length(stat)*5, height=6)
    }
  }

  
  out = data.frame(Celltype=rownames(my), my)
  colnames(out) = c('Celltype', colnames(my))
  write.table(out, paste(prefix, 'region_celltype_gene.xls', sep='_'), quote = F, sep="\t", row.names=F)
  
  #添加几列信息
  miares$Enrichment= -log10(miares$pvalue)
  miares$Depletion= -log10(1-miares$pvalue)
  miares$final_value= ifelse(miares$Enrichment > miares$Depletion,miares$Enrichment,miares$Depletion)
  miares$final_class= ifelse(miares$Enrichment > miares$Depletion,"Enrichment","Depletion")
  miares$final_value[miares$final_value >= pvalue_log10_neg] = pvalue_log10_neg
  
  #画图部分
  labeldf1=as.data.frame(table(region_specific$region))
  colnames(labeldf1)=c("region","gene_num")
  labeldf1$region=as.character(labeldf1$region)
  labeldf1$region_gene_num=paste0(labeldf1$region," (",labeldf1$gene_num," genes)")
  labeldf1=labeldf1%>%arrange(region)
  
  labeldf2=as.data.frame(table(celltype_specific$celltype))
  colnames(labeldf2)=c("celltype","gene_num")
  labeldf2$celltype=as.character(labeldf2$celltype)
  labeldf2$celltype_gene_num=paste0(labeldf2$celltype," (",labeldf2$gene_num,")")
  labeldf2=labeldf2%>%arrange(celltype)
  
  miares$region=factor(miares$region,levels = sort(unique(miares$region)))
  miares$term=factor(miares$term,levels = sort(unique(miares$term)))
  miares.Enrichment=miares%>%filter(final_class == "Enrichment")
  miares.Depletion=miares%>%filter(final_class == "Depletion")
  
  library(ggnewscale)
  library(colorspace)
  pb=ggplot()+
    geom_tile(data = miares.Enrichment,mapping = aes(x=region,y=term,fill=final_value),color="black",size=0.5)+
    scale_fill_gradientn("Enrichment",colours = brewer.pal(9, "Reds"))+
    new_scale_fill() +
    geom_tile(data = miares.Depletion,mapping = aes(x=region,y=term,fill=final_value),color="black",size=0.5)+
    scale_fill_gradientn("Depletion",colours = brewer.pal(9, "Blues"))+
    scale_x_discrete(expand = c(0,0))+
    scale_y_discrete(expand = c(0,0),breaks=labeldf2$celltype,labels=labeldf2$celltype_gene_num)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text.y.left = element_text(size = axis.text.y.left.size,color = "black"),
      axis.ticks.length.y.left = unit(axis.ticks.y.left.length,"cm"),
      axis.text.x.bottom = element_blank(),
      axis.ticks.x.bottom = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = legend.title.size,vjust = 1)
    )
  
  updf=as.data.frame(levels(miares$region))
  colnames(updf)="region"
  updf$region=factor(updf$region,levels = updf$region)
  q4 <- sequential_hcl(length(updf$region), palette = "BluGrn")
  pa=ggplot(data = updf,aes(x=region,y=0,fill=region))+
    geom_tile()+
    scale_x_discrete(expand = c(0,0))+
    scale_fill_manual(values = q4,breaks = labeldf1$region,labels=labeldf1$region_gene_num)+
    scale_y_continuous(expand = c(0,0),breaks = 0,labels = "Cell types (no. of genes)")+
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = legend.text.size_uppanel),
      legend.direction = "vertical",
      
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x.bottom = element_blank(),
      axis.text.y.left = element_text(size = axis.text.y.left_uppanel,color = "black")
    )+
    guides(fill = guide_legend(override.aes = list(size=10), ncol = 2))
  
  library(patchwork)
  celltype_num=length(unique(celltype_specific$celltype))
  pa / pb + plot_layout(heights = c(1,celltype_num))
  ggsave(paste(prefix, 'MIA.pdf', sep='_'),width = plotwidth,height = plotheight,units = "cm")
  print('MIA分析运行完成')
  #返回有用的矩阵
  return(miares)
}

# 读取config.ini 文件，获取sample list 和Assay
prefix<-opt$prefix
outdir<- opt$outdir
mkdirs(outdir,prefix)
setwd(paste(outdir,prefix,sep='/'))
ini<-opt$configini
ini.list <- read.config(file = ini)


#space
space_rds<-readRDS(opt$spacerds)
#RNA
rna_rds<-readRDS(opt$RNArds)

anchors_ims<-as.numeric(ini.list$Para$anchors_ims)
celltype<-ini.list$Para$RNAcelltype
seurat_clusters<-ini.list$Para$spacecluster

#检查celltype是否存在rds中
if(is.null(rna_rds@meta.data[,celltype])){
	print(paste("celltype变量在",opt$RNArds,"文件中不存在，流程退出",sep="_"))
	q();
}
#检查seurat_clusters是否存在spacerds中
if(is.null(space_rds@meta.data[,seurat_clusters])){
	print(paste("seurat_clusters变量在",opt$spacerds,"文件中不存在，流程退出",sep="_"))
	q();
}
#检查两者的基因名字是否有高于1000的重复
gene_space<-row.names(space_rds@assays$Spatial@counts)
gene_RNA<-row.names(rna_rds@assays$RNA@counts)
common_gene<-intersect(gene_space,gene_RNA)
if(length(common_gene)>as.numeric(ini.list$Para$common_gene)){
	print("RNA和空间数据集中的gene重复超过预期值，继续分析")
}else{
	print(paste("RNA和空间数据集中的gene重复数为",length(common_gene),"流程退出",sep=""))
	q();
}

#输出每种细胞类型的细胞数量
cell_num<-as.data.frame(table(rna_rds@meta.data[,celltype]))
names(cell_num)<-c("celltype","count")
cell_num<-cell_num[order(cell_num$count),]
write.table(cell_num,paste(prefix,"celltype_count.xls",sep="_"),sep="\t",quote=F,row.names=F)

### 1.单细胞转录组流程 #####################################################
### 添加注释信息
if (is.null(rna_rds@reductions$umap)){
	print("RNA数据集中不存在umap, 需要运行以下步骤")
	rna_rds <- ScaleData(rna_rds, verbose = FALSE)
	rna_rds <- RunPCA(rna_rds, npcs = as.numeric(ini.list$Para$integration_pca_dims), verbose = FALSE)
	rna_rds <- RunUMAP(rna_rds, reduction = "pca",dims = 1:as.numeric(ini.list$Para$RunUMAP_dim))
}else{
	print("RNA数据集中存在umap，不再重新分析")
}

p=DimPlot(rna_rds, group.by = celltype, reduction = "umap", pt.size = 1, label = T, repel = T, label.size = 6)
ggsave(p, file=paste(prefix, 'celltype_umap.pdf', sep='_'), width = 12, height = 8)

### 找差异基因
# 控制三个阈值：logfc.threshold p_val_adj d
Idents(rna_rds)=celltype
maintype_marker=FindAllMarkers(rna_rds, logfc.threshold = ini.list$Para$rna_marker_gene_logfc.threshold, only.pos = T)  #0.5
maintype_marker=maintype_marker%>%filter(p_val_adj < ini.list$Para$rna_marker_gene_padj)  #1e-05
maintype_marker$d=maintype_marker$pct.1 - maintype_marker$pct.2
maintype_marker=maintype_marker%>%filter(d > ini.list$Para$rna_pct.1_pct.2)  #0.2
maintype_marker=maintype_marker%>%arrange(cluster,desc(avg_log2FC))
maintype_marker=as.data.frame(maintype_marker)

### 2.空间转录组流程 ###########################################################
### 整合坐标、region
#（此处把坐标类比单细胞转录组分析中的注释信息）
Idents(space_rds)=seurat_clusters
p1=DimPlot(space_rds, group.by = seurat_clusters, pt.size = 1, label = T, repel = T, label.size = 6)
p2=SpatialDimPlot(space_rds)+theme(legend.position = "right")
p3<-plot_grid(p1, p2)
ggsave(p3, file=paste(prefix, 'space_umap.pdf', sep='_'), width = 20, height = 8)

### 找region特异基因

region_marker=FindAllMarkers(space_rds,logfc.threshold = ini.list$Para$space_marker_gene_logfc.threshold, only.pos = T)  #0
region_marker=region_marker%>%filter(p_val_adj < ini.list$Para$space_marker_gene_padj)  #0.1
region_marker$d=region_marker$pct.1 - region_marker$pct.2
region_marker=region_marker%>%filter(d > ini.list$Para$space_pct.1_pct.2)  #0.05
region_marker=region_marker%>%arrange(cluster,desc(avg_log2FC))
region_marker=as.data.frame(region_marker)

# 说明：
# 1. 上述找两个DEG数据框的方法不唯一，阈值也不唯一
# 2. 第二个DEG数据框也可以是cluster的marker
# 3. MIA分析模式在单细胞和空间转录组场景都可以应用，空转场景是看细胞亚群的富集程度，单细胞场景是做细胞亚群注释


### 画图看看
library(RColorBrewer)
library(scales)
region=unique(space_rds@meta.data[,seurat_clusters])
region_num=length(region)

### 3.MIA分析 ##################################################################
region_specific=region_marker[,c("cluster","gene")]
colnames(region_specific)[1]="region"

celltype_specific=maintype_marker[,c("cluster","gene")]
colnames(celltype_specific)[1]="celltype"

N=length(union(rownames(rna_rds),rownames(space_rds)))

miares=syMIA(region_specific,celltype_specific,N)
write.table(miares, paste(prefix, 'mia_result.xls', sep='_'), sep='\t', quote=F, row.names=F)

tmp = miares[,c(1,2,3)]
tmp = tmp[order(tmp$region, tmp$term),]
data = tmp$pvalue
p_data = matrix(data = data, nrow = length(unique(tmp$term)), ncol = length(unique(tmp$region)), dimnames = list(unique(tmp$term), unique(tmp$region)))
p_data = as.data.frame(p_data)
p=pheatmap::pheatmap(p_data, cluster_row = FALSE,cluster_col=FALSE,color = colorRampPalette(c("white","#87CEFA"))(50))
ggsave(p, file=paste(prefix, 'pheatmap.pdf', sep='_'), width=8, height=5)

library("dplyr")
region_term = tmp %>% group_by(region) %>% slice(which.min(pvalue))
region_term2 = data.frame(region=region_term$region, term=region_term$term)
write.table(region_term2, paste(prefix, 'region_celltype.xls', sep='_'), sep='\t', quote=F, row.names=F)

labers=region_term2[match(as.numeric(as.character(space_rds@active.ident)),region_term2[,1]),2]
space_rds$term = labers
saveRDS(space_rds, file = paste(prefix, 'space_anno.rds', sep='_'))
print('全部分析完成')

#空间辅助注释UMAP
p2 <- DimPlot(space_rds, reduction = "umap", label = F,group.by = "term",label.size = 4,pt.size = 0,raster=FALSE)
Idents(space_rds) = space_rds$term
p3 <- SpatialDimPlot(space_rds)+theme(legend.position = "right")
p4 <- plot_grid(p2, p3)
ggsave(p4, file=paste(prefix, 'space_anno_umap.pdf', sep='_'), width = 20, height = 8)

eachcell <- SpatialDimPlot(space_rds, cells.highlight = CellsByIdentities(object = space_rds, idents = unique(Idents(space_rds))),facet.highlight = TRUE)
ggsave(eachcell, file=paste(prefix, "SpatialDimPlot_Pit1_eachcell.pdf", sep='_'), width=10, height=8)
