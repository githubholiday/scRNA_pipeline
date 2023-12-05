library(dplyr)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(Seurat)
library(configr)
library(RColorBrewer)
library(scales)
library(getopt)

library('getopt')
para<- matrix(c(
	'help',	'h',	0,	"logical",
	'infile',	'i',	1,	"character",
	'outpre',	'o',	1,	"character"
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
	--infile	i	character	the file of MIA enrichment and group column [forced]
	--outpre	o	character	the prefix of outfile[forced]
	\n")
	q(status=1)
}


#绘图的参数
axis.text.y.left.size=14      #主panel纵轴的文本大小
axis.ticks.y.left.length=0.2  #主panel纵轴的刻度线长度
legend.title.size=14          #主panel的legend的title大小
legend.text.size_uppanel=14   #上层panel的legend的文本大小
axis.text.y.left_uppanel=14
plotwidth=18
plotheight=25


miares<- read.csv(opt$infile,sep="\t",header=T)
miares <- data.frame(miares)

miares$region=factor(miares$region,levels = sort(unique(miares$region)))
miares$term=factor(miares$term,levels = sort(unique(miares$term)))
miares$group=factor(miares$group,levels = sort(unique(miares$group)))
miares.Enrichment=miares%>%filter(final_class == "Enrichment")
miares.Depletion=miares%>%filter(final_class == "Depletion")

#miares %>% arrange(group)%>% mutate(order=row_number(),order_factor=factor(order,labels=group)) -> data_sort1 
#miares <- data_sort1
library(ggnewscale)
library(colorspace)
miares.Enrichment$region=factor(miares.Enrichment$region,levels=unique(miares.Enrichment$region))
miares.Depletion$region=factor(miares.Depletion$region,levels=unique(miares.Depletion$region))
pb=ggplot()+
    geom_tile(data = miares.Enrichment,mapping = aes(x=region,y=term,fill=final_value),color="black",size=0.5)+
    scale_fill_gradientn("Enrichment",colours = brewer.pal(9, "Reds"))+
    new_scale_fill() +
    geom_tile(data = miares.Depletion,mapping = aes(x=region,y=term,fill=final_value),color="black",size=0.5)+
    scale_fill_gradientn("Depletion",colours = brewer.pal(9, "Blues"))+
    scale_x_discrete(expand = c(0,0),breaks=miares$region,labels=miares$region)+
    scale_y_discrete(expand = c(0,0),breaks=miares$term,labels=miares$term)+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text.y.left = element_text(size = axis.text.y.left.size,color = "black"),
      axis.ticks.length.y.left = unit(axis.ticks.y.left.length,"cm"),
      axis.text.x.bottom = element_text(size = axis.text.y.left.size,color = "black"),
      axis.ticks.x.bottom = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = legend.title.size,vjust = 1)
    )
    #pdf("/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/yangzhang/test/tmp/bar/LH_plot1.pdf")
    #print(pb)
    #dev.off()
  
  #updf=as.data.frame(levels(miares$region))
  #colnames(updf)="region"
  #updf$region=factor(updf$region,levels = updf$region)
  #print(  updf$region)
  #main_type <-  updf$region
  updf=as.data.frame(miares[,c("group","region")])
  updf=updf[!duplicated(updf$region),]
#  colnames(updf)="group"
  updf$group=factor(updf$group)
  updf$region=factor(updf$region,levels=updf$region)
  print(updf$region)
  main_type <-  updf$group
  print(main_type)
  #q4<- hcl_palettes("sequential (single-hue)", n = length(updf$region), plot = TRUE)
  q4 <- sequential_hcl(length(main_type), palette = "Red-Yellow")
  print(q4)
  pa=ggplot(data = updf,aes(x=region,y=0,fill=group))+
    geom_tile()+
    scale_x_discrete(expand = c(0,0))+
    #scale_fill_manual(values = q4,breaks = labeldf1$region,labels=labeldf1$region_gene_num)+
    scale_fill_manual(values = q4,breaks =main_type,labels=main_type)+
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
    #pdf("/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/yangzhang/test/tmp/bar/LH_plot2.pdf")
    #print(pa)
    #dev.off()
  library(patchwork)
  celltype_num=length(unique(miares$term))
  pa / pb + plot_layout(heights = c(1,celltype_num))
  ggsave(opt$outpre,width = plotwidth,height = plotheight,units = "cm")
  print('MIA分析运行完成')