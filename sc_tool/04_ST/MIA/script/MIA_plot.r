library(dplyr)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(Seurat)
library(configr)
library(RColorBrewer)
library(scales)
library(ggnewscale)
library(colorspace)

axis.text.y.left.size=14      #主panel纵轴的文本大小
axis.ticks.y.left.length=0.2  #主panel纵轴的刻度线长度
legend.title.size=14          #主panel的legend的title大小
legend.text.size_uppanel=14   #上层panel的legend的文本大小
axis.text.y.left_uppanel=14
plotwidth=18
plotheight=25


miares<- read.csv("/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/03_Spatial_tr/02_Analysis_xifen/result/Annotation/MIA/test/MIA.xls",header=TRUE,sep="\t")
miares <- data.frame(miares)
miares$region=factor(miares$region,levels = sort(unique(miares$region)))
miares$term=factor(miares$term,levels = sort(unique(miares$term)))
miares$group=factor(miares$group,levels = sort(unique(miares$group)))
miares.Enrichment=miares%>%filter(final_class == "Enrichment")
print(miares.Enrichment)
miares.Depletion=miares%>%filter(final_class == "Depletion")
miares %>% arrange(group)%>% mutate(order=row_number(),order_factor=factor(order,labels=region)) -> data_sort1
miares <- data_sort1  

miares.Enrichment%>% arrange(group)%>% mutate(order=row_number(),order_factor=factor(order,labels=region)) -> data_sort2
miares.Enrichment <- data_sort2

miares.Depletion%>% arrange(group)%>% mutate(order=row_number(),order_factor=factor(order,labels=region)) -> data_sort3
miares.Depletion <- data_sort3

updf=as.data.frame(miares[,c("group","region")])
updf=updf[!duplicated(updf$region),]
updf %>% arrange(group)%>% mutate(order=row_number(),order_factor=factor(order,labels=region)) -> data_sort1
updf <- data_sort1
updf=updf[!duplicated(updf$region),]



pb=ggplot()+
    geom_tile(data = miares.Enrichment,mapping = aes(x=factor(miares.Enrichment$region,levels=unique(miares.Enrichment$region)),y=term,fill=final_value),color="black",size=0.5)+
    scale_fill_gradientn("Enrichment",colours = brewer.pal(9, "Reds"))+
    new_scale_fill() +
    geom_tile(data = miares.Depletion,mapping = aes(x=region,y=term,fill=final_value),color="black",size=0.5)+
    scale_fill_gradientn("Depletion",colours = brewer.pal(9, "Blues"))+
    scale_x_discrete(expand = c(0,0),breaks=updf$region,labels=updf$region)+
    #scale_y_discrete(expand = c(0,0),breaks=labeldf2$celltype,labels=labeldf2$celltype_gene_num)+
    scale_y_discrete(expand = c(0,0),breaks=miares$term,labels=miares$term)+
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
    pdf("/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/03_Spatial_tr/02_Analysis_xifen/result/Annotation/MIA/test/LH_plot1.pdf")
    print(pb)
    dev.off()
  
  #updf=as.data.frame(levels(miares$region))
  #colnames(updf)="region"
  #updf$region=factor(updf$region,levels = updf$region)
  #print(  updf$region)
  #main_type <-  updf$region

  #updf=as.data.frame(levels(miares$group))
  
  #colnames(updf)="group"
  updf$group=factor(updf$group)
  print(  updf$group)
  main_type <-  updf$group
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
  pdf("/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/03_Spatial_tr/02_Analysis_xifen/result/Annotation/MIA/test/LH_plot2.pdf")
  print(pa)
    dev.off()
  library(patchwork)
  celltype_num=length(unique(miares$term))
  pa / pb + plot_layout(heights = c(1,celltype_num))
  prefix = "/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/03_Spatial_tr/02_Analysis_xifen/result/Annotation/MIA/test/LH"
  ggsave(paste(prefix, 'MIA.pdf', sep='_'),width = plotwidth,height = plotheight,units = "cm")
  print('MIA分析运行完成')