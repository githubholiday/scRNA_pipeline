@@@@project_info
MainMenu:引言
SubMenu:方法介绍
P:#,;#基于 MIA（Multimodal intersection analysis）的算法可以将RNA数据中的细胞类型与空间转录组中的位置对应起来。该算法首先利用FindAllMarkers函数查找两个数据集中的高可变marker基因。利用超几何分布检验来评估两者数据之间的交集情况是否存在显著富集enrichment或缺失deleption（P<0.05）从而得到每个cluster可能的细胞类型。

MainMenu:项目信息
@@@@spaceinfo
SubMenu:空间数据信息
P:#,;#下图为空间数据集中的 cluster 的UMAP图和在空间切片上的展示图：
Image:upload/*Space/*/*space_umap.png,400,4,空间数据展示图
P:#,;#通过此图可以明确空间数据集中的 spot 分群信息，以及 cluster 在切片上的空间位置分布。
P:#,;#空间数据展示图下载链接：
Excel:upload/*Space/*/*space_umap.p*,,,空间数据展示图下载链接：
@@@@RNAinfo
SubMenu:RNA数据信息
P:#,;#下图为RNA数据集细胞类型可视化展示UMAP图：
Image:upload/*RNA/*/*celltype_umap.png,400,4,细胞类型UMAP展示图
P:#,;#通过此图，可以了解RNA数据集中所包含的所有细胞类型，不同的细胞类型在UMAP图上没有重叠，说明细胞类型的特异性较好。
P:#,;#细胞类型UMAP展示图下载链接：
Excel:upload/*RNA/*/*celltype_umap.p*,,,细胞类型UMAP展示图下载链接：
P:#,;#除此之外，统计在RNA数据集中不同细胞类型的数量，可以了解不同细胞类型在数据集中的分布情况：
Table:upload/*RNA/*/*celltype_count.xls,,8,850,,0,细胞类型数量统计表
PRE:
（1）celltype：细胞类型；
（2）count：对应细胞类型的细胞数目。
PRE
P:#,;#细胞类型数量统计表下载链接：
Excel:upload/*RNA/*/*celltype_count.xls,,,细胞类型数量统计表下载链接：
@@@@MIA
MainMenu:联合分析
SubMenu:marker基因交集
P:#,;#统计RNA和空间数据集中的marker基因共有的情况，如下表所示：
Table:upload/*_MIA/*/*_region_celltype_gene.xls,,8,850,,0,共有marker基因数量统计表
PRE:
（1）celltype：细胞类型；
（2）第二列-：空间数据集cluster和对应细胞类型共有marker基因数量。
PRE
P:#,;#共有marker基因数量统计表下载链接：
Excel:upload/*_MIA/*/*_region_celltype_gene.xls,,,共有marker基因数量统计表下载链接：
P:#,;#各细胞类型同空间数据集共有的 top3 marker基因在空间转录组中的表达情况如下：
Image:upload/*_MIA/*/*top3_gene_FeaturePlot.png,400,4,细胞类型top3 marker基因展示图
P:#,;#细胞类型top3 marker基因展示图下载链接：
Excel:upload/*_MIA/*/*top3_gene_FeaturePlot.p*,,,细胞类型top3 marker基因展示图下载链接：
SubMenu: cluster预测结果
P:#,;#使用MIA方法推断特定细胞类型在一个给定空间数据集cluster的富集情况。得到每个cluster预测出的细胞类型，如下方示例表格展示同一个空间上的cluster与各种细胞类型之间的富集程度或者是缺失程度。
Table:upload/*_MIA/*/example_mia_result.xls,,8,850,,0,MIA结果示例表
PRE:
（1）region：空间转录组中的cluster编号；
（2）term：RNA转录组中的细胞类型；
（3）pvalue：P值；
（4）Enrichment：富集值；
（5）Depletion：缺失值；
（6）final_value：富集值和缺失值中的最大值；
（7）final_class：该细胞类型在该cluster中的状态，是富集还是缺失。
PRE
P:#,;#注：Depletion=-log10(1-P)  ，Enrichment= -log10(P)。P值越小，Depletion越小, Enrichment值越大。
P:#,;#MIA结果表下载链接：
Excel:upload/*_MIA/*/*_mia_result.xls,,,MIA结果表下载链接：

P:#,;#利用heatmap图展示不同cluster在各个细胞类型中的P值情况：
Image:upload/*_MIA/*/*pheatmap.png,400,4,cluster同细胞类型显著性热图
P:#,;#横轴表示不同cluster，纵轴表示不同细胞类型的P值。
P:#,;#cluster同细胞类型显著性热图下载链接：
Excel:upload/*_MIA/*/*pheatmap.p*,,,cluster同细胞类型显著性热图下载链接：

P:#,;#将enrichment 和Depletion两个值在热图中展示，可以看出每个cluster最可能的细胞类型:
Image:upload/*_MIA/*/*_MIA.png,400,4,MIA结果展示图
P:#,;#注：括号里面的genes数值为高可变基因数量。
P:#,;#MIA结果展示图下载链接：
Excel:upload/*_MIA/*/*_MIA.p*,,,MIA结果展示图下载链接：
P:#,;#统计出空间数据中每个cluster预测出的可能细胞类型信息，如下表所示：
Table:upload/*_MIA/*/*_region_celltype.xls,,8,850,,0,cluster预测细胞类型结果表
PRE:
（1）region：空间转录组中的cluster编号；
（2）term：预测的细胞类型。
PRE
P:#,;#cluster预测细胞类型结果表下载链接：
Excel:upload/*_MIA/*/*_region_celltype.xls,,,cluster预测细胞类型结果表下载链接：
P:#,;#空间数据集映射RNA数据之后的细胞类型的UMAP图：
Image:upload/*_MIA/*/*_space_anno_umap.png,400,4,空间注释UMAP结果展示图
P:#,;#不同的颜色代表不同的细胞类型，通过此图可以看出空间数据集中存在的细胞类型种类。
P:#,;#空间注释UMAP图下载链接：
Excel:upload/*_MIA/*/*_space_anno_umap.p*,,,空间注释UMAP图下载链接：
P:#,;#分别展示各个细胞类型在空间组织区域的分布情况。如下图所示：
Image:upload/*_MIA/*/*_SpatialDimPlot_Pit1_eachcell.png,400,4,不同细胞类型分布展示图
P:#,;#展示细胞类型在空间数据集的位置分布情况，其中红色高亮为对应细胞类型区域，灰色为其他区域。
P:#,;#不同细胞类型分布展示图下载链接：
Excel:upload/*_MIA/*/*_SpatialDimPlot_Pit1_eachcell.p*,,,不同细胞类型分布展示图下载链接：
@@@@appendix
MainMenu:附录
SubMenu:参考文件
P:#,;#Moncada R , Barkley D , Wagner F , et al. Integrating microarray-based spatial transcriptomics and single-cell RNA-seq reveals tissue architecture in pancreatic ductal adenocarcinomas[J]. Nature Biotechnology.
SubMenu: 软件参数
P:#,;#本次分析所使用的软件见下表：
Table:upload/*public-picture/software.xls,,5,950,,0,软件列表
PRE:
（1）Number：序号；
（2）software：软件名称或包名称；
（3）version：软件版本或包版本；
（4）function：软件作用；
（5）parameters：参数说明。
PRE
P:#,;#分析软件列表下载链接：
Excel:upload/*public-picture/software.xls,,,结果下载链接：
SubMenu: 结果目录
ShowDir:$(REPORT_DIR)

