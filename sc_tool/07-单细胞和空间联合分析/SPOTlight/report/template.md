@@@@project_info
MainMenu: 引言
SubMenu: 软件原理
P:#,;#SPOTlight将空间转录组与注释后scRNA-seq数据进行整合，是用于估计每个spot中每种细胞类型的比例的算法。该算法以种子型非负矩阵分解（NMF）回归为中心，使用细胞类型标记基因和非负最小二乘法（NNLS）进行初始化，识别出每个spot中的细胞类型和比例。随后对空间转录组中的每个spot进行细胞类型的解析。
Image:upload/*public-picture/image.png,400,1,原理图

@@@@spaceinfo
MainMenu: 空间数据信息
SubMenu: 样本聚类信息展示
P:#,;#下图为空间数据集中的cluster的UMAP图和在空间切片上的展示图：
Image:upload/*space/*spatial_cluster_umap.png,400,4,空间数据集的UMAP图和空间切片图
P:#,;#通过此图可以明确空间数据集中的spot的分群信息，以及cluster在切片上的空间位置关系。
P:#,;#空间数据集的UMAP图下载链接：
Excel:upload/*space/*spatial_cluster_umap.p*,,,空间数据集的UMAP图下载链接：

@@@@RNAinfo
MainMenu: RNA数据信息
SubMenu: 细胞类型信息展示
P:#,;#首先利用umap图展示RNA数据集中包含的细胞类型：
Image:upload/*RNA/*singlecell_cluster_umap.png,400,1,RNA数据集中的细胞类型UMAP图
P:#,;#通过此图，可以了解RNA数据集中所包含的所有的细胞类型。不同的细胞类型在UMAP图上没有重叠，说明细胞类型的特异性较好。
P:#,;#RNA数据集中细胞类型的UMAP图下载链接：
Excel:upload/*RNA/*singlecell_cluster_umap.p*,,,RNA数据集中细胞类型的UMAP图下载链接：

P:#,;#下方表格为统计在RNA数据集中不同细胞类型的数量，可以了解不同细胞类型在数据集中的分布情况：
Table:upload/*RNA/*celltype_count.xls,,,1000,,0,单细胞的细胞类型统计表
PRE:
（1）celltype：细胞类型；
（2）count：细胞数量。
PRE
P:#,;#细胞类型统计表下载链接：
Excel:upload/*RNA/*celltype_count.xls,400,1,单细胞的细胞类型统计表下载链接：

@@@@Integrate
MainMenu: 联合分析结果
SubMenu: spot中的细胞类型占比统计
P:#,;#得到每个spot的细胞类型占比统计，如下表：
Table:upload/*Integrate/*all_spot_celltype_example.xls,,,1000,,0,每个spot的细胞类型占比统计表
PRE:
（1）spot：barcode ID；
（2）第二行-：不同的细胞类型在该spot中的占比。
PRE
P:#,;#通过此表格，可以看出每个spot预测到的细胞类型及其各自的比例。
P:#,;#spot的细胞类型占比统计结果表格下载链接：
Excel:upload/*Integrate/*all_spot_celltype.xls,400,,每个spot的细胞类型占比统计表下载链接：


P:#,;#首先将每个spot的细胞类型及比例在空间切片图上进行展示:
Image:upload/*Integrate/*SpatialScatterpie.png,400,1,每个spot的细胞类型及比例在空间切片图上的展示图
P:#,;#每一个spot中用饼图展示各个细胞类型的占比情况。不同颜色代表不同的细胞类型，其中占比低于0.1的细胞类型将不会展示在饼图中。
P:#,;#每个spot的细胞类型及比例在空间切片图上的展示图下载链接：
Excel:upload/*Integrate/*SpatialScatterpie.p*,400,1,每个spot的细胞类型及比例在空间切片图上的展示图下载链接：


P:#,;#然后在空间切片图上分别展示各个细胞类型占比:
Image:upload/*Integrate/*single_celltype*.png,400,1,空间切片图上各个细胞类型占比
P:#,;#左图为空间转录组中的聚类结果；右图为SPOTlight分析的结果，颜色代表该细胞类型在该spot中的占比，颜色越红，代表占比越大。
P:#,;#空间切片图上各个细胞类型占比展示图下载链接：
Excel:upload/*Integrate/*single_celltype*.p*,400,1,空间切片图上各个细胞类型占比展示图下载链接：

P:#,;#同时对空间转录组中不同的cluster分别展示spot中细胞类型占比统计结果：
Image:upload/*Integrate/*SpatialScatterpie_*.png,400,1,各个cluster中细胞类型占比饼图
P:#,;#各个cluster中细胞类型占比饼图下载链接：
Excel:upload/*Integrate/*SpatialScatterpie_*.p*,400,1,各个cluster中细胞类型占比饼图下载链接：

P:#,;#最后将每种细胞类型在不同cluster中占比情况以箱式图的形式进行展示，如下：
Image:upload/*Integrate/*boxplot*.png,400,1,不同cluster中细胞类型占比箱式图
P:#,;#一个细胞绘制一个图；横坐标是空间的spot cluster分类。纵坐标是spot中的该细胞类型所占的百分比。
P:#,;#不同cluster中细胞类型占比箱式图下载链接：
Excel:upload/*Integrate/*boxplot*.p*,400,1,不同cluster中细胞类型占比箱式图下载链接：

SubMenu: 不同细胞类型在空间分布上的相似性
P:#,;#根据每个spot中每种细胞类型所占的比例来计算各细胞类型的相似性（pearson），值为正，说明两种细胞类型在空间上的分布趋势是正相关的，越接近1，分布趋势越相似；值为负，代表负相关，两种细胞类型在空间上的分布趋势是相反的，越接近-1，分布趋势相反。
Image:upload/*Integrate/*CorrelationMatrix.png,400,1,不同细胞类型在空间分布上的相似性
P:#,;#不同细胞类型在空间分布上的相似性图下载链接：
Excel:upload/*Integrate/*CorrelationMatrix.p*,400,1,不同细胞类型在空间分布上的相似性图下载链接：

SubMenu: 细胞类型间的空间共定位
P:#,;#细胞类型之间的关联性越强，我们在同一个点中发现它们的频率就越高。因此，基于每个spot中的不同细胞类型的有无，利用get_spatial_interaction_graph函数统计所有spot中两两细胞类型共同出现的次数。以下网络图展示了不同细胞类型在空间上的共定位情况。
Image:upload/*Integrate/*Interactions.png,400,1,细胞类型间的空间共定位
P:#,;#线的宽度表示的是两种细胞在空间上的相近性，线越粗，表示两个细胞在空间上越相近。
P:#,;#细胞类型间的空间共定位图下载链接：
Excel:upload/*Integrate/*Interactions.p*,400,1,细胞类型间的空间共定位图下载链接：

@@@@appendix
MainMenu: 附录

SubMenu: 参考文献
P:#,;# Elosua M , Nieto P , Mereu E , et al. SPOTlight:Seeded NMF regression to Deconvolute Spatial Transcriptomics Spots with Single-Cell Transcriptomes[J]. Nucleic Acids Research, 2021.

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
