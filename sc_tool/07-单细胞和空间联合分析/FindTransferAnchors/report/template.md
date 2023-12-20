@@@@project_info
MainMenu: 引言
SubMenu: 软件原理
P:#,;#采用seurat工具将RNA转录组的细胞类型信息映射到空间数据集上。一个关键步骤是在两个数据集之间识别锚点anchors。 这些锚点是处于同一生物学状态的跨数据集细胞对。 分别使用 FindTransferAnchors 函数计算用于锚点的识别，并得到的锚点信息传递给TransferData函数，得到每个spot的最可能的一种细胞类型预测结果。
Image:upload/*public-picture/image.png,400,1,原理图

@@@@spaceinfo
MainMenu: 空间数据信息
SubMenu: 样本聚类信息展示
P:#,;#下图为空间数据集中的cluster的UMAP图和在空间切片上的展示图：
Image:upload/*space/*cluster_umap.png,400,1,空间数据集的UMAP图和空间切片图。
P:#,;#空间数据集的UMAP图下载链接：
Excel:upload/*space/*cluster_umap.p*,,,空间数据库的UMAP图下载链接：
P:#,;#通过此图可以明确空间数据集中的spot的分群信息，以及cluster在切片上的空间位置关系。
@@@@RNAinfo
MainMenu: RNA数据信息
SubMenu: 细胞类型信息展示
P:#,;#首先利用umap图展示RNA数据集中包含的细胞类型：
Image:upload/*RNA/*celltype_umap.png,400,1,RNA数据集中的细胞类型UMAP图。
P:#,;#RNA数据集中的细胞类型的UMAP图下载链接：
Excel:upload/*RNA/*celltype_umap.p*,,,空间数据库的UMAP图下载链接：
P:#,;#通过此图，可以了解RNA数据集中所包含的所有的细胞类型、不同的细胞类型在UMAP图上没有重叠，说明细胞类型的特异性较好。
P:#,;#下方表格为统计在RNA数据集中不同细胞类型的数量，可以了解不同细胞类型在数据集中的分布情况：
Table:upload/*RNA/*celltype_count.xls,,,1000,,0,细胞类型统计表
PRE:
（1）celltype：细胞类型；
（2）count：细胞数量。
PRE

P:#,;#细胞类型统计表下载链接：
Excel:upload/*RNA/*celltype_count.xls,400,,细胞类型统计表下载链接：

@@@@Integrate
MainMenu: 联合分析结果
SubMenu: spot预测结果
P:#,;#使用 FindTransferAnchors 和TransferData函数，来整合RNA数据集和单细胞转录组和空间转录组信息，得到每个spot预测的细胞类型，如下方示例表格展示：
Table:upload/*Integrate/*predictions_example.xls,,,1000,,0,细胞类型统计表
PRE:
（1）spot：spot的编号；
（2）predicted.id：该spot预测的细胞类型；
（3）prediction.score.*：该spot在不同细胞类型中的预测分数；
（4）prediction.score.max：最大的预测分数。
PRE
P:#,;#spot细胞类型预测结果表格下载链接：
Excel:upload/*Integrate/*predictions.xls,400,,spot细胞类型预测结果表格下载链接：
P:#,;#通过此表格，可以看出每个spot预测到的细胞类型，预测分数越高，说明属于某种细胞类型的可能性越高。

P:#,;#另外，统计出了空间数据中细胞类型预测的数量，可以看出空间数据集中的细胞类型的组成：
Table:upload/*Integrate/*Spot_cell_count.xls,,,1000,,0,spot的细胞类型统计表
PRE:
（1）celltype：细胞类型；
（2）count：spot数量。
PRE
P:#,;#spot的细胞类型统计表下载链接：
Excel:upload/*Integrate/*Spot_cell_count.xls,400,,spot的细胞类型统计表下载链接：
P:#,;#空间数据集注释后的细胞类型的UMAP图：
Image:upload/*Integrate/*predicted.id_umap1.png,400,1,空间数据集中的细胞类型UMAP图。
P:#,;#不同的颜色代表不同的细胞类型，通过此图可以看出空间数据集中存在的细胞类型种类。
P:#,;#空间数据集中的细胞类型UMAP图下载链接：
Excel:upload/*Integrate/*predicted.id_umap1.p*,,,空间数据集中的细胞类型UMAP图下载链接：

P:#,;#所有的细胞类型在空间转录组切片中的展示如下：
Image:upload/*Integrate/*SpatialDimPlot.png,400,1,空间转录组切片的细胞类型展示图。
P:#,;#不同的颜色代表不同的细胞类型。
P:#,;#空间转录组切片的细胞类型展示图下载链接：
Excel:upload/*Integrate/*SpatialDimPlot.p*,,,空间转录组切片的细胞类型展示图下载链接：
P:#,;#同时，还可以分别展示各个细胞类型在空间组织区域的分布情况。如下图所示：
Image:upload/*Integrate/*SpatialDimPlot_Pit1_eachcell.png,400,1,各个细胞类型在空间组织区域的分布图。
P:#,;#每一张子图代表一种细胞类型。红色代表此spot是该细胞类型。
P:#,;#各个细胞类型在空间组织区域的分布图下载链接：
Excel:upload/*Integrate/*SpatialDimPlot_Pit1_eachcell.p*,,,各个细胞类型在空间组织区域的分布图下载链接:
SubMenu: 不同cluster中细胞类型统计
P:#,;#下图为空间数据集的每个cluster中不同细胞类型的占比统计图：
Image:upload/*Integrate/*heatmap_cluster.png,400,1,每个cluster中不同细胞类型的占比统计图。
P:#,;#横轴为cluster。纵轴为细胞类型，数值为在cluster中的细胞类型占比。数值越大，颜色越接近黄色，表示该cluster中该细胞类型较多。
P:#,;#每个cluster中不同细胞类型的占比统计图下载链接：
Excel:upload/*Integrate/*heatmap_cluster.p*,,,每个cluster中不同细胞类型的占比统计图下载链接：
P:#,;#细胞类型在不同cluster中的数量统计：
Table:upload/*Integrate/*celltype_count_space_cluster.xls,,,1000,,0,细胞类型在不同cluster中的数量统计表
PRE:
（1）celltype：细胞类型；
（2）cluster：空间数据集中的cluster编号；
（3）count：空间数据集该cluster中确认为该细胞类型的spot数量。
PRE
P:#,;#细胞类型在不同cluster中的数量统计表下载链接：
Excel:upload/*Integrate/*celltype_count_space_cluster.xls,,细胞类型在不同cluster中的数量统计表下载链接：
@@@@appendix
MainMenu: 附录
SubMenu: 参考文献
P:#,;#Stuart T , Butler A , Hoffman P , et al. Comprehensive Integration of Single-Cell Data[J]. Cell, 2019, 177(7):1888-1902.e21.
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

