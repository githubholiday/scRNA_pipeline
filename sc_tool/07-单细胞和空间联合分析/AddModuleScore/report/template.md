@@@@project_info
MainMenu: 引言
SubMenu: 软件原理
P:#,;#基于 AddModuleScore函数可以将RNA数据中的细胞类型与空间转录组中的位置对应起来。该函数的主要功能就是为基因集来进行打分。AddModuleScore函数的原理是从RNA数据集中的每个细胞类型获得marker基因集，计算此基因集在不同的空间的组织区域进行富集评分。通过富集评分把单细胞数据映射到空间转录组上，富集评分就可以看出这个细胞类型在每一个spot的富集情况。

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
SubMenu: 富集评分
P:#,;#将每一种细胞类型的top10marker基因作为一组基因集, 利用AddModuleScore()函数分析该基因集在每个空间数据集中cluster的score值。score值表示的是背景基因的平均值在于找每个基因的所在的bin，在该bin内随机抽取其他基因作为背景，最后所有的目标基因算一个平均值，所有的背景基因算一个平均值，两者相减就是该基因集的score值。score值越高，表明代表这个细胞类型在每一个spot中的比例越大。
P:#,;#将score值用小提琴图展示：
Image:upload/*Integrate/*celltype_Vlnplot.png,400,,小提琴图
P:#,;#横轴为空间cluster编号，纵轴为计算得到的score值。每个点代表一个spot。
P:#,;#小提琴图下载链接：
Excel:upload/*Integrate/*celltype_Vlnplot.p*,,,小提琴图下载链接：
P:#,;#每种细胞类型在spot中的score值用Featureplot图展示。如下:
Image:upload/*Integrate/*celltype_FeaturePlot.png,400,,Featureplot图
P:#,;#每个点表示一个spot，颜色越紫，表示该细胞类型在spot中的富集程度越高。
P:#,;#Featureplot图下载链接：
Excel:upload/*Integrate/*celltype_FeaturePlot.p*,,,Featureplot图下载链接：
P:#,;#每个spot中不同细胞类型的score值可通过查询以下的表格中的数值：
Table:upload/*Integrate/*celltype_AddModuleScore_example.xls,,,1000,,0,细胞类型score值统计表
PRE:
（1）spot: spot编号；
（2）orig.ident：该spot所属的样本类型；
（3）seurat_clusters：该spot所在的空间数据集中的cluster编号；
（4）第4行-：不同细胞类型的score值。
PRE
P:#,;#细胞类型score值统计表完整版的下载链接：
Excel:upload/*Integrate/*celltype_AddModuleScore.xls,400,,细胞类型score统计表完整版的下载链接：
SubMenu: 细胞类型鉴定结果
P:#,;#基于每个spot在不同细胞类型中的score值，我们选取score值最高的细胞作为spot的细胞类型。
P:#,;#统计出空间数据中细胞类型预测的数量，可以看出空间数据集中的细胞类型的组成：
Table:upload/*Integrate/*celltype_count_space.xls,,,1000,,0,空间数据中细胞类型预测统计表
PRE:
（1）celltype：细胞类型；
（2）count：spot的数量。
PRE
P:#,;#空间数据中细胞类型预测统计表下载链接：
Excel:upload/*Integrate/*celltype_count_space.xls,,,空间数据中细胞类型预测统计表下载链接：
P:#,;#在空间数据集的不同cluster中，分别统计细胞类型的组成：
Table:upload/*Integrate/*celltype_count_space_cluster.xls,,,1000,,0,空间数据集的不同cluster中的细胞组成表
PRE:
（1）celltype：细胞类型；
（2）cluster：空间数据集中的cluster编号；
（3）count：空间数据集该cluster中确认为该细胞类型的spot的数量。
PRE
P:#,;#空间数据集的不同cluster中的细胞组成表下载链接：
Excel:upload/*Integrate/*celltype_count_space_cluster.xls,,,空间数据集的不同cluster中的细胞组成表下载链接：
P:#,;#利用百分比堆积图展示不同空间cluster中的各个细胞的占比情况：
Image:upload/*Integrate/*celltype_count_space_cluster.png,400,1,细胞占比图
P:#,;#不同颜色代表不同的细胞类型；横轴为空间数据集中的cluster编号。
P:#,;#细胞占比图下载链接：
Excel:upload/*Integrate/*celltype_count_space_cluster.p*,,,细胞占比图下载链接：
SubMenu: 细胞类型空间可视化
P:#,;#空间数据集用RNA数据注释之后的细胞类型的UMAP图：
Image:upload/*Integrate/*celltype_space_umap.png,400,1,空间数据集的细胞类型umap图
P:#,;#不同颜色代表不同的细胞类型；每个点代表一个spot编号。通过此图可以看出空间数据集中存在的细胞类型种类。
P:#,;#空间数据集的细胞类型umap图下载链接：
Excel:upload/*Integrate/*celltype_space_umap.p*,,,空间数据集的细胞类型umap图下载链接：
P:#,;#所有的细胞类型在空间转录组切片中的展示如下：
Image:upload/*Integrate/*SpatialDimPlot.png,400,1,细胞类型在空间转录组切片展示图
P:#,;#不同颜色代表不同的细胞类型；通过此图可以看出不同区域的细胞存在情况。
P:#,;#细胞类型在空间转录组切片展示图下载链接：
Excel:upload/*Integrate/*SpatialDimPlot.p*,,,细胞类型在空间转录组切片展示图下载链接：
P:#,;#分别展示各个细胞类型在空间组织区域的分布情况。如下图所示：
Image:upload/*Integrate/*SpatialDimPlot_Pit1_eachcell.png,400,1,各个细胞类型在空间切片展示图
P:#,;#每一张图代表一种细胞类型在空间切片的存在情况。红色代表此spot中存在该细胞，灰色代表不存在该细胞。
P:#,;#各个细胞类型在空间切片展示图下载链接：
Excel:upload/*Integrate/*SpatialDimPlot_Pit1_eachcell.p*,,,各个细胞类型在空间切片展示图下载链接：

@@@@appendix
MainMenu: 附录
SubMenu: 参考文献
P:#,;#Ji A L , Rubin A J , Thrane K , et al. Multimodal Analysis of Composition and Spatial Architecture in Human Squamous Cell Carcinoma (vol 182, pg 497, 2020)[J]. Cell, 2020(6):182.
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













