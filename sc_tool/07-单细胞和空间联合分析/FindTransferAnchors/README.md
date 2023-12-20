### 模块： MK_FindTransferAnchors

*模块功能：
*模块版本：v1.0.0
*邮箱： mengli@genome.cn

### 使用示例及参数说明：
Usage:
make -f MK_FindTransferAnchors config= outdir= rnards= spacerds= prefix= info= spot_cell
参数说明：
	config: [文件|可选]  模块配置文件，和软件相关参数，默认为$(makefile_dir)/config/config.txt
	outdir: [路径|必需]  输出路径
	rnards: [文件|必需]  带有细胞注释结果的scRNA-seq rds 文件，需要在assay下有RNA的数据槽
	spacerds: [文件|必需]  空间转录组的 rds 文件
	prefix: [字符|必需]  文件前缀名称
	info: [文件|必需]  分析的参数配置文件

target说明：
	spot_cell：通过MK_FindTransferAnchors函数注释空间转录组spot
	report：整理目录并生成report

### 运行环境、软件及数据库：
	北京 sge 集群
	软件：R /annoroad/data1/bioinfo/PMO/yaomengcheng/bk_Anaconda3/envs/monocle3/bin/Rscript
	输入文件大小：200M
	申请CPU：6
	申请内存：80G
	实际CPU：4
	实际内存：40G
	运行时长：5min
### 输入文件示例:
见test/input/
./
|-- config.ini #分析参数的配置文件
|-- normal_celltype_SCT.rds #RNA数据集的rds文件
`-- NP183765_0.8.rds #空间数据集的rds文件

### 输出文件示例
见test/output/
|-- readme.doc
`-- sample
    |-- sample_celltype_count_space_cluster.xls
    |-- sample_heatmap_cluster.pdf
    |-- sample_heatmap_cluster.png
    |-- sample_heatmap_spot_metadata.xls
    |-- sample_predicted.id_umap1.pdf
    |-- sample_predicted.id_umap1.png
    |-- sample_predictions.xls
    |-- sample_SpatialDimPlot.pdf
    |-- sample_SpatialDimPlot_Pit1_eachcell.pdf
    |-- sample_SpatialDimPlot_Pit1_eachcell.png
    |-- sample_SpatialDimPlot.png
    `-- sample_Spot_cell_count.xls

文件解释说明：
1：sample_celltype_count_space_cluster.xls：细胞类型在不同cluster中的数量统计
（1）celltype：细胞类型；
（2）cluster：cluster编号；
（3）count：细胞数量。
2：sample_heatmap_cluster.p*：空间数据集中的每个cluster中的不同细胞类型的占比统计图
横轴为cluster，纵轴为细胞类型，数值为在cluster中的细胞类型占比。数值越大，颜色越接近黄色，表示该cluster中该细胞类型较多。

3：sample_heatmap_spot_metadata.xls：每个spot的metadata数据信息
（1）spot：spot的编号；
（2）seurat_clusters：spot的cluster编号；
（3）predicted.id：该spot预测的细胞类型；
（4）prediction.score.*：该spot在不同细胞类型中的预测分数；
（5）prediction.score.max：最大的预测分数。

4：sample_predicted.id_umap1.p*：空间数据集整合RNA数据后的细胞类型的UMAP图
不同的颜色代表不同的细胞类型。

5：sample_predictions.xls：每个spot的各个细胞类型score值
（1）spot：spot的编号；
（2）seurat_clusters：spot的cluster编号；
（3）predicted.id：该spot预测的细胞类型；
（4）prediction.score.*：该spot在不同细胞类型中的预测分数；
（5）prediction.score.max：最大的预测分数。

6：sample_SpatialDimPlot.p*：所有的细胞类型在空间转录组中的切片图展示。
每个点代表一个spot。不同的颜色代表不同的细胞类型。

7：sample_SpatialDimPlot_Pit1_eachcell.p*：每个细胞类型分别在空间转录组中的切片图展示。

8：sample_Spot_cell_count.xls：空间数据中细胞类型预测的数量
（1）cell：细胞类型；
（2）spot count：spot数量。

### 注意事项
RNAcelltype = celltype #RNA数据集的细胞类型的标签
spacecluster = seurat_clusters #空间数据集的要分析的分类标签
marker_gene_min.pct = 0.1 #FindAllMarkers函数中的min.pct数值
marker_gene_logfc.threshold = 0.25 #FindAllMarkers函数中的logfc.threshold数值
marker_gene_padj = 0.05 #marker基因的筛选p_val_adj小于该数值的基因
marker_gene_num = 10 #绘图展示前10的marker 基因
anchors_ims = 20 #FindTransferAnchors函数中的dims变量最大值
normalization.method = SCT #FindTransferAnchors函数中标准化的方法
common_gene = 1000 #RNA数据集和空间数据集中的共有的最少细胞数

