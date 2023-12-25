### 模块： MK_MIA

*模块功能：单细胞转录组和空间转录组联合分析-MIA分析
*模块版本：v1.0.0
*邮箱： leiguo@genome.cn

### 使用示例及参数说明：
Usage:
	make -f MK_MIA config= outdir= rnards= spacerds= prefix= info= MIA
参数说明：
	config: [文件|可选]  模块配置文件，和软件相关参数，默认为$(makefile_dir)/config/config.txt
	outdir: [路径|必需]  输出路径
	rnards: [文件|必需]  带有细胞注释结果的scRNA-seq rds 文件
	spacerds: [文件|必需]  空间转录组的 rds 文件
	prefix: [字符|必需]  文件前缀名称
	info: [文件|必需]  分析的参数配置文件

target说明：
	"MIA：MIA富集分析"
	"report：整理目录并生成report"


### 运行环境、软件及数据库：
	北京 SEG 集群
	软件：R /annoroad/data1/bioinfo/PMO/yaomengcheng/bk_Anaconda3/envs/monocle3/bin/Rscript
	版本：4.0.3

	申请CPU：8
	申请内存：20G
	运行时长：50m

### 输入文件示例
见test/input/
./
|-- config.ini #分析参数的配置文件
|-- normal_celltype_SCT.rds #RNA数据集的rds文件
`-- NP183765_0.8.rds #空间数据集的rds文件


### 输出文件示例
见test/output/
sample/
|-- sample_MIA.pdf
|-- sample_mia_result.xls
|-- sample_pheatmap.pdf
|-- sample_region_celltype_gene.xls
|-- sample_region_celltype.xls
|-- sample_RNA.celltype_umap.pdf
|-- sample_space_anno.rds
|-- sample_space_anno_umap.pdf
|-- sample_space.region.pdf
|-- sample_SpatialDimPlot_Pit1_eachcell.pdf
`-- SpatialDimPlot.pdf


1. sample_space_umap.pdf
空间数据集中的 cluster 的UMAP图（左图） 和 在空间切片上的展示图（右图）。

2. sample_celltype_count.xls
（1）celltype：RNA数据集的细胞类型；
（2）count：对应细胞类型的细胞数目。

3. sample_celltype_umap.pdf
RNA数据集细胞类型可视化展示UMAP图

4. sample_region_celltype_gene.xls
（1）celltype：RNA数据集的细胞类型；
（2）第二列-：空间数据集cluster和对应细胞类型共有marker基因数量。

5. *_top3_gene_FeaturePlot.pdf
RNA数据集的各细胞类型同空间数据集共有的 top3 marker基因在空间转录组中的表达分布图

6. sample_mia_result.xls
（1）region：空间转录组中的的cluster编号；
（2）term：RNA转录组中的细胞类型；
（3）pvalue：P值；
（4）Enrichment：富集值；
（5）Depletion：缺失值；
（6）final_value：富集值和缺失值中的最大值；
（7）final_class：该细胞类型在该cluster中的状态，是富集还是缺失。

7. sample_pheatmap.pdf
展示不同cluster在各个细胞类型中的P值情况，横轴表示不同cluster，纵轴表示对应不同细胞类型的P值。

8. sample_MIA.pdf
将enrichment 和Depletion两个值在热图中展示，可以看出每个cluster最可能的细胞类型，括号里面的genes数值为高可变基因数量。

9. sample_region_celltype.xls
（1）region：空间转录组中的的cluster编号；
（2）term：预测的细胞类型。

10. sample_space_anno_umap.pdf
使用MIA预测的细胞类型，空间数据集中的细胞类型的UMAP图（左图） 和 在空间切片上的展示图（右图）。

11. sample_SpatialDimPlot_Pit1_eachcell.pdf
每个细胞类型分别在空间转录组中的切片图展示。


### 注意事项
[Para]
RNAcelltype = celltype        #RNA注释信息列
spacecluster = seurat_clusters   #空间cluster列
rna_marker_gene_logfc.threshold = 0.25  #RNA marker基因筛选指标，fc > 0.25
rna_marker_gene_padj = 0.05   #RNA marker基因筛选指标，矫正后的p < 0.05
rna_pct.1_pct.2=0.2   #RNA marker基因筛选指标，pct1 - pct2 差值 > 0.2
space_marker_gene_logfc.threshold = 0  #空间 marker基因筛选指标
space_marker_gene_padj = 0.1 #空间 marker基因筛选指标
space_pct.1_pct.2=0.05  #空间 marker基因筛选指标
common_gene=1000  #RNA和空间相同基因的数量
