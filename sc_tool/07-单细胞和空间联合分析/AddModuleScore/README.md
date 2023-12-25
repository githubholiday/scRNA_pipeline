### 模块： MK_AddModuleScore

*模块功能：
*模块版本：v1.0.0
*邮箱： mengli@genome.cn

### 使用示例及参数说明：
Usage:
	"make -f MK_AddModuleScore config= outdir= rnards= spacerds= prefix= info= spot_cell"
参数说明：
	"config: [文件|可选]  模块配置文件，和软件相关参数，默认为/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/limeng/16_Space_RNA/seurat/AddModuleScore/MK_AddModuleScore//config/config.txt "
	"outdir: [路径|必需]  输出路径 "
	"rnards: [文件|必需]  带有细胞注释结果的scRNA-seq rds 文件 "
	"spacerds: [文件|必需]  空间转录组的 rds 文件 "
	"prefix: [字符|必需]  文件前缀名称 "
	"info: [文件|必需]  分析的参数配置文件 "
target说明：
	"spot_cell：通过AddModuleScore函数注释空间转录组spot"
	"report：整理目录并生成report"

### 运行环境、软件及数据库：
	北京 sge 集群
	软件：R /annoroad/data1/bioinfo/PMO/yaomengcheng/bk_Anaconda3/envs/monocle3/bin/Rscript
	输入文件大小：200M
	申请CPU：6
	申请内存：80G
	实际CPU：4
	实际内存：40G
	运行时长：5min

### 输入文件示例
见test/input/

./
|-- config.ini #分析参数的配置文件
|-- normal_celltype_SCT.rds #RNA数据集的rds文件
`-- NP183765_0.8.rds #空间数据集的rds文件

### 输出文件示例
见test/output/
./
├── readme.doc 
└── sample
    ├── sample_all_markers.xls
    ├── sample_celltype_AddModuleScore.xls
    ├── sample_celltype_count_space_cluster.p*
    ├── sample_celltype_count_space_cluster.xls
    ├── sample_celltype_count_space.xls
    ├── sample_celltype_count.xls
    ├── sample_celltype_FeaturePlot.p*
    ├── sample_celltype_space_umap.p*
    ├── sample_celltype_umap.p*
    ├── sample_celltype_Vlnplot.p*
    ├── sample_cluster_umap.p*
    ├── sample_SpatialDimPlot.p*
    ├── sample_SpatialDimPlot_Pit1_eachcell.p*
    └── sample_spot_cell.rds

文件解释说明：
sample_all_markers.xls：RNA数据集中每种细胞的P< 0.05 的marker基因信息。
sample_celltype_AddModuleScore.xls：每个spot的AddModuleScore打分值统计。
sample_celltype_count_space_cluster.p*：空间数据集中每个cluster中的各个细胞类型的占比图。横轴为空间数据集的cluster，纵轴为各个细胞类型的占比，不同颜色代表不同的细胞类型。
sample_celltype_count_space_cluster.xls：各个细胞类型在空间数据集不同cluster中的数量。
sample_celltype_count_space.xls：各个细胞类型在空间数据集中的数量。
sample_celltype_count.xls：RNA数据集中的各个细胞类型的数量。
sample_celltype_FeaturePlot.p*：每种细胞类型在spot中的score值用Featureplot图展示。Score值越大，颜色越深。
sample_celltype_space_umap.p*：空间数据集中不同细胞类型的umap 图，不同颜色代表不同的细胞。
sample_celltype_umap.p*：RNA数据集中不同细胞类型的umap 图，不同颜色代表不同的细胞。
sample_SpatialDimPlot.p*：空间数据集的细胞类型结果的切片图。
sample_SpatialDimPlot_Pit1_eachcell.p*：空间数据集的每种细胞类型结果的切片图。
sample_spot_cell.rds：加入细胞类型信息的空间数据集的rds文件。

### 注意事项
RNAcelltype = celltype #RNA数据集的细胞类型的标签
spacecluster = seurat_clusters #空间数据集的要分析的分类标签
marker_gene_min.pct = 0.1 #FindAllMarkers函数中的min.pct数值
marker_gene_logfc.threshold = 0.25 #FindAllMarkers函数中的logfc.threshold数值
marker_gene_padj = 0.05 #marker基因的筛选p_val_adj小于该数值的基因
marker_gene_num = 10 #绘图展示前10的marker 基因
common_gene = 1000 #RNA数据集和空间数据集中的共有的最少细胞数

