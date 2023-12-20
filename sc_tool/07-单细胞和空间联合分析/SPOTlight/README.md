### 模块： mk_SPOTlight

*模块功能：使用SPOTlight软件对单细胞10x转录组和空间转录组进行联合分析。
*模块版本：v1.0.0
*邮箱： yangzhang@genome.cn

### 使用示例及参数说明：

Usage:
	 make -f mk_SPOTlight scdata= stdata= prefix= configini= outdir= scriptdir= config= SPOTlight
参数说明：
	 config: [文件|可选]  模块配置文件，和软件相关参数，默认为.//config/config.txt 
	 scdata: [文件|必需]  单细胞10X的rds，需要带有注释
	 stdata: [文件|必需]  空间转录组的rds 
	 prefix: [字符|必需]  输出结果的前缀 
	 configini: [文件|必需]  参数配置文件 
	 outdir: [路径|必需]  分析结果输出路径 
	 scriptdir: [路径|可选]  source的脚本所在路径，默认为script路径 
Usage:
	 make -f mk_SPOTlight outdir= config= report
参数说明：
	 config: [文件|可选]  模块配置文件，和软件相关参数，默认为.//config/config.txt 
	 outdir: [路径|必需]  分析结果输出路径，会在该路径下生成report目录，report目录下即为报告结果

### 输入文件示例
见test/input/
.
├── scdata.rds       单细胞10X转录组的rds，带有注释结果
└── stdata.rds       空间转录组的rds

### 运行环境及软件：
	北京238 R4.2.3（SPOTlight，scater，scran，Seurat，dplyr，NMF，ggsci，configr，ggplot2 ，scatterpie，spatialexperiment ，singlecellexperiment ，ggcorrplot ，arrangements）

### 资源消耗及运行时长
	单细胞rds：233M；空转rds：167M
	1CPU, 20G , 15min

### 输出文件示例
.
└── *
    ├── *_all_spot_celltype_example.xls
    ├── *_all_spot_celltype.xls
    ├── *_boxplot_B_cell.pdf
    ├── *_boxplot_Epithelial.pdf
    ├── *_boxplot_Myeloid.pdf
    ├── *_boxplot_Pit1_Epithelial.pdf
    ├── *_boxplot_Stromal.pdf
    ├── *_boxplot_T_cell.pdf
    ├── *_celltype_count.xls
    ├── *_CorrelationMatrix.pdf
    ├── *_CorrelationMatrix.xls
    ├── *_Interactions.pdf
    ├── *_res.rds
    ├── *_singlecell_cluster_umap.pdf
    ├── *_single_celltype_B_cell.pdf
    ├── *_single_celltype_Epithelial.pdf
    ├── *_single_celltype_Myeloid.pdf
    ├── *_single_celltype_Pit1_Epithelial.pdf
    ├── *_single_celltype_Stromal.pdf
    ├── *_single_celltype_T_cell.pdf
    ├── *_spatial_cluster_umap.pdf
    ├── *_SpatialScatterpie_1.pdf
    ├── *_SpatialScatterpie_2.pdf
    ├── *_SpatialScatterpie_3.pdf
    ├── *_SpatialScatterpie_4.pdf
    ├── *_SpatialScatterpie_5.pdf
    ├── *_SpatialScatterpie_6.pdf
    ├── *_SpatialScatterpie_7.pdf
    ├── *_SpatialScatterpie_8.pdf
    └── *_SpatialScatterpie.pdf

主要结果文件说明：
见report/README_SPOTlight.doc

### 注意事项
config.ini中各参数说明：
cellnumber = 100 从单细胞注释的每种细胞类型中选择多少个细胞作为代表细胞。选择这些代表细胞提取特征，作为该细胞类型的特点，用于解卷积。一般选择100个足够，如果某种细胞类型不足100，则全部采用。
genenumber = 3000 从单细胞注释的结果中选择多少个高可变基因，一般选择3000。
RNAcelltype = celltype 含有注释的单细胞rds中，注释信息所在的数据槽名称。
spacecluster = seurat_clusters 空间转录组中含有聚类信息的数据槽名称。
repeat_gene_num = 3000 空间转录组和10X单细胞中含有的共同基因的数目，不低于此数据，才会进行后续的分析。
mgs_col = mean.AUC 解卷积时选择的打分列名称。
mgs_col_threshold = 0.8 解卷积时选择的打分列的数据阈值，超过此值才会保留下来。
npcs = 20 单细胞rds中如果没有umap的结果，则会重新运行RunPCA，参数为此参数。
umap_dim = 20 单细胞rds中如果没有umap的结果，则会在RunPCA后面RunUMAP，参数为此参数。
pie_size = 2.5 各细胞类型占比饼图中每个spot的大小。
