Monocle3_multi_mk:
*模块功能：拟时间分析，支持使用单样本的rds文件或多样本合并后的rds文件，同时支持10XRNA转录组和空间转录组。
*模块版本： V1
*邮箱： mengli@genome.cn

###使用示例：
	make -f Monocle3_multi_mk difffile= outdir= prefix= configini= Monocle3

###软件：
	R：
		路径：/annoroad/data1/bioinfo/PMO/yaomengcheng/bk_Anaconda3/envs/monocle3/bin/Rscript
		版本：R version 4.0.3

###运行环境：
	北京的sge集群

###输入参数：
	difffile： 输入10x比较分析的rds文件，10x多样本合并后的rds文件：Integrating/*/2_clusters/*_immune_combined.rds
	outdir：   输出目录
	configini: Monocle3的部分参数文件,包含seurat_clusters、cores、resolution参数描述，示例：/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/limeng/04_10xRNA/Monocle3_all/script/config.ini
	prefix：   输出文件前缀名称
### 资源消耗
	requests.cpu = 10
	requests.memory = 50G

### 运行时长
	单样本：4h左右

### 输入文件示例
/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/RNA/10XGenomics/XS05KF2022120037/XS05KF2022120037-01/xijinguan/Analysis-test/Analysis/DiffExpr/KIB1/KIB1_all_analysis.rds

### 输出文件
见test/single_rds/result/
sample
├── 1_trajectory
│   ├── sample_cell_type.fine.partition.pdf
│   ├── sample_learn_graph_cells.no.use_partition.pdf
│   ├── sample_learn_graph_cells.pdf
│   └── sample_order_cells_cells.no.use_partition.pdf
├── 2_track_gene
│   ├── sample_Genes_Jitterplot.pdf
│   ├── sample_Trajectory_DEG_genes.example.xls
│   └── sample_Trajectory_DEG_genes.xls
├── 3_cogene_model
│   ├── sample_Genes_Module.pdf
│   └── sample_Genes_Module.xls
└── readme.doc



### 输出文件说明
1_trajectory：
1.1*_cell_type.fine.partition.p*
首先对seurat的结果，使用 monocle3 进行分群（partition），该步骤对其他结果没有影响。图中点代表细胞，左图不同颜色代表seurat的聚类结果，右图不同颜色代表monocle3的分群（partition）结果。
1.2*_learn_graph_cells.no.use_partition.p*
使用partition，进行轨迹学习，使用learn_graph（）函数。partition就是将所有的细胞亚群连成一个轨迹还是自己根据上步中partition的结果，每个partition单独生成轨迹。图为seurat分群和monocle3的分群轨迹图（不进行partition分群）。
1.3*_learn_graph_cells.pdf
不使用partition，进行轨迹学习，即自己根据上步中partition的结果，每个partition单独生成轨迹。轨迹图中不同样色为seurat分群和monocle3的分群轨迹图。
1.4*_order_cells_cells.no.use_partition.pdf
默认选择了一个根节点，该图为发育拟时间值的分布。一般颜色最深的为发育的起始位点。左图代表partition的发育时间图，右图代表no partition的发育时间图。

2_track_gene：
同时，通过对时间轨迹中所有差异基因的表达模式进行研究，找到具有相似趋势的基因然后分组以查看它们的共同点。使用graph_test分析获得莫兰指数（morans_I）,其值在-1至1之间，0代表此基因没有空间表达效应，1代表此基因在空间距离相近的细胞中的表达值高度相似，挑选的top10的基因进行展示。
2.1*_Trajectory_DEG_genes.*.xls
差异基因表格展示。
2.2*_Genes_Jitterplot.png
Top10的差异基因表达趋势图。

3_cogene_model：
Monocle3同时可以进行寻找基因共表达关系分析，使用Louvain community analysis计算,有些模块高度特定于其他分区，而其他模块则可能跨多个分区共享。
3.1*_Genes_Module.p*
为细胞共表达基因模块热图，图中横坐标为seurat聚类，纵坐标为模块名称。
3.1*_Genes_Module.xls
为细胞共表达基因模块表格，表格详细情况如下。

表格文件说明：
a)	*_Trajectory_DEG_genes.*.xls
[1]	gene short_name      具有差异的基因的名称。
[2]	p_val                 基因差异的p值，值越小越具有差异。
[3]	morans_test_statistic 莫兰指数（morans_I）的T检验结果。
[4]	morans_I              莫兰指数（morans_I）,其值在-1至1之间，0代表此基因没有空间表达效应，1代表此基因在空间距离相近细胞中的表达值高度相似。
[5]	status               基因的分析时状态，一般为ok。
[6]	q_val                 p值的校正值，也可认为是FDR，该值越小差异越显著。

b)	*_Genes_Module.xls
[1]	id                基因名称。
[2]	module            模块名。
[3]	supermodule       super模块名。
[4]	dim_1             Louvain community analysis计算的维度值1。
[5]	dim_2             Louvain community analysis计算的维度值2。



如有其他疑问，可以参考以下链接：
https://cole-trapnell-lab.github.io/monocle3/docs

### 注意事项

