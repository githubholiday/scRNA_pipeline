Stat_marker_mk:
*模块功能：判别cell_marker.txt文件中marker是否有表达；统计marker在cluster中的表达情况;统计出每个cluster可能的细胞类型。兼容10x RNA 转录组和10x 空间转录组的rds文件。
*模块版本： V1
*邮箱： mengli@genome.cn

###使用示例：
	make -f Stat_marker_mk outdir= rds= marker= prefix= configini= force= group= stat_marker

###软件：
	R：
		路径：/annoroad/data1/bioinfo/PMO/yaomengcheng/bk_Anaconda3/envs/monocle3/bin/Rscript
		版本：R version 4.0.3
	python3:
		路径：/annoroad/share/software/install/Python-3.3.2/bin/python3
		版本：Python-3.3.2

###运行环境：
	北京的sge集群

###输入参数：
	outdir：输出目录【必须】
	rds：整合后的rds文件【必须】
	marker：cell_marker.txt【必须】
	prefix：文件前缀，建议为组合名，不允许有空格！【必须】
	configini：每个比较组的config.ini文件【必须】
	group: 组合名字，必须是：处理组/对照组 【必须】
	force: 如果部分marker基因没有表达，是否继续执行后续分析，yes:继续分析；no:不继续分析【必须】

### 资源消耗
	requests.cpu = 2
	requests.memory = 15G
随着cluster数量的增加而增加
### 运行时长
	1h左右

### 输入文件示例
/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/limeng/04_10xRNA/Cellident/test/test_Stat_marker/Stat_marker.sh

### 输出文件
见../../test/test_Stat_marker/result

prefix
├── cell_ident.xls
├── cell_marker_notexist.xls
├── cluster_cell_marker.xls
├── cluster_marker.xls
├── prefix_top_dotplot.pdf
├── prefix_top_dotplot.png
├── prefix_VlnPlot.pdf
├── prefix_VlnPlot.png
└── readme.doc

### 输出文件说明
cell_marker_notexist.xls: 客户提供的marker不存在列表
(1)CellType：细胞类型
(2)markers not exist：不存在的marker 基因
 包含大小写匹配

cell_ident.xls： cluster预测细胞类型表
(1)Cluster：细胞分群
(2)CellType：该细胞分群可能含有的细胞类型

cluster_marker.xls： marker在cluster中的表达情况
(1)p_val：差异富集p-value；
(2)avg_log2FC：差异表达倍数log2的平均值；大于0的话是up基因，小于0的话是down 基因。
(3)pct.1：在该类中有表达的细胞比例；
(4)pct.2：除该类外其它类中有表达的细胞比例；
(5)p_val_adj：校正后的p-value；
(5)cluster：该基因所在的聚类；
(6)gene：基因名称；

cluster_cell_marker.xls  marker基因在不同cluster中的表达情况

*top_dotplot.p*: 客户提供的marker基因在不同cluster中的dotplot图；
点的大小编码一个cluster中表达特定基因的阳性细胞的百分比，而颜色反映cluster中特定基因的平均表达水平。

*_VlnPlot.p*：客户提供的marker基因在不同cluster中的小提琴图；
横轴为不同的cluster，纵轴为是不同的marker基因的表达量。


### 注意事项
