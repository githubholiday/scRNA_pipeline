# 10X Integration
* 模块功能：整合分析后进行聚类，marker，注释
* 模块版本： v0.0.1
* 作者：leiguo
* 邮箱：leiguo@genome.cn

### 软件环境
* R: v4.0.3
* Seurat: v4.1.1

### 资源消耗
跟据cluster个数，查找marker基因时间不同。
* 25G,  1cpu,  ~ 120m

### 使用方法
make -f this_make rds_dir= cmp_name= integrating_dir= config_file= Integrating
target：
* Integrating：聚类，marker

参数：
* rds_dir: 多样本合并后的rds文件 【必选|路径】；
* cmp_name: 比较组名称，示例A_VS_B 【必选|字符串】；
* integrating_dir：合并rds输出目录 【必选|路径】；
* config_file：比较组及相关参数配置文件 【必选|路径】。

make -f this_make infile= species= result_dir= cmp_name= celltype anno_plot
target：
* celltype：人/小鼠 scibet注释
* anno_plot：scibet注释画图

参数：
* infile：聚类后的rds 【必选|路径】；
* species：物种，human / mouse 【必选|字符串】；
* result_dir：输出目录 【必选|路径】；
* cmp_name: 比较组名称，示例A_VS_B 【必选|字符串】；


### 输入文件示例
input/
config.ini
1_QC/P_VS_H.rds  #按比较组合并，且整合分析后rds

config.ini示例：
[sample]
sample1 = P1/P2/P3/H1/H2/H3   #样本名称
sample2 = P/P/P/H/H/H         #样本对应分组名称

[cmp]
cmp1 = P/H   #比较组

[Para]
object_list_normalization.method = LogNormalize   #标准化方法
object_list_scale.factor = 10000
object_list_nfeatures_findvariablefeatures = 2000
object_list_findvariablefeatures_method = vst     #整合分析方法
qc_pca_plot_w_h = 12,8                            #图片长宽
sct = no
integration_cca_dims = 20    #维度数
integration_pca_dims = 20    #维度数
integration_runpca_npcs = 30
reduction_dims_num = 20     
reduction_resolution = 0.8   #聚类分辨率
reduction_w_h = 24,8         #聚类图长宽
marker_gene_min.pct = 0.1    
marker_gene_logfc.threshold = 0.25   #marker logfc值
marker_gene_test.use = wilcox        #marker检验方法
de_gene_logfc.threshold = 0.25       #差异logfc值 
de_gene_test.use = wilcox            #差异检验方法
de_gene_min.pct = 0.1


### 输出文件示例
output/
1_QC/
	P_VS_H_pca.qc.pdf  #pca降维质控
2_clusters
	AverageExpression.xls  #基因平均表达量
	AverageExpression_heatmap.pdf  #平均表达量热图
	stas.celltype.csv	#各样本cluster细胞数
	*cluster_groups*	#分组聚类图
	*cluster_samples*	#样本聚类图
	P_VS_H_cellType.stats.pdf  #各样本cluster细胞数占比图
3_marker
	P_VS_H_all.markers.csv  #比较组所有marker基因结果
	P_VS_H_all_clusters.pdf  #marker基因在所有cluster中表达小提琴图
	P_VS_H_top_dotplot.pdf   #marker基因top点状图
	P_VS_H_top10_marker_heatmap.pdf
	P_VS_H_top1markers_FeaturePlot.pdf
	P_VS_H_umap_cluster_anno.pdf
CellIdent
	P_VS_H_scibet_stat.xls #scibet注释细胞类型结果
	P_VS_H_marker_heatmap.pdf  #注释marker基因热图
	P_VS_H_celltype_anno_stat.xls
	P_VS_H_celltype_predict.csv
	P_VS_H_query_model_heatmap.pdf
    
### 输出文件
1. 处理_对照_pca.qc.pdf 
此图为降维聚类质控结果图，一共5幅图，第一幅图为多个样品合并后，PCA降维后，前15pc主成分的热图；第二幅图为多个样品降维后，在pca散点图中的分布；第三幅图为前两个pc中，各个基因的贡献度；第四幅为pc；第五幅图，x轴为pca主成分数，纵坐标为标准偏差，一般会选择在拐点处的pc数目进行下游分析。
2. 处理_对照_umap_cluster_samples.pdf
降维聚类后的umap结果图，左图为map图中不同的颜色代表不同的样品，右图中不同颜色代表不同的亚群。
3.cell.tsne.csv
列名的意义：
	orig.ident	细胞barcode
	tSNE_1	TSNE降维分析后的x轴坐标
	tSNE_2	TSNE降维分析后的y轴坐标
	sample	sample样品来源
	cluster	细胞所属亚群编号
降维后tsne图的坐标轴结果文件，可以用excel表格打开。
4.cell.umap.csv
降维后umap图的坐标轴结果文件，可以用excel表格打开。
列名的意义：
     orig.ident  细胞barcode
     UMAP_1  UMAP降维分析后的x轴坐标
     UMAP_2  UMAP降维分析后的y轴坐标
     sample  sample样品来源
     cluster 细胞所属亚群编号
5.处理_对照_cellType.stats.pdf
每个样品中不同亚群的分布统计结果图
6.stas.celltype.csv
每个亚群中不同样品的细胞数分布结果文件，可以用excel表格打开。
####marker
1.处理_对照.markers.csv
所有亚群中marker基因结果，一般用于细胞亚群命名用。
2.处理_对照_all_clusters.pdf
按照avg_logFC排序后，每个亚群中top1基因的小提琴图，这里只显示了top1，如果想要更多的基因图，可以以后单独画。
3.处理_对照_top1markers_FeaturePlot.pdf
按照avg_logFC排序后，每个亚群中top1基因的表达分布图。红色代表高表达，灰色代表低表达。
4.marker_violin/  top1基因单个小题图
	4.1基因名称_clusters.pdf
	某个基因的单个小提琴图，基因在不同亚群中的表达分布。
5.处理_对照.all.markers.add_significant.csv
所有亚群中marker基因的显著性结果文件，其实就是在【处理_对照.markers.csv】文件中增加了是否显著性标签：yes or no
6.处理_对照.all.markers.anno.xls
所有亚群中marker基因的各个数据库的注释结果，可以查看某个基因的功能或者代谢通路
7.处理_对照_cluster_anno
	7.1 处理_对照_cluster聚类编号.anno.xls
	按照每个亚群对其marker基因进行注释结果文件。
8.处理_对照.all.markers.csv  文件每列说明
第一列表头为空       这一列为基因名称，
p_val                第七列亚群（第七列为亚群编号）的细胞与除了此亚群所有细胞比较的，统计学检验的p值。
avg_logFC            第七列亚群（第七列为亚群编号）的细胞与除了此亚群所有细胞平均表达量差异倍数的自然对数值，其实这里进行了一些数据处理。
pct.1                第七列亚群（第七列为亚群编号）的细胞中表达该基因的细胞数占总细胞数的比例。
pct.2                除了第七列亚群（第七列为亚群编号）的细胞以外其他所有细胞中表达该基因的细胞数占所有其他细胞数的比例。
p_val_adj            第七列亚群（第七列为亚群编号）的细胞与除了此亚群所有细胞比较的，统计学检验的p值的校正值。
cluster              比较的细胞亚群编号
gene                 比较分析的基因名称，一般与第一列一致，但是也有不一致的情况，比如该参考基因组中，基因symbol有重复的，软件将会自动的在基因名称后面添加一个.数字，数字是根据出现次数增加的。

###scibet
scibet结果说明:(目前仅支持对人和小鼠进行细胞注释)
输出结果一共有三个文件：
*_celltype_predict.csv 
*_marker_heatmap.pdf/png 
*_query_model_heatmap.pdf/png

文件内容说明：

1. *_celltype_predict.csv，此文件为细胞类型预测结果文件（每个细胞），一共三列：
第一列：列名为空，一般为序列号，没有什么其他意义
第二列：Barcode，细胞的barcode，区分细胞的。
第三列：seurat_clusters，为seurat输出结果的亚群名称
第四列：prd_celltype 预测的细胞类型。

2. *_marker_heatmap.pdf
挑选出前100个特异基因，用z-score气泡图来展示这些基因在所有亚群中的表达情况。
 
3. *_query_model_heatmap.pdf
细胞亚群预测结果热图，每一列为seurat输出结果的亚群名称（需要鉴定的细胞类型），每一行为细胞亚群的名称，颜色接近黄色，表示属于该细胞类型的可能性越高。

4. *scibet_stat.xls ,scibet亚群统计的结果
第一列Cluster，为聚类类别。
其余列为细胞数量最多的前三种细胞类型及对应的数量（括号中为该类细胞数占该cluster细胞数的比例）。
若小于三种则展示所有类型。

5. *celltype_anno_stat_example.xls/*celltype_anno_stat.xls 细胞亚群注释整合表
基于seurat聚类的准确性和scibet的亚群自动鉴定，再结合积累的marker基因的信息，将seurat的聚类结果、scibet的注释结果以及给定的细胞对应marker信息整合成一个表格，整合的方法是选取scibet中每个cluster注释上细胞数最多的celltype作为该cluster的细胞类型。
第一列Cluster,为聚类类别。
第二列Anno 为注释后的细胞类别，
第三列Expected_Marker 为预期细胞类别对应的marker基因数，(如果未提供预期细胞类别文件，则该列均为0)
第四列Covered_Marker 为覆盖到的marker基因数（如果未提供预期细胞类别文件，则该列均为0），
第五列Ratio为覆盖到的基因占预期细胞类别maker基因的比例。

5. *tsne_cluster_anno.p*/*umap_cluster_anno.p*
注释后的细胞tsne图谱/umap图谱
细胞亚群注释后tsne的图谱,左边为各细胞在样本中的分布展示，右边为各细胞在细胞类群上的分布

如果是差异比较组分析会有以下两个文件：
6. *cellType.stats.p* 样本中不同聚类细胞占比图
据聚类分析的结果，可以研究每个样本中不同聚类的细胞占比情况，以下对每个样本进行了聚群细胞占比统计，并绘制了聚类占比的堆叠柱形图。

7. *cellType.clusterstats.p* 聚类中不同样本细胞占比图
同时各个不同的聚类中，其样本的占比情况也可以绘制成占比堆叠图，以清晰的展示不同样本或者不同聚类的占比情况。

### 注意事项
无

