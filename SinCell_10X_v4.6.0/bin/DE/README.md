# 10X DE
* 模块功能：按比较组对cluster或celltype进行差异分析
* 模块版本： v0.0.1
* 作者：leiguo
* 邮箱：leiguo@genome.cn

### 软件环境
* R: v4.0.3
* Seurat: v4.1.1

### 资源消耗
* 30G,  1cpu,  ~80m

### 使用方法
make -f this_make rds= combine= result_dir= config_file= label= DE
target：
* DE: 按比较组进行差异分析

参数：
* rds: 合并rds路径 【必选|路径】；
* combine: 组合名称，多个样本的合并的名称 【必选|字符串】；
* result_dir：差异分析输出目录，以比较组创建目录 【必选|路径】；
* config_file：比较组及相关参数配置文件 【必选|路径】；
* label: 进行差异分析的列，seurat_clusters或Cell_type,默认seurat_clusters,如果指定Cell_type，避免特殊字符目录将以Celltype1、2进行命名 【选填|字符串】。

make -f this_make combine= cmp_name= result_dir= species= GENE= Annotation_diff
target：
* Annotation_diff: 差异注释及统计, pval=0.05, qval=0.05, lgfc=0.25

参数：
* combine: 组合名称，多个样本的合并的名称 【必选|字符串】；
* cmp_name: 比较组名称，示例A_VS_B 【必选|字符串】；
* result_dir：差异分析输出目录，以比较组创建目录 【必选|路径】；
* species：注释参考基因组版本号 【必选|字符串】；
* GENE：差异基因对应gene名称 【必选|路径】。

### 输入文件示例
input/
config.ini
P_VS_H_immune_combined.rds  #多样本整合分析rds

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
P1_VS_H/1_P1_VS_H/
1_P1_VS_H_diff_gene.csv    #差异基因结果
1_P1_VS_H_FeaturePlot.pdf  #差异基因在各样本中表达
1_P1_VS_H_signif_Plot.pdf  #差异基因在比较组中表达箱线图
1_P1_VS_H_VlnPlot.pdf      #差异基因在各cluster中表达小提琴图
1_P1_VS_H_diff_gene_symbol.xls  #标注显著性差异结果文件
1_P1_VS_H_diff_gene_symbol.anno.xls  #差异基因注释文件
    
### 输出文件
1. 群编号_处理_VS_对照_diff_gene.csv  两组比较的差异结果，两组是否有差异，处理条件基因表达的影响。
第一列            第一列表头为空，第一列为差异分析的基因。
p_val              两组比较的统计学分析的p值。
avg_logFC          两组比较的平均表达量的差异倍数的自然对数值。
pct.1              基因在第一个比较组表达的细胞数占比（一般为处理组）。
pct.2              基因在第二个比较组表达的细胞数占比（一般为对照组）。
p_val_adj          两组差异分析的p值的校正值。
2. 亚群编号_处理_VS_对照_FeaturePlot.pdf
某个差异基因在UMAP降维上表达的分布结果，一般对【亚群编号_处理_VS_对照_diff_gene.csv】文件按照avg_logFC进行排序，选择top10的基因进行展示。
3. 亚群编号_处理_VS_对照_signif_Plot.pdf
两组差异的top10基因通过箱线图来展示是否有差异，数值为两组比较的p值，不过这里的p值可能与文件的p值不一致，因为这个差异分析的p值是通过第三方R包signif单独分析的，与我们文件：亚群编号_处理_VS_对照_diff_gene.csv结果采用方法不一致。
4. 亚群编号_处理_VS_对照_VlnPlot.pdf
两组差异的top10基因在所有亚群中分布表达分布的小提琴图。看此图的时候，可能需要注意的是我们要看具体的某个亚群。
按比较组合并，且整合分析后rds，用于后续分析
5. 亚群编号_处理_VS_对照_diff_gene_symbol.xls
Gene_ID: 基因ID；
Gene_Symbol：基因名称
Up/Down: 上下调
Significant：是否显著

### 注意事项
无

