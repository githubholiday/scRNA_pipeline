# 10X filter_cell
* 模块功能：使用Seurat软件对cellranger结果进行细胞过滤
* 模块版本： v0.0.1
* 作者：leiguo
* 邮箱：leiguo@genome.cn

### 软件环境
* R: v4.0.3
* Seurat: v4.1.1

### 资源消耗
* 4G,  1cpu,  ~3m

### 使用方法
make -f mk_seurat_filter dataDir= output_dir= sample_id= species_tag= Filter_cell
target：
* Filter_cell：Seurat过滤细胞
* 过滤指标：
1. nFeature_RNA > 200    #细胞表达基因数
2. nFeature_RNA < 10000
3. percent.mt < 20%
4. percent.HB < 5%
同时满足以上指标的细胞将被保留。

参数：
* dataDir：cellranger稀疏矩阵结果目录【必选|路径】；
* sample_id: 样本名称【必选|字符串】；
* species_tag：物种ID，仅 人:9606 小鼠:10090过滤红细胞【必选|数字】;
* output_dir：输出目录 【必选|路径】。
* mincell: 基因表达的最小细胞数，小于该值的基因将被过滤 【选填|默认：3】；
* mingene：细胞包含最小基因数，小于该值的细胞将被过滤 【选填|默认：200】；
* maxgene：细胞包含最大基因数，大于该值的细胞将被过滤 【选填|默认：10000】；
* mito：细胞中线粒体基因表达占比，大于该值的细胞将被过滤 【选填|默认：20】；
* hb：细胞中红细胞基因表达占比，大于该值的细胞将被过滤 【选填|默认：5】。


### 输入文件示例
input/
barcodes.tsv
genes.tsv
matrix.mtx

### 输出文件示例
output/
H1_filter_cell.csv            #过滤后细胞统计信息
H1_filter_UMI.csv             #过滤后细胞UMI表达数据
H1_filter_stat_cell.csv       #各指标过滤细胞数
H1_filter_cell.rds            #过滤细胞保存rds，后续使用
H1_Gene_UMI_mito_percent.pdf  #过滤后各统计细胞分布图
    
### 输出文件
1. sample_filter_cell.csv
样品中细胞统计信息
（1）细胞barcode；
（2）orig.ident：样品名称。
（3）nCount_RNA：细胞中表达的UMI数目；
（4）nFeature_RNA：细胞中表达的基因数目；
（5）percent.HB：红细胞相关基因占比，非人和小鼠为0；
（6）percent.mt：线粒体相关基因占比。
2.样品_filter_UMI.csv
各细胞中基因的UMI表达量，每一行代表一条基因，每一列代表一个细胞。
由于文件较大，建议不要直接用文本编辑器打开；如需要查看感兴趣的基因及细胞，通过基因列表和barcode信息提取相应的数据。
3.样品_Gene_UMI_mito_percent.p*
样品表达量结果展示，小提琴图的每个点代表一个细胞，从左至右分别为基因表达数目、unique UMI总数和线粒体基因比例的结果。
4.sample_Cell_filter_stat.xls
（1）Sample：样本名称；
（2）total_cell：总的有效细胞数；
（3）remaining_cell；过滤后保留的有效细胞数（过滤掉以下类型的细胞：表达基因数量过高过低的细胞，线粒体基因比例大于20%的细胞，红细胞）；
（4）low_nFeature：表达基因数低于200的细胞数（大概率是将死细胞或者细胞碎片）；
（5）high_nFeature：表达基因数高于10000的细胞数（大概率是双细胞或者多细胞）；
（6）high_percent.mt：线粒体基因比例高于20%的细胞数（线粒体比例高，除了特殊细胞类型外，可能是异常状态的细胞）；
（7）high_HB：红细胞基因比例高于5%的细胞数；
（8）all_filtered_cell：所有被过滤掉的细胞数量：表达基因数量过高或过低的细胞，线粒体基因比例大于20%的细胞，红细胞。
5.sample_filter_cell.rds
过滤后保存rds文件用于后续分析使用。


### 注意事项
无

