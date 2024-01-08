# 10X release_cellranger
* 模块功能：10x cellranger等分析结果快速整理释放
* 模块版本： v0.0.1
* 作者：leiguo
* 邮箱：leiguo@genome.cn

### 软件环境
* Python3: 3.3.2

### 资源消耗
* 5G,  1cpu,  ~10m, 跟据文件大小生成md5时间不同

### 使用方法
make -f mk_release indir= release_dir= release_conf= project_id= release= Release
target：
* Release：整理释放目录

参数：
* indir: 10x分析目录Analysis【必须|路径】
* release_dir: 输出整理目录【必须|路径】
* project_id：项目号，用于整理目录和交付【必须|项目号】
* release_conf：连接配置文件，参考upload整理方法【可选|路径】
* release: 是否交付华为云，默认不交付 【可选|yes or no】


### 输入文件示例
10x标准分析Analysis/
CellRanger_Count
Integrating
Seurat_Filter


release_conf示例：(参考upload整理方法)https://doc.weixin.qq.com/doc/w3_AOgA8gZSADsEogZTNyTQZSGRrALhS?scode=AN4ATQeCAAsDMs1cMqAE8A2AZPADc
示例：INDIR/CellRanger_Count/*/outs/cloupe.cloupe     OUTDIR/release/*1/      link    2
第1列：输入文件，一定要具体到文件，哪怕用通配符表示，比如BIN/appendix/* 表示  appendix下面所有的文件。 保留字符INDIR，BIN，分别 对应参数里的-i和-b参数
第2列： 表示输出位置，可以是文件，可以是路径，如果是路径，则必须以“/"结尾。
第3列： 处理过程， 目前只有copy 和 link 2种，表示拷贝文件，或者软连接文件
第4列： 表示文件重要程度 
	0	文件不存在，不影响该模块生成
	2	文件不存在， 脚本断掉，不会生成upload


### 输出文件示例
output/
├── Project_id_cellranger_result
│   ├── EEQ1
│   │   ├── cloupe.cloupe
│   │   ├── EEQ1.rds
│   │   ├── filtered_gene_bc_matrices.h5
│   │   ├── filter_matrix
│   │   │   ├── barcodes.tsv.gz
│   │   │   ├── features.tsv.gz
│   │   │   └── matrix.mtx.gz
│   │   ├── possorted_genome_bam.bam
│   │   ├── possorted_genome_bam.bam.bai
│   │   ├── raw_feature_bc_matrix.h5
│   │   └── raw_matrix
│   │       ├── barcodes.tsv.gz
│   │       ├── features.tsv.gz
│   │       └── matrix.mtx.gz
│   └── SQ_VS_EQ
│       └── SQ_VS_EQ.rds
├── file_not_exist.info #整理目录日志文件
├── remove_moudle.info
├── upload.log
└── warning_moudle.info

### 输出文件
该目录下的文件为CellRanger软件分析的结果文件，可直接用于seurat分析
1. 样品/filter_matrix/matrix.mtx.gz
用坐标格式存储的基因表达量结果，用三列信息来表示一个细胞中一条基因的unique UMI总数。如用来表示“AAACCTGAGGTGATAT-1”所标识的细胞中基因“ENSG00000198727”表达量为32的一行中：
（1）第1列：“3366”代表Gene在genes.tsv中的序号，对应第3366行记录的基因“ENSG00000198727”；
（2）第2列：“1”代表Cell在barcodes.tsv中的序号，对应第1行Cell Barcode“AAACCTGAGGTGATAT-1”标识的细胞；
（3）第3列：“32”代表unique UMI总数，即基因表达量为32。
由于文件较大，建议不要直接用文本编辑器打开；如需要查看感兴趣的基因及细胞，通过基因列表和barcode信息提取相应的数据。
2. 样品/filter_matrix/features.tsv.gz
记录基因ID和基因名的文件,基因所在的行数就是样品.matrix.mtx中的第一列
3. 样品/filter_matrix/barcodes.tsv.gz
记录cell barcode标记的文件，细胞barcode所在的行数就是样品.matrix.mtx中的第二列
4. 样品/raw_matrix
该目录为未过滤矩阵结果，文件参考filter_matrix介绍
5. 样品/filtered_gene_bc_matrices.h5, 样品/raw_feature_bc_matrix.h5
过滤后和原始的H5结果文件
6. 样品/cloupe.cloupe
CellRanger分析结果文件，可以导入Loupe™ Cell Browser查看。
7. 样品/possorted_genome_bam.bam(.bai)
比对bam文件及索引
8. 样品/样品.rds
seurat读取过滤后矩阵，并过滤低质量细胞后转化rds。
9. 比较组/比较组.rds
多样本整合分析后，包含聚类、marker等信息rds。


### 注意事项
无

