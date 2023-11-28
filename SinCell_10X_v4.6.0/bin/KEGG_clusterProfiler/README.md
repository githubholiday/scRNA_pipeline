mk_Correct:

*模块功能：基于clusterProfiler进行KEGG富集分析
*模块版本： V1
*负责人：任雪
*邮箱： xueren@genome.cn

### 使用示例：
	标准版（需联网或添加KEGG.db数据库），当前SGE不可用：
	make -f KEGG_clusterProfiler.mk genelist= head=[T or F]  species= outdir= prefix= KEGG_clusterProfiler
	非标版：不需要联网，可通用
	make -f KEGG_clusterProfiler.mk genelist= head=[T or F]  species=no outdir= prefix= term2gene= KEGG_clusterProfiler

### 软件：
	软件：clusterProfiler R包
	路径：export LD_LIBRARY_PATH=/opt/glibc-2.14/lib:$$LD_LIBRARY_PATH && /annoroad/data1/bioinfo/PMO/suyalei/software/Anaconda/minconda3/envs/r_4/bin/Rscript
	版本：v4.6.0

### 运行环境：
	北京SGE集群
	说明：可以在其他集群运行，必须有clusterProfiler的环境

### 输入参数：
genelist： [file | must]  用于富集分析的文件，第1列是基因list，富集分析时会对一列去重
head: 	  [String | choice |default T]  genelist是否有表头, 默认T, 有表头，且基因列表头名称必须为Gene,如果有up和Down的区分，对应表头必须为Up/Down
species:   [String | choice ] 当用标准库时需要提供拉丁名， 具体可参考/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/FTP/Database/ClusterProfile/OrgDb.list  第3列的写法，当不在标准库里或者不提供时，则按非标处理 
outdir:   [String | must ]  输出路径
prefix:   [String | must ] 输出文件前缀，一般是样本名或组名等
term2gene:   [String | choice ] 参考基因组建库时ko.list 文件，第列是基因名称，第2列是KO号， 当用非标建库时必填
title:   [String | choice ] 绘图时表头文件，默认和prefix一致
number:   [String | choice ] 绘图时用到的最显著的条目，默认30
padjust:   [String | choice |default 0.05 ] 富集分析时的显著阈值，默认0.05
maxgene:   [String | choice |default 10000 ] 用于富集分析的最大基因数，默认10000
network:   [String | choice |default T] KEGG的富集分析是否联网，默认连网T, 如果为F，则调用local KEGG.db
config:   [String | default ] 软件等配置文件，默认config/config.txt

###新增
category:[String | default ] 分析时的所属大类：all, plant，animal,archaea,bacteria,fungi,fungus,protists

### 数据库配置参数：
ORGDB:   [String | default ] clusterProfiler 里面的标准物种库
KO2MAP:   [String | default ] KEGG数据库更新时 KO与map的对应关系，第一列是ko号，第2列是以|分割的map编号
MAP2NAME:   [String | default ] KEGG数据库更新时 map文件，第1列是map编号，第2列是map名称,用物种所属大类本身的列表

MAP2NAME=/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/FTP/Database/KEGG/20221108/data/pathway_$(category).list
KO2MAP=/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/FTP/Database/KEGG/20221108/data/ko2map/ko2map.xls
ORGDB=/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/FTP/Database/ClusterProfile/OrgDb.list

### 使用说明：
1） 该模块支持用clusterProfiler标准库分析， species提供的必须在标准库内（目前支持 homo_sapiens），见填写说明。 此时可以不提供term2gene
2） 当使用标准库时， 如果环境联网则直接联网分析， 如果不联网，设置network=F,采用软件安装的KEGG.db 分析（目前还不支持）
3） 当species不在标准库中时 则采用非标准模式分析，必须提供term2gene。


### 资源消耗
	基因列表文件：669行
	申请CPU：1
	申请内存：20G
	实际内存：目前5G不够

### 运行时长
	单组合运行时长：13min

### 输入文件
genelist：
必须保证第1列是基因信息，如果有表头，第1列表头必须是Gene，如果有Up和Down列，该列对应的表头必须是 Up/Down。在满足基因需求的基础上支持多种类型输入
1）无表头的，1列基因list文件
2）有表头的多列文件， 可以包含Up/Down信息，则report文件里会统计up和down的基因数目
3）无表头的多列文件

term2gene：
参考基因组建库时ko.list 文件，第列是基因名称，第2列是KO号， 当用非标建库时必填
ENSMUSG00000082683	K11094
ENSMUSG00000054256	K14411


### 输入文件示例
genelist=test_local/APAP_control.sig.xls
term2gene=/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/Public/Database/GenomeDatabase/animal/Mus_musculus/Mus_musculus.GRCm38.91/RNA/annotation/ko/ko.list

### 输出文件
主要输出文件为：

|-- prefix.genelist 同输入的genelist，输入无表头时会添加上表头，用于生成report文件。   [中间文件]
|-- prefix.kegg.barplot.pdf   富集分析柱状图  [交付文件]
|-- prefix.kegg.barplot.png   富集分析柱状图  [交付文件]
|-- prefix.kegg.clusterProfiler.result  clusterProfiler kegg 富集分析原始结果  [中间文件] 
|-- prefix.kegg.clusterProfiler.result.id  clusterProfiler kegg 富集分析时基因ID信息 [中间文件]
|-- prefix.kegg.dotplot.pdf   富集分析气泡图  [交付文件]
|-- prefix.kegg.dotplot.png   富集分析气泡图  [交付文件]
|-- prefix.kegg.enrichment.xls  显著富集条目输出文件（默认 padjust<0.05）[交付文件]
|-- prefix.kegg.report.xls    富集分析report文件，包括不显著的结果 [交付文件]
|-- prefix.map2gene.list  非标准模式下，生成的基因和map对应关系文件，用于clusterProfiler分析 [中间文件]
`-- kegg.example.report.xls  结果示例文件 [交付文件]


### 输出文件示例
.test_local
|-- APAP_control.genelist
|-- APAP_control.kegg.barplot.pdf
|-- APAP_control.kegg.barplot.png
|-- APAP_control.kegg.clusterProfiler.result
|-- APAP_control.kegg.clusterProfiler.result.id
|-- APAP_control.kegg.dotplot.pdf
|-- APAP_control.kegg.dotplot.png
|-- APAP_control.kegg.enrichment.xls
|-- APAP_control.kegg.report.xls
|-- APAP_control.map2gene.list
`-- kegg.example.report.xls


其中，*enrichment.xls各列信息如下：
(1) ID : 通路ID
(2) Description：通路名称
(3) GeneRatio：富集到该通路的基因数目/用于富集分析的基因数目;
(4) BgRatio: 富集到该通路的的背景基因数目/用于富集分析时的背景基因数目;
(5) pvalue：检验后的p值；
(6) p.adjust：BH方法校正后的p值；
(7) qvalue：检验后的q值；
(8) geneID：富集到该通过的基因ID
(9) Count：富集到该通过的基因个数
(10) Significant：该map是否显著富集，yes，为显著；no，为不显著。


*report.xls各列信息如下：
（1）Map：kegg通路编号
（2）Name：kegg通路名称
（3）Count1：富集到该通路的基因数目；
（4）Count2：用于富集分析的基因数目；
（5）Count3：富集到该通路的的背景基因数目；
（6）Count4：用于富集分析时的背景基因数目；
（7）pval：检验后的p值；
（8）p.adjust：BH方法校正后的p值；
（8）qval：检验后的q值；
（9）*Gene：富集到该通路上的基因
（10）*Count：富集到该通路上的基因数目
（11）Links：该map的数据库链接；
（12）Result：该map是否显著富集，yes，为显著；no，为不显著。

### 注意事项
1）目前仅对显著的30个进行绘图
2） 当没有显著富集结果时，则不会产生图片。
3)  只有当genelist 里面存在一列 名叫 Up/Down 时，report中才会存在Up_Gene 和 Up_Count,以及Down_Gene 和 Down_Count
4)  不同形式的输入文件genelist， 结果只有Up和Down 是否存在的区分，没有其他区别


## 个性化分析[不放置到报告中]
### 模块功能
联网进行pathview 绘图（基于数据库中通路同，把关心的基因的相应数据以颜色的形式填充到图上）

### 使用示例
make -f KEGG_clusterProfiler.mk geneinfo= maplist=  species= kegg_dir= map_dir= prefix= column= ORGDB= Map_pathview

### 软件：
	软件：clusterProfiler，pahtview R包
	路径：export LD_LIBRARY_PATH=/opt/glibc-2.14/lib:$$LD_LIBRARY_PATH && /annoroad/data1/bioinfo/PMO/suyalei/software/Anaconda/minconda3/envs/r_4/bin/Rscript
	版本：v4.6.0

### 运行环境：
	203联网环境下，必须有clusterProfiler的环境

### 参数介绍：
geneinfo: 必须有表头，第一列是基因，其他列有绘图信息的文件，类似de.report.xls, 绘图信息必须是数字（比如up是1，down是-1，或者FC值等）,需要提供绝对路径
column: [String | choice ] geneinfo中第几列用于map绘图分析，默认是第1列, 只有基因没有绘图信息时，默认信息都是1，相当于给基因标注
maplist: [String | must ] 必须有表头，第一列是mapID，会对该文件内的mapID进行绘图，一般是KEGG_clusterProfiler 产生的*enrichment.xls 文件，需要提供绝对路径
map_dir: [String | must ] 输出路径， 必须进行该路径，进行绘图 
KEGG_MAP_DIR: [String | choice ] kegg数据库各物种通路图，如果该路径下没有指定map的png和xml文件， 则需要联网下载，如果存在则不需要下载 ,默认统一存放（/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/FTP/Database/ClusterProfile/KEGG/），也可以自己指定，需要提供绝对路径

species:   [String | must ] 用标准库时需要提供拉丁名， 具体可参考/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/FTP/Database/ClusterProfile/OrgDb.list 第3列的写法，当不在标准库里则不能进行map绘制
prefix:   [String | must ] 输出文件前缀，一般是样本名或组名等

### 输入示例：
geneinfo=test_map/APAP_geneinfo.xls
maplist =test_map/APAP_control.kegg.enrichment.xls

### 输出文件：
.
├── prefix.gene.id.xls 基因ID对应的数据库ID和基因名称
└── speID*.prefix.png 


### 输出文件示例
.
├── APAP_control.gene.id.list
├── mmu04068.APAP_control.png
└── mmu04668.APAP_control.png

### 注释事项
（1）当KEGG_MAP_DIR 里面没有mapID的图时必须联网才能下载分析，可以自己指定，联网下载。
（2）物种必须是OrgDB里面的物种，否则缺少物种对应信息
（3）geneinfo 和 maplist，以及KEGG_MAP_DIR 必须是绝对路径，因为生成的结果要放置到当前路径下
（4）绘图信息必须是数字， 比如up和down的信息转为1和-1等
（5）分析中有些通路会因为各种原因绘制失败，请知
