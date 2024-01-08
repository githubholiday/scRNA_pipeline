mk_Correct:

*模块功能：基于clusterProfiler进行GO富集分析
*模块版本： V1
*负责人：任雪
*邮箱： xueren@genome.cn

### 使用示例：
	标准版：
	make -f GO_clusterProfiler.mk genelist= head=[T or F]  species= outdir= prefix= GO_clusterProfiler
	非标版：不需要联网，可通用
	make -f GO_clusterProfiler.mk genelist= head=[T or F]  species=no outdir= prefix= term2gene= GO_clusterProfiler

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
term2gene:   [String | choice ] 参考基因组建库时go.list 文件，第列是基因名称，第2-N列是GO号， 当用非标建库时必填
title:   [String | choice ] 绘图时表头文件，默认和prefix一致
number:   [String | choice ] 绘图时用到的每个类别最显著的条目，默认10， 3类总共计划绘图30条。
padjust:   [String | choice |default 0.05 ] 富集分析时的显著阈值，默认0.05
maxgene:   [String | choice |default 10000 ] 用于富集分析的最大基因数，默认10000
config:   [String | default ] 软件等配置文件，默认config/config.txt

### 数据库配置参数：
GO2NAME:   [String | default ] GO数据库更新时 GO描述文件，3列，第1列GO编号，第2列描述，第3列ONTOLOGY
ORGDB:   [String | default ] clusterProfiler 里面的标准物种库

GO2NAME=/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/FTP/Database/GO/20221108/data/go.class.clusterprofile
ORGDB=/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/FTP/Database/ClusterProfile/OrgDb.list

### 使用说明：
1） 该模块支持用clusterProfiler标准库分析， species提供的必须在标准库内（目前支持 homo_sapiens,mus_musculus），见填写说明。 此时可以不提供term2gene
2） 当species不在标准库中时 则采用非标准模式分析，必须提供term2gene。
3)  如果明确使用非标库，speceies=no 即可。


### 资源消耗
	基因列表文件：669行
	申请CPU：1
	申请内存：20G
	实际内存：15G

### 运行时长
	单组合运行时长：13min

### 输入文件
genelist：
必须保证第1列是基因信息，如果有表头，第1列表头必须是Gene，如果有Up和Down列，该列对应的表头必须是 Up/Down。在满足基因需求的基础上支持多种类型输入
1）无表头的，1列基因list文件
2）有表头的多列文件，第1列表头必须是Gene，可以包含Up/Down信息，则report文件里会统计up和down的基因数目
3）无表头的多列文件

term2gene：
参考基因组建库时go.list 文件，第列是基因名称，第2-N列是GO号， 当用非标建库时必填
ENSMUSG00000082683	GO:0071013	GO:0001650
ENSMUSG00000054256	GO:0005737


### 输入文件示例
genelist=test_std/APAP_control.sig.xls
term2gene=/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/Public/Database/GenomeDatabase/animal/Mus_musculus/Mus_musculus.GRCm38.91/RNA/annotation/go/go.list

### 输出文件
主要输出文件为：

|-- prefix.genelist 同输入的genelist，输入无表头时会添加上表头，用于生成report文件。   [中间文件]
|-- prefix.go.barplot.pdf   富集分析柱状图  [交付文件]
|-- prefix.go.barplot.png   富集分析柱状图  [交付文件]
|-- prefix.go.clusterProfiler.result  clusterProfiler go 富集分析原始结果  [中间文件] 
|-- prefix.go.dotplot.pdf   富集分析气泡图  [交付文件]
|-- prefix.go.dotplot.png   富集分析气泡图  [交付文件]
|-- prefix.go.enrichment.xls  显著富集条目输出文件（默认 padjust<0.05）[交付文件，报告不展示]
|-- prefix.go.report.xls    富集分析report文件，包括不显著的结果 [交付文件]
|-- prefix.term2gene.list  非标准模式下，生成的基因和GO term对应关系文件，用于clusterProfiler分析 [中间文件]
`-- go.example.report.xls  结果示例文件 [交付文件]


### 输出文件示例
.test_local
|-- APAP_control.genelist
|-- APAP_control.go.barplot.pdf
|-- APAP_control.go.barplot.png
|-- APAP_control.go.clusterProfiler.result
|-- APAP_control.go.dotplot.pdf
|-- APAP_control.go.dotplot.png
|-- APAP_control.go.enrichment.xls
|-- APAP_control.go.report.xls
|-- APAP_control.term2gene.list
`-- go.example.report.xls


其中，*enrichment.xls各列信息如下：
(1) ONTOLOGY : GO Term的分类
(2) ID : GO Term的ID
(3) Description：GO Term的描述
(4) GeneRatio：富集到该Term的基因数目/用于富集分析的基因数目;
(5) BgRatio: 富集到该Term的的背景基因数目/用于富集分析时的背景基因数目;
(6) pvalue：检验后的p值；
(7) p.adjust：BH方法校正后的p值；
(8) qvalue：检验后的q值；
(9) geneID：富集到该Term的基因ID；
(10) Count：富集到该Term的基因个数；
(11) Significant：该Term是否显著富集，yes，为显著；no，为不显著。


*report.xls各列信息如下：
（1）ID：GO Term的ID
（2）Ontology：该Term 所属分类
（3）Description：GO Term的描述
（4）Count1：富集到该Term的基因数目；
（5）Count2：用于富集分析的基因数目；
（6）Count3：富集到该Term的的背景基因数目；
（7）Count4：用于富集分析时的背景基因数目；
（8）pval：检验后的p值；
（9）p.adjust：BH方法校正后的p值；
（10）qval：检验后的q值；
（11）*Gene：富集到该Term上的基因
（12）*Count：富集到该Term上的基因数目
（13）Links：该GO Term的数据库链接；
（14）Result：该Term是否显著富集，yes，为显著；no，为不显著。

*barplot.png GO富集条形图
纵坐标表示GO条目，横坐标表示富集到该条目的基因数量，颜色表示padjust，颜色越红表示越显著。

*dotplot.png GO富集气泡图
纵坐标表示GO条目，横坐标表示富集到该条目的基因数量占总基因的比例，颜色表示padjust，颜色越红表示越显著；气泡大小表示富集到该条目的基因数量，气泡越大表示基因数量越多。


### 注意事项
1） 目前仅对每个类型最显著的10个进行绘图，一共30个
2） 当没有显著富集结果时，则不会产生图片。
3)  只有当genelist 里面存在一列 名叫 Up/Down 时，report中才会存在Up_Gene 和 Up_Count,以及Down_Gene 和 Down_Count
4)  不同形式的输入文件genelist， 结果只有Up和Down 是否存在的区分，没有其他区别。
5） 3个分类的结果是一起展示，如果某个分类的结果不存在，则不展示
