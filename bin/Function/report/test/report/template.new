@@@@clusterProfiler
MainMenu:富集分析
@@@@GO_clusterProfiler
SubMenu:GO富集分析
P:#,;#基因本体（Gene Ontology，GO）是一个在生物信息学领域中广泛使用的本体，是基因功能国际标准分类体系，提供了一套动态更新的标准词汇表来描述生物体中基因和基因产物的属性，可以挖掘出所研究的生物学问题相关的生物学过程。GO分为三个Ontology，分别是：分子功能（Molecular Function，MF）、细胞组分（Cellular Component，CC）和生物过程（Biological Process，BP）。可以通过GO富集分析， 确认候选基因中是否有显著富集到特定GO条目上的基因组合，进一步研究候选基因在不同的GO层面与性状，疾病的关系，深入探索生物学分子机制。
P:#,;#采用ClusterProfiler（T Wu et.al, 2021）计算目标基因中显著富集的GO条目。通过GO功能显著性富集分析能确定候选基因行使的主要生物学功能。 设定padjust<0.05 为显著性阈值。候选基因GO统计结果示例见下表：
Table:upload/Integrating_analysis/*/diff_gene/*/cluster*/GO/*.go.example.report.xls,,,650,,0,候选基因GO统计结果示例表
PRE:
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
PRE
Excel:upload/Integrating_analysis/*/diff_gene/*/cluster*/GO/*.go.report.xls,,,候选基因GO统计结果表下载

@@@@GO_enrich
SubMenu:GO富集分析
P:#,;#选取每个类别最显著的10个GO条目（如果不足10则用该类别全部条目）用条形图展示GO富集分析结果如下：
Table:upload/Integrating_analysis/*/diff_gene/*/cluster*/GO/*.barplot.png,,,650,,0,GO富集分析条形图
Excel:upload/Integrating_analysis/*/diff_gene/*/cluster*/GO/*.barplot.p*,,,GO富集分析条形图下载
PRE:
纵坐标表示GO条目，横坐标表示富集到该条目的基因数量，颜色表示padjust，颜色越红表示越显著。
PRE
P:#,;#选取每个类别最显著的10个GO条目（如果不足10则用该类别全部条目）用气泡图展示GO富集分析结果如下：
Table:upload/Integrating_analysis/*/diff_gene/*/cluster*/GO/*.dotplot.png,,,650,,0,GO富集分析气泡图
Excel:upload/Integrating_analysis/*/diff_gene/*/cluster*/GO/*.dotplot.p*,,,GO富集分析气泡图下载
PRE:
纵坐标表示GO条目，横坐标表示富集到该条目的基因数量占总基因的比例，颜色表示padjust，颜色越红表示越显著；气泡大小表示富集到该条目的基因数量，气泡越大表示基因数量越多。
PRE

@@@@KEGG_clusterProfiler
SubMenu:KEGG富集分析
P:#,;#KEGG（Kyoto Encyclopedia of Genes and Genomes，京都基因与基因组百科全书）是基因组破译方面的数据库。在给出染色体中一套完整基因的情况下，它可以对蛋白质交互（互动）网络在各种各样的细胞活动过程起的作用做出预测。KEGG的PATHWAY数据库整合当前在分子互动网络（比如通路、联合体）的知识，GENES/SSDB/KO数据库提供关于在基因组计划中发现的基因和蛋白质的相关知识，COMPOUND/GLYCAN/REACTION数据库提供生化复合物及反应方面的知识。
其中基因数据库（GENES Database）含有所有已知的完整基因组和不完整基因组。有细菌、蓝藻、真核生物等生物体的基因序列，如人、小鼠、果蝇、拟南芥等等；通路数据库（PATHWAY Database）储存了基因功能的相关信息，通过图形来表示细胞内的生物学过程，例如代谢、膜运输、信号传导和细胞的生长周期；配体数据库（LIGAND Database）包括了细胞内的化学复合物、酶分子和酶反应的信息。
P:#,;#采用ClusterProfiler（T Wu et.al, 2021）计算目标基因中显著富集的map通路。通过KEGG功能显著性富集分析能确定候选基因行使的主要生物学功能。 设定padjust<0.05 为显著性阈值。候选基因KEGG通路富集统计结果示例见下表：
Table:upload/Integrating_analysis/*/diff_gene/*/cluster*/KEGG/*.kegg.example.report.xls,,,650,,0,候选基因KEGG统计结果示例表
PRE:
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
PRE
Excel:upload/Integrating_analysis/*/diff_gene/*/cluster*/KEGG/*.kegg.report.xls,,,候选基因KEGG统计结果表下载

@@@@KEGG_enrich
SubMenu:KEGG富集分析
P:#,;#选取最显著的30个pathway通路（如果不足30则用全部通路）用条形图展示KEGG富集分析结果如下：
Table:upload/Integrating_analysis/*/diff_gene/*/cluster*/KEGG/*.barplot.png,,,650,,0,KEGG富集分析条形图
Excel:upload/Integrating_analysis/*/diff_gene/*/cluster*/KEGG/*.barplot.p*,,,KEGG富集分析条形图下载
PRE:
纵坐标表示通路名称，横坐标表示富集到该通路的基因数量，颜色表示padjust，颜色越红表示越显著。
PRE
P:#,;#选取最显著的30个pathway通路（如果不足30则用全部通路）用气泡图展示KEGG富集分析结果如下：
Table:upload/Integrating_analysis/*/diff_gene/*/cluster*/KEGG/*.dotplot.png,,,650,,0,KEGG富集分析气泡图
Excel:upload/Integrating_analysis/*/diff_gene/*/cluster*/KEGG/*.dotplot.p*,,,KEGG富集分析气泡图下载
PRE:
纵坐标表示通路名称，横坐标表示富集到该通路的基因数量占总基因的比例，颜色表示padjust，颜色越红表示越显著；气泡大小表示富集到该通路的基因数量，气泡越大表示基因数量越多。
PRE
