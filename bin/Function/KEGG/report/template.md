@@@@clusterProfiler
## 富集分析
@@@@KEGG_clusterProfiler
### KEGG富集分析
KEGG（Kyoto Encyclopedia of Genes and Genomes，京都基因与基因组百科全书）是基因组破译方面的数据库。在给出染色体中一套完整基因的情况下，它可以对蛋白质交互（互动）网络在各种各样的细胞活动过程起的作用做出预测。KEGG的PATHWAY数据库整合当前在分子互动网络（比如通路、联合体）的知识，GENES/SSDB/KO数据库提供关于在基因组计划中发现的基因和蛋白质的相关知识，COMPOUND/GLYCAN/REACTION数据库提供生化复合物及反应方面的知识。
其中基因数据库（GENES Database）含有所有已知的完整基因组和不完整基因组。有细菌、蓝藻、真核生物等生物体的基因序列，如人、小鼠、果蝇、拟南芥等等；通路数据库（PATHWAY Database）储存了基因功能的相关信息，通过图形来表示细胞内的生物学过程，例如代谢、膜运输、信号传导和细胞的生长周期；配体数据库（LIGAND Database）包括了细胞内的化学复合物、酶分子和酶反应的信息。
    采用ClusterProfiler（T Wu et.al, 2021）计算目标基因中显著富集的map通路。通过KEGG功能显著性富集分析能确定候选基因行使的主要生物学功能。 设定padjust<0.05 为显著性阈值。候选基因KEGG通路富集统计结果示例见下表：
![富集分析表]{{table_kegg_report}}
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

KEGG功能富集分析得到的通路注释结果下载链接：
![富集分析表]{{download_kegg_report}}

@@@@KEGG_enrich
选取最显著的30个pathway通路（如果不足30则用全部通路）用条形图展示KEGG富集分析结果如下：
![富集分析图]{{image_kegg_bar}}
纵坐标表示通路名称，横坐标表示富集到该通路的基因数量，颜色表示padjust，颜色越红表示越显著。
选取最显著的30个pathway通路（如果不足30则用全部通路）用气泡图展示KEGG富集分析结果如下：
![富集分析图]{{image_kegg_dot}}
纵坐标表示通路名称，横坐标表示富集到该通路的基因数量占总基因的比例，颜色表示padjust，颜色越红表示越显著；气泡大小表示富集到该通路的基因数量，气泡越大表示基因数量越多。

### 参考文献
T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141
