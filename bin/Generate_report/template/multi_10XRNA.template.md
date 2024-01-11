@@@@project_info
# 项目信息
## 基本思想
10x Genomics的Chromium Controller是基于微流控技术的单细胞分选系统，针对单个样本可以获得5,00-10,000的单细胞油包水液滴，随后配套试剂可用于制备单细胞转录组文库，最后通过测序实现对高通量单细胞数据进行细胞群体的分类以及细胞群体间基因表达的差异分析。10X Genomics单细胞表达分析可应用的细胞类型众多，如肿瘤细胞、免疫细胞、干细胞、神经细胞、生殖细胞、胚胎细胞等，是研究肿瘤异质性、免疫细胞群体功能和胚胎发育过程的优秀方法。     



## 实验流程
### 文库制备
首先获得单细胞悬浮液，并对其进行活性检测，细胞活性>85%且单分散性好杂质含量低为宜，随后进行微流控通道油包水液滴文库的制备及上机操作。为保证单个细胞的有效捕获效率，一般要求细胞浓度尽量控制在700-1200cell/uL范围内。文库制备过程中，单个Gel Bead与单个细胞被单独的油包水液滴包裹形成GEM，每个Gel Bead上有独特的Barcode和UMI序列以及用于起始逆转录反应的Poly-dT引物序列;随后，在该GEM反应体系中，细胞发生破碎裂解，释放出mRNA并与Poly-dT引物序列在逆转录酶作用下起始逆转录反应，生成Full cDNA。接着cDNA进行扩增并构建文库，单细胞分选和逆转录过程示意图如下：     



![文库制备流程图]{{image_count_1}}
### 上机测序
文库构建完成后，先使用Qubit3.0进行初步定量，稀释文库至1ng/uL，随后使用Agilent 2100对文库的insert size进行检测，insert size符合预期后，使用StepOnePlus Real-Time PCR System 荧光定量PCR仪进行Q-PCR，对文库的有效浓度进行准确定量（文库有效浓度指标为不低于10nM），以保证文库质量。     



质量合格的文库用$(PLATFROM)平台进行测序。测序策略为$(SEQ)。     



## 信息分析流程
测序所得原始下机序列（Raw reads），首先我们根据10x转录组文库结构进行数据的截取，获得包含barcode和UMI序列信息的Read1和目标序列Read2；然后利用10x官方软件CellRanger 进行数据的分析处理。     



安诺优达10x Genomics单细胞转录组测序信息分析流程主要分为四部分：CellRanger数据分析、Seurat细胞分群和差异分析、基于Seurat结果的功能注释分析和拟时间分析。     



（1）CellRanger数据分析：包含测序数据质控、质量评估和基因表达定量、细胞鉴定等分析 ；    


（2）Seurat细胞分群和差异分析：基于CellRanger获得的细胞基因表达矩阵，进行细胞再筛选、分群和差异分析等；    


（3）功能注释：主要对与找到的marker基因和cluster间差异基因进行功能注释和富集分析；    


（4）拟时间分析：通过构建细胞间的变化轨迹重塑细胞随着时间的变化过程。    



信息分析技术路线如下：     



![信息分析流程图]{{image_count_2}}
注：分析结果以实际产出为准，其中细胞类型注释分析仅对样本为人或者小鼠的进行。     



@@@@basci_QC
# 基本质控
## Rawdata统计
本项目中样本经CellRanger分析过程中的主要统计指标，如下表所示：     



![CellRanger分析统计表]{{table_count_1}}
（1）Sample：样本名称；    


（2）Estimated_Number_of_Cells：有效细胞数，有效细胞的判断，见下文附录；    


（3）Mean_Reads_per_Cell：每个细胞的平均reads数；    


（4）Median_Genes_per_Cell：每个细胞中基因数目的中位数；    


（5）Number_of_Reads：比对到参考基因组的reads数；    


（6）Valid_Barcodes：带有正确标记的reads百分比，每个标记对应到每个细胞；    


（7）Sequencing_Saturation：测序数据饱和度；    


（8）Reads_Mapped_to_Genome：比对到参考基因组的序列数百分比；    


（9）Reads_Mapped_Confidently_to_Transcriptome：比对到参考转录组的reads百分比；    


（10）Fraction_Reads_in_Cells：比对到参考基因组且来源于高质量细胞的reads百分比。    



CellRanger分析统计表下载链接：     



![CellRanger分析统计结果下载]{{download_count_1}}
CellRanger分析报告下载链接：     



![CellRanger分析报告下载]{{download_count_2}}
## Filterdata统计
CellRanger分析后，Seurat选用目前普遍设定的阈值：     



（1）过滤基因表达数目低于200，高于10000的细胞；     



（2）过滤线粒体基因（以MT开头的基因，不区分大小写）比例大于20%的细胞；     



（3）除此之外，对于人或者小鼠的样本还进行红细胞的过滤，过滤红细胞基因表达大于5%的细胞。     



统计过滤细胞后样品中每个细胞中的nCount_RNA（number of UMI），nFeature_RNA（number of Gene）和线粒体基因占nFeature_RNA的比例，结果如下图所示：     



![表达量统计图]{{image_count_3}}
注：图中小提琴图的每个点代表一个细胞，从左至右分别为基因表达数目nFeature_RNA、unique UMI总数nCount_RNA和线粒体基因所占比例percent.mt，没有线粒体注释信息的物种则线粒体基因占比为0。     



过滤细胞后样品中基因表达量统计结果下载链接：     



![过滤细胞后样品中基因表达量统计结果下载链接]{{download_count_3}}
Seurat过滤细胞数量统计表如下所示：     



![Seurat过滤细胞数量统计表]{{table_count_2}}
（1）sample：样本名称；    


（2）total_cell：总的有效细胞数；    


（3）high_quality_cell：过滤后保留的有效高质量细胞数（过滤掉以下类型的细胞：表达基因数量过高过低的细胞，线粒体基因比例大于20%的细胞，红细胞基因比例大于5%的细胞）；    


（4）low_nFeature：表达基因数低于200的细胞数（大概率是死细胞或者细胞碎片）；    


（5）high_nFeature：表达基因数高于10000的细胞数（大概率是双细胞或者多细胞）；    


（6）high_MT：线粒体基因比例高于20%的细胞数（线粒体比例高，除了特殊细胞类型外，可能是异常状态的细胞）；    


（7）high_HB：红细胞基因比例高于5%的细胞；    


（8）all_filtered_cell：所有被过滤掉的细胞数量：表达基因数量过高或过低的细胞，线粒体基因比例大于20%的细胞，红细胞基因比例大于5%的细胞。    



Seurat过滤细胞数量统计表下载：     



![Seurat分析后细胞数量统计结果下载]{{download_count_4}}
@@@@integrated_analysis
# 聚类分析
@@@@qc
## 多样本合并聚类
将不同的样品用不同的颜色在UMAP上进行展示，与聚类图比较可以看出每个样品细胞在聚类中的分布情况。本展示的聚类分辨率（resolution）为0.6。数据集越大，细胞类群越复杂，需要更大的resolution，获得细胞的类群越多。     



聚类结果使用UMAP1和UMAP2进行展示，如下图：     



![UMAP聚类结果图]{{image_count_4}}
注：左边图展示了不同样本分布情况，不同的颜色表示不同的样本；右边图展示了样本合并后聚类的情况，不同颜色表示不同的cluster。     



UMAP聚类结果图下载链接：     



![UMAP聚类结果图下载链接]{{download_count_5}}
@@@@cluster_stat
## 聚类细胞占比统计
为了进一步研究样本在不同聚类中的细胞的差异，故绘制占比堆叠图，清晰地展示样本在不同聚类中的细胞占比情况。     



样本在不同聚类中的细胞占比统计表如下：     



![样本在不同聚类中的细胞占比统计表]{{table_count_3}}
（1）Cluster：细胞聚类号；    


（2）第2列-最后：样本名，每个单元格表示样本在不同细胞类型中的细胞数（细胞占比）。    



样本在不同聚类中的细胞占比统计表下载链接：     



![样本在不同聚类中的细胞占比统计表下载链接]{{download_count_6}}
样本在不同聚类中的细胞占比堆叠柱形图如下：     



![样本在聚类中细胞占比柱形图]{{image_count_5}}
横轴为聚类，纵轴为每个聚类中各样本的细胞占有的比例。     



样本的细胞占比堆叠柱形图下载链接：     



![样本的细胞占比柱形图下载链接]{{download_count_7}}
@@@@marker_gene
## Marker基因分析
聚类分析得到了不同的聚类，为了进一步研究每个聚类的特征，鉴定出每个聚类所特有的基因，即显著差异的基因或者称之为marker基因。     



marker基因分析过程如下（以cluster1为例）：     



（1）计算cluster1中某个基因的平均表达量；     



（2）计算所有cluster中除了cluster1的对应基因的平均表达量；     



（3）采用wilcox统计方法进行差异分析，得到差异基因。     



marker基因鉴定的阈值为：     



（1）min.pct为0.10；     



（2）avg_log2FC为0.25；     



（3）avg_log2FC>0为上调（Up）, avg_log2FC<0为下调（Down）。     



统计每个细胞cluster中的marker基因，并进行功能注释，结果示例如下：     



![marker基因功能注释示例表]{{table_count_4}}
（1）Gene_ID：基因ID；    


（2）gene_name：基因名；    


（3）p_val：差异表达分析的p-value；    


（4）avg_log2FC：差异表达倍数的log值；    


（5）pct.1：基因在Cluster中样本1有表达的细胞比例；    


（6）pct.2：基因在Cluster中样本2有表达的细胞比例；    


（7）p_val_adj：校正后的p-value；    


（8）cluster：聚类号；    


（9）Up/Down：上调还是下调表达，Up上调，Down为下调；    


（10）Significant：是否为显著性差异；    


（11）NR:Seq-id：基因同NR数据库的最优比对结果；    


（12）NR:Score：基因同NR数据库的比对得分；    


（13）NR:Evalue：基因同NR数据库的比对Evalue值；    


（14）NR:Description：NR数据库中该基因的功能描述；    


（14）NT:Seq-id：基因同NT数据库的最优比对结果；    


（16）NT:Score：基因同NT数据库的比对得分；    


（17）NT:Evalue：基因同NT数据库的比对Evalue值；    


（18）NT:Description：NT数据库中该基因的功能描述；    


（19）Uniprot:UniProtKB-AC：基因同Uniprot数据库的最优比对结果；    


（20）Uniprot:Score：基因同Uniprot数据库的比对得分；    


（21）Uniprot:Evalue：基因同Uniprot数据库的比对Evalue值；    


（22）Uniprot:Description：Uniprot数据库中该基因的功能描述；    


（23）COG:gene：比对上的COG数据库中的基因名；    


（24）COG:Score：与COG数据库的比对得分；    


（25）COG:Eval：与COG数据库的比对Evalue值；    


（26）COG:num：比对上的COG数据库中的基因ID；    


（27）Pfam:pfam_ID：比对上的蛋白家族Pfam的基因ID；    


（28）Pfam:pfam_Name：比对上的蛋白家族Pfam的基因名；    


（29）Pfam:pfam_Description：比对上的蛋白家族Pfam的功能描述；    


（30）GO:biological_process：注释到的描述生物进程的GO Term；    


（31）GO:cellular_component：注释到的描述细胞组分的GO Term；    


（32）GO:molecular_function：注释到的描述分子功能的GO Term；    


（33）KEGG:KO：注释到的KEGG中的ID；    


（34）KEGG:Description：KEGG中的功能描述。    



marker基因功能注释表下载链接：     



![marker基因功能注释下载链接]{{download_count_8}}
选择各cluster的前10的marker基因绘制表达量热图，如下：     



![marker基因表达量热图]{{image_count_6}}
注：横轴为不同的细胞聚类，纵轴为每个细胞聚类的前10个marker基因。黄色表示高表达，紫色表示低表达。     



marker基因表达量热图下载链接：     



![marker基因在不同细胞表达的分布图]{{download_count_9}}
为了更为清晰的展示marker基因在各个聚类的表达情况，选取每个cluster的top1的marker基因绘制小提琴图和气泡图进行展示。     



marker基因表达量小提琴图如下：     



![marker基因的小提琴图]{{image_count_7}}
注：横轴为不同的聚类，纵轴为marker基因的表达量。通过此图，可以直观的比较marker基因在不同cluster中的表达量变化情况。     



marker基因的小提琴图下载链接：     



![marker基因的小提琴图下载链接]{{download_count_10}}
marker基因表达量气泡图如下：     



![marker基因气泡图]{{image_count_8}}
注：横轴为不同的cluster，纵轴为marker基因的表达量。点的颜色表示表达量高低，从蓝到红表示表达量从低到高，即越红表示表达量越高。点的大小表示某cluster中有该基因表达的细胞占比，点越大，说明该基因表达的细胞占比越高。     



marker基因气泡图下载链接：     



![marker基因气泡图下载链接]{{download_count_11}}
@@@@celltype
# 细胞类型注释
@@@@heatmap
## 亚群预测热图
本分析利用软件scibet对比较分析中的细胞类型进行鉴定，scibet软件可以对人和鼠的细胞类型进行大类的注释，提供一定参考，如果需要更加精细的注释可以提供基因marker进行手动注释或者使用其他软件辅助注释。     



细胞亚群预测结果热图如下：     



![细胞亚群预测结果热图]{{image_count_9}}
注：横轴为聚类号，纵轴为预测的细胞类型名。颜色从紫色到黄色，表示属于该聚类预测为该细胞类型的可能性逐渐增大。     



细胞亚群预测结果热图下载链接：     



![细胞亚群预测结果热图下载链接]{{download_count_12}}
@@@@celltype_stat
## 亚群注释结果统计
基于Seurat聚类的准确性和scibet的亚群自动注释，再结合积累的marker基因的信息，将Seurat的聚类结果、scibet的注释结果以及给定的细胞对应marker信息整合成一个表格，整合的方法是选取scibet中每个cluster注释上细胞数最多的celltype作为该cluster的细胞类型。如果其中细胞类型与给定的markerlist类型一致，则判断该marker列表与cluster covered marker的个数及具体的基因名称，并在表格中展示；如果celltype不一致，则展示对应cluster的top10的marker基因。     



scibet亚群统计表:     



![scibet亚群统计表]{{table_count_5}}
（1）CLUSTER：聚类号；    


（2）CellTypeN：细胞数量最多的前三种细胞类型；    


（3）CellTypeN_Count：细胞类型对应的细胞数量（括号中为该类细胞数占该cluster细胞数的比例）。    



scibet亚群统计表下载链接：     



![scibet亚群统计表下载链接]{{download_count_13}}
@@@@umap
## 亚群注释UMAP图谱
基于以上亚群注释统计的CellType1类绘制细胞UMAP图谱:     



![注释后的细胞UMAP图谱]{{image_count_10}}
注：细胞亚群注释后UMAP的图谱,左边为各细胞在样本中的分布展示，右边为各细胞在细胞类群上的分布。     



细胞亚群注释后图谱下载链接：     



![细胞亚群注释后图谱下载链接]{{download_count_14}}
@@@@celltype_cluster_stat
## 聚类细胞占比统计
为了进一步研究样本在不同细胞类型中的差异，故绘制占比堆叠图，清晰地展示样本在不同细胞类型中的的细胞占比情况。     



样本在不同细胞类型中的细胞占比统计表如下：     



![样本在不同细胞类型中的细胞占比统计表]{{table_count_6}}
（1）celltype：细胞类型名；    


（2）第2列-最后：样本名，每个单元格表示样本在不同细胞类型中的细胞数（细胞占比）。    



样本在不同细胞类型中的细胞占比统计表下载链接：     



![样本在不同细胞类型中的细胞占比统计表下载链接]{{download_count_15}}
样本在不同细胞类型中的细胞占比柱形图如下：     



![样本在细胞类型中的细胞占比柱形图]{{image_count_11}}
横轴为细胞类型，纵轴为每个细胞类型中各样本的细胞占有的比例。     



样本在细胞类型中的细胞占比柱形图下载链接：     



![样本在细胞类型中的细胞占比柱形图下载链接]{{download_count_16}}
@@@@de_analysis
# 差异分析
@@@@de_cluster
## 比较组细胞聚类统计
研究各比较组的细胞在聚类中的差异，将各比较组的细胞用不同的颜色在UMAP上进行展示。     



各比较组的UMAP聚类结果，如下图：     



![各比较组UMAP聚类结果图]{{image_count_12}}
各比较组UMAP聚类结果图下载链接：     



![合并组UMAP聚类结果图下载链接]{{download_count_17}}
各比较组在不同聚类中的细胞占比统计表如下：     



![比较组在不同聚类中的细胞占比统计表]{{table_count_7}}
（1）Cluster：细胞聚类号；    


（2）第2列-最后：样本名，每个单元格表示不同比较组在不同细胞类型中的细胞数（细胞占比）    



比较组在不同聚类中的细胞占比统计表下载链接：     



![合并组在不同聚类中的细胞占比统计表下载链接]{{download_count_18}}
比较组在不同聚类中的细胞占比堆叠柱形图如下：     



![合并组在聚类中细胞占比柱形图]{{image_count_13}}
横轴为聚类，纵轴为每个聚类中各比较组的细胞占有的比例。     



比较组的细胞占比堆叠柱形图下载链接：     



![合并组的细胞占比柱形图下载链接]{{download_count_19}}
@@@@de_celltype
## 亚群注释比较组统计
为了研究各比较组在不同细胞类型中的差异，统计不同细胞类型中的的各比较组的细胞占比。     



比较组在不同细胞类型中的细胞占比统计表如下：     



![比较组在不同细胞类型中的细胞占比统计表]{{table_count_8}}
（1）celltype：细胞类型名；    


（2）第2列-最后：样本名，每个单元格表示比较组在不同细胞类型中的细胞数（细胞占比）。    



比较组在不同细胞类型中的细胞占比统计表下载链接：     



![不同比较组组中的细胞占比统计表下载链接]{{download_count_20}}
比较组在不同细胞类型中的细胞占比堆叠柱形图如下：     



![比较组在细胞类型中的细胞占比柱形图]{{image_count_14}}
横轴为细胞类型，纵轴为每个类型中各比较组的细胞占有的比例。     



比较组在不同细胞类型中的细胞占比堆叠柱形图下载链接：     



![合并组的细胞占比柱形图下载链接]{{download_count_21}}
@@@@de_gene
## 差异基因分析
对同一聚类中的不同比较组进行差异分析，但是由于某些特异性的聚类可能只有一个样品的细胞，因此并不是所有的聚类都有结果。     



差异分析基因鉴定的阈值为：     



（1）min.pct为0.10；     



（2）|logfc|为0.25；     



（3）p_value<=0.05，q_value <=0.05；     



（4）logfc>0.25为上调（Up），logfc <0.25为下调（Down）；     



不同比较组在各聚类中的差异基因数量统计表，文件名为比较组名，统计表示例如下：     



![各聚类中不同比较组的差异基因数量统计表]{{table_count_9}}
（1）Cluster：聚类；    


（2）Up_gene：比较组在该聚类中上调差异基因数；    


（3）Down_gene：比较组在该聚类中下调差异基因数；    


（4）Total_gene：比较组在该聚类中总的差异基因数量。    



不同比较组在各聚类中的差异基因数量统计表下载链接：     



![各聚类中差异基因统计结果下载]{{download_count_22}}
对不同比较组在各聚类中的差异基因进行功能注释。结果示例如下：     



![差异基因功能注释示例表]{{table_count_10}}
（1）Gene_ID：基因ID；    


（2）gene_name：基因名；    


（3）p_val：差异表达分析的p-value；    


（4）avg_log2FC：差异表达倍数的log值；    


（5）pct.1：基因在Cluster中处理组有表达的细胞比例；    


（6）pct.2：基因在Cluster中对照组有表达的细胞比例；    


（7）p_val_adj：校正后的p-value；    


（8）Up/Down：上调还是下调表达，Up上调，Down为下调；    


（9）Significant：是否为显著性差异；    


（10）NR:Seq-id：基因同NR数据库的最优比对结果；    


（11）NR:Score：基因同NR数据库的比对得分；    


（12）NR:Evalue：基因同NR数据库的比对Evalue值；    


（13）NR:Description：NR数据库中该基因的功能描述；    


（14）NT:Seq-id：基因同NT数据库的最优比对结果；    


（15）NT:Score：基因同NT数据库的比对得分；    


（16）NT:Evalue：基因同NT数据库的比对Evalue值；    


（17）NT:Description：NT数据库中该基因的功能描述；    


（18）Uniprot:UniProtKB-AC：基因同Uniprot数据库的最优比对结果；    


（19）Uniprot:Score：基因同Uniprot数据库的比对得分；    


（20）Uniprot:Evalue：基因同Uniprot数据库的比对Evalue值；    


（21）Uniprot:Description：Uniprot数据库中该基因的功能描述；    


（22）COG:gene：比对上的COG数据库中的基因名；    


（23）COG:Score：与COG数据库的比对得分；    


（24）COG:Eval：与COG数据库的比对Evalue值；    


（25）COG:num：比对上的COG数据库中的基因ID；    


（26）Pfam:pfam_ID：比对上的蛋白家族Pfam的基因ID；    


（27）Pfam:pfam_Name：比对上的蛋白家族Pfam的基因名；    


（28）Pfam:pfam_Description：比对上的蛋白家族Pfam的功能描述；    


（29）GO:biological_process：注释到的描述生物进程的GO Term；    


（30）GO:cellular_component：注释到的描述细胞组分的GO Term；    


（31）GO:molecular_function：注释到的描述分子功能的GO Term；    


（32）KEGG:KO：注释到的KEGG中的ID；    


（33）KEGG:Description：KEGG中的功能描述。    



不同比较组在各聚类中的差异基因功能注释表下载链接：     



![同一聚类中不同样品的差异分析注释的结果下载链接]{{download_count_23}}
不同比较组在各聚类中的差异基因表达量进行展示，跟据avg_log2FC绝对值从大到小排序，取每个比较组在各聚类中的top10的差异基因表达量绘制小提琴图、箱线图。     



差异基因表达量小提琴图，如下：     



![比较组的差异基因表达量小提琴图]{{image_count_15}}
注：每个比较组有单独的小提琴图，不同的颜色代表不同的组。横轴为聚类，纵轴为基因表达量。     



不同比较组的差异基因表达量小提琴图结果下载链接：     



![同比较组的小提琴图结果下载链接]{{download_count_24}}
利用箱线图来展示两比较组之间的差异基因是否存在显著性差异，结果展示如下：     



![两组样本差异显著性箱线图]{{image_count_16}}
注：每个小图表示一个基因，横轴为不同比较组，纵轴为该基因的表达量，p值为两组差异显著性评分，一般认为小于0.05为具有显著性差异。     



注：这里的p值可能与文件的p值不一致，因为这个差异分析的p值是通过第三方R包signif单独分析的。     



不同比较组差异基因箱线图下载链接：     



![不同比较组差异基因箱线图下载链接]{{download_count_25}}
@@@@GO_clusterProfiler
## GO富集分析
基因本体（Gene Ontology，GO）是一个在生物信息学领域中广泛使用的本体，是基因功能国际标准分类体系，提供了一套动态更新的标准词汇表来描述生物体中基因和基因产物的属性，可以挖掘出所研究的生物学问题相关的生物学过程。GO分为三个Ontology，分别是：分子功能（Molecular Function，MF）、细胞组分（Cellular Component，CC）和生物过程（Biological Process，BP）。可以通过GO富集分析， 确认候选基因中是否有显著富集到特定GO条目上的基因组合，进一步研究候选基因在不同的GO层面与性状，疾病的关系，深入探索生物学分子机制。     



采用ClusterProfiler（T Wu et.al, 2021）计算目标基因中显著富集的GO条目。通过GO功能显著性富集分析能确定候选基因行使的主要生物学功能。 设定padjust<0.05 为显著性阈值。     



![差异基因GO统计结果示例表]{{table_count_11}}
（1）ID：GO Term的ID；    


（2）Ontology：该Term 所属分类；    


（3）Description：GO Term的描述；    


（4）Count1：富集到该Term的基因数目；    


（5）Count2：用于富集分析的基因数目；    


（6）Count3：富集到该Term的的背景基因数目；    


（7）Count4：用于富集分析时的背景基因数目；    


（8）pval：检验后的p值；    


（9）p.adjust：BH方法校正后的p值；    


（10）qval：检验后的q值；    


（11）*Gene：富集到该Term上的基因；    


（12）*Count：富集到该Term上的基因数目；    


（13）Links：该GO Term的数据库链接；    


（14）Result：该Term是否显著富集，yes，为显著；no，为不显著。    



差异基因GO统计结果表下载链接：     



![差异基因GO统计结果表下载链接]{{download_count_26}}
@@@@GO_enrich
选取每个类别最显著的10个GO条目（如果不足10则用该类别全部条目）用条形图展示GO富集分析结果如下：     



![GO富集分析条形图]{{image_count_17}}
GO富集分析条形图下载链接     



![GO富集分析条形图下载链接：]{{download_count_27}}
纵轴为GO条目，横轴为富集到该条目的基因数量，颜色表示padjust，颜色越红表示越显著。    



选取每个类别最显著的10个GO条目（如果不足10则用该类别全部条目）用气泡图展示GO富集分析结果如下：     



![GO富集分析气泡图]{{image_count_18}}
GO富集分析气泡图下载链接：     



![GO富集分析气泡图下载链接：]{{download_count_28}}
纵轴为GO条目，横轴为富集到该条目的基因数量占总基因的比例，颜色表示padjust，颜色越红表示越显著；气泡大小表示富集到该条目的基因数量，气泡越大表示基因数量越多。    



@@@@KEGG_clusterProfiler
## KEGG富集分析
KEGG（Kyoto Encyclopedia of Genes and Genomes，京都基因与基因组百科全书）是基因组破译方面的数据库。在给出染色体中一套完整基因的情况下，它可以对蛋白质交互（互动）网络在各种各样的细胞活动过程起的作用做出预测。KEGG的PATHWAY数据库整合当前在分子互动网络（比如通路、联合体）的知识，GENES/SSDB/KO数据库提供关于在基因组计划中发现的基因和蛋白质的相关知识，COMPOUND/GLYCAN/REACTION数据库提供生化复合物及反应方面的知识。     



其中基因数据库（GENES Database）含有所有已知的完整基因组和不完整基因组。有细菌、蓝藻、真核生物等生物体的基因序列，如人、小鼠、果蝇、拟南芥等等；通路数据库（PATHWAY Database）储存了基因功能的相关信息，通过图形来表示细胞内的生物学过程，例如代谢、膜运输、信号传导和细胞的生长周期；配体数据库（LIGAND Database）包括了细胞内的化学复合物、酶分子和酶反应的信息。     



采用ClusterProfiler（T Wu et.al, 2021）计算目标基因中显著富集的map通路。通过KEGG功能显著性富集分析能确定候选基因行使的主要生物学功能。 设定padjust<0.05 为显著性阈值。     



差异基因KEGG通路富集分析结果示例如下表：     



![差异基因KEGG统计结果示例表]{{table_count_12}}
（1）Map：kegg通路编号；    


（2）Name：kegg通路名称；    


（3）Count1：富集到该通路的基因数目；    


（4）Count2：用于富集分析的基因数目；    


（5）Count3：富集到该通路的的背景基因数目；    


（6）Count4：用于富集分析时的背景基因数目；    


（7）pval：检验后的p值；    


（8）p.adjust：BH方法校正后的p值；    


（8）qval：检验后的q值；    


（9）*Gene：富集到该通路上的基因；    


（10）*Count：富集到该通路上的基因数目；    


（11）Links：该map的数据库链接；    


（12）Result：该map是否显著富集，yes，为显著；no，为不显著。    



差异基因KEGG统计结果下载链接：     



![差异基因KEGG统计结果下载链接：]{{download_count_29}}
@@@@KEGG_enrich
选取最显著的30个pathway通路（如果不足30则用全部通路）用条形图展示KEGG富集分析结果如下：     



![KEGG富集分析条形图]{{image_count_19}}
纵轴为通路名称，横轴为富集到该通路的基因数量，颜色表示padjust，颜色越红表示越显著。    



KEGG富集分析条形图下载链接：     



![KEGG富集分析条形图下载链接：]{{download_count_30}}
选取最显著的30个pathway通路（如果不足30则用全部通路）用气泡图展示KEGG富集分析结果如下：     



![KEGG富集分析气泡图]{{image_count_20}}
纵轴为通路名称，横轴为富集到该通路的基因数量，颜色表示padjust，颜色越红表示越显著。    



KEGG富集分析气泡图下载链接：     



![KEGG富集分析气泡图下载链接：]{{download_count_31}}
@@@@wikipathway
## WikiPathway富集分析
wikiPathway数据库是由科学家维护的生物学通路数据库，包含脊椎动物、无脊椎动物、植物和微生物，且更新速度非常快。本项目利用基于R语言的clusterprofiler包进行wikiPathway富集和GSEA分析。默认仅对人和小鼠进行该部分的分析。这里取富集分析top10生物学通路进行展示。     



wikiPathway富集分析柱形图：     



![wikiPathway富集分析柱形图]{{image_count_21}}
注：富集分析柱形图。横轴为富集到该term上基因个数，纵轴为富集的term，图例为富集显著性打分值校正的P值。     



wikiPathway富集分析气泡图：     



![wikiPathway富集分析气泡图]{{image_count_22}}
注：富集分析气泡图。横轴为富集到该term上的基因比例，纵轴为富集的term，图例为基因count和p.adjust值。     



wikiPathway富集分析基因网络图：     



![wikiPathway富集分析网络图]{{image_count_23}}
注：富集分析基因网络图，中间节点为富集的term，周围节点为基因，节点大小由相关的基因个数决定。     



Gene Set Enrichment Analysis（GSEA）使用预定义的基因集（通常来自功能注释或先前实验的结果），将基因按照在两类样本中的差异表达定的基因集合是否在这个排序表的顶端或者底端富集。基因集富集分析检测基因集合而不是单个基因的表达变化，因此可以包含这些细微的表达变化。     



对两组样本间差异基因进行wikiPathway生物学通路GSEA分析，同样取top10个生物学通路进行展示。     



GSEA分析结果图：     



![GSEA分析结果图]{{image_count_24}}
注：GSEA分析结果图展示了差异基因在基因集（本项目指wikiPathway）上的富集情况。结果图分为如下三部分，第一部分为基因Enrichment Score的折线图，横轴为该基因集下的每个基因，纵轴为对应的Running ES，折线的峰值就是这个基因集的Enrichemnt score。第二部分为hit，用线条标记位于该基因集下的基因。第三部分为所有基因的rank值分布图, 默认采用Signal2Noise算法。     



GSEA分析结果图下载链接：     



![GSEA分析结果图下载链接]{{download_count_32}}
@@@@monocle3
# 拟时间分析
## 拟时间分析图谱
单细胞基因表达研究使人们能够在复杂的生物过程和高度异源的细胞群中描述转录调控过程。每个细胞都是某个转录调控过程的快照。拟时间分析使用算法来学习每个细胞所经历的基因表达变化，作为动态生物过程的一部分。通过了解基因表达变化的整体“轨迹”，就可以将每个细胞放置在轨迹中的适当位置，从而得到整个细胞群落内部亚群之间的排布方式。除了构建单细胞轨迹之外，拟时间分析还能够进行差异表达分析来揭示重要的基因和细胞。Monocle 软件是cole-trapnell 实验室研发的使用最广泛拟时间分析工具，号称单细胞分析三剑客之一。 Monocle 利用机器学习的方法（Reversed Graph Embedding）来对细胞进行排序，这个方法能够可靠而准确地解决复杂的生物学过程。 基于 Monocle3，可以使用 rds 文件在网页版中进行轨迹调整。另外其可以使用 Seurat 的 UMAP 结果进行后续分析，故本分析基于此进行。     



下图展示不同细胞类型分群结果和 monocle3 的分群结果图:     



![拟时间分布图]{{image_count_25}}
注：左图代表不同细胞类型的聚类结果，右图代表 monocle3 的轨迹分群结果，其中每个点代表一个细胞，每个颜色代表不同的类群。     



拟时间分群图下载链接：     



![拟时间分群图下载链接：]{{download_count_33}}
轨迹学习时可以选择使用分区（Partition）与否，使用 learn_graph()函数。分区指是否将所有的细胞亚群连成一个轨迹还是根据上步中分群结果，每个亚群单独生成轨迹图，以下展示两种情况下的轨迹图。     



下图为不进行分区（partition）分群时的轨迹图，不同的颜色表示monocle3分群和细胞分群情况：     



![分群轨迹图]{{image_count_26}}
图中每个点代表一个细胞，不同的颜色代表不同的分群，左图代表 monocle3 的聚类轨迹图，右图代表 Seurat 的不同细胞类型的聚类轨迹图。     



不使用partition模式轨迹图下载链接：     



![分群轨迹图下载链接：]{{download_count_34}}
下图为进行分区（partition）时的轨迹图，不同的颜色表示 monocle3 分群和 Seurat 的细胞分群情况：     



![分群轨迹图]{{image_count_27}}
图中每个点代表一个细胞，不同的颜色代表不同的分群，左图代表 monocle3 的聚类轨迹图，右图代表 Seurat 的细胞分群聚类轨迹图。     



进行分区（partition）模式轨迹图下载链接：     



![分群轨迹图下载链接：]{{download_count_35}}
如果之前没有明显的发育关系，在构建发育轨迹的时候，就尽量不要构建一个发育轨迹，建议使用分区（partition），如果是细胞有发育关系，就可以构建一个发育轨迹，就可以不使用分区（partition）。     



为了展示分析，默认选择了一个群作为起始位置，可尽量根据细胞亚群的发育关系，后期再确定真实的发育起始位置。     



下图展示了细胞轨迹从起始到终止的时间图：     



![分群轨迹图]{{image_count_28}}
图中每个点代表一个细胞，不同的颜色代表不同的分群，左图代表分区的轨迹时间图，右图代表不分区的轨迹时间图。     



轨迹时间图下载链接：     



![轨迹时间图下载链接：]{{download_count_36}}
## 共表达基因模块热图
Monocle3同时可以进行寻找基因共表达关系分析，有些模块高度特定于其他分区，而其他模块则可能跨多个分区共享。     



下图为细胞共表达基因模块热图     



![基因共表达模块热图]{{image_count_29}}
图中横轴为Seurat聚类，纵轴为模块名称。     



基因共表达模块热图下载链接：     



![基因共表达模块热图：]{{download_count_37}}
基因共表达模块表格下载链接：     



![基因共表达模块表格下载]{{download_count_38}}
# 附录
## 软件与方法说明
### 有效细胞判定说明
有效细胞的判定主要基于UMI数目分布来确定的。其评估方法如下：     



（1）根据ChromiumTM Single Cell 3’/5’Solution 的细胞捕获率，预估细胞数量N；     



（2）根据每个barcode的UMI数量降序排列，计算前N个barcode的UMI数量的99分位数，记为Y；     



（3）满足条件(UMI数量 > Y × 10%)的barcode即为有效细胞。     



### 软件版本说明
本次分析所使用的软件见下表：     



![软件列表]{{table_count_13}}
（1）No：序号；    


（2）software：软件名称或包名称；    


（3）version：软件版本或包版本；    


（4）function：软件作用    


（5）offical website：软件或包的官网。    



### 英文版方法说明
CellRanger software was applied to demultiplex the Illumina BCL output into FASTQ files. The Cell Ranger count was then applied to each FASTQ file to align reads to the reference genome and generate barcode and unique molecular identifier counts. We followed the Seurat integrated analysis and comparative analysis workflows to do all scRNA-seq analyses (Stuart et al., 2019). For quality control and filtering out low quality cells, only cells expressing more than 200 genes (defined as genes detected in at least 3 cells) and fewer than 20% mitochondrial genes were selected.All single cells passed quality control for further batch correction and unbiased clustering.     



The datasets were integrated based on ‘anchors’ identified between datasets (nfeatures=2000, normalization.method=‘SCT’) before performing linear dimensional reduction by principal-component (PC) analysis. The top 30 PCs were included in a UMAP dimensionality reduction. After obtaining the top 30 PCs, we computed the shared nearest-neighbor graph and we identified clusters in the network using the Louvain algorithm with the resolution . The UMAP (Uniform Manifold Approximation and Projection) method was used for visualization of unsupervised clustering. Differential gene expression or marker gene was determined by the ‘findMarkers’ function with the default Wilcoxon’s rank-sum test either as one versus the rest or as a direct comparison with parameters min.pct=0.1 and logfc.threshold=0.25. Cell cluster identities were determined using scibet software. We used the clusterProfiler package for differential expression gene GO and KEGG pathway annotations and enrichment analysis.     



英文版详细method下载链接：     



![Method文档下载。]{{download_count_39}}
## 常见售后问题解答
FAQ文档下载链接：     



![FAQ文档下载]{{download_count_40}}
## 参考文献
Abdi H, Williams L J. Principal component analysis[J]. Wiley interdisciplinary reviews: computational statistics, 2010, 2(4): 433-459.     



Butler A, Hoffman P, Smibert P, et al. Integrating single-cell transcriptomic data across different conditions, technologies, and species[J]. Nature biotechnology, 2018, 36(5): 411.     



Camp J G, Sekine K, Gerber T, et al. Multilineage communication regulates human liver bud development from pluripotency[J]. Nature, 2017, 546(7659): 533.     



Ding C, He X. K-means clustering via principal component analysis[C]//Proceedings of the twenty-first international conference on Machine learning. ACM, 2004: 29.     



Dobin A, Davis C A, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner[J]. Bioinformatics, 2013, 29(1): 15-21.     



Maaten L, Hinton G. Visualizing data using t-SNE[J]. Journal of Machine Learning Research, 2008, 9(Nov): 2579-2605.     



Macosko E Z, Basu A, Satija R, et al. Highly Parallel Genome-wide Expression Profiling of Individual Cells Using Nanoliter Droplets.[J]. Cell, 2015, 161(5):1202.     



Hinton G, Roweis S. Stochastic Neighbor Embedding[J]. Advances in Neural Information Processing Systems, 2010, 41(4):833--840.     



Shannon, P. et al. Cytoscape: A software environment for integrated models of biomolecular interaction networks. Genome Research 13, 2498-2504 (2003).     



Van der Maaten, L.J.P. Accelerating t-SNE using Tree-Based Algorithms. Journal of Machine Learning Research 15, 3221-3245 (2014).     



Zhong S, Zhang S, Fan X, et al. A single-cell RNA-seq survey of the developmental landscape of the human prefrontal cortex[J]. Nature, 2018, 555(7697): 524.     



Wang Y ,  Wang R ,  Zhang S , et al. iTALK: an R Package to Characterize and Illustrate Intercellular Communication.  2019.     



Cabello-Aguilar S , Alame M , Kon-Sun-Tack F , et al. SingleCellSignalR: inference of intercellular networks from single-cell transcriptomics[J]. Nucleic Acids Research, 2020.     



Patel AP, Tirosh I, Trombetta JJ, et al. Single-cell RNA-seq highlights intratumoral heterogeneity in primary glioblastoma. Science.2014;344(6190):1396-1401. doi:10.1126/science.1254257.     



La Manno G, Soldatov R, Zeisel A, et al.RNA velocity of single cells.Nature. 2018;560(7719):494-498. doi:10.1038/s41586-018-0414-6.     



Gao R , Bai S , Ying C H , et al. Delineating copy number and clonal substructure in human tumors from single-cell transcriptomes[J]. Nature Biotechnology, 2021:1-10.     



     



Pipeline Version:$(VERSION)     
