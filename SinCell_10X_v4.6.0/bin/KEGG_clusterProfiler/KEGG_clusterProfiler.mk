Bin=$(shell dirname $(abspath $(firstword $(MAKEFILE_LIST))))/
config=$(Bin)/config/config.txt
Scr=$(Bin)/script/
include $(config)

HELP:
	@echo "进行KEGG富集分析:"
	@echo "标准版（需联网或添加KEGG.db数据库），当前SGE不可用:"
	@echo "make -f KEGG_clusterProfiler.mk genelist= head=[T or F]  species= outdir= prefix=  title= KEGG_clusterProfiler"
	@echo "非标准版："
	@echo "make -f KEGG_clusterProfiler.mk genelist= head=[T or F]  species=no outdir= prefix= term2gene= title= KEGG_clusterProfiler"
	@echo "参数说明："
	@echo "genelist： [file | must]  用于富集分析的文件，第1列是基因list，富集分析时会对一列去重"
	@echo "head: 	  [String | choice |default T]  genelist是否有表头, 默认T, 有表头，且基因列表头名称必须为Gene,如果有up和Down的区分，对应表头必须为Up/Down" 
	@echo "species:   [String | choice ] 当用标准库时需要提供拉丁名， 具体可参考$(Scr)/OrgDb.list 第3列的写法，当不在标准库里或者不提供时，则按非标处理 "
	@echo "outdir:   [String | must ]  输出路径"
	@echo "prefix:   [String | must ] 输出文件前缀，一般是样本名或组名等"
	@echo "term2gene:   [String | choice ] 参考基因组建库时ko.list 文件，第列是基因名称，第2列是KO号， 当用非标建库时必填"
	@echo "title:   [String | choice ] 绘图时表头文件，默认和prefix一致"
	@echo "number:   [String | choice ] 绘图时用到的最显著的条目，默认30"
	@echo "padjust:   [String | choice |default 0.05 ] 富集分析时的显著阈值，默认0.05"
	@echo "maxgene:   [String | choice |default 10000 ] 用于富集分析的最大基因数，默认10000"
	@echo "network:   [String | choice |default T] KEGG的富集分析是否联网，默认连网T, 如果为F，则调用local KEGG.db"
	@echo "配置参数："
	@echo "ORGDB:   [String | default ] clusterProfiler 里面的标准物种库"
	@echo "KO2MAP:   [String | default ] KEGG数据库更新时 KO与map的对应关系，第一列是ko号，第2列是以|分割的map编号"
	@echo "MAP2NAME:   [String | default ] KEGG数据库更新时 map文件，第1列是map编号，第2列是map名称, 需要指定到物种信息"
	@echo "category:[String | default ] 分析时的所属大类：all, plant，animal,archaea,bacteria,fungi,fungus,protists"
	@echo "config:   [String | default ] 软件等配置文件，默认config/config.txt"
	@echo "  "
	@echo "联网进行pathview 绘图："
	@echo "make -f KEGG_clusterProfiler.mk geneinfo= maplist=  species= kegg_dir= map_dir= prefix= column= ORGDB= Map_pathview"
	@echo "geneinfo: 必须有表头，第一列是基因，其他列有绘图信息的文件，类似de.report.xls, 绘图信息必须是数字（比如up是1，down是-1，或者FC值等）,需要提供绝对路径"
	@echo "column: [String | choice ] geneinfo中第几列用于map绘图分析，默认是第2列"
	@echo "maplist: [String | must ] 必须有表头，第一列是mapID，会对该文件内的mapID进行绘图，一般是KEGG_clusterProfiler 产生的*enrichment.xls 文件，需要提供绝对路径"
	@echo "map_dir: [String | must ] 输出路径， 必须进行该路径，进行绘图 "
	@echo "KEGG_MAP_DIR: [String | choice ] kegg数据库各物种通路图，如果该路径下没有指定map的png和xml文件， 则需要联网下载，如果存在则不需要下载 ,默认统一存放，也可以自己指定，需要提供绝对路径"
	@echo "species:   [String | must ] 用标准库时需要提供拉丁名， 具体可参考$(Scr)/OrgDb.list 第3列的写法，当不在标准库里则不能进行map绘制"
	@echo "prefix:   [String | must ] 输出文件前缀，一般是样本名或组名等"

title =$(prefix)
head=T
network=T
number=30
padjust=0.05
maxgene=10000
category=all

.PHONY:KEGG_clusterProfiler
KEGG_clusterProfiler:
	@echo "KEGG_clusterProfiler starts at "`date`
	mkdir -p ${outdir}
	[ ! -e "$(term2gene)" ] || $(PYTHON3) $(Scr)/get_kegg.py -gk $(term2gene) -km $(KO2MAP) -o $(outdir)/$(prefix).map2gene.list
	$(RSCRIPT) $(Scr)/KEGG.R -i $(genelist) -s $(species) -o $(outdir)/$(prefix).kegg.clusterProfiler.result -t $(head) -g $(outdir)/$(prefix).map2gene.list -n $(MAP2NAME) -p $(padjust) -m $(maxgene) -e $(network) -c $(ORGDB)
	grep -E "Significant|yes" $(outdir)/$(prefix).kegg.clusterProfiler.result >$(outdir)/$(prefix).kegg.enrichment.xls
	@echo "KEGG_clusterProfiler draw plot start"
	$(RSCRIPT) $(Scr)/enrichment_plot.R -i $(outdir)/$(prefix).kegg.enrichment.xls -t $(title) -p $(outdir)/$(prefix).kegg -n $(number) -f F -y "Description"
	if [ "$(head)" = "F" ] ; then \
		awk 'BEGIN{print "Gene"}{print $$1}' $(genelist) >$(outdir)/$(prefix).genelist ;\
	else \
		cat $(genelist) >$(outdir)/$(prefix).genelist ;\
	fi 
	$(GLIBC) $(Scr)/generate_kegg -d $(outdir)/$(prefix).genelist -k $(outdir)/$(prefix).kegg.clusterProfiler.result -c $(outdir)/$(prefix).kegg.clusterProfiler.result.id -o $(outdir)/$(prefix).kegg.report.xls 
	$(PYTHON3) $(Scr)/anno.py -i $(outdir)/$(prefix).kegg.report.xls -o $(outdir)/kegg.example.report.xls
	@echo "###### KEGG clusterProfile  ends at "`date`

column=1
.PHONY:Map_pathview
Map_pathview:
	@echo "map pathview starts at "`date`
	mkdir -p $(map_dir)
	cd $(map_dir) && $(RSCRIPT) $(Scr)/pathway.R -p $(prefix) -g $(geneinfo) -s $(species) -m $(maplist) -c $(ORGDB) -k $(KEGG_MAP_DIR) -n $(column)
	@echo "map pathview starts at "`date`
