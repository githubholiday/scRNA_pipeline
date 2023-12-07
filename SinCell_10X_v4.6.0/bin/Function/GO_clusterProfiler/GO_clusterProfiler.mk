Bin=$(shell dirname $(abspath $(firstword $(MAKEFILE_LIST))))/
ifdef config
	include $(config)
else
	include $(Bin)/config/config.txt
endif
Scr=$(Bin)/script/
include $(config)

HELP:
	@echo "进行GO富集分析"
	@echo "标准库分析："
	@echo "make -f GO_clusterProfiler.mk  genelist= head=[T or F]  species= outdir= prefix= title= GO_clusterProfiler"
	@echo "非标准库分析："
	@echo "make -f GO_clusterProfiler.mk genelist= head=[T or F]  species=no outdir= prefix= term2gene= title= GO_clusterProfiler"
	@echo "参数说明："
	@echo "genelist： [file | must]  用于富集分析的文件，第1列是基因list，富集分析时会对一列去重"
	@echo "head: 	  [String | choice |default T]  genelist是否有表头, 默认T, 有表头，且基因列表头名称必须为Gene,如果有up和Down的区分，对应表头必须为Up/Down" 
	@echo "species:   [String | choice ] 当用标准库时需要提供拉丁名， 具体可参考$(Scr)/OrgDb.list 第3列的写法，当不在标准库里或者不提供时，则按非标处理 "
	@echo "outdir:   [String | must ]  输出路径"
	@echo "prefix:   [String | must ] 输出文件前缀，一般是样本名或组名等"
	@echo "term2gene:   [String | choice ] 参考基因组建库时go.list 文件，第列是基因名称，第2-N列是GO编号， 当用非标建库时必填"
	@echo "title:   [String | choice ] 绘图时表头文件，默认和prefix一致"
	@echo "number:   [String | choice ] 绘图时用到的每个类别最显著的条目，默认10"
	@echo "padjust:   [String | choice |default 0.05 ] 富集分析时的显著阈值，默认0.05"
	@echo "maxgene:   [String | choice |default 10000 ] 用于富集分析的最大基因数，默认10000"
	@echo "配置参数："
	@echo "GO2NAME:   [String | default ] GO数据库更新时 GO描述文件，3列，第1列GO编号，第2列描述，第3列ONTOLOGY"
	@echo "ORGDB:   [String | default ] clusterProfiler 标准库配置文件"
	@echo "config:   [String | default ] 软件等配置文件，默认config/config.txt"
	

title =$(prefix)
head=T
number=10
padjust=0.05
maxgene=10000
.PHONY:GO_clusterProfiler
GO_clusterProfiler:
	@echo "GO_clusterProfiler starts at "`date`
	mkdir -p ${outdir}
	[ ! -e "$(term2gene)" ] || $(PYTHON3) $(Scr)/get_go.py -i $(term2gene) -o $(outdir)/$(prefix).term2gene.list
	$(RSCRIPT) $(Scr)/GO.R -i $(genelist) -s $(species) -o $(outdir)/$(prefix).go.clusterProfiler.result -t $(head) -g $(outdir)/$(prefix).term2gene.list -n $(GO2NAME) -p $(padjust) -m $(maxgene) -c $(ORGDB)
	grep -E "Significant|yes" $(outdir)/$(prefix).go.clusterProfiler.result >$(outdir)/$(prefix).go.enrichment.xls
	@echo "GO_clusterProfiler draw plot start"
	$(RSCRIPT) $(Scr)/enrichment_plot.R -i $(outdir)/$(prefix).go.enrichment.xls -t $(title) -p $(outdir)/$(prefix).go -n $(number) -f T -y "GO Term"
	if [ "$(head)" = "F" ] ; then \
		awk 'BEGIN{print "Gene"}{print $$1}' $(genelist) >$(outdir)/$(prefix).genelist ;\
	else \
		cat $(genelist) >$(outdir)/$(prefix).genelist ;\
	fi 
	$(GLIBC) $(Scr)/generate_go -d $(outdir)/$(prefix).genelist -g $(outdir)/$(prefix).go.clusterProfiler.result -o $(outdir)/$(prefix).go.report.xls
	$(PYTHON3) $(Scr)/anno.py -i $(outdir)/$(prefix).go.report.xls -o $(outdir)/go.example.report.xls
	@echo "###### GO clusterProfile  ends at "`date`
