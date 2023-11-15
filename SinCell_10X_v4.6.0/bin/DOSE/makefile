tmpdir = $(dir  $(abspath $(firstword $(MAKEFILE_LIST))))
include $(tmpdir)/../../software/software.txt

bindir = $(tmpdir)
species?=hsa
pcutoff = 1
qcutoff = 1
tag?=DOSE
de_dir=$(outdir)/Integrating/$(combine)/Diff_gene/$(cmp_name)
de_file=$(de_dir)/*/*diff_gene_symbol.xls
de_suffix=_diff_gene_symbol.xls
gsea_dir=$(de_dir)
ifeq ($(species_type),10090)
        species = mmu
else ifeq ($(species_type),9606)
        species = hsa
else ifeq ($(species_type),9913)
        species = cbr
else ifeq ($(species_type),10116)
        species = rat
else
        species = N
endif

Help:
	@echo "### enrich analysis  ###"
	@echo "该模块主要运行clusterProfile 和 DOSE"
	@echo "  By huiping lu"
	@echo Description:
	@echo USAGE:	
		@echo "make -f this.makefile outdir= difffile= header= species= prefix= upcutoff= downcutoff= pcutoff= qcutoff= EnrichGSE "
		@echo "make -f this.makefile outdir= difffile= header= prefix= upcutoff= downcutoff= EnrichGSE"	
		@echo "make -f this.makefile outdir= seq= project_id= project_name= tag= Arrange Web_Report"
		@echo "默认分析物种为人，默认不进行富集显著性筛选"
		@echo 
.PHONY:EnrichGSE
EnrichGSE:
	@echo Enrich analysis start at `date`
	mkdir -p $(gsea_dir)
	if [ $(species) = NO ];then\
		echo "No analysis for cell type !!!";\
	else \
		for clst in ` ls $(de_file) ` ; do name=` basename $$clst | sed 's/$(de_suffix)//' ` ;\
			if [ -f $$clst ] ;\
			then \
				cut -f1,4 $$clst|sed 's/^gene://' >$$clst.tmp	;\
				$(Rscript35) $(tmpdir)/clusterProfile.R -g $$clst.tmp -t $(header) -s $(species) -f $$name -u $(upcutoff) -d $(downcutoff) -p $(pcutoff) -q $(qcutoff) -o $(gsea_dir) -b $(bindir) ;\
			else \
				echo No file $$clst ;\
			fi ;\
		done ;\
	fi ;\
	rm -rf $(de_file).tmp
	for i in `ls $(gsea_dir)/*/WikiPathway/*.pdf` ;do name=`echo $$i |sed 's/.pdf/.png/'`; $(CONVERT) $$i $$name; done 
	@echo Enrich analysis finish at `date`

.PHONY:Arrange
Arrange:
	@echo Arrange result start at `date`
	[ -d $(outdir)/report/upload/ ] && rm -rf $(outdir)/report/upload/ || echo result_dir clear
	mkdir -p $(outdir)/report/upload
	cp $(tmpdir)/DOSE.template $(outdir)/report/
	cp -rL $(tmpdir)/Readme/* $(outdir)/report/upload/
	ln -snf $(outdir)/Disease/DO/* $(outdir)/report/upload/Disease/DO/
	ln -snf $(outdir)/Disease/DGN/* $(outdir)/report/upload/Disease/DGN/
	ln -snf $(outdir)/Pathway/* $(outdir)/report/upload/wikiPathway/
	ln -snf $(outdir)/Reactome/* $(outdir)/report/upload/Reactome/
	$(PYTHON3) $(tmpdir)/anno.py -i $(outdir)/report/upload/Disease/DO/*enrich*.xls -o $(outdir)/report/upload/Disease/DO/enrichDO.example.xls
	$(PYTHON3) $(tmpdir)/anno.py -i $(outdir)/report/upload/Disease/DO/*GSEA*.xls -o $(outdir)/report/upload/Disease/DO/GSEADO.example.xls
	$(PYTHON3) $(tmpdir)/anno.py -i $(outdir)/report/upload/Disease/DGN/*enrich*.xls -o $(outdir)/report/upload/Disease/DGN/enrichDGN.example.xls
	$(PYTHON3) $(tmpdir)/anno.py -i $(outdir)/report/upload/Disease/DGN/*GSEA*.xls -o $(outdir)/report/upload/Disease/DGN/GSEADGN.example.xls 
	$(PYTHON3) $(tmpdir)/anno.py -i $(outdir)/report/upload/wikiPathway/*enrich*.xls -o $(outdir)/report/upload/wikiPathway/enrichPathway.example.xls
	$(PYTHON3) $(tmpdir)/anno.py -i $(outdir)/report/upload/wikiPathway/*GSEA*.xls -o $(outdir)/report/upload/wikiPathway/GSEAPathway.example.xls
	$(PYTHON3) $(tmpdir)/anno.py -i $(outdir)/report/upload/Reactome/*enrich*.xls -o $(outdir)/report/upload/Reactome/enrichReactome.example.xls
	$(PYTHON3) $(tmpdir)/anno.py -i $(outdir)/report/upload/Reactome/*GSEA*.xls -o $(outdir)/report/upload/Reactome/GSEAReactome.example.xls 
	@echo Arrange result end at `date`
.PHONY:Web_Report
Web_Report:
	@echo report start at `date`
	mkdir -p $(outdir)/report/upload
	echo -e "SEQ:$(seq)\nPROJECT_ID:$(project_id)\nPROJECT_NAME:$(project_name)\nREPORT_DIR:$(outdir)/report/upload\n" > $(outdir)/report/report.conf
	echo make -f  $(Public_webreport)/makefile report_dir=$(outdir)/report projectType=$(tag) Web_Report  > $(outdir)/report/work.sh
	make -f  $(Public_webreport)/makefile report_dir=$(outdir)/report projectType=$(tag) Web_Report
