tmpdir=$(dir  $(abspath $(firstword $(MAKEFILE_LIST))))
BIN=$(tmpdir)../
include $(BIN)../software/software.txt


project_id=$(shell echo $(project) | awk -F'_' '{print $$1}')

QC_dir=$(outdir)/QC
check_email=$(BIN)/QC/check.email.txt
filter_cell=$(outdir)/Seurat_Filter/*/*_filter_stat_cell.csv

HELP:
	@echo this script is used to do QC for 10X Genomics
	@echo usage:
	@echo "the makefile outdir= project= sample_list= QC_Basic"
	@echo "the makefile outdir= Filter_Stat"
	@echo 

.PHONY:Filter_Stat
Filter_Stat:
	echo Filter_Cell_Stat start at `date`
	mkdir -p $(outdir)/QC/filter
	$(PYTHON3) $(BIN)/QC/filter_stat_summary.py -c $(filter_cell) -o $(outdir)/QC/filter/Cell_filter_stat.xls
	echo Filter_Cell_Stat end at `date`

.PHONY:QC_Basic
QC_Basic:
	echo QC start at `date`
	mkdir -p $(QC_dir)/Basic
	if [ -s $(QC_dir)/Basic/sample.list ] ;\
	then \
	rm -f $(QC_dir)/Basic/sample.list ;\
	fi
	echo 'for name in `ls $(outdir)/CellRanger_Count/*/outs/metrics_summary.csv`;do name1=$${name%%/outs/metrics_summary.csv} ;file=$${name1#*CellRanger_Count/};echo $$file >>$(QC_dir)/Basic/sample.list;ln -snf $(outdir)/CellRanger_Count/$$file/outs/metrics_summary.csv $(QC_dir)/Basic/$$file.metrics_summary.csv.tmp;ln -snf $(outdir)/CellRanger_Count/$$file/stat/$$file\_filter_cells_stat.csv $(QC_dir)/Basic/$$file\_filter_cells_stat.csv ;paste -d , $(QC_dir)/Basic/$$file.metrics_summary.csv.tmp $(QC_dir)/Basic/$$file\_filter_cells_stat.csv >$(QC_dir)/Basic/$$file.metrics_summary.csv;done' >$(QC_dir)/Basic/rename_metrics.sh 
	sh $(QC_dir)/Basic/rename_metrics.sh
	$(PYTHON3) $(BIN)/QC/merge.metrics_summary.py -i $(QC_dir)/Basic/*metrics_summary.csv -o $(QC_dir)/Basic/analysis_summary.xls.tmp
	$(PERL) $(BIN)/QC/get_presummary.pl $(QC_dir)/Basic/analysis_summary.xls.tmp $(QC_dir)/Basic/analysis_summary.xls $(QC_dir)/Basic/analysis_summary_short.xls
	#$(PYTHON3) $(BIN)/QC/qc.py -i $(QC_dir) -t $(BIN)/QC/10XGenomics_template.xls -n $(BIN)/QC/10XGenomics_qc.list -p $(project_id) -o $(QC_dir)/Basic/summary -e yes -d $(QC_dir) -m ~/.email/.email.txt -s $(QC_dir)/Basic/sample.list
	echo QC end at `date`

QC_Check:
	@echo "analysis qc check begin"
	#$(PYTHON3) $(BIN)/QC/qc.py -i $(QC_dir) -t $(BIN)/QC/10XGenomics_template.xls -n $(BIN)/QC/10XGenomics_qc.list -p $(project_id) -pn $(project_name) -o $(QC_dir)/Total -e yes -m $(check_email) -s 0 -d $(QC_dir) --check -f
	@echo "analysis qc check end"

QC_Total:
	@echo "All of QC is begin!"
	mkdir -p $(QC_dir)/Total
	#$(PYTHON3) $(BIN)/QC/qc.py -i $(QC_dir) -t $(BIN)/QC/10XGenomics_template.xls -n 0 -p $(project_id) -s 0 -o $(QC_dir)/Total -e yes -m ~/.email/.email.txt -d $(QC_dir) -f 
	@echo "All of QC is finished!"

QC_path_release:
	if [ -s $(QC_dir)/Total_check.log ] ;\
	then \
		if [ `cat $(QC_dir)/Total_check.log | grep 指控结果 | cut -f 2 -d '='` == '质控通过' ] ; \
		then \
			echo "质控通过" ; \
		else \
			echo "质控不通过" ; \
		fi ; \
	fi
	make -f $(BIN)/release_cellranger/mk_release indir=$(outdir) release_dir=$(outdir)/release_cellranger project_id=$(project_id) Release
