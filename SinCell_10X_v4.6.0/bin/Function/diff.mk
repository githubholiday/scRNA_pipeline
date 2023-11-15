tmpdir = $(dir  $(abspath $(firstword $(MAKEFILE_LIST))))
BIN = $(tmpdir)/..

ifdef config
	include $(config)
else
	include $(tmpdir)/config/config.txt
endif
include $(BIN)/../species/$(species).txt

de_dir=$(indir)/Integrating/$(compare)/5_diff_gene_condition
 
#
shell_dir=$(outdir)/shell/Function_shell/
go_shell_dir=$(outdir)/shell/Function_shell/GO_shell/
go_shell=$(go_shell_dir)/go.sh
kegg_shell_dir=$(outdir)/shell/Function_shell/KEGG_shell/
kegg_shell=$(kegg_shell_dir)/kegg.sh

#投递任务的最大任务数
maxjob=50

Annotation_diff:
	echo Annotation for DIFF start
	mkdir -p $(de_dir)
	for clst in `ls $(DE)`;do name=` basename $$clst | sed 's/\.csv//' ` ;\
		dir=`dirname $$clst` ;\
		$(PERL) $(BIN)/Function/anno/add_sig_DE.pl -result $$clst -tsv $(GENE) -p $(pval) -q $(qval) -lgfc $(lgfc) -out $$dir/$$name\_symbol.xls; \
		$(PYTHON3) $(tmpdir)/add_annotation.py -i $$dir/$$name\_symbol.xls -ar $(REF_ANNOTATION) -g $(GTF) -out $$dir/$$name\_symbol.anno.xls -c 0 0 ;done
	$(PYTHON3) $(BIN)/Function/anno/DEnumber_stat_summary.py -d $(de_file) -s $(compare) -o $(de_dir)/$(compare).de_stat.xls
	echo Annotation for DIFF end	

GO:
	echo "#############GO analsyis start at "`date`
	[ -d $(go_shell_dir) ] || mkdir -p $(go_shell_dir)
	$(PYTHON3) $(BIN)/Function/generate_shell.py -i $(de_file) -o $(go_shell) -a go -s $(species) -c $(compare) 
	echo '$(QSUB_SGE) --queue sci.q --lines 1 --maxjob $(maxjob) --resource "vf=20G -l p=2" --maxcycle 1 --quota 1000G --jobidstart 0 $(go_shell) && echo GO analysis finished ' > $(go_shell_dir)/go_qsub.sh
	$(QSUB_SGE) --queue sci.q --lines 1 --maxjob $(maxjob) --resource "vf=20G -l p=2" --maxcycle 1 --quota 1000G --jobidstart 0 $(go_shell) && echo GO anaysis finished
	echo "#############GO analsyis end at "`date`

KEGG:
	echo "#############KEGG analsyis start at "`date`
	[ -d $(kegg_shell_dir) ] || mkdir -p $(kegg_shell_dir)
	$(PYTHON3) $(BIN)/Function/generate_shell.py -i $(de_file) -o $(kegg_shell) -a kegg -s $(species) -c $(compare) -ca $(category)
	echo '$(QSUB_SGE) --queue sci.q --lines 1 --maxjob $(maxjob) --resource "vf=20G -l p=2" --maxcycle 1 --quota 1000G --jobidstart 0 $(kegg_shell) && echo KEGG analysis finished '> $(kegg_shell_dir)/kegg_qsub.sh
	$(QSUB_SGE) --queue sci.q --lines 1 --maxjob $(maxjob) --resource "vf=20G -l p=2" --maxcycle 1 --quota 1000G --jobidstart 0 $(kegg_shell) && echo KEGG analysis finished
	echo "#############KEGG analysis end at" `date`
