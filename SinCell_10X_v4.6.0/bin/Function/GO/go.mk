#species=hg19
tmpdir=$(dir  $(abspath $(firstword $(MAKEFILE_LIST))))
include $(tmpdir)/../../../software/software.txt
include $(tmpdir)/../../../species/$(species).txt


BIN=$(tmpdir)..
comma:=,
cmp=$(subst $(comma),_,$(compare))
go_dir = $(outdir)/GO/
# GO_Candidate
#de_report = $(indir)/cmp/$(tool)/$(cmp)/$(cmp).report.xls
go_candidate=$(go_dir)/go_enrich.list


# GO_Up_Down
de_list=$(go_dir)/de.list
up_list = $(go_dir)/up.list
down_list = $(go_dir)/down.list
out_prefix = $(go_dir)/$(cmp)
go_path = $(GO_Base)/go.path
go_alias = $(GO_Base)/go.alias
go_class = $(GO_Base)/go.class
blastx_dir = $(dir  $(abspath $(firstword $(GO_annotate))))
blastx_gene = $(blastx_dir)/../gtf/uniprot/gene.uniq.blastx
uniport_anno = $(Uniport_Base)/sprot.anno
go_list = $(go_dir)/go.list

HELP:
	@echo usage:
	@echo help:
	@echo 'GO Analysis'

GO_Candidate:
	echo go candidate start at `date`
	[ -d $(go_dir) ] || mkdir -p $(go_dir)
	grep -v -w 'no' $(de_report) > $(go_candidate)
	sed -i 's/Gene_ID/Gene/g' $(go_candidate)
	sed -i 's#\tUp\t#\tup\t#g' $(go_candidate)
	sed -i 's#\tDown\t#\tdown\t#g' $(go_candidate)
	echo go candidate end at `date`

GO_Enrich:
	echo go enrich start at `date`
	[ -d $(go_dir) ] || mkdir -p $(go_dir)
	if [ `wc -l $(go_candidate) | awk '{print $$1}'` -gt 1 ];\
		then \
			make -f $(BIN)/GO_clusterProfiler/GO_clusterProfiler.mk genelist=species=no term2gene=$(GO_annotate) prefix=$(cmp) outdir=$(go_dir) genelist=$(go_candidate) config=$(BIN)/config/config.txt GO_clusterProfiler ;\
		else \
			echo "no de gene , do not analysisGO_Enrich " ;\
		fi;
	echo go enrich end at `date`

GO_Up_Down:
	echo go up_down plot start at `date`
	mkdir -p $(go_dir)
	grep -w 'yes' $(de_report) | cut -f1 | sort -u > $(de_list)
	grep -w 'yes' $(de_report) | grep -w 'up' | cut -f1 | sort -u > $(up_list)
	grep -w 'yes' $(de_report) | grep -w 'down' | cut -f1 | sort -u > $(down_list)
	$(PYTHON3) $(tmpdir)/go_with_gene.py -c $(go_class) -p $(go_path) -s $(go_alias) -g $(GO_annotate) -i $(de_list) -o $(out_prefix) >$(go_dir)/go.warning
	$(PYTHON3) $(tmpdir)/go_up_down.py -p $(go_path) -i $(out_prefix).go_with_gene -u $(up_list) -d $(down_list) -o $(out_prefix).Up_Down
	$(PYTHON3) $(tmpdir)/UpDownGenes_GoEnrichment.py -f $(out_prefix).Up_Down.xls -o $(go_dir) -s $(cmp)
	$(TRANSPOSE) -t $(out_prefix).Up_Down.xls > $(out_prefix).Up_Down.xls.t
	cat $(out_prefix).Up_Down.xls.t |sed /Percent/d |awk 'BEGIN{FS="\t";OFS="\t"}{if($$1=="Ontology")a=$$0;else{print $$0;}}END{print a;}' > $(out_prefix).Up_Down.xls.tmp
	$(PYTHON3) $(BIN)/plot_html/Bar_plot_go.py -i $(out_prefix).Up_Down.xls.tmp -t $(cmp) -y "Number of Genes" -x "" -o $(out_prefix).Up_Down.html -m group -W 1200 -H 600
	rm $(out_prefix).Up_Down.xls.t $(out_prefix).Up_Down.xls.tmp
	echo got up_down plot end at `date`
