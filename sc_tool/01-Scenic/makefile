makefile_dir=$(dir $(firstword $(MAKEFILE_LIST)))
makefile_name=$(notdir $(firstword $(MAKEFILE_LIST)))
script=$(makefile_dir)/

R=/annoroad/data1/software/bin/miniconda/envs/Monocle/bin/Rscript
PYTHON3=/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/miniconda/envs/pyscenic/bin/python3

database=/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/ScRNA/01_20230901/Scenic/database/

HELP:
	@echo Description:进行Scenic分析
	@echo Usage: make [target] [parameters]
	@echo "Target:"
	@echo -e "get_expression: 从rda文件中提取表达矩阵并转化为loom文件"
	@echo -e " make -f infile=rda文件 outdir=输出目录 sample=样本名 object=rda文件中的样本名 get_expression"
	@echo -e "\nScenic_grn: 进行基因调控网络分析"
	@echo -e " make -f loom=loom文件 outdir=输出目录 sample=样本名  Scenic_grn"
	@echo -e "\nScenic_ctx: 进行转录因子分析"
	@echo -e " make -f loom=loom文件 outdir=输出目录 sample=样本名  Scenic_ctx"
	@echo -e "\nScenic_aucell: 进行细胞分析"
	@echo -e " make -f loom=loom文件 outdir=输出目录 sample=样本名  Scenic_aucell"

get_expression:
	cho "############### get_expression start at `date` ###############"
	mkdir -p $(outdir)/$(sample)
	$(R) $(script)/prepare.r -i $(infile) -o $(outdir)/$(sample)/$(sample).csv  -s $(object)
	$(PYTHON3) $(script)/csv2loom.py -i $(outdir)/$(sample)/$(sample).csv  -o $(outdir)/$(sample)/$(sample).loom
	cho "############### get_expression end at `date` ###############"

loom=$(outdir)/$(sample)/$(sample).loom
Scenic_grn:
	mkdir -p $(outdir)/$(sample)
	echo "############### pyscenic grn start at `date` ###############"
	source /annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/miniconda/bin/activate  pyscenic && pyscenic grn --num_workers 20 --output $(outdir)/${sample}/$(sample).adj.sample.tsv --method grnboost2 $(loom) $(database)/hs_hgnc_tfs.txt 1>>$(outdir)/${sample}/log.o.txt 2>>$(outdir)/${sample}/log.e.txt
	echo "############### pyscenic grn end at `date` ###############"

loom=$(outdir)/$(sample)/$(sample).loom
Scenic_ctx:
	echo "############### pyscenic ctx start at `date` ###############"
	mkdir -p $(outdir)/$(sample)
	source /annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/miniconda/bin/activate  pyscenic && pyscenic ctx $(outdir)/${sample}/$(sample).adj.sample.tsv $(database)/hg19-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather  --annotations_fname $(database)/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname $(loom) --mode "dask_multiprocessing" --output $(outdir)/${sample}/$(sample).regulons.csv --num_workers 3 --mask_dropouts 1>>$(outdir)/${sample}/log.o.txt 2>>$(outdir)/${sample}/log.e.txt
	echo "############### pyscenic ctx end at `date` ###############"

loom=$(outdir)/$(sample)/$(sample).loom
Scenic_aucell:
	echo "############### pyscenic aucell start at `date` ###############"
	mkdir -p $(outdir)/$(sample)
	source /annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/miniconda/bin/activate  pyscenic && pyscenic aucell $(loom) $(outdir)/${sample}/$(sample).regulons.csv --output $(outdir)/$(sample)/$(sample).aucell.csv --num_workers 3 1>>$(outdir)/${sample}/log.o.txt 2>>$(outdir)/${sample}/log.e.txt
	echo "############### pyscenic aucell end at `date` ###############"

Pheatmap:
	echo "############### pyscenic Pheatmap start at `date` ###############"
	mkdir -p $(outdir)/$(sample)
	$(R) $(script)/get_cell_type.r -r $(rds) -o $(outdir)/$(sample)/$(sample).meta.csv
	$(PYTHON3) $(script)/group.py -m $(outdir)/$(sample)/$(sample).meta.csv -o $(outdir)/$(sample)/$(sample).meta.group.csv
	$(PYTHON3) $(script)/heatmap_pre.py -a $(outdir)/$(sample)/$(sample).aucell.csv -g $(outdir)/$(sample)/$(sample).meta.group.csv -o $(outdir)/$(sample)/$(sample).heatmap.xls
	$(R) $(script)/replicates_pheatmap.R -i $(outdir)/$(sample)/$(sample).heatmap.xls -o $(outdir)/$(sample)/$(sample).heatmap.pdf -f $(outdir)/$(sample)/$(sample).heatmap_PlotOrder.xls -a $(outdir)/$(sample)/$(sample).meta.group.csv -l "blue,white,red" -s column -r TRUE -c FALSE -T 60 -t 40 -C FALSE -R TRUE -W 0.3 -H 10 -F 4 -B FALSE -p 0 -G 0


