DIR=$(dir $(abspath $(firstword $(MAKEFILE_LIST))))
Bin=$(DIR)/script
ifeq ($(strip $(config)),)
	Bconfig=$(DIR)/config/config.txt
else
	Bconfig=$(config)
endif
include $(Bconfig)

origin_refdir=$(outdir)/refdata-cellranger
ifdef $(species)
	ref_dir=$(origin_refdir)/refdata-cellranger_$(species)
endif

fq_dir=$(indir)/Analysis/$(fq_id)/filter/$(lib)
input_matrix_dir = $(indir)/$(sample_id)/
input_h5_matrix = $(indir)/$(sample_id)/*h5 


count_dir=$(outdir)/CellRanger_Count

outs=$(count_dir)/$(sample_id)/outs
## 过滤后的路径
filter_dir=$(outs)/filtered_gene_bc_matrices*

## 过滤后的h5文件
h5_matrix=$(outs)/filtered_gene_bc_matrices_h5.h5
## 过滤后的三件套路径
dataDir=$(filter_dir)/*/
## 原始的三件套路径
rawdataDir=$(outs)/raw_feature_bc_matrix/

barcodes=$(dataDir)/barcodes.tsv
genes=$(dataDir)/genes.tsv
matrix=$(dataDir)/matrix.mtx

down_dir=$(outs)/analysis

re_matrix=$(cluster_dir)/info/$(sample_id).filtered_gene_bc_matrices_h5.h5

cpu=8
mem=40
ifdef force_cell
	forcecell = --force-cells=$(force_cell)
endif
count_para = --id=$(sample_id) --fastqs=$(fq_dir) --sample=$(fq_id) --transcriptome=$(ref_dir) --localmem=$(mem) --localcores=$(cpu) $(forcecell) --include-introns true
aggr_list=$(count_dir)/$(sample_id)_list.csv
aggr_para = --id=$(sample_id) --csv=$(aggr_list) --normalize=mapped --localmem=$(mem) --localcores=$(cpu)
clst_para = --localmem=$(mem) --localcores=$(cpu) --id=$(sample_id)

Help:
	'''
	Count_Matrix 用CellRanger得到分析结果，表达量matrixs等
	'''


Count_Matrix:
	[ -d $(count_dir)/$(sample_id) ] && rm -rf $(count_dir)/$(sample_id) || echo fresh analysis
	mkdir -p $(count_dir)
	if [ "$(data_type)" = "fastq" ] &&  [ -d $(fq_dir) ] ;\
	then \
		cd $(count_dir) && $(CellRanger) count $(count_para) ;\
		[ -d $(outs)/filtered_feature_bc_matrix/ ] && mv $(outs)/filtered_feature_bc_matrix/ $(outs)/filtered_gene_bc_matrices && mkdir -p $(outs)/filtered_gene_bc_matrices/ref/ &&  mv $(outs)/filtered_gene_bc_matrices/*gz  $(outs)/filtered_gene_bc_matrices/ref/  ;\
		[ -f $(outs)/filtered_feature_bc_matrix.h5 ] && mv  $(outs)/filtered_feature_bc_matrix.h5 $(outs)/filtered_gene_bc_matrices.h5 ;\
	else \
		echo ${fq_dir} is not found ;\
	fi

	if [ "$(data_type)" = "matrix" ] &&  [ -d $(input_matrix_dir) ] ;\
	then \
		[ -d $(outs)/raw_feature_bc_matrix/ ] && rm -rf $(outs)/raw_feature_bc_matrix/ ;\
		mkdir -p $(outs) ;\
		ln -snf $(input_matrix_dir) $(outs)/raw_feature_bc_matrix ;\
	else \
		echo ${input_matrix_dir} is not found ;\
	fi

	if [ "$(data_type)" = "h5" ] &&  [ -f $(input_h5_matrix) ] ;\
	then \
		[ -d $(outs)/raw_feature_bc_matrix/ ] && rm -rf $(outs)/raw_feature_bc_matrix/ ;\
		mkdir -p $(outs)/raw_feature_bc_matrix/ ;\
		ln -snf $(input_h5_matrix) $(outs)/raw_feature_bc_matrix ;\
	else \
		echo ${input_h5_matrix} is not found ;\
	fi

	

# Stat target generates statistics and analysis for scRNA-seq data.
# It creates a directory for storing the statistics and analysis results.
# If the data type is "fastq", it runs the seurat_stat.R and seurat_analysis_statcell.R scripts on the dataDir and rawdataDir respectively.
# If the data type is "matrix", it runs the seurat_stat.R and seurat_analysis_statcell.R scripts on the input_matrix_dir and rawdataDir respectively.
# If the data type is "h5", it runs the seurat_stat.R and seurat_analysis_statcell.R scripts on the input_h5_matrix and rawdataDir respectively.
# It also converts the generated PDF files to PNG format.
Stat:
	mkdir -p $(count_dir)/$(sample_id)/stat
	if [ $(data_type) = "fastq" ] ;\
	then \
		$(RSCRIPT) $(Bin)/seurat_stat.R -d $(dataDir) -s $(sample_id) -o $(count_dir)/$(sample_id)/stat  ;\
		$(RSCRIPT) $(Bin)/seurat_analysis_statcell.R -d $(rawdataDir) --sample $(sample_id) -o $(count_dir)/$(sample_id)/stat --species $(species) ;\
		for i in `ls $(count_dir)/$(sample_id)/stat/*.pdf` ;do name=`echo $$i |sed 's/.pdf/.png/'`; $(CONVERT) $$i $$name; done ;\
	else \
		echo ${fq_dir} is not found ;\
	fi

	if [ $(data_type) = "matrix" ] ;\
	then \
		$(RSCRIPT) $(Bin)/seurat_stat.R -d $(input_matrix_dir) -s $(sample_id) -o $(count_dir)/$(sample_id)/stat ;\
		$(RSCRIPT) $(Bin)/seurat_analysis_statcell.R -d $(rawdataDir) --sample $(sample_id) -o $(count_dir)/$(sample_id)/stat --species $(species) ;\
		for i in `ls $(count_dir)/$(sample_id)/stat/*.pdf` ;do name=`echo $$i |sed 's/.pdf/.png/'`; $(CONVERT) $$i $$name; done ;\
	else \
		echo ${input_matrix_dir} is not found ;\
	fi

	if [ $(data_type) = "h5" ] ;\
	then \
		$(RSCRIPT) $(Bin)/seurat_stat.R -d $(input_h5_matrix) -s $(sample_id) -o $(count_dir)/$(sample_id)/stat -t h5 ;\
		$(RSCRIPT) $(Bin)/seurat_analysis_statcell.R -d $(rawdataDir) --sample $(sample_id) -o $(count_dir)/$(sample_id)/stat --species $(species) ;\
		for i in `ls $(count_dir)/$(sample_id)/stat/*.pdf` ;do name=`echo $$i |sed 's/.pdf/.png/'`; $(CONVERT) $$i $$name; done ;\
	else \
		echo ${input_h5_matrix} is not found ;\
	fi


Aggregation:	
	[ -d $(count_dir)/$(sample_id) ] && rm -rf $(count_dir)/$(sample_id) || echo fresh analysis
	mkdir -p $(count_dir)
	$(PYTHON3) $(tmpdir)/prepare_aggr.py -i $(count_dir) -s $(sample_id) -l $(group_list) -o $(count_dir)/$(sample_id)_list.csv 
	cd $(count_dir) && export PATH=$(CellRanger):$$PATH && $(CellRanger)/cellranger aggr $(aggr_para)
	[ -d $(outs)/filtered_feature_bc_matrix/ ] && mv $(outs)/filtered_feature_bc_matrix/ $(outs)/filtered_gene_bc_matrices && mkdir -p $(outs)/filtered_gene_bc_matrices/ref/ &&  mv $(outs)/filtered_gene_bc_matrices/*gz  $(dataDir)/
	[ -f $(dataDir)/barcodes.tsv.gz ] && gunzip $(dataDir)/barcodes.tsv.gz && gunzip $(dataDir)/features.tsv.gz && gunzip $(dataDir)/matrix.mtx.gz 
	[ -f $(dataDir)/features.tsv ] && mv $(dataDir)/features.tsv $(outs)/ && cut -f  1-2 $(outs)//features.tsv >$(outs)/filtered_gene_bc_matrices/ref//genes.tsv
	[ -f $(outs)/filtered_feature_bc_matrix.h5 ] && mv  $(outs)/filtered_feature_bc_matrix.h5 $(outs)/filtered_gene_bc_matrices_h5.h5

Aggr_plot:
	mkdir -p $(count_dir)/$(sample_id)/stat
	[ -f $(outs)/aggregation.csv ] && mv $(outs)/aggregation.csv $(outs)/aggregation_csv.csv
	$(PYTHON3) $(tmpdir)/get_cluster.py -i $(outs)/analysis/tsne/2_components/projection.csv -c $(outs)/analysis/clustering/graphclust/clusters.csv -o $(count_dir)/$(sample_id)/stat/$(sample_id)_graphcluster.xls -l $(outs)/aggregation_csv.csv
	$(RSCRIPT3) $(tmpdir)/scatter_ggplot2.r -i $(count_dir)/$(sample_id)/stat/$(sample_id)_graphcluster.xls -o $(count_dir)/$(sample_id)/stat/$(sample_id)
	for i in `ls $(count_dir)/$(sample_id)/stat/*.pdf` ;do name=`echo $$i |sed 's/.pdf/.png/'`; $(CONVERT) $$i $$name; done

Gene_Cluster:
	mkdir -p $(cluster_dir)/info
	ln -snf $(h5_matrix) $(re_matrix)
	$(PYTHON3) $(tmpdir)/tsv2csv.py -i $(barcodes) -o $(cluster_dir)/info/$(sample_id).barcodes.csv -t 'Barcode' 
	$(PYTHON3) $(tmpdir)/tsv2csv.py -i $(genes) -o $(cluster_dir)/info/$(sample_id).genes.csv -t 'Gene,Alias'
	cd $(cluster_dir) && export PATH=$(CellRanger):$$PATH && $(CellRanger)/cellranger reanalyze $(clst_para) --matrix $(re_matrix) --barcodes=$(cluster_dir)/info/$(sample_id).barcodes.csv --genes=$(cluster_dir)/info/$(sample_id).genes.csv

Create_Ref:
	echo Reference preparation for for 10x Genomics Cell Ranger started 
	[ -d $(origin_refdir) ] && rm -rf $(origin_refdir) || echo new origin_refdir
	mkdir -p $(origin_refdir)
	ln -snf $(ref_gtf) $(origin_refdir)
	ln -snf $(ref_fa) $(origin_refdir)
	#export LD_LIBRARY_PATH=$(BCL_LIB) && export PATH=$(BCL_PATH):$$PATH && 
	$(CellRanger)/cellranger mkgtf $(origin_refdir)/$(notdir $(ref_gtf)) $(origin_refdir)/new_$(notdir $(ref_gtf))
	cd $(origin_refdir) && export PATH=$(CellRanger):$$PATH && $(CellRanger)/cellranger mkref --nthreads=$(cpu) --genome=$(species) --fasta=$(origin_refdir)/$(notdir $(ref_fa)) --genes=$(origin_refdir)/new_$(notdir $(ref_gtf))
	mv $(origin_refdir)/$(species) $(origin_refdir)/refdata-cellranger-$(species)
	echo [Processed] Reference preparation for for 10x Genomics Cell Ranger Finished
