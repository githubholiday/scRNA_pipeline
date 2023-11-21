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
count_dir=$(outdir)/CellRanger_Count
outs=$(count_dir)/$(sample_id)/outs
filter_dir=$(count_dir)/$(sample_id)/outs/filtered_gene_bc_matrices*
h5_matrix=$(outs)/filtered_gene_bc_matrices_h5.h5
dataDir=$(filter_dir)/*/
rawdataDir=$(count_dir)/$(sample_id)/outs/raw_feature_bc_matrix/
barcodes=$(dataDir)/barcodes.tsv
genes=$(dataDir)/genes.tsv
matrix=$(dataDir)/matrix.mtx
down_dir=$(count_dir)/$(sample_id)/outs/analysis

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


rawdataDir=$(indir)/$(sample_id)/outs/raw_feature_bc_matrix/
filterDir=$(indir)/$(sample_id)/outs/filtered_gene_bc_matrices/ref/

Help:
	@echo -e "Count_Matrix 用CellRanger得到分析结果，表达量matrixs等"
	@echo -e "make -f makefile count_dir= sample_id= fq_dir= count_para= outs= fq_id= ref_dir=  Count_Matrix"
	@echo -e "\t Parameters:\n"
	@echo -e "\t count_dir: CellRanger的输出目录"
	@echo -e "\t sample_id: 最终的样本名称"
	@echo -e "\t fq_id: fastq文件的ID"
	@echo -e "\t fq_dir: fastq文件目录,目录下为该样本的fastq.gz文件"
	@echo -e "\t count_para: CellRanger count的参数"
	@echo -e "\t outs: CellRanger count的输出目录,一般为count_dir/sample_id/outs"
	@echo -e "\t ref_dir:参考基因组路径"

	@echo -e "Stat 用seurat分析结果，得到分析结果的统计信息"
	@echo -e "make -f makefile filterDir= rawDir= sample_id= outdir= species= Stat"
	@echo -e "\t Parameters:\n"
	@echo -e "\t filterDir: CellRanger count的输出目录,一般为 indir/sample_id/outs/"filtered_gene_bc_matrices"
	@echo -e "\t rawDir: CellRanger count的输出目录,一般为 indir/sample_id/outs/raw_feature_bc_matrix"
	@echo -e "\t species: 物种名称"
	@echo -e "\t sample_id: 最终的样本名称"
	@echo -e "\t outdir: 输出目录"


Count_Matrix:
	[ -d $(count_dir)/$(sample_id) ] && rm -rf $(count_dir)/$(sample_id) || echo fresh analysis
	mkdir -p $(count_dir)
	if [ -d $(fq_dir) ] ;\
	then \
		cd $(count_dir) && $(CellRanger) count $(count_para) ;\
	else \
		echo ${fq_dir} is not found ;\
	fi
	[ -d $(count_dir)/$(sample_id)/outs/filtered_feature_bc_matrix/ ] && mv $(count_dir)/$(sample_id)/outs/filtered_feature_bc_matrix/ $(count_dir)/$(sample_id)/outs/filtered_gene_bc_matrices && mkdir -p $(count_dir)/$(sample_id)/outs/filtered_gene_bc_matrices/ref/ &&  mv $(count_dir)/$(sample_id)/outs/filtered_gene_bc_matrices/*gz  $(count_dir)/$(sample_id)/outs/filtered_gene_bc_matrices/ref/ 
	[ -f $(outs)/filtered_feature_bc_matrix.h5 ] && mv  $(outs)/filtered_feature_bc_matrix.h5 $(outs)/filtered_gene_bc_matrices.h5


Stat:
	mkdir -p $(outdir)
	$(RSCRIPT) $(Bin)/seurat_stat.R -d $(filterDir) -s $(sample_id) -o$(outdir)
	$(RSCRIPT) $(Bin)/seurat_analysis_statcell.R -d $(rawDir) --sample $(sample_id) -o $(outdir) --species $(species)
	for i in `ls $(outdir)/*.pdf` ;do name=`echo $$i |sed 's/.pdf/.png/'`; $(CONVERT) $$i $$name; done
