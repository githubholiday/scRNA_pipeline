
makefile_dir=$(dir $(firstword $(MAKEFILE_LIST)))
makefile_name=$(notdir $(firstword $(MAKEFILE_LIST)))
script=$(makefile_dir)/script

Rscript=/annoroad/data1/software/bin/miniconda/envs/Monocle/bin/Rscript
PYTHON3=/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/miniconda/envs/pyscenic/bin/python3

database=/annoroad/data1/bioinfo/PROJECT/RD/Cooperation/RD_Group/tuchengfang/Work/12_Project/ScRNA/01_20230901/Scenic/database/
csvtk=/annoroad/data1/software/bin/miniconda/envs/TGS_16s_pbtools/bin/csvtk


HELP:
	@echo Description:进行Scenic分析
	@echo Usage: make [target] [parameters]
	@echo "Target:"
	@echo -e "GSVA: 使用GSVA包进行GSVA分析"
	@echo -e " make -f infile=rda文件 outdir=输出目录 prefix=输出文件前缀 GSVA"
	@echo -e "\nget_cell_type: 获取rds文件中细胞和cluster以及细胞类型的对应关系"
	@echo -e " make -f rds=样本的rds文件 outdir=输出目录 prefix=输出文件前缀  get_cell_type"
	@echo -e "\nmerge_file: 将GSVA输出文件和细胞与cluster的对应关系表合并到一起，能确定转录因子和细胞的对应的关系，并得到cluster和细胞类型与转录因子的对应关系值（取cluster中每个细胞的富集值均值作为该cluster的富集值)"
	@echo -e " make -f cell_type_file=get_cell_type的输出文件 outdir=输出目录 prefix=输出文件前缀 merge_file"
	@echo -e "\ngsva_heatmap: 根据merge_file的输出文件，进行热图绘制，分别绘制cluster和细胞类型的图"
	@echo -e " make -f infile=输入文件 outdir=输出目录 prefix=输出文件前缀  gsva_heatmap"
	@echo -e "infile:第一行为转录因子，第一列为cluster或细胞类型，其余行列为对应的富集分数值"


GSVA:
	cho "############### GSVA start at `date` ###############"
    mkdir -p $(outdir)
    $(Rscript) $(script_dir)/gsva_3rna.r -r $(rda) -d $(database) -o $(outdir)/$(prefix).gsva.csv
	cho "############### GSVA end at `date` ###############"

get_cell_type:
	cho "############### get_cell_type start at `date` ###############"
	mkdir -p $(outdir)
	$(Rscript) $(script_dir)/get_cell_type.r -r $(rds) -o $(outdir)/$(prefix).cell_type.csv
	cho "############### get_cell_type end at `date` ###############"

merge_file:
	echo "############### merge_file start at `date` ###############"
	mkdir -p $(outdir)
	$(csvtk) -t join $(cell_type_file) $(gsva_out) > $(outdir)/$(prefix).gsva.merge.csv 
	$(PYTHON3) $(script_dir)/mean_gsva.py -i $(outdir)/$(prefix).gsva.merge.csv -o $(outdir) -p $(prefix)
	echo "############### merge_file end at `date` ###############"

gsva_heatmap:
	echo "############### gsva_heatmap start at `date` ###############"
	mkdir -p $(outdir)
	$(Rscript) $(script_dir)/gsva_pheatmap.r -i $(infile) -o $(outdir)/$(prefix).pdf
	$(CONVERT) -density 300 $(outdir)/$(prefix).pdf $(outdir)/$(prefix).png
	echo "############### gsva_heatmap end at `date` ###############"


