tmpdir=$(dir  $(abspath $(firstword $(MAKEFILE_LIST))))/
Rscript=/annoroad/data1/bioinfo/PMO/yaomengcheng/Anaconda3/envs/cellassign/bin/Rscript
python3=/annoroad/share/software/install/Python-3.3.2/bin/python3
RSCRIPT_monocle3=/annoroad/data1/bioinfo/PMO/yaomengcheng/Anaconda3/bin/Rscript
CONVERT=/usr/bin/convert
result_dir=$(outdir)

ifeq ($(multi),no)
	sample_type = single 
else
	sample_type = multi
endif

HELP:
		@echo 程序功能：
		@echo 单细胞转录组进行亚群自动化鉴定，通过输入【seurat的rds】文件进行细胞亚群自动化鉴定。
		@echo
		@echo -e 'marker：           是否进行marker基因分析，找到差异最大的marker基因，默认不分析'
		@echo -e 'rds：              result for analysis by seurat3.0 ,目前仅支持seurat3.0以后版本的rds结果文件，rds的meta矩阵中必须有orig.ident列,cluster列选用active.ident[forced]'
		@echo -e 'db：               参数数据集，一共是7个，其中有5个是人的参考数据集HumanPrimaryCellAtlasData,BlueprintEncodeData,DatabaseImmuneCellExpressionData,NovershternHematopoieticData,MonacoImmuneData, 小鼠的两个ImmGenData,MouseRNAseqData。'
		@echo -e 'species：          物种，目前为人和小鼠，可以选择物种，然后会有默认的数据集。'
		@echo -e 'outdir：           输出目录，结果输出目录。'
		@echo -e 'prefix：           结果输出文件的前缀名称。'
		@echo -e 'multi：            是否单样本或者多样本。'

		@echo
		@echo 使用方法1：
		@echo make -f makefile rds= db= species=  outdir= prefix= multi=no singler
		@echo
		@echo 程序更新：
		@echo v1.0.6 2020-09-07 by 'Yao Mengcheng' [mengchengyao\@genome.cn];
		@echo 2022-03-17 by 'zhaoqian' [qianzhao\@genome.cn];
singler:
		echo `date "+%Y%m%d %H:%M:%S"` "- singleR - INFO - 开始获取输入文件"
		[ -d $(result_dir) ] && echo "singler outdir is exist" || mkdir -p $(result_dir)
		$(Rscript) $(tmpdir)/singleR_v1.0.6.r -r $(rds) -d $(db) -s $(species) -o $(result_dir) --prefix $(prefix)
		$(python3) $(tmpdir)/stat_singler_cluster.py -c $(result_dir)/result/$(prefix)_singleR_celltype.xls -o $(result_dir)/result/
		if [ $(sample_type) == "single" ] ;then \
			$(RSCRIPT_monocle3) $(tmpdir)/Rename_singsample_cluster.R -r $(rds) -a $(result_dir)/result/singleR_celltype_stat_top1.xls  -o $(result_dir)/result/ --prefix $(prefix) ;\
		else \
			$(RSCRIPT_monocle3) $(tmpdir)/Rename_cluster.R -r $(rds) -a $(result_dir)/result/singleR_celltype_stat_top1.xls -o $(result_dir)/result/ --prefix $(prefix)  ;\
		fi
		for i in `ls $(result_dir)/result/*.pdf` ; do name=`echo $$i |sed 's/.pdf/.png/'`; $(CONVERT) $$i $$name; done
		cp $(tmpdir)/singleR结果说明.pdf $(result_dir)/
		echo `date "+%Y%m%d %H:%M:%S"` "- Finished singleR - INFO - 成功获取输入文件"

test:
	ls $(indir)/*_all_UMI.csv |while read i ;\
		do \
		sample_name=`basename $$i _all_UMI.csv` ;\
		echo $$sample_name ;\
	done


