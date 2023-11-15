BIN=$(dir  $(abspath $(firstword $(MAKEFILE_LIST))))
include $(BIN)../../software/software.txt

CONFIG=$(INDIR)/../config.ini
OUTDIR="/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/FTP/golden_data/RD_All_Pipeline_Test_Example/安诺人类单细胞数据库v20210826/2_Anno_Base_Analysis/"

HELP:
	@echo Description: 此脚本用于自动收集10X单细胞转录组中的人类样品基本分析结果到指定路径
	@echo Usage :
	@echo make -f makefile CONFIG= INDIR= OUTDIR= Analysis_Project_Collection
	@echo
	@echo [参数说明]    :
	@echo CONFIG        : [必选] 10X 转录组标准分析生成的流程参数配置文件（config.ini）
	@echo INDIR         : [必选] 10X 转录组标准分析存放shell、CellIdent等结果文件夹的目录
	@echo OUTDIR        : [可选] 人单细胞转录组数据库存放路径，makefile上面已经指定
	@echo
	@echo [Version]:
	@echo v1.0.0 2021-09-27 by 'Simon Lee' [huayunli\@genome.cn];



Analysis_Project_Collection:
	echo "Analysis project collection start at `date`"
	$(PYTHON3) $(BIN)/extract_config.py -c $(CONFIG) -i $(INDIR) -o $(OUTDIR)
	echo "Analysis project collection finished at `date`"
