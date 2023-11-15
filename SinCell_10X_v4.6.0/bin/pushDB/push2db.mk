BIN_DIR=$(dir $(abspath $(firstword $(MAKEFILE_LIST))))/
include $(BIN_DIR)../../software/software.txt

tmpfile=$(input).pushdb.tmp
project_id=$(shell echo $(project) | awk -F'_' '{print $$1}')
TSV2SQL=$(BIN_DIR)/mariadb_BJ.py
HELP:
	@echo Description:
	@echo -e "\t" 将流程监控指标插入数据库
	@echo
	@echo Usage：
	@echo -e "\t" make -f makefile config= input= project= table_name= target= INSERT
	@echo -e "\t" make -f makefile config= input= project= table_name= target= TRANSPOSE_INSERT
	@echo -e "\t" make -f makefile config= input= project= table_name= target= UPDATE
	@echo -e "\t" make -f makefile config= input= project= table_name= target= TRANSPOSE_UPDATE
	@echo
	@echo Parameters:
	@echo -e "\t" config: config配置文件，包含数据库格式及相关插入信息
	@echo -e "\t" input: 输入文件
	@echo -e "\t" project: 子项目编号
	@echo -e "\t" table_name: 要操作的表名，一个流程对应一个固定的表名
	@echo -e "\t" target: config配置文件的target
	@echo
	@echo 程序更新:
	@echo v1.0.0 2019-08-28 by 'zhang yue' [yuezhang\@genome.cn];

INSERT:
	$(PYTHON3) $(TSV2SQL) -c $(config) -i $(tmpfile) -p $(project) -t $(table_name) -k $(target) -e

TRANSPOSE_INSERT:
	$(TRANSPOSE) -t $(input) >$(tmpfile)
	$(PYTHON3) $(TSV2SQL) -c $(config) -i $(tmpfile) -p $(project_id) -t $(table_name) -k $(target) -e
	rm $(tmpfile)

UPDATE:
	$(PYTHON3) $(TSV2SQL) -c $(config) -i $(input) -p $(project) -t $(table_name) -k $(target) -u

TRANSPOSE_UPDATE:
	$(TRANSPOSE) -t $(input) >$(tmpfile)
	$(PYTHON3) $(TSV2SQL) -c $(config) -i $(tmpfile) -p $(project) -t $(table_name) -k $(target) -u
	rm $(tmpfile)

GENOME:
	$(ICONV) -f GBK -t UTF-8 -o $(tmpfile) $(input)
	$(PYTHON3) $(BIN_DIR)/genome_size.py -i $(tmpfile) -g $(genome) -o $(tmpfile)
	$(PYTHON3) $(TSV2SQL) -c $(config) -i $(tmpfile) -p $(project) -t $(table_name) -k $(target) -u
	rm $(tmpfile)
