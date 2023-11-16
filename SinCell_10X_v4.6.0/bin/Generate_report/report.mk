tmpdir=$(dir $(abspath $(firstword $(MAKEFILE_LIST))))
Bin=$(tmpdir)/script
ifeq ($(strip $(config)),)
	Bconfig=$(tmpdir)/config/config.txt
else
	Bconfig=$(config)
endif
include $(Bconfig)

#Single_Report,Multi_Report
report_dir=$(outdir)/report
template_dir=$(tmpdir)/template
version=v$(shell echo $(pipeline) | awk -F'_v' '{print $$2}' | awk -F'/' '{print $$1}')


Help:
	@echo -e "功能:10X单细胞转录组报告生成模块"
	@echo -e "Target说明:"
	@echo -e "\t Report:适用于10X单细胞转录组的报告生成"
	@echo -e "\t Prepare:生成流程报告需要的report.conf文件"
	@echo -e "Prepare使用说明"
	@echo -e "make -f report.mk report_dir= seq= project_id= project_name= version= platform=  Prepare"
	@echo -e "\t[参数说明]"
	@echo -e "\treport_dir:结题报告路径，一般为outdir/report"
	@echo -e "\tseq:测序类型，形如PE150"
	@echo -e "\tproject_id:子项目编号"
	@echo -e "\tproject_name:项目名称，用于命名报告名称"
	@echo -e "\tversion:软件版本，用于报告中版本展示"
	@echo -e "\tplatform:测序平台，一般为MGI 或 Illumina"
	@echo -e "Report 使用说明"
	@echo -e "make -f report.mk indir= report_dir= template_file= upload_conf=  Report"
	@echo -e "\t[参数说明]"
	@echo -e "\treport_dir：结题报告生成路径，一般为outdir/report"
	@echo -e "\tindir：结果文件所在路径，一般为Analysis目录"
	@echo -e "\ttemplate_file：结题报告模板文件，可由template_dir/template_type.template组合成，template_dir为流程报告脚本路径，即该mk的路径，template_type为报告模板的前缀，一般为Single_10XRNA 或 Multi_10XRNA"
	@echo -e "\tupload_conf：整理upload目录的config文件，默认为template_dir/upload.conf,template_dir为该mk所在目录"


Clean_Upload:
	for file in `find $(upload_dir) -type l` ;\
	do \
		stat -L $$file > /dev/null 2> /dev/null ;\
		if [ $$? -gt 0 ] ;\
		then \
			rm $$file ;\
		fi; \
	done


.PHONY:Prepare
Prepare:
	[ -d $(report_dir) ] && echo $(report_dir) dir exist || mkdir -p $(report_dir)
	echo -e "SEQ:$(seq)\nPROJECT_ID:$(project_id)\nPROJECT_NAME:$(project_name)\nREPORT_DIR:$(report_dir)/upload\nVERSION:$(version)\nPLATFROM:$(platform)" > $(report_dir)/report.conf


template_file=$(template_dir)/$(template_type).template
upload_conf=$(template_dir)/$(template_type).upload.conf
no_tag=public-picture
lims_conf=$(tmpdir)/config/lims.ini
.PHONY:Report
Report:
	echo generate web report start at `date`
	$(PYTHON3) $(get_upload) -i $(indir) -o $(report_dir) -t $(template_file) -c $(upload_conf) -d $(no_tag) -b $(tmpdir) -ot $(report_dir)/report.template -n
	ssh 192.168.1.3 $(PYTHON3) $(Report) -i $(report_dir)/report.template -c $(report_dir)/report.conf -u admin -t cloud
	echo ssh 192.168.1.3 $(PYTHON3) $(Report) -i $(report_dir)/report.template -c $(report_dir)/report.conf -u admin -t cloud > $(report_dir)/report.sh
	$(PYTHON3) $(pipline_stat) -p $(project_id) -d $(report_dir) -c $(lims_conf)
	echo generate web report end at `date`


.PHONY:Convert
Convert:
	if [ -s $(CONVERT) ] ; \
	then \
		for i in `find $(indir) -name *.pdf` ; do name=`echo $$i |sed 's/.pdf/.png/'`; $(CONVERT) $$i $$name; done ; \
	else \
		echo "$(CONVERT) May not exist on the task running node " ; \
		exit 1 ; \
	fi