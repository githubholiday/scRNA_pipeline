# !/annoroad/data1/software/bin/miniconda/envs/Transcriptome_pipeline_2023_Python3.5.5/bin/python
"""
For SinCell_10X Flow-line Automatic configuration 
rewrite based on Transcriptome_flow-line_autoConfig.py
-in, --indir, 输入目录，需要包含Filter info Analysis-test子目录, 必须参数
-t, --type, 流程类型，配置流水线用，默认10XGenomics
-c, --category, 物种类型配置文件，跟据参考基因组获取物种类型，animal，plant等类型，默认/annoroad/data1/bioinfo/PMO/Public/database/RNA/database/species_category.ini
-ref, --refer_dir, cellranger参考基因组路径及前缀，默认/nas/data1/bioinfo/PROJECT/Commercial/Cooperation/Database/RNA/SinCell_10X/CellRanger_reference/refdata-cellranger-
-cf, --pipeline_config, 流程配置文件，包含job_config路径，默认脚本同级目录config_pipeline.ini
-py, --python3, python软件路径，用于投递脚本生成，默认/annoroad/data1/software/bin/miniconda/envs/python3_base/bin/python3
-pg, --pipeline_generate, 生成流程的脚本路径，默认/annoroad/data1/software/bin/pipeline_generate/bin/current/pipeline_generate.py
-pipd, --pipelineDir, 流程bin路径，默认为当前脚本所在目录
-fc， --forcecell, cellranger强制细胞数，默认为空
-pub, --publicDir, 10x单细胞转录组流程路径，用于识别模块，默认/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/Public/Pipeline/Stable/RNA/SinCell_10X/current/bin

"""
import argparse
import os
import sys
import re
import pandas as pd
import configparser
bin = os.path.abspath(os.path.dirname(__file__))
sys.path.append(bin + '/lib')
from PipMethod import myconf,generateShell,mkdir
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")


__author__='Yao Mengcheng'
__mail__= 'mengchengyao@genome.cn'
__modifier__= 'Yao Mengcheng'
__date__= '20191029'
'''
这个脚本之所以用之前流水线配置脚本的基础函数，主要是为了以后方便与流水线进行兼容，便有将之前的流水线差异分析升级到目前的差异分析类型。
by 姚盟成 20191029
'''
config_diff=configparser.ConfigParser(allow_no_value=True)

def raw_config(type):
	'''
差异分析的配置文件中的配置信息
	'''
	config_diff.add_section("sample")
	config_diff.add_section("cmp")
	config_diff.add_section("Para")
	#增加section 
	#config_diff.set("Para", "object_list_min.cells", '3')#增加option 
	if type=="10090":
		config_diff.set("Para", "object_list_mitoName", 'mt')
		config_diff.set("Para", "object_list_HB", 'mouse_hb')
	elif type=="9606":
		config_diff.set("Para", "object_list_mitoName", 'MT')
		config_diff.set("Para", "object_list_HB", 'human_hb')
	else:
		config_diff.set("Para", "object_list_mitoName", 'mt')
		config_diff.set("Para", "object_list_HB", '')
	#config_diff.set("Para", "object_list_mt.percent", '20')
	#config_diff.set("Para", "object_list_min_nFeature_RNA", '200')
	#config_diff.set("Para", "object_list_max_nFeature_RNA", '10000')
	config_diff.set("Para", "object_list_normalization.method", 'LogNormalize')
	config_diff.set("Para", "object_list_scale.factor", '10000')
	config_diff.set("Para", "object_list_nfeatures_FindVariableFeatures", '2000')
	config_diff.set("Para", "object_list_FindVariableFeatures_method", 'vst')
	config_diff.set("Para", "qc_pca_plot_w_h", '12,8')
	config_diff.set("Para", "sct", 'no')
	#config_diff.set("Para", "reduction_diff", 'yes')
	#config_diff.set("Para", "diff_reduction_resolution", '0.1,0.25,0.5,0.8,1,1.2')
	config_diff.set("Para", "integration_cca_dims", '20')
	config_diff.set("Para", "integration_pca_dims", '20')
	config_diff.set("Para", "integration_RunPCA_npcs", '30')
	config_diff.set("Para", "reduction_dims_num", '20')
	config_diff.set("Para", "reduction_resolution", '0.8')
	config_diff.set("Para", "reduction_w_h", '24,8')
	config_diff.set("Para", "marker_gene_min.pct", '0.1')
	config_diff.set("Para", "marker_gene_logfc.threshold", '0.25')
	config_diff.set("Para", "marker_gene_test.use", 'wilcox')
	config_diff.set("Para", "DE_gene_logfc.threshold", '0.25')
	config_diff.set("Para", "DE_gene_test.use", 'wilcox')
	config_diff.set("Para", "DE_gene_min.pct", '0.1')
	config_diff.set("Para", "seurat_clusters", 'seurat_clusters')
	config_diff.set("Para", "seurat_title", 'clusters')
	config_diff.set("Para", "cores", '6')
	config_diff.set("Para", "resolution", '1e-3')
	config_diff.set("Para", "maxcellnum", '20000')

class readinfo:
	Genome_version = None
	ppi_species = None
	pipeclass = None
	info_df = pd.DataFrame()
	de_cmp_df = pd.DataFrame()
	deg_cmp_df = pd.DataFrame()
	sample_order={}
	sample_group={}
	cmp_group={}
	def __init__(self, info_file, outdir,config,ppidf, ref_prefix,pipeclass):
		self.info_file = info_file
		self.outdir = outdir
		self.config = config
		self.ppidf = ppidf
		self.ref_prefix = ref_prefix
		self.pipeclass = pipeclass
		self.read_info()
		self.read_cmp()
		self.read_genome_version()
		
	def read_info(self):
		df = pd.read_excel(self.info_file, sheet_name=1)
		sample_is, sample_ie = 0, 0
		for i,row in df.iterrows():
			if row[0] == "样品名称": 
				sample_is = i  
				break

		df = pd.read_excel(self.info_file, sheet_name=0,header=sample_is+1)
		for i,row in df.iterrows():
			#if row["样品名称"] == "差异比较组合" or pd.isnull(row["样品名称"]) or pd.isnull(row["样品编号"]) or pd.isnull(row["物种拉丁名"]) or pd.isnull(row["样品描述"]) or pd.isnull(row["结题报告中样品名称"]) or pd.isnull(row["分组"]):
			if row["样品名称"] == "差异比较组合" or pd.isnull(row["样品名称"]) or pd.isnull(row["样品编号"]) or pd.isnull(row["物种拉丁名"]) or pd.isnull(row["结题报告中样品名称"]) or pd.isnull(row["分组"]):
				break
			sample_ie = i
		info_df = df[["样品名称","样品编号","物种拉丁名","样品描述","结题报告中样品名称","分组","组织部位","细胞类型","有何特殊处理"]].ix[0:sample_ie,:]

		self.info_df = info_df
		sample_group={}
		for i,row in info_df.iterrows():
			self.config.set("sample",'\t'.join(row.astype('str')))
			for one_group in row.astype('str')[5].split('/'):
				if one_group not in sample_group:
					sample_group[one_group]=[row.astype('str')[4]]
				else:
					sample_group[one_group].append(row.astype('str')[4])
			#self.config.set("cmp",row.astype('str')[4]) #将cmp中的单个样本给去掉了
		

		self.sample_group = sample_group
		self.config.set("Para","Para_num",str(info_df.shape[0]))

		sample_dup = info_df["样品名称"][info_df["样品名称"].duplicated(keep='first')]
		report_sample_dup = info_df["结题报告中样品名称"][info_df["结题报告中样品名称"].duplicated(keep='first')]
		if sample_dup.shape[0] != 0 or report_sample_dup.shape[0] != 0:
			if sample_dup.shape[0] != 0:
				sys.stderr.write ("以下 样品名称  有重复,请修正:\n")
				sys.stderr.write (" ".join(set(sample_dup)))
			if report_sample_dup.shape[0] != 0:
				sys.stderr.write ("\n以下 结题报告中样品名称  有重复,请修正:\n")
				sys.stderr.write (" ".join(set(report_sample_dup))+"\n")
			sys.exit(1)

		pipe_info = os.path.join(self.outdir,"info.txt")
		info_df.to_csv(pipe_info,sep='\t', header=None,index=None,encoding="utf-8")
		self.config.set("Para",'Para_info',pipe_info)
		sample_list_file = os.path.join(self.outdir,"sample.list")
		info_df[['结题报告中样品名称']].to_csv(sample_list_file,sep='\t', header=None,index=None,encoding="utf-8")
		self.config.set("Para",'Para_list',sample_list_file)

	def read_cmp(self):

		df = pd.read_excel(self.info_file, sheet_name=0)
		cmp_is, cmp_ie = 0, 0
		for i,row in df.iterrows():
			if row[0] == "合并比较组合": 
				cmp_is = i  
				break


		df = pd.read_excel(self.info_file, sheet_name=0,header=cmp_is+2)

		for i,row in df[["合并组合","组合名","比较组样品",'比较组']].iterrows():
			if pd.isnull(row["合并组合"]) or pd.isnull(row["组合名"]) or pd.isnull(row["比较组样品"]) or pd.isnull(row["比较组"]):
				break
			cmp_ie = i
		cmp_df = df[["合并组合","组合名","比较组样品","比较组"]].ix[0:cmp_ie,1:]
		cmp_df.dropna(how='any',inplace=True)

		cmp_group = {}
		sample_order = {}
		for i,row in cmp_df.iterrows():
			if row.isnull().any() : continue
			#print(i,row["比较组"])
			self.config.set("group",'\t'.join(row.astype('str'))+'\t'+row["比较组"].replace('/','_VS_'))
			self.config.set("cmp",row.astype('str')[0])
			if row.astype('str')[0] not in cmp_group:
				cmp_group[row.astype('str')[0]]=[row.astype('str')[2],]
			else:
				cmp_group[row.astype('str')[0]].append(row.astype('str')[2])
			if row.astype('str')[0] not in sample_order:
				sample_order[row.astype('str')[0]]=row.astype('str')[1]

			else:
				continue
		self.cmp_group = cmp_group
		self.sample_order =sample_order
		pipe_cmp = os.path.join(self.outdir,"combine.txt")
		if not len(cmp_df.index) == 0 :
			cmp_df.to_csv(pipe_cmp,sep='\t', header=None,index=None,encoding="utf-8")
			self.pipeclass = '1'
			sample_list_file = os.path.join(self.outdir,"sample.list")
			cmp_df[['组合名']].to_csv(sample_list_file,sep='\t', mode='a',header=None,index=None,encoding="utf-8")

		self.config.set("Para",'Para_cmbFile',pipe_cmp)

	def read_genome_version(self):

		df = pd.read_excel(self.info_file, sheet_name=0,index_col = 0)
		df = df.T
		self.Ref_version = df.ix[0,'CellRanger参考基因组']
		self.config.set("Para",'Para_ref',self.ref_prefix + self.Ref_version)
		print(df.ix[0,'CellRanger参考基因组'])
		self.Genome_version = df.ix[0,'注释基因组版本选择']
		self.config.set("Para",'Para_species',self.Genome_version)
		ppi_species = df.ix[0,'蛋白互作参考物种']
		ppi_df = self.ppidf[self.ppidf['STRING_name_compact'] == ppi_species]
		for i,row in ppi_df.iterrows():
			self.ppi_species = str(row['#taxon_id'])
		self.config.set("Para",'Para_ppi_species',self.ppi_species)

class read_filter_config_feedback:
	def __init__(self, filter_config_file, feedback_file,config,filtertype,seq_platform):
		self.sub_project_id = None
		self.finished = 'no'
		self.seq_platform= seq_platform
		self.filter_config_file = filter_config_file
		self.feedback_file = feedback_file
		self.filtertype = filtertype
		self.config = config
		self.filter_config()
		self.feedback_config()

	def filter_config(self):
		config_filter = myconf()
		config_filter.readfp(open(self.filter_config_file,encoding="utf-8"))
		projectName =''

		seq_type=''
		if self.filtertype=='k8s':
			self.sub_project_id = config_filter.get("Para","Para_project")
			projectName = config_filter.get("Para","Para_project_name")
			try:
				seq_type = config_filter.get("Para","Para_required_seq")
			except:
				seq_type = "PE150"
			if self.seq_platform :
				if not self.seq_platform == "BGI" and self.seq_platform == 'Illumina':
					print("-seq_platform 只接受  BGI  或者 Illumina  字符，请确认")
					sys.exit(0)
				else :
					print("已经通过外部参数seq_platform给定测序平台为 {0}".format(self.seq_platform))
			else :
				if 'Para_sequencer' in config_filter['Para']:
					self.seq_platform = config_filter.get("Para","Para_sequencer")
					if self.seq_platform == "MGI":
						self.seq_platform = "BGI"
				else :
					print("外部参数-seq_platform 以及 Filter_Result/shell/prepare/config.ini文件中无Para_sequencer，故测序平台默认为Illumina")
					self.seq_platform = 'Illumina'
		else:
			self.sub_project_id = config_filter.get("CONTROL","project")
			projectName = config_filter.get("CONTROL","project_name")
			try:
				seq_type = config_filter.get("CONTROL","need_seq")
			except:
				seq_type = "PE150"
			self.seq_platform = 'Illumina'

		#self.sub_project_id = config_filter.get("CONTROL","project")
		#projectName = config_filter.get("CONTROL","project_name")
		reobj = re.compile('任务单.*')
		projectName = projectName.lstrip(self.sub_project_id)

		projectName = re.sub(reobj,'',projectName)
		projectName = projectName.rstrip('测序')
		projectName = projectName.rstrip('建库')
		projectName = projectName.rstrip('过滤')
		projectName += '结题报告'
		#seq_type = config_filter.get("CONTROL","seq_type")
		SE = seq_type[:2]
		self.config.set("Para","Para_projectName",projectName)
		self.config.set("Para","Para_seq",seq_type)
		self.config.set("Para","Para_END",SE)
		self.config.set("Para","Para_platform",self.seq_platform)

	def feedback_config(self):
		config_feedback = myconf()
		if not os.path.isfile(self.feedback_file):
			return

		config_feedback.readfp(open(self.feedback_file,encoding="utf-8"))
		if config_feedback.has_section("DATA") and config_feedback.has_option("DATA","finish"):
			self.finished =  config_feedback['DATA']["finish"]
		else:
			return


def default_Para(config,indir):
	[config.add_section(i) for i in ['sample','group','cmp','Para']]
	config.set("Para","Para_clean",os.path.join(indir,"Filter","Filter_Result"))
	config.set("Para","Para_lib","clean")
	config.set("Para", "Para_mincells", '3')#增加option 
	config.set("Para", "Para_mt", '20')
	config.set("Para", "Para_hb", '5')
	config.set("Para", "Para_min_nFeature", '200')
	config.set("Para", "Para_max_nFeature", '10000')


class generate_pipeline_qsub:
	def __init__(self, python3, pipeline_generate,pipe_type,config_file,sub_project_id,outdir,pipeline_config_file,pipelineDir, pipeclass):
		self.work_shell = None
		self.pipeclass = pipeclass
		self.python3 = python3
		self.pipeline_generate = pipeline_generate
		self.pipe_type = pipe_type
		self.config_file = config_file
		self.sub_project_id = sub_project_id
		self.pipeline_config_file = pipeline_config_file
		self.pipelineDir = pipelineDir
		self.outdir = outdir
		self.generate_pipeline()
		self.pipetype = None

	def generate_pipeline(self):
		work_shell = os.path.join(self.outdir,'{0}_{1}_qsub_sge.sh'.format(self.pipe_type,self.sub_project_id))
		self.work_shell = work_shell
		self.pipeline_configi_file = os.path.abspath(self.pipeline_config_file)

		if str(self.pipeclass) == '1' and self.pipe_type == '10XGenomics' :
			pipetype = 'SinCell_10X'
			self.pipetype = pipetype
		else:
			self.pipetype = self.pipe_type

		mkdir(['{0}/Analysis'.format(self.outdir)])
		content = "{0} {1} -i {2} -o {3}/pipeline && \\\n".format(self.python3,self.pipeline_generate,self.pipeline_config_file,self.outdir)
		content += "sleep 30s && \\\n"
		content += "{0} {1}/pipeline/pipeline.py -i {2} -j {3}_{5} -b {4} -o {1}/Analysis -name {3}_{5} -r".format(self.python3,self.outdir,self.config_file,self.sub_project_id,os.path.dirname(self.pipelineDir),self.pipe_type)
		generateShell(work_shell,content)

def concession(indir):
	concession = False
	concession_conf = myconf()
	concession_file = os.path.join(indir,"concession.ini")
	if not os.path.isfile(concession_file):
		return concession

	concession_conf.readfp(open(concession_file,encoding="utf-8"))
	if concession_conf.has_section("sample") and concession_conf.has_option("sample","all_concession") and concession_conf['sample']["all_concession"]=="yes":
		return True
	else:
		return concession
def _mkdir(dir):
	if os.path.exists(dir):
		os.system("rm -rf {0}".format(dir))
		os.system("mkdir -p {0}".format(dir))
	else:
		os.system("mkdir -p {0}".format(dir))
	return 
def _mkdir2(dir):
	if os.path.exists(dir):
		pass
	else:
		os.system("mkdir -p {0}".format(dir))
	return 0

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter,
		epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\n'.format(__author__,__mail__,__date__))
	parser.add_argument('-in','--indir',help='which dir have Filter info Analysis concession',dest='indir',type=str,default=os.path.abspath(os.getcwd()))
	parser.add_argument('-t','--type',help='pipeline type, default is 10XGenomics',dest='type',type=str, default='10XGenomics')   
	parser.add_argument('-c','--category',help='species category',dest='category',type=str,default=os.path.join('/annoroad/data1/bioinfo/PMO/Public/database/RNA/database/species_category.ini'))
	parser.add_argument('-ref','--refer_dir',help='CellRanger reference dir',default='/nas/data1/bioinfo/PROJECT/Commercial/Cooperation/Database/RNA/SinCell_10X/CellRanger_reference/refdata-cellranger-')

	parser.add_argument('-cf','--pipeline_config',help='pipeline_config.ini',dest='pipeline_config',type=str,default=os.path.join(bin,'config_pipeline.ini'))
	parser.add_argument('-py','--python3',help='python3 path',dest='python3',type=str,default='/annoroad/data1/software/bin/miniconda/envs/python3_base/bin/python3')
	parser.add_argument('-pg','--pipeline_generate',help='path of pipeline_generate.py',dest='pipeline_generate',type=str,default='/annoroad/data1/software/bin/pipeline_generate/bin/current/pipeline_generate.py')
	parser.add_argument('-pipd','--pipelineDir',help='pipeline_bin',dest='pipelineDir',type=str,default=bin)
	parser.add_argument('-fc','--forcecell',help='forcecell',dest='forcecell',type=str)
	parser.add_argument('-cm','--clustermarker',help='clustermarker',dest='clustermarker',type=str)
	parser.add_argument('-seq_platform','--seq_platform',help='sequence platform(BGI or Illumina)',dest='seq_platform',required=False)
	args=parser.parse_args()
	config = myconf()
	config_category = myconf()
	config_category.readfp(open(args.category,encoding="utf-8"))
	"判断输出目录是否建立"
	if not os.path.isdir(args.indir):
		sys.stderr.write("Your dir is not make!\n{0}\n".format(args.indir))
		sys.exit(1)

	outdir = os.path.join(args.indir,"Analysis-238")
	config_file = os.path.join(outdir,"config.ini")
	Path(config_file).touch()
	Path(config_file).unlink()
	
	config_piptype = myconf()
	config_piptype.readfp(open(args.pipeline_config, encoding="utf-8"))
	ppi_species = config_piptype.get('config', 'ppi_species')
	ppidf = pd.read_table(ppi_species, header=0,index_col=None,encoding="utf-8")

	"设置默认参数"
	default_Para(config,args.indir)

	"设置从过滤目录获得的参数"
	feedback_file = os.path.join(args.indir,"Filter","Filter_Result","email","email_feedback.txt")
	filter_config_file = os.path.join(args.indir,"Filter","Filter_Result","config.ini")
	filter_k8s_config_file = os.path.join(args.indir,"Filter","Filter_Result/shell/prepare","config.ini")
	seq_platform=''
	if args.seq_platform:
		seq_platform=args.seq_platform
	if os.path.isfile(filter_k8s_config_file):
		filtertype='k8s'
		filter_feedback = read_filter_config_feedback(filter_k8s_config_file,feedback_file,config,filtertype,seq_platform)
	elif os.path.isfile(filter_config_file) :
		filtertype='normal'# filter_config_file, feedback_file,config
		filter_feedback = read_filter_config_feedback(filter_config_file,feedback_file,config,filtertype,seq_platform)
	else:
		sys.stderr.write("Filter文件夹中无config.ini,请核实数据\n")
		sys.exit(1)
	
	'''
	if not os.path.isfile(feedback_file) or not os.path.isfile(filter_config_file):
		sys.stderr.write("filter is not ok\n")
		sys.exit(1)
	'''

	"设置从subproject_info.xls获得的参数"
	subprojectID_info_file = os.path.join(args.indir,'info',"{0}_info.xls".format(filter_feedback.sub_project_id))
	subprojectID_info_file2 = os.path.join(args.indir,'info',"{0}_info.xlsx".format(filter_feedback.sub_project_id))

	if os.path.isfile(subprojectID_info_file):
		pass
	elif os.path.isfile(subprojectID_info_file2):
		subprojectID_info_file = subprojectID_info_file2
	else:
		sys.stderr.write("信息搜集表还没准备好,请及时上传\n")
		sys.exit(1)
	print(subprojectID_info_file)
	print(args.refer_dir)
	pipeclass=None
	print(subprojectID_info_file)
	info = readinfo(subprojectID_info_file,outdir,config,ppidf, args.refer_dir,pipeclass)
	if args.forcecell:
		config.set("Para","Para_forcecell",args.forcecell)
	else:
		config.set("Para","Para_forcecell",'')
	if args.clustermarker:
		config.set("Para","Para_cellmarker_anno",args.clustermarker)
	else:
		config.set("Para","Para_cellmarker_anno",'no')
	if info.pipeclass == '1'  and info.info_df.shape[0] > 1:
		config.set("Para","Para_ProjectType",'multiSample')
		config.set("Para","Para_project","{0}_{1}".format(filter_feedback.sub_project_id,'multiSample'))
	else:
		config.set("Para","Para_ProjectType",'singleSample')
		config.set("Para","Para_project","{0}_{1}".format(filter_feedback.sub_project_id,'singleSample'))
	config.set("Para","Para_project_id",filter_feedback.sub_project_id)
	config.set("Para","Para_category",config_category.get("category",info.Genome_version))
	config.write(open(config_file, "w+"))
	
	sample_group=info.sample_group
	cmp_group=info.cmp_group
	sample_order=info.sample_order
	
	raw_config(config.get("Para",'Para_ppi_species'))
	
	#print(sample_group, cmp_group, sample_order)
	if len(cmp_group) > 0:
		for one_cmp in cmp_group:
			tmp_config=config_diff
			group=[]
			i=1
			for ggroup in cmp_group[one_cmp]:
				tmp_config.set('cmp','cmp'+str(i),ggroup)
				i+=1
				for ggg in ggroup.split('/'):
					group.append(ggg)#list(sorted(set(cmp_group[one_cmp])))
			group=list(sorted(set(group)))
			sample1=[]
			sample2=[]
			for gg in group:
				sample1+=sample_group[gg]
				sample2+=[gg]*len(sample_group[gg])
			tmpe_sample_order_list=sample_order[one_cmp].strip(' ').strip("/").split('/')
			index=[sample1.index(i) for  i in tmpe_sample_order_list ]
			sample2_order=[ sample2[i] for i in index]
			for i in tmpe_sample_order_list:
				if tmpe_sample_order_list.count(i) > 1:
					print('{} 样本在比较组{}重复出现了{}次，请注意，核实修改信搜后重新投递'.format(i, one_cmp, tmpe_sample_order_list.count(i)))
					sys.exit(1)
			tmp_config.set('sample','sample1','/'.join(tmpe_sample_order_list))
			tmp_config.set('sample','sample2','/'.join(sample2_order))
			outdir_diff='{outdir}/Analysis/Integrating/{cmp}'.format( outdir=outdir, cmp=one_cmp)
			_mkdir2(outdir_diff)
			config_diff_file=outdir_diff+'/config.ini'
			tmp_config.write(open(config_diff_file,'w'))
	else:
		for g in sample_group:
			for s in sample_group[g]:
				tmp_config=config_diff
				sample1=s
				sample2=g
				tmp_config.set('sample','sample1',str(sample1))
				tmp_config.set('sample','sample2',str(sample2))
				outdir_diff='{outdir}/Analysis/Integrating/{cmp}'.format( outdir=outdir, cmp=sample1)
				_mkdir2(outdir_diff)
				config_diff_file=outdir_diff+'/config.ini'
				tmp_config.write(open(config_diff_file,'w'))




	print(info.pipeclass)
	if info.pipeclass == '1' and info.info_df.shape[0] > 1:
		config_pipeline_file = config_piptype.get('config','multiSample')
	else:
		config_pipeline_file = config_piptype.get('config','singleSample')
	config_pipeline_file = config_pipeline_file.replace('BIN',bin)
	pip_qsub = generate_pipeline_qsub(args.python3, args.pipeline_generate,args.type,config_file,filter_feedback.sub_project_id,outdir,config_pipeline_file,args.pipelineDir,info.pipeclass)

if __name__=="__main__": 
	main()
