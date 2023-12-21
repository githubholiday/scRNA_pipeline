# !/annoroad/data1/software/bin/miniconda/envs/Transcriptome_pipeline_2023_Python3.5.5/bin/python
"""
For SinCell_10X Flow-line Automatic configuration 
rewrite based on Transcriptome_flow-line_autoConfig.py
-in, --indir, 输入目录，需要包含Filter info Analysis-test子目录, 必须参数
-t, --type, 流程类型，配置流水线用，默认10XGenomics
-c, --category, 物种类型配置文件，跟据参考基因组获取物种类型，animal，plant等类型，默认/annoroad/data1/bioinfo/PMO/Public/database/RNA/database/species_category.ini
-ref, --refer_dir, cellranger参考基因组路径及前缀，默认/nas/data1/bioinfo/PROJECT/Commercial/Cooperation/Database/RNA/SinCell_10X/CellRanger_reference/refdata-cellranger-
-cf, --pipeline_config, 流程配置文件，包含job_config路径，默认脚本同级目录config_pipeline.ini
-py, --python3, python软件路径，用于投递脚本生成，默认 python3
-pg, --pipeline_generate, 生成流程的脚本路径，默认 /annoroad/data1/software/bin/pipeline_generate/bin/current/pipeline_generate.py
-pipd, --pipelineDir, 流程bin路径，默认为当前脚本所在目录
-fc， --forcecell, cellranger强制细胞数，默认为空
-pub, --publicDir, 10x单细胞转录组流程路径，用于识别模块，默认/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/Public/Pipeline/Stable/RNA/SinCell_10X/current/bin

"""
import argparse
import os
import sys
import configparser
import datetime
import glob
import json
bindir = os.path.abspath(os.path.dirname(__file__))
filename=os.path.basename(__file__)

__author__='Tu chengfang'
__mail__= 'hh@qq.com'
__modifier__= 'Holiday T'
__date__= '20191029'
'''
这个脚本之所以用之前流水线配置脚本的基础函数，主要是为了以后方便与流水线进行兼容，便有将之前的流水线差异分析升级到目前的差异分析类型。
by 姚盟成 20191029
import json
infile = open("A",'r')
json_dict = json.load(infile)
subProjectID = json_dict["sub_project_id"] #调用某个值

'''

class Log():
	def __init__( self, filename, funcname = '' ):
		self.filename = filename 
		self.funcname = funcname
	def format( self, level, message ) :
		date_now = datetime.datetime.now().strftime('%Y%m%d %H:%M:%S')
		formatter = ''
		if self.funcname == '' :
			formatter = '\n{0} - {1} - {2} - {3} \n'.format( date_now, self.filename, level, message )
		else :
			
			formatter = '\n{0} - {1} - {2} -  {3} - {4}\n'.format( date_now, self.filename, self.funcname, level, message )
		return formatter
	def info( self, message ):
		formatter = self.format( 'INFO', message )
		sys.stdout.write( formatter )
	def debug( self, message ) :
		formatter = self.format( 'DEBUG', message )
		sys.stdout.write( formatter )
	def warning( self, message ) :
		formatter = self.format( 'WARNING', message )
		sys.stdout.write( formatter )
	def error( self, message ) :
		formatter = self.format( 'ERROR', message )
		sys.stderr.write( formatter )
	def critical( self, message ) :
		formatter = self.format( 'CRITICAL', message )
		sys.stderr.write( formatter )
		
config_diff=configparser.ConfigParser(allow_no_value=True)
def my_run(cmd):
	if os.system( cmd) == 0 :
		my_log.info("执行成功：{0}".format( cmd ))
	else :
		my_log.info("执行失败，退出:{0}".format( cmd ))
		sys.exit(1)


def read_info_file(info_file, info_json, info_conf, table2json_script ):
	cmd = '{0} -c {1} -x {2} -j {3}'.format( table2json_script, info_conf, info_file, info_json)
	my_run( cmd )
	
def get_info_file( info_dir ):
	'''
    获取信息收集表文件，如果无或者超过1个，都会报错退出
    '''
	info_file = glob.glob( info_dir+'/*')
	if len(info_file) == 0 :
		my_log.error("{0} 目录下没有信息收集表文件，退出".format( info_dir))
		sys.exit(1)
	elif len(info_file) > 1:
		my_log.error("{0} 目录下有多个信息收集表文件，退出".format( info_dir))
		sys.exit(1)
	else:
		return info_file[0]

class myconf(configparser.ConfigParser):
	def __init__(self, defaults=None):
		configparser.ConfigParser.__init__(self, defaults=None, allow_no_value=True)

	def optionxform(self, optionstr):
		return optionstr

def read_config( config_file ):
	'''
	读取流程配置文件，获取软件等信息
	'''
	config_dic = {}
	with open(config_file, 'r') as infile:
		for line in infile:
			if line.startswith('#'):continue
			if line.strip() == '' : continue
			tmp = line.rstrip().split('=',1)
			key = tmp[0]
			value = tmp[1]
			if key in config_dic :
				my_log.error("{0} 在config文件中重复".format(key))
				sys.exit(1)	
			else :
				config_dic[key] = value
	return config_dic

def my_mkdir( dir_list ):
	for each_dir in dir_list :
		if not os.path.exists(each_dir) :
			os.makedirs( each_dir )

class Pipe_Info():
	def __init__( self, info_conf, analysis_dir, config_dic, filter_dir ):
		self.info_conf = info_conf
		self.analysis_dir = analysis_dir
		self.config_dic = config_dic
		self.filter_dir = filter_dir
		self.info_file = '{0}/prepare/info.txt'.format( self.analysis_dir )
		self.sample_list = '{0}/prepare/sample.list'.format( self.analysis_dir )
		self.config_file = '{0}/prepare/config.ini'.format( self.analysis_dir )
		self.cmp_file = '{0}/prepare/combine.txt'.format( self.analysis_dir )
		self.analysis_type = 'singleSample'
		self.config = myconf()
		self.load_json()
		self.get_group_file()
		self.sample_num = self.generate_pb_conf()
	
	def load_json( self ):
		json_file = open( self.info_conf, 'r')
		self.json_dic = json.load( json_file )
		self.sub_project_id = self.json_dic['sub_project_id']
		self.project_name = self.json_dic['project_name']
		self.cellranger_ref = self.json_dic['cellranger_ref']
		self.ref = self.json_dic['ref']
	#ok
	def default_para( self ):
		label_list = ['sample', 'group','cmp','Para','Combine']
		[self.config.add_section(i) for i in label_list]
		
	def get_ref_dir(self):
		'''
		获取cellranger的参考基因组路径
		'''
		ref_dir = self.config_dic['cellranger_ref_dir']
		cellranger_ref = '{0}/{1}'.format(ref_dir, self.cellranger_ref)
		if os.path.exists( cellranger_ref):
			self.config.set('Para','Para_ref',cellranger_ref)
		else:
			my_log.error("{0} 路径下无{1}".format( ref_dir, cellranger_ref))
			sys.exit(1)
	def get_ref_type(self):
		'''
		从species_category.ini文件中获取物种是动物、植物、真菌……
		'''
		ref_type = 'animal'
		with open(self.config_dic['species_category'], 'r') as infile:
			for line in infile:
				if line.startswith(self.ref):
					tmp = line.rstrip().split('=')
					ref_type = tmp[1].strip(' ')
		self.config.set('Para','Para_species',self.ref)
		self.config.set('Para','Para_category',ref_type)
		
			
		

	def config_sample( self ):
		'''
        获取config.ini文件中的sample（样本）信息，同时将样本信息输出到sample.list文件中，便于后期报告生成
        '''
		sample_info_dict = self.json_dic['samples']
		sample_list_file = open( self.sample_list, 'w')
		info_file = open( self.info_file, 'w')
		sample_list = []
		for sample in sample_info_dict :
			sample_info = sample_info_dict[sample][0]#[0]
			sample_list_file.write(sample+'\n')
			if sample_info[5] == '0':sample_info[5] = ''
			info_file.write( '\t'.join(sample_info)+'\n')
			config_sample_str = '\t'.join(sample_info)
			self.config.set('sample',config_sample_str)
			sample_list.append(sample)
		combine_value = "Combine\t"+'/'.join(sample_list)
		self.config.set('Combine',combine_value)
	###ok					
	def config_cmp( self) :
		'''
        获取config.ini文件中的cmp(比较组)信息
        '''
		cmp_info_list = self.json_dic['group_pair']
		if len(cmp_info_list) > 0 :
			self.analysis_type = 'multiSample'
		cmp_file = open(self.cmp_file,'w')
		for each_cmp in cmp_info_list:
			cmp_str = '_VS_'.join(each_cmp)
			cmp_file.write(cmp_str+'\n')
			self.config.set('cmp',cmp_str)
			self.config_group(each_cmp, cmp_str)
			
			
	def config_group( self, group_list, cmp_str ):
		#比较组中的样本,组名\t组中样本(sample1/sample2)
		group_info = self.json_dic['group_sample']
		for group in group_list :
			sample_list = []
			if group in group_info:
				sample = group_info[group]
				sample_list.append((i[0] for i in sample))
				value = group+"\t"+"\\".join(sample_list)
				self.config.set('group',value)
			

	def config_para( self ):
		'''
        获取config.ini文件中的cmp(比较组)信息,可以根据流程需求随时添加
        '''
		pipe_config = os.path.abspath("{0}/../../config/config.txt".format(bindir))
		self.config.set('Para','Para_project_id',self.sub_project_id)
		self.config.set('Para','Para_project_name',self.project_name)
		self.config.set('Para','Para_outdir','{0}/Analysis'.format( self.analysis_dir))
		self.config.set('Para','Para_list',self.sample_list)
		self.config.set('Para','Para_info',self.info_file)
		self.config.set('Para','Para_group_file',self.group_file)
		self.config.set('Para','Para_lib','clean')
		self.config.set('Para','Para_mt','20')
		self.config.set('Para','Para_hb','5')
		self.config.set('Para','Para_min_nFeature',"200")
		self.config.set('Para','Para_max_nFeature',"10000")
		self.config.set('Para','Para_all',"Combine")
		self.config.set('Para','Para_ProjectType',self.analysis_type)
		self.get_ref_dir()
		self.get_ref_type()
		
		
	def config_write( self ):
		self.default_para()
		self.config_sample()
		self.config_cmp()
		self.config_para()
		self.config.write( open(self.config_file, 'w'))
		my_log.info("config文件输出完成:{0}".format( self.config_file))

	def get_job_config(self):
		'''
		根据分析类型获取job_config文件
	    '''
		self.job_config = '{0}/job_config/{1}_job_config.txt'.format(bindir, self.analysis_type)
		if not os.path.exists( self.job_config ):
			my_log.error("文件不存在:{0},请确认分析类型是否正确：[16s,18s,ITs]，退出".format(  self.job_config  ))
			sys.exit(1)
			
	def generate_work_shell(self, run ):
		'''
        生成投递脚本
        '''
		self.get_job_config()
			
		python3 = self.config_dic['PYTHON3']
		generate_pipeline = self.config_dic['generate_pipeline']
		work_cmd = '{python3} {generate_pipeline} -i {self.job_config} -o {self.analysis_dir}/prepare/pipeline &&'.format(python3=python3, generate_pipeline=generate_pipeline,self=self )
		work_cmd += '\n{python3} {self.analysis_dir}/prepare/pipeline/pipeline.py -i {self.config_file} -j {self.sub_project_id}_10X -b {bindir}/../ -o {self.analysis_dir}/Analysis -name {self.sub_project_id}_10X -r\n'.format(python3=python3, self=self, bindir=bindir)
		
		work_shell = '{0}/{1}_qsub_sge.sh'.format( self.analysis_dir, self.sub_project_id)
		with open( work_shell, 'w') as outfile:
			outfile.write("#!/bin/bash\n")
			outfile.write(work_cmd)
			run_cmd = 'nohup /bin/bash {0} &'.format( work_shell )
			if run :
				my_run(run_cmd)
			else:
				my_log.info("请手动投递脚本:{0}".format( work_shell ))

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter,
		epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\n'.format(__author__,__mail__,__date__))
	parser.add_argument('-in','--indir',help='which dir have Filter info Analysis concession',dest='indir',type=str,required=True) 
	parser.add_argument('-o','--outdir',help='outdir of analysis',dest='outdir')
	parser.add_argument('-c','--config',help='config',dest='config',default="{0}/../../software/software.txt".format(bindir))
	args=parser.parse_args()
	
	indir = args.indir
	if args.outdir :
		outdir = args.outdir
	else:
		outdir = args.indir+"Analysis"

	
	analysis_dir = '{0}/Analysis/'.format( outdir )
	prepare_dir = '{0}/prepare'.format( analysis_dir )
	my_mkdir( [analysis_dir, prepare_dir])
    
    #获取下机数据路径
	filter_dir = '{0}/Filter/Filter_Result/'.format( indir )
	if not os.path.exists( filter_dir ):
		my_log.error("{0} 不存在，请准备好下机数据".format(filter_dir))
		sys.exit(1)
    # 读取配置文件
	config_dic = read_config(args.config)

	"设置从subproject_info.xls获得的参数"
	info_dir = '{0}/info'.format( args.indir )
	info_file = get_info_file(info_dir)
	info_json = '{0}/info.json'.format( prepare_dir)
	info_conf = '{0}/tab2json/config.json'.format( bindir)
	table2json_script = '{0}/tab2json/table2json'.format(bindir)
	read_info_file( info_file, info_json, info_conf, table2json_script )
	
    #
	
	
	my_pipe = Pipe_Info(info_json, analysis_dir, config_dic, filter_dir)
	my_pipe.config_write()
	my_pipe.generate_work_shell( args.run)

if __name__=="__main__": 
	my_log = Log(filename)
	main()
