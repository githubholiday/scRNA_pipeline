#! python3
"""
rewrite based on Transcriptome_flow-line_autoConfig.py
-in, --indir, 输入目录，需要包含Filter info Analysis，info子目录, 必须参数
-o,--outdir, 输出路径，如果不给定的，默认为indir
-c,--config, 流程配置文件，默认为 ../../software/software.txt
--ppi, 物种的taxid号，主要用于流程中去除线粒体等特殊处理
-r,--run 是否投递，给定该参数则投递，不给定则不投递。

"""
import argparse
import os
import sys
import configparser
import datetime
import glob
import json
import time
bindir = os.path.abspath(os.path.dirname(__file__))
filename=os.path.basename(__file__)

__author__='Tu chengfang'
__mail__= 'hh@qq.com'
__modifier__= 'Holiday T'
__date__= '20231221'

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
			time.sleep(0.5)

class Pipe_Info():
	def __init__( self, info_conf, analysis_dir, config_dic, filter_dir, ppi_species,pipe_config_file ):
		self.info_conf = info_conf
		self.analysis_dir = analysis_dir
		self.result_dir = '{0}/Analysis'.format( self.analysis_dir)
		self.config_dic = config_dic
		self.filter_dir = filter_dir
		self.ppi = ppi_species
		self.info_file = '{0}/prepare/info.txt'.format( self.analysis_dir )
		self.sample_list = '{0}/prepare/sample.list'.format( self.analysis_dir )
		self.config_file = '{0}/prepare/config.ini'.format( self.analysis_dir )
		self.cmp_file = '{0}/prepare/combine.txt'.format( self.analysis_dir )
		self.analysis_type = 'singleSample'
		self.pipe_config_file = pipe_config_file
		self.config = myconf()
		self.load_json()
		
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
		cellranger_ref_dir = '{0}{1}'.format(ref_dir, self.cellranger_ref)
		if os.path.exists( cellranger_ref_dir):
			self.config.set('Para','Para_ref',cellranger_ref_dir)
		else:
			my_log.error("{0} 路径下无 {1}".format( ref_dir, self.cellranger_ref ))
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
		sample信息为："样品名称","样品编号","结题报告中样品名称","分组","强制细胞数","样品描述","组织部位/细胞类型","有何特殊处理"
		
		获取所有样本的列表
        '''
		sample_info_dict = self.json_dic['samples']
		sample_list_file = open( self.sample_list, 'w')
		info_file = open( self.info_file, 'w')
		self.all_sample_list = []
		for sample in sample_info_dict :
			sample_info = sample_info_dict[sample][0]#[0]
			sample_report_name = sample_info[1]
			#将样本信息写入到sample.list中
			sample_list_file.write(sample+'\n')
			#如果强制细胞数写成0，则默认为空
			group = sample_info[2]
			force_cell = sample_info[5]
			if force_cell == '0':
				sample_info[5] = ''
			#将样本和组的对应关系存下来
			#sample信息中奖样本名称加上
			config_sample_info = [sample]+sample_info
			#写入到info.txt文件中
			config_sample_str = '\t'.join(config_sample_info)
			info_file.write( config_sample_str +'\n')
			#配置到config文件中农
			self.config.set('sample',config_sample_str)
			self.all_sample_list.append(sample_report_name)
		combine_value = "Combine\t"+'/'.join(self.all_sample_list)
		self.config.set('Combine',combine_value)
		#获取Combine目录下的config.ini文件
		self.combine_config(self.all_sample_list )
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
		'''
		功能：获取config文件中group的信息
		group_list:[A,B]
		cmp_str:A_VS_B
		'''
		#group_sample ：group:[A,A,A]
		group_info = self.json_dic['group_sample']
		cmp_sample_list = []
		cmp_group_list = []
		for group in group_list :
			sample_list = []
			if group in group_info:
				sample_list = group_info[group]
				for sample in sample_list :
					cmp_sample_list.append(sample)
					cmp_group_list.append(group)
				value = cmp_str+"\t"+"/".join(sample_list)
				self.config.set('group',value)
		#获取每个比较组目录下的config.ini文件
		self.cmp_config(cmp_str, cmp_sample_list, cmp_group_list, group_list)

	def config_para( self ):
		'''
        获取config.ini文件中的cmp(比较组)信息,可以根据流程需求随时添加
        '''
		pipe_config = os.path.abspath("{0}/../../config/config.txt".format(bindir))
		self.config.set('Para','Para_project_id',self.sub_project_id)
		self.config.set('Para','Para_project_name',self.project_name)
		self.config.set('Para','Para_outdir','{0}/Analysis'.format( self.analysis_dir))
		self.config.set('Para','Para_clean',self.filter_dir)
		self.config.set('Para','Para_list',self.sample_list)
		self.config.set('Para','Para_info',self.info_file)
		self.config.set('Para','Para_cmbFile',self.cmp_file)
		self.config.set('Para','Para_lib','clean')
		self.config.set('Para','Para_mt','20')
		self.config.set('Para','Para_hb','5')
		self.config.set('Para','Para_mincells','3')
		self.config.set('Para','Para_min_nFeature',"200")
		self.config.set('Para','Para_max_nFeature',"10000")
		self.config.set('Para','Para_all',"Combine")
		self.config.set('Para','Para_ProjectType',self.analysis_type)
		self.get_ref_dir()
		self.get_ref_type()
		self.config.set('Para','Para_ppi_species',self.ppi)
		self.config.set('Para','Para_config', self.pipe_config_file)
		
	
	def integrating_config(self, config):
		'''
		合并比较分析目录下的config文件,其公共部分
		'''
		if self.ppi == '9606':
			config.set('Para','object_list_mitoname','MT')
			config.set('Para','object_list_hb','human_hb')
		elif self.ppi == '10090':
			config.set('Para','object_list_mitoname','mt')
			config.set('Para','object_list_hb','mouse_hb')
		else:
			config.set('Para','object_list_mitoname','mt')
			config.set('Para','object_list_hb','')
		config.set('Para','object_list_normalization.method','LogNormalize')
		config.set('Para','object_list_scale.factor.method','10000')
		config.set('Para','object_list_nfeatures_findvariablefeatures','2000')
		config.set('Para','object_list_findvariablefeatures_method','vst')
		config.set('Para','qc_pca_plot_w_h','12,8')
		config.set('Para','sct','no')
		config.set('Para','integration_cca_dims','20')
		config.set('Para','integration_pca_dims','20')
		config.set('Para','integration_runpca_npcs','30')
		config.set('Para','reduction_dims_num','20')
		config.set('Para','reduction_resolution','0.6')
		config.set('Para','reduction_w_h','24,8')
		config.set('Para','marker_gene_min.pct','0.1')
		config.set('Para','marker_gene_logfc.threshold','0.25')
		config.set('Para','marker_gene_test.use','wilcox')
		config.set('Para','de_gene_logfc.threshold','wilcox')
		config.set('Para','de_gene_min.pct','0.1')
		config.set('Para','seurat_clusters','seurat_clusters')
		config.set('Para','seurat_title.pct','clusters')
		config.set('Para','cores','6')
		config.set('Para','resolution','1e-3')
		config.set('Para','maxcellnum','20000')
		
	def combine_config(self , all_sample_list ):
		'''
		获取Combine目录下的config配置文件
		'''
		combine_dir = '{0}/Integrating/Combine'.format( self.result_dir )
		my_mkdir([combine_dir])
		config_file = '{0}/config.ini'.format(combine_dir)
		combine_config = myconf()
		label_list = ['sample','Para']
		[combine_config.add_section(i) for i in label_list]
		combine_config.set('sample','sample1','/'.join(all_sample_list))
		self.integrating_config(combine_config)
		combine_config.write( open(config_file, 'w'))
	
	def cmp_config( self, cmp_str, sample_list, group_list, cmp_list ):
		'''
		cmp_str:A_VS_B
		sample_list:该比较组中所有样本的列表
		group_list:该比较组中所有样本对应的组名列表
		cmp_list:[A,B]
		'''
		cmp_dir = '{0}/Integrating/{1}'.format( self.result_dir,cmp_str )
		my_mkdir([cmp_dir])
		config_file = '{0}/config.ini'.format(cmp_dir)
		cmp_config = myconf()
		label_list = ['sample', 'cmp', 'Para']
		[cmp_config.add_section(i) for i in label_list]
		cmp_config.set('sample','sample1','/'.join(sample_list))
		cmp_config.set('sample','sample2','/'.join(group_list))
		cmp_config.set('cmp','cmp1','/'.join(cmp_list))
		self.integrating_config(cmp_config)
		cmp_config.write( open(config_file, 'w'))

	def pipe_config_write( self ):
		'''
		流程运行配置config文件生成
		'''
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
		self.job_config = '{0}/job_config/{1}.job_config.txt'.format(bindir, self.analysis_type)
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
	parser.add_argument('--ppi',help='ppi code 9606-human, 10090-mm',required=True)
	parser.add_argument('-r','--run',help='run or not',action='store_true')
	args=parser.parse_args()
	
	indir = args.indir
	if args.outdir :
		outdir = args.outdir
	else:
		outdir = args.indir

	
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

	my_pipe = Pipe_Info(info_json, analysis_dir, config_dic, filter_dir, args.ppi, args.config)
	my_pipe.pipe_config_write()
	my_pipe.generate_work_shell( args.run)

if __name__=="__main__": 
	my_log = Log(filename)
	main()
