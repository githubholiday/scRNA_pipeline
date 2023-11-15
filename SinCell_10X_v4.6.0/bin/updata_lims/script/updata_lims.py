#! -*- coding : utf-8 -*-
'''
1. 更新lims数据库
'''

import os, sys
import configparser
import argparse
import glob
import time
import json
bindir = os.path.abspath(os.path.dirname(__file__))
sys.path.append('{0}/lib'.format(bindir))
from readconfig import Config_Parser
import Lims_SQL

__author__ = 'leiguo'
__mail__ = 'leiguo@genome.cn'

def my_log( level, message ) :
    logging.basicConfig(level = logging.INFO,format = '%(asctime)s - %(filename)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)
    if level == 'info' :
        return logger.info( message )
    if level == 'warning' :
        return logger.warning( message )
    if level == 'debug' :
        return logger.debug( message )
    if level == 'error' :
        return logger.error( message )

def check_file_exists( *file_list ) :
    for file in file_list :
        if os.path.exists( file ) :
            my_log( 'info', 'file : {0}'.format( file ) )
        else :
            my_log( 'error', 'file is not exists : {0}'.format( file ) )

def make_dir( dir ) :
    if not os.path.exists( dir ) :
        try :
            os.makedirs( dir )
            time.sleep(1)
            my_log( 'info', 'mkdir {0} sucessful!'.format( dir) )
        except :
            my_log( 'error', 'mkdir {0} failed!'.format( dir) )
    else :
        my_log( 'info', '{0} is exist'.format( dir ) )

def myrun( cmd ) :
    if os.system( cmd ) == 0 :
        my_log( 'info', '{0} run sucessfully !'.format( cmd ) )
    else :
        my_log( 'error', '{0} run failed !'.format( cmd ) )

class Update_database():
    def __init__(self, config, sci_bioinfo, value):
        self.Lims = Lims_SQL.LIMS(config)
        self.sci = Lims_SQL.LIMS(sci_bioinfo)
        self.value = value
        self.project = value['project_code']
        self.sample = value['sample']

    #update
    def update_db(self):
        select_condition = [('project_code', self.project),('sample', self.sample)]
        species_infos = self.Lims.select('tb_info_sampleinfo', ['species'], [('project_code', self.project),('sample_name', self.sample)])
        species = ''
        if len(species_infos): 
            species = species_infos[0][0]
        self.value['species'] = species
        infos = self.sci.select('10xgenomics', ['id'], select_condition)
        now_time = int(time.time()*1000)

        value_list = []
        for k in self.value:
            tmp = (k, self.value[k])
            value_list.append(tmp)

        if len(infos) > 0:
            value_list.append(('update_time', now_time))
            self.sci.update('10xgenomics', value_list, select_condition)
        else:
            value_list.append(('creat_time', now_time))
            value_list.append(('update_time', now_time))
            col_list = [i[0] for i in value_list]
            insert_value_list = [i[1] for i in value_list]
            print(col_list, insert_value_list)
            self.sci.insert(  table_name='10xgenomics', col_list=col_list, value_list=insert_value_list )

        print("插入数据库完成")

    def close_db(self):
        self.Lims.close()
        self.sci.close()

def get_stat(file, sample):
    value = {}
    with open(file, 'r') as fi:
        for index, line in enumerate(fi):
            info = line.strip().split('\t')
            if index == 0:
                sample_index = info.index(sample)
            else:
                value[info[0]] = info[sample_index]
    return value


def trans_name(value, index_name):
    new_value = {}
    for k in value:
        if k in index_name:
            new_value[index_name[k]] = value[k]
    return new_value


def get_version(path):
    if not os.path.exists(path):
        print('{} 不存在'.format(path))
        sys.exit()
    cmd = '{} --version'.format(path)
    r = os.popen(cmd)
    info = r.read()
    version = info.strip().split(' ')[-1]
    print('软件版本为{}'.format(version))
    return version

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog='author:\t{0}\nmail:\t{1}\n'.format(__author__, __mail__))
    parser.add_argument('-c', '--config', help='project config file', dest='config', required=True)
    parser.add_argument('-s', '--sample', help='sample name', dest='sample', required=True)
    parser.add_argument('-t', '--type', help='10X or C4', dest='type', required=True)
    parser.add_argument('-a', '--stat', help='cellranger/C4 stat file', dest='stat', required=True)
    parser.add_argument('-f', '--filter', help='filter cell file', dest='filter', required=True)
    parser.add_argument('-d', '--doublet', help='doublet cell file', dest='doublet', required=True)
    parser.add_argument('-w', '--software', help='cellranger/C4 software path', dest='software', required=True)
    parser.add_argument('-j', '--json_config', help='json config of index name', dest='json_config', default='{}/json/index.json'.format(bindir))
    args = parser.parse_args()
    
    config = Config_Parser(args.config)
    sample_list = config.return_title_block('sample')
    for i in sample_list:
        tmp = i.split('\t')
        if tmp[4] == args.sample:
            info = tmp
            break
    else:
        print('config.ini中无{}样本信息，请确认'.format(args.sample))
        sys.exit(1)
    
    lims_config = '{}/lims.ini'.format(bindir)
    sci_bioinfo = '{}/sci_bioinfo.ini'.format(bindir)
    index_config = args.json_config
    with open(index_config, 'r') as fi:
        index = json.load(fi)
        index_name = index[args.type]

    
    
    value = {} #插入数据库信息 K -- V
    
    value['software_version'] = get_version(args.software)
    value['project_code'] = config.all_block('Para', 'Para_project_id')
    value['analysis_ref'] = config.all_block('Para', 'Para_ref')
    value['type'] = args.type
    value['sample'] = info[0]
    value['sample_description'] = info[3]
    value['sample_name'] = info[4]
    value['tissue_from_info'] = info[6]
    value['celltype_from_info'] = info[7]
    value['treatment_from_info'] = info[8]
    
    value.update(get_stat(args.stat, args.sample))
    value.update(get_stat(args.filter, args.sample))
    value.update(get_stat(args.doublet, args.sample))

    new_value = trans_name(value, index_name)
    #######
    Update_data = Update_database(lims_config, sci_bioinfo, new_value)
    Update_data.update_db()
    close_database = Update_data.close_db()
    
if __name__ == "__main__":
    main()
