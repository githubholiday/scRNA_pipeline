'''
根据给定的输入文件列表，生成投递任务的shell脚本
'''

#! /usr/bin/env python3
import argparse
import time
import sys
import re
import os
bindir = os.path.abspath(os.path.dirname(__file__))

__author__ = 'Liu Huiling'
__mail__ = 'huilingliu@genome.cn'
__doc__ = 'the description of program'
'''
The up and down gene count in second_level_go 
'''
pat1=re.compile('^s+$')

def info_from_infile( infile):
	'''
	infile:string,输入文件路径
	'''
	infile_path = os.path.abspath(infile)
	infile_dir_path = os.path.dirname(infile_path)
	infile_name = os.path.basename(infile_path)
	#cmp_name = infile_name.split("_gene_symbol.anno.xls")[0]
	cmp_name = infile_name.split('_diff_gene_symbol.xls')[0]
	return infile_dir_path, infile_name, cmp_name

def para_check(args, map_dic):
	check_status = 0 

	if len(args.input) == 0 :
		print("输入文件不存在，请检查，退出~")
		check_status += 1
	
	if args.action not in map_dic :
		print("-a must be go or kegg，please retry")
		check_status += 1
	if check_status > 0 :
		print("check your parameters and try again,bye bye~")
		sys.exit(1)

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-i','--input',help='input file',dest='input',nargs='+',required=True)
	parser.add_argument('-o','--outdir',help='outdir of the output',dest='output',required=True)
	parser.add_argument('-a','--action',help='go or kegg',dest='action',required=True)
	parser.add_argument('-s','--species',help='species',dest='species')
	parser.add_argument('-ca','--category',help='category name',dest='category')
	args=parser.parse_args()
	map_dic = {'go':"GO_shell/diff_go.sh",'kegg':"KEGG_shell/diff_kegg.sh"}
	para_check( args, map_dic)
	file_num =0
	with open( args.output, 'w') as outfile:
		for infile in args.input:
			file_num += 1
			infile_dir_path, infile_name, cmp_name = info_from_infile( infile )
			if args.action == 'go':
				go_cmd = 'make -f {BIN}/GO/go.mk de_report={infile} go_dir={infile_dir_path}/GO cmp={cmp_name} species={species} GO_Candidate GO_Enrich'.format(BIN=bindir, infile=infile, infile_dir_path=infile_dir_path, cmp_name=cmp_name, species=args.species)
				outfile.write(go_cmd+'\n')
			
			elif args.action == 'kegg' :
				kegg_cmd = 'make -f {BIN}/KEGG/kegg.mk de_report={infile} kegg_dir={infile_dir_path}/KEGG cmp={cmp_name} species={species} category={category} KEGG_Candidate KEGG_Enrich'.format(BIN=bindir, infile=infile, infile_dir_path=infile_dir_path, cmp_name=cmp_name, species=args.species, category=args.category)
				outfile.write(kegg_cmd+'\n')
			else:
				print("-a must be go or kegg")
	print("共处理 {0} 个文件".format(file_num))




if __name__ == '__main__':
	main()

