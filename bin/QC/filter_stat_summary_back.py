#! /usr/bin/env python3
import argparse
import time
import sys
import re
import os
import csv
import pandas as pd
bindir = os.path.abspath(os.path.dirname(__file__))

__author__ = ''
__mail__ = ''
__doc__ = 'the description of program'

pat1=re.compile('^s+$')
'''
统计过滤后的细胞数，基因中位数
'''
def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-c','--cell',help='cell file',dest='cell', required=True,nargs='+')
	parser.add_argument('-o','--output',help='output file',dest='output',required=True)
	args=parser.parse_args()

	info ={}
	for cellfile in args.cell:
		if not os.path.exists(cellfile): continue
		name = os.path.basename(cellfile).split('_')[0]
		df = pd.read_csv(cellfile,index_col = None, header = 0)
		print(df.shape[0])
		if name in info.keys():
			info[name].update({'Number of Cells':df.shape[0],'Median Genes per Cell':int(df['nFeature_RNA'].median())})
		else:
			info.update({name:{'Number of Cells':df.shape[0],'Median Genes per Cell':int(df['nFeature_RNA'].median())}})
	infodf=pd.DataFrame(info)
	infodf=infodf.fillna(0)
	infodf.T.to_csv(args.output, sep='\t',index_label='Sample',columns=['Number of Cells','Median Genes per Cell'])

if __name__ == '__main__':
	main()
