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
合并过滤后的细胞数,及筛选后的细胞
'''
def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-c','--cell',help='cell stat file',dest='cell', required=True,nargs='+')
	parser.add_argument('-o','--output',help='output file',dest='output',required=True)
	args=parser.parse_args()

	result=pd.DataFrame()
	for cellfile in args.cell:
		if not os.path.exists(cellfile): continue
		name = os.path.basename(cellfile).split('_')[0]
		df = pd.read_csv(cellfile,index_col = None, header = 0)
		result=result.append(df)
	infodf=result.T
	infodf=infodf.fillna(0)
	print(infodf)
	infodf.to_csv(args.output, sep='\t',header=0)

if __name__ == '__main__':
	main()
