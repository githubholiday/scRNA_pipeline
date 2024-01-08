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
统计比较组差异分析的基因数
'''
def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-d','--diff',help='diff file',dest='diff', required=True,nargs='+')
	parser.add_argument('-s','--sample',help='sample',dest='sample', required=True)
	parser.add_argument('-o','--output',help='output file',dest='output',required=True)
	args=parser.parse_args()

	info ={}
	for difffile in args.diff:
		if not os.path.exists(difffile): continue
		name = os.path.basename(difffile).split('_')[0]
		df = pd.read_csv(difffile,index_col = None,sep='\t', header = 0)
		down_df=df[(df['Significant'] =='yes') & (df['Up/Down']=='Down')]
		up_df=df[(df['Significant']=='yes') & (df['Up/Down']=='Up')]
		total=df[df['Significant']=='yes']
		if name in info.keys():
			info[name].update({'Up_gene':up_df.shape[0],'Down_gene':down_df.shape[0],'Total_gene':total.shape[0]})
		else:
			info.update({name:{'Up_gene':up_df.shape[0],'Down_gene':down_df.shape[0],'Total_gene':total.shape[0]}})
	infodf=pd.DataFrame(info)
	infodf=infodf.fillna(0)
	infodf=infodf.T
	infodf.index=infodf.index.astype(int)
	infodf=infodf.sort_index()
	print(infodf)
	
	infodf.to_csv(args.output, sep='\t',index_label='Cluster',columns=['Up_gene','Down_gene','Total_gene'])

if __name__ == '__main__':
	main()
