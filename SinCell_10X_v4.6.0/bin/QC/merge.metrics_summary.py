#! /usr/bin/env python3
import argparse
import time
import sys
import re
import os
import csv
import pandas as pd
bindir = os.path.abspath(os.path.dirname(__file__))

__author__ = 'Liu Huiling'
__mail__ = 'huilingliu@genome.cn'
__doc__ = 'the description of program'

pat1=re.compile('^s+$')
def read_csv(file):
	name = os.path.basename(file).split('.')[0]
	df = pd.read_csv(file, index_col = None, header = 0)
	df.index = [name]
	return df

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-i','--input',help='input file',dest='input', required=True,nargs='+')
	parser.add_argument('-o','--output',help='output file',dest='output',required=True)
	args=parser.parse_args()

	#result = pd.DataFrame.empty
	df_sum = pd.DataFrame()
	flag = 0
	for file in args.input:
		if not os.path.exists(file): continue
		name = os.path.basename(file).split('.|_')[0]
		df = read_csv(file)
		if flag == 0 :
			df_sum = df
			flag = 1
		else:
			df_sum = pd.concat([df_sum, df], axis=0)
	df_sum.T.to_csv(args.output, sep='\t',index_label='Sample')

if __name__ == '__main__':
	main()
