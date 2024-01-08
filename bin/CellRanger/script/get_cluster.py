#! /usr/bin/env python3
import argparse
import time
import sys
import re
import os
import pandas as pd
bindir = os.path.abspath(os.path.dirname(__file__))

__author__ = 'Liu Huiling'
__mail__ = 'huilingliu@genome.cn'
__doc__ = 'the description of program'

pat1=re.compile('^s+$')

def read_label(file):
	label = {}
	with open (file,'r') as f:
		for i,info in enumerate(f):
			if i ==0 : continue
			label[str(i)] = info.split(',')[0]
	f.close()
	return label
def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-i','--input',help='input file',dest='input',required=True)
	parser.add_argument('-c','--cluster',help='cluster file',dest='cluster',required=True)
	parser.add_argument('-l','--label',help='aggregation_csv.csv',dest='label')
	parser.add_argument('-o','--output',help='output file',dest='output',required=True)
	args=parser.parse_args()

	tsne = pd.read_csv(args.input, index_col=0, header=0)
	cluster = pd.read_csv(args.cluster, index_col=0, header=0)
	result = pd.concat([tsne, cluster], axis=1, join='outer')
	dict = {}
	for i in result['Cluster'].unique():
		dict[i] = 'Cluster' + str(i)
	result['Cluster'] = result['Cluster'].astype('object').replace(dict)
	result['group'] = result.index.str.get(-1)
	#print(result['Cluster'].head())
	if args.label:
		label = read_label(args.label)
		#print(label.keys())
		#print(label.values())
		result['group'] = result['group'].replace(label)
	
	#print(result.head())
	result.to_csv(args.output, sep='\t',header=True, index=True)

if __name__ == '__main__':
	main()
