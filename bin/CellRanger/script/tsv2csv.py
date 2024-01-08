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

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-i','--input',help='input file',dest='input',type=open,required=True)
	parser.add_argument('-t','--title',help='header of output file',dest='title',default='Name')
	parser.add_argument('-o','--output',help='output file',dest='output',required=True)
	args=parser.parse_args()

	df = pd.read_table(args.input, index_col=None, header=None)
	df.columns = [args.title.split(',')]
	df.to_csv(args.output,index=False)
	#['0'.format(args.title)])

if __name__ == '__main__':
	main()
