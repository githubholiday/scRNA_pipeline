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

pat1=re.compile('^s+$')

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-i','--input',help='input dir',dest='input',required=True)
	parser.add_argument('-o','--output',help='output file',dest='output',type=argparse.FileType('w'),required=True)
	parser.add_argument('-s','--sample',help='combine sample name',dest='sample',required=True)
	parser.add_argument('-l','--list',help='lists seperated by /',dest='list',required=True)
	args=parser.parse_args()

	args.output.write('library_id,molecule_h5\n')
	for s in args.list.split('/'):
		args.output.write('{0},{1}\n'.format(s,os.path.join(args.input,s,'outs/molecule_info.h5')))

if __name__ == '__main__':
	main()
