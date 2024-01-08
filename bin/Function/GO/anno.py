#! /usr/bin/env python3
import argparse
import sys
import os
import re
bindir = os.path.abspath(os.path.dirname(__file__))

__author__='Liu Tao'
__mail__= 'taoliu@annoroad.com'

pat1=re.compile('^\s+$')

def parseId(name):
	return name.split(r'|')[1]

def readAnno(f_file, col, type):
	r_dict = {}
	header = type
	for line in f_file:
		if line.startswith('#') or re.search(pat1,line):continue
		tmp=line.rstrip().split('\t')
		if line.startswith('UniProtKB-AC') : 
			if col == None : header = line.rstrip()
			continue
		if col == None:
			r_dict[tmp[0]] = tmp
		else:
			if tmp[col] == '' :
				continue
			else:
				r_dict[tmp[0]] = [i.replace(" ",'') for i in tmp[col].split(';')]
	return r_dict,header


def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--input',help='input file',dest='input',type=open,required=True)
	parser.add_argument('-a','--anno',help='anno file',dest='anno',type=open,required=True)
	parser.add_argument('-o','--output',help='output file',dest='output',type=argparse.FileType('w'),required=True)
	parser.add_argument('-c','--col',help='col file',dest='col',type=int)
	parser.add_argument('-t','--type',help='annot type',dest='type')
	args=parser.parse_args()

	anno , header = readAnno(args.anno, args.col,args.type)

	if args.col == None:
		args.output.write('{0}\t{1}\t{2}\t{3}\n'.format("GeneID","Score","Evalue",header))
	for line in args.input:
		if line.startswith('#') or re.search(pat1,line):continue
		tmp=line.rstrip().split('\t')
		u_id = parseId(tmp[1])
		if u_id in anno:
			if args.col == None:
				args.output.write('{0}\t{1}\t{2}\t{3}\n'.format(tmp[0],tmp[11],tmp[10],"\t".join(anno[u_id])))
			else:
				args.output.write('{0}\t{1}\n'.format(tmp[0],"\t".join(anno[u_id])))


if __name__ == '__main__':
	main()
