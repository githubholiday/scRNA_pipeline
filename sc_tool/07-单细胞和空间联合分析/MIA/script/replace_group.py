#! /usr/bin/env python3
import argparse
import sys
import os
import re
bindir = os.path.abspath(os.path.dirname(__file__))
import time

__author__='Su Lin'
__mail__= 'linsu@annoroad.com'
__doc__='the decription of program'

pat1=re.compile('^\s+$')

def get_cluster_group(celltype_file):
	cluster_group_dict = {}
	with open(celltype_file, 'r') as infile:
		for line in infile:
			if line.startswith('cluster'):	continue
			tmp = line.rstrip('\n').split('\t')
			cluster_group_dict[tmp[0]] = tmp[1]
	return cluster_group_dict


def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--mia enrichment file',help='input file',dest='input',required=True)
	parser.add_argument('-c','--celltype',help='cluster and celltype',dest='celltype',required=True)
	parser.add_argument('-o','--outfile',help='outfile',dest='outfile',required=True)
	args=parser.parse_args()
	cluster_group_dict = get_cluster_group(args.celltype)
	with open(args.input,'r') as f , open(args.outfile, 'w') as out:
		for line in f:
			if line.startswith('region'):	
				head = line.rstrip('\n').split('\t')
				out.write('\t'.join(head)+'\tgroup\n')
				continue
			tmp = line.rstrip('\n').split('\t')
			cluster = line[0]
			if cluster in cluster_group_dict:
				group = cluster_group_dict[cluster]
			else:
				group = "na"
			tmp.append(group)
			out.write('\t'.join(tmp)+'\n')
if __name__ == '__main__':
	print('Start : ',time.ctime())
	main()
	print('End : ',time.ctime())
