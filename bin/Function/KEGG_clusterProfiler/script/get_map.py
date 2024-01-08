#!/usr/bin/env python3
"""
do something
"""
import argparse
import os
import sys
import re
import logging
bin = os.path.abspath(os.path.dirname(__file__))
__author__='Ren Xue'
__mail__= 'xueren@genome.cn'
__date__= 'Sat Jan 28 09:39:11 2023'

def main():
	parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter,epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\n'.format(__author__,__mail__,__date__))
	parser.add_argument('-des','--describe',help='ko.txt',dest='ko',type=str,required=True)
	parser.add_argument('-km','--ko2map',help='ko.list',dest='map',type=str,required=True)
	parser.add_argument('-o','--outfile',help='outfile name',dest='outfile',type=str,required=True)
	args=parser.parse_args()
	out=open(args.outfile,"w")
	ko={}
	with open(args.ko,"r") as file:
		for line in file:
			tmp = line.rstrip("\n").split("\t")
			ko[tmp[0]] = tmp[1]
	with open(args.map,"r") as file:
		for line in file:
			tmp = line.rstrip("\n").split("\t")
			tmp2 = tmp[1].split("|")
			if not tmp[0] in ko:continue
			for j in range(len(tmp2)):
				des = ko[tmp[0]]
				out.write("{0}\t{1}\n".format(tmp2[j],des))
	out.close()


if __name__=="__main__":
	main()
