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
__date__= 'Mon 31 Aug 2020 05:40:57 PM CST'



def main():
	parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter,epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\n'.format(__author__,__mail__,__date__))
	parser.add_argument('-f','--file',help='input name',dest='input',type=str,required=True)
	parser.add_argument('-i','--indir',help='indir',dest='indir',type=str,default=".")
	args=parser.parse_args()
	if os.path.exists(args.input):
		inf=open(args.input,"r")
		for line in inf:
			tmp=line.rstrip("\n").split("\t")
			path=args.indir + "/" + tmp[-1].rstrip("/")
			os.system("mkdir -p {0}".format(path))
			for i in range(len(tmp)-1):
				if "tar.gz" in tmp[i]:
					os.system("tar -zxvf {0} -C {1}/ ".format(tmp[i],path))
				else:
					os.system("cp -rf {0} {1}/".format(tmp[i],path))
	else:
		print("input file do not exist\n")

if __name__=="__main__":
	main()
