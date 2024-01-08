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
	parser.add_argument('-i','--infile',help='infile name',dest='infile',type=str,required=True)
	parser.add_argument('-o','--outfile',help='outfile name',dest='outfile',type=str,required=True)
	args=parser.parse_args()
	out=open(args.outfile,"w")
	with open(args.infile,"r") as file:
		for line in file:
			tmp = line.rstrip("\n").split("\t")
			for j in range(1,len(tmp)):
				out.write("{0}\t{1}\n".format(tmp[j],tmp[0]))
	out.close()


if __name__=="__main__":
	main()
