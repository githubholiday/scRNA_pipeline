import os, sys
import argparse

os.getcwd()
os.listdir(os.getcwd()) 

import loompy as lp;
import numpy as np;
import scanpy as sc;
__author__ =  "Tuchengfnag"
__mail__="chengfangtu@genome.cn"

def main():
	parser=argparse.ArgumentParser(
			description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='Author:\t{0}\nE-mail:\t{1}\n'.format(__author__,__mail__)
			)
	parser.add_argument('-i','--infile',dest='infile',help='file of matrix',required=True)
	parser.add_argument('-o','--output',dest='output',help='output of loom file',required=True)
	args=parser.parse_args()
	
	x=sc.read_csv(args.infile)
	row_attrs = {"Gene": np.array(x.var_names),}
	col_attrs = {"CellID": np.array(x.obs_names)}
	lp.create(args.output,x.X.transpose(),row_attrs,col_attrs)
if __name__=='__main__':
	main()