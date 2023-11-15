#! /usr/bin/env python3
import argparse
import sys
import re
import pandas as pd
import os
import logging
bindir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(bindir)
from myplot import *

__author__='zhang yue'
__mail__= 'yuezhang@genome.cn'
__date__= '2018/11/20'

pat1=re.compile('^\s+$')

LOG = os.path.basename(__file__)

def my_log( level, message ) :
	logging.basicConfig(level = logging.INFO,format = '%(asctime)s - %(filename)s - %(levelname)s - %(message)s')
	logger = logging.getLogger(__name__)
	if level == 'info' :
		return logger.info( message )
	if level == 'warning' :
		return logger.warning( message )
	if level == 'debug' :
		return logger.debug( message )
	if level == 'error' :
		return logger.error( message )

def check_file_exists( *file_list ) :
	for file in file_list :
		if os.path.exists( file ) :
			my_log( 'info', 'file : {0}'.format( file ) )
		else :
			my_log( 'error', 'file is not exists : {0}'.format( file ) )
			sys.exit(1)

def make_dir( dir ) :
	try :
		os.makedirs( dir )
		time.sleep(1)
		my_log( 'info', 'mkdir {0} sucessful!'.format( dir) )
	except :
		my_log( 'error', 'mkdir {0} failed!'.format( dir) )

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--input',help='input file',dest='input',required=True)
	parser.add_argument('-o','--outfile',help='output file',dest='outfile',required=True)
	parser.add_argument('-x','--datax',help='x data',dest='datax')
	parser.add_argument('-y','--datay',help='y data, if more than 1, separated by ","',dest='datay')
	parser.add_argument('-t','--title',help='title',dest='title',default='FPKM Distribution')
	parser.add_argument('-ylab','--ylabel',help='ylabel',dest='ylabel',default='FPKM Value')
	parser.add_argument('-xlab','--xlabel',help='xlabel',dest='xlabel',default='Samples')
	args=parser.parse_args()

	check_file_exists( args.input )
	outdir = os.path.dirname( args.outfile )
	if not os.path.exists(outdir):
		my_log( 'info', '{0}  directory is not exists, i creat it'.format(outdir) )
		make_dir( outdir )

	matrix = pd.read_table(args.input, header=0, index_col=0, encoding='utf-8')
	tt = myPlot( matrix )
	args.datay = args.datay.split(',') if args.datay else []
	tt.box_plot( x=args.datax, y= args.datay,xlabel=args.xlabel, ylabel=args.ylabel, title=args.title, filename=args.outfile )

if __name__ == '__main__':
	main()
