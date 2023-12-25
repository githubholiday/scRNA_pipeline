#! /usr/bin/env python3
import argparse
import sys
import os
import re
import datetime
import pandas as pd
import glob

bindir = os.path.abspath(os.path.dirname(__file__))
filename=os.path.basename(__file__)

__author__='tu chengfang '
__mail__= 'chengfangtu@genome.cn'


class Log():
	def __init__( self, filename, funcname = '' ):
		self.filename = filename 
		self.funcname = funcname
	def format( self, level, message ) :
		date_now = datetime.datetime.now().strftime('%Y%m%d %H:%M:%S')
		formatter = ''
		if self.funcname == '' :
			formatter = '\n{0} - {1} - {2} - {3} \n'.format( date_now, self.filename, level, message )
		else :
			
			formatter = '\n{0} - {1} - {2} -  {3} - {4}\n'.format( date_now, self.filename, self.funcname, level, message )
		return formatter
	def info( self, message ):
		formatter = self.format( 'INFO', message )
		sys.stdout.write( formatter )
	def debug( self, message ) :
		formatter = self.format( 'DEBUG', message )
		sys.stdout.write( formatter )
	def warning( self, message ) :
		formatter = self.format( 'WARNING', message )
		sys.stdout.write( formatter )
	def error( self, message ) :
		formatter = self.format( 'ERROR', message )
		sys.stderr.write( formatter )
	def critical( self, message ) :
		formatter = self.format( 'CRITICAL', message )
		sys.stderr.write( formatter )


def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--input',help='input file',dest='input',required=True)
	parser.add_argument('-t','--type',help='type',dest='type',required=True)
	parser.add_argument('-n','--name',help='source',dest='name',required=True)
	parser.add_argument('-o','--output',help='output',dest='output',required=True)
	args=parser.parse_args()
	
	sample_name = os.path.basename(args.input).split('_net')[0]
	with open( args.input, 'r') as infile, open(args.output, 'w') as outfile :
		selected_index = 0
		out_list = [0,1,4,5,7,8,9,10]
		for line_index, line in enumerate(infile) :
			tmp = line.rstrip().split('\t')
			if line_index == 0 :
				for i_index, value in enumerate( tmp ) :
					if value == args.type :
						selected_index = i_index
				out_head = [tmp[i] for i in out_list]
				header = [sample_name]+out_head
				outfile.write( '\t'.join(header) + '\n' )
			else :
				if tmp[selected_index] == args.name :
					value = ["1"] + [tmp[i] for i in out_list]
					outfile.write( '\t'.join(value) + '\n' )
				else:
					continue
					
if __name__ == '__main__':
	main()