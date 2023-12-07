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

def AnnoCount(line_list):
	num_dash=line_list.count('--')
	num_dot=line_list.count('.')
	return (num_dash+num_dot)

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--input',help='input file',dest='input',required=True,nargs='+')
	parser.add_argument('-o','--output',help='output file',dest='output',type=argparse.FileType('w'),required=True)
	parser.add_argument('--yes','--yes',help='if yes in line',action='store_true')
	args=parser.parse_args()

	out = []
	anno_count = 100
	anno_line = ''
	for file in args.input:
		for index,line in enumerate(open(file,'r')):
			if index == 0:
				lines = line.rstrip('\n').split('\t')
				lines[0] = upperFirstWord(lines[0])
				lines = [ i for i in lines if i !="" ]
				out.append(lines)
				continue
			if index == 1:
				anno_line = line
				continue
			lines = line.rstrip('\n').split('\t')
			count = AnnoCount(lines)
			if count < anno_count:
				anno_count = count
				if args.yes:
					if 'yes' in lines:
						anno_line = line
				else:
					anno_line = line
		out.append(anno_line.rstrip('\n').split('\t'))
		my_out=list(map(list,zip(*out)))
		for lines in my_out:
			if len(lines[1])>50:
				lines[1] = lines[1][0:50]+'...'
			args.output.writelines("%s\n" %('\t'.join(lines)))
		break

def upperFirstWord(inStr):
	return "%s" % (inStr[:1].upper() + inStr[1:])

if __name__ == '__main__':
	print('Start : ',time.ctime())
	main()
	print('End : ',time.ctime())
