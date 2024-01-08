#! /usr/bin/env python
'''
Description:
	this file is used to get overlap for 2 file, and output intercouse line
'''
import argparse
import sys
import os

__author__='Liu Tao'
__mail__='taoliu@annoroad.com'
def Overlap(dict1,dict2):
	out_dict={}
	for i in dict1:
		if i in dict2:
			#out_dict[i]="{0}\t{1}".format(dict1[i],dict2[i])
			out_dict[i]="{0}".format('\t'.join(dict2[i]))
	return out_dict

def store_record(id,record,args,count,line,q_value):
	tmp=line.split('\t')
	if not id in record:
		record[id]={}
	if args.col[count] == 'all':
		record[id][count] = line
	else:
		try:
			a_col = int(args.col[count])
			record[id][count] = tmp[a_col]
			if float(tmp[a_col]) < q_value:
				return 1
			else:
				return 0 
		except :
			print(args.col[count],'should be int')
			sys.exit()
	#return record

def main():
	parser=argparse.ArgumentParser(
			description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='Author:\t{0}\nE-mail:\t{1}\n'.format(__author__,__mail__)
			)
	parser.add_argument('-f','--files',dest='files',help='input more than ones file',required=True,nargs='+')
	parser.add_argument('-head','--header',dest='header',help='header line number [0]',default=[0],type=int,nargs ='+')
	parser.add_argument('-i','--index',dest='index',help='which col is index [0] ',default=[0],type=int,nargs='*')
	parser.add_argument('-c','--col',dest='col',help='which col should be store [all] ',nargs='*',required=True)
	parser.add_argument('-n','--name',dest='name',help='names ',nargs='*',type=str,default=[])
	parser.add_argument('-s','--sep',dest='sep',help='separte symbol in id column',default='')
	parser.add_argument('-q','--qvalue',dest='q_value',help='q value [ 0.5] ',default=0.05,type=float)
	parser.add_argument('-t','--type',dest='type',help='BP,CC,MF,kegg',required=True)
	args=parser.parse_args()
	
	#if len(args.files)<2 :
	#	print('please input more than one files',file=sys.stderr)
	#	sys.exit()
	#feature=['_CC','_BP','_MF','_kegg_']
	if args.type == 'kegg':
		feature = args.type+'_'
	feature = '_'+args.type
	count=0
	record={}
	col_count=0
	items = []
	for file in args.files:
		tag=0
		try :
			args.index[count]
		except IndexError:
			args.index.append(args.index[count-1])

		try:
			args.col[count]
		except IndexError:
			args.col.append(args.col[count-1])

		try:
			args.header[count]
		except IndexError:
			args.header.append(args.header[count-1])

		try :
			args.name[count]
		except IndexError:
			basename = os.path.basename(file)
			name = basename.rstrip('kegg.report.xls')
			#[name] = [name.split(feature)[0] for i in name ]
			name = name.split(feature)[0]
			#name = '_'.join(basename.split('_')[1:3])
			args.name.append(name)

		with open(file,'r') as f_file:
			for line in f_file:
				tag+=1
				if args.header[count] >= tag :continue
				line=line.rstrip()
				tmp=line.split('\t')
				if not col_count : col_count =len(line)
				ids=tmp[args.index[count]]
				ids = ids.replace(r"'",'')
				ids = ids.replace(r'"','')
				if args.sep:
					ids = ids.split(args.sep)
					for id in ids:
						if not id in items:items.append(id)
						store_record(id,record,args,count,line,args.q_value)
				else:
					if store_record(ids,record,args,count,line,args.q_value) :
						if not ids in items:
							items.append(ids)
		count+=1
	
	print("name\t"+"\t".join(args.name))
	if 'name' in record:
		output='name\t'
		for file_num in range(len(args.files)):
			if file_num in record['name']:
				output+=record['name'][file_num]+'\t'
			else:
				output+='miss\t'
		output=output.rstrip()
		print(output)

	for id in items:
		if id == 'name' : continue
		output=id+'\t'
		for file_num in range(len(args.files)):
			if file_num in record[id]:
				output+=record[id][file_num]+'\t'
			else:
				output+='1\t'
		output=output.rstrip()
		print(output)

if __name__=='__main__':
	main()
