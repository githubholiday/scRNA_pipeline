#! /usr/bin/env python3
import argparse
import time
import sys
import re
import os
import pandas as pd
import numpy as np
bindir = os.path.abspath(os.path.dirname(__file__))

__author__ = 'Liu Huiling'
__mail__ = 'huilingliu@genome.cn'
__doc__ = 'the description of program'
'''
The up and down gene count in second_level_go 
'''
pat1=re.compile('^s+$')

def Read_Gene(file):
	gene_set =set()
	with open(file,'r') as f:
		for line in f:
			gene_set.add(line.rstrip())
	return gene_set

def Read_Go(file): # second_level_go
	df = pd.read_table(file, low_memory=False, index_col=0, header=None)
	result = df.iloc[:,1:2].drop_duplicates()
	result.columns = ['Accession']
	return result

def Up_Down_Count(infile,up_list,down_list):
	df = pd.read_table(infile, index_col=False, header=0)
	#background = df.loc[:,['Accession','Background_Gene']]
	de = df.loc[:,['Accession','Candidate_Gene']]
	'''
	back_count = df.loc[:,['Accession','Background_Count']]
	back_dict = {}
	for i,row in background.iterrows():
		tmp = row.values.tolist()
		back_dict[tmp[0]] = [tmp[0],len(tmp[1].split('|'))]
	back_count = pd.DataFrame.from_dict(back_dict, orient='index')
	back_count.columns = ['Accession','Background_Count']
	'''

	dict = {}
	for i,row in de.iterrows():
		tmp = row.values.tolist()
		#print(tmp)
		count,genes = Gene_Count(tmp[1],up_list,down_list)
		dict[tmp[0]] = [tmp[0],genes['de'],count['de'],genes['up'],count['up'],genes['down'],count['down']]
	up_down_count = pd.DataFrame.from_dict(dict, orient='index')
	up_down_count.columns = ['Accession','DE_Gene','DE_Count','Up_Gene','Up_Count','Down_Gene','Down_Count']

	result = pd.merge(df, up_down_count, on='Accession')
	result = result.drop(['Candidate_Gene','Candidate_Count'], axis=1)
	return result

def Gene_Count(genes,up_list,down_list):
	up_set = set()
	down_set = set()
	count = {}
	gene_set = {}
	for g in genes.split('|'):
		if g in up_list:
			up_set.add(g)
		elif g in down_list:
			down_set.add(g)
		else:
			print('unexpected error!')
			exit()
	count['up'] = len(up_set)
	count['down'] = len(down_set)
	count['de'] = len(up_set|down_set)
	gene_set['up'] = '|'.join(up_set)
	gene_set['down'] = '|'.join(down_set)
	gene_set['de'] = genes
	for key in count:
		if count[key] == 0:
			gene_set[key] = '--' 
	return count,gene_set

def Percentage(gene_count,second_level_go,up_count,down_count):
	up_down_count = gene_count.loc[:,['Ontology','Term_name','Accession','Up_Count','Down_Count']]
	
	result = pd.merge(second_level_go, up_down_count, on='Accession')
	up_rate = result['Up_Count']/up_count*100 
	result.insert(4,'Up_Percent',up_rate.apply(lambda x: '{0:.2f}'.format(x)))
	down_rate = result['Down_Count']/down_count*100
	result['Down_Percent'] = down_rate.apply(lambda x: '{0:.2f}'.format(x))
	result.drop('Accession', axis=1, inplace=True)
	result = result.replace('nan',0)
	print('The annotated second_level_go number: {0}'.format(len(result.index)))
	return result

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-i','--input',help='input file',dest='input',required=True)
	parser.add_argument('-o','--output',help='output file',dest='output',required=True)
	parser.add_argument('-u','--up',help='up gene list file',dest='up')
	parser.add_argument('-d','--down',help='down gene list file',dest='down')
	parser.add_argument('-p','--path',help='path file',dest='path')
	args=parser.parse_args()

	up_set = Read_Gene(args.up)
	print('The total Up gene number: {0}'.format(len(up_set)))
	down_set = Read_Gene(args.down)
	print('The total Down gene number: {0}'.format(len(down_set)))

	second_level_go = Read_Go(args.path)
	print('The total second_level_go number: {0}'.format(len(second_level_go.index)))
	gene_count = Up_Down_Count(args.input,up_set,down_set)
	gene_count.to_csv(args.output, sep='\t', index=False, header =True)

	percentage = Percentage(gene_count, second_level_go, len(up_set), len(down_set))
	percentage.to_csv(args.output+'.xls', sep='\t', index=False, header =True,na_rep=0)

if __name__ == '__main__':
	main()
