#! /usr/bin/env python3
import argparse
import sys
import time
import os
import re
import pandas as pd
__author__='Liu Huiling'
__mail__= 'huilingliu@genome.cn'
pat1=re.compile('^\s+$')
'''
This program can be used to summarize the genes and calculate the gene number for each go term, 
	based on all expressed genes and differencial expressed genes.

Gene_in_GO: transform "gene in GOs" into "GO contain genes" 
	rules: replace the gos in go.list with go exsists in term
Extend_Go + Extend_Gene: extend the "GO contain genes" based on go.path file
	rules: search from leaf to root
Candidate_in_Go: fileter "GO contain genes" with candidate genes
	rules: "GO contain genes" && genes include in candidate
Merge_GO: merge go that have union keys in Extend_Go result for background and candidate genes
	rules: only go have candidate genes are maintained, 
			calculate the gene number at the same time
Statistic_in_Go: output the m, M, n, N for each go which include candidate genes
Extract_Class: get go term for each Ontology, 'cellular_component','biological_process','molecular_function'
	rules: 

Output :
	prefix.go_with_gene.type
	prefix.count.type
	while type must be one of ['cellular_component','biological_process','molecular_function']

Example:
	python3 this.py -p go.path -s go.alias -c go.class -g go.list -i de.list -o out_prefix
'''
def Read_Alias(alias_file): #read the alias for go
	dict = {}
	with open(alias_file,'r') as f:
		for i,line in enumerate(f):
			tmp = line.rstrip().split('\t')
			for go in tmp[1:]:
				dict[go] = tmp[0]
	f.close()
	return dict

def Read_gene(list_file): #find all gene in a seperate go
	gene_set = set()
	with open(list_file,'r') as f:
		for i,line in enumerate(f):
			if line.startswith('Gene'):continue
			gene_set.add(line.rstrip('\n'))
	f.close()
	return gene_set

def Gene_in_GO(go_file,alias_dict): #find all gene in a seperate go
	warning = set()
	with open(go_file,'r') as f:
		dict = {}
		total_set = set()
		for i,line in enumerate(f):
			tmp = line.rstrip().split('\t')
			total_set.add(tmp[0])
			for go in tmp[1:]:
				if go in alias_dict:
					warning.add("{0}\treplaced by\t{1}".format(go,alias_dict[go]))
					go = alias_dict[go] #replace the alias with go exists in term.txt 
				dict.setdefault(go,set()).add(tmp[0])
	f.close()
	print ('\n'.join(warning))
	return dict,total_set

def Extend_Go(go_with_gene,go_path):
	#go_path = pd.read_table(path_file, low_memory=False, header=None, index_col=0) # col0 = ontology,col1= ontology accession
	#go_path = go_path.drop(go_path.loc[:,[1]], axis=1) #drop col1= ontology accession
	for i,row in go_path.iterrows(): # i, onthology; row, go in path
		row = row.dropna(axis=0,how='all')
		tmp_list = row.values.tolist()
		tmp_list.reverse() #from 0 to len(tmp_list), leaf to root
		while len(tmp_list) > 1:
			if not tmp_list[0] in  go_with_gene:
				tmp_list.remove(tmp_list[0]) #remove leaf which have no gene-annotation
				continue
			parent = tmp_list[1] #parent for leaf
			if not parent in go_with_gene : go_with_gene.setdefault(parent,set())
			go_with_gene[parent].update(go_with_gene[tmp_list[0]])
			tmp_list.remove(tmp_list[0]) #remove leaf
	return go_with_gene

def Candidate_in_Go(go_with_gene,gene_set): ##get candidate gene in a seperate go
	#gene_set = Read_gene(list_file)
	go_with_candidate = {}
	part_set = set()

	for go in go_with_gene:
		candidate = go_with_gene[go] & gene_set # the candidate means gene list have go annotation
		if len(candidate) == 0:continue
		part_set |= candidate
		go_with_candidate[go] = candidate
	return go_with_candidate,part_set

def Merge_Go(background,candidate): #merge go in with both background and candidate genes
	union_keys = set(background.keys()) & set(candidate.keys())
	print('The union go term number is {0}'.format(len(union_keys)))
	union_dict = {}
	count_dict = {}
	for key in union_keys:
		#union_dict[key] =['|'.join(candidate[key]),'|'.join(background[key])]
		union_dict[key] =['|'.join(background[key]),len(background[key]),'|'.join(candidate[key]),len(candidate[key])]
		count_dict[key] =[len(candidate[key]),len(background[key])]
	return union_dict,count_dict

def Statistic_in_Go(count_dict,n,N): #count m and M
	result = pd.DataFrame.from_dict(count_dict, orient='index')
	result['n'] = n
	result['N'] = N
	result.columns = ['m','M','n','N']
	result = result.astype('int')
	return result

def Extract_Class(class_file,info,group_name):
	Ontology = ['cellular_component','biological_process','molecular_function']
	
	info['Accession'] = info.index
	df = pd.read_table(class_file, header = 0, index_col = None)
	group_df = df.groupby(df.Ontology,as_index=False)
	for name,group in group_df:
		if not name in Ontology :
			print('Illegal Ontology in class file !')
			exit()
		if name != group_name : continue
		result = pd.merge(group.loc[:,['Ontology','Term_name','Accession']], info, on='Accession', how ='inner', copy =False)
	return result

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--input',help='input candidate gene list',dest='input',required=True)
	parser.add_argument('-o','--output',help='output prefix',dest='output',required=True)
	parser.add_argument('-g','--go',help='go list in database file',dest='go')
	parser.add_argument('-p','--path',help='path file',dest='path',required=True)
	parser.add_argument('-s','--synonym',help='go alias file',dest='synonym',required=True)
	parser.add_argument('-c','--classify',help='go class file',dest='classify',required=True)
	args=parser.parse_args()
	
	Ontology = ['cellular_component','biological_process','molecular_function']
	# get go alias info
	print("Start Reading {0}".format(args.synonym))
	go_alias_dict = Read_Alias(args.synonym)
	print("Finish Reading {0}".format( args.synonym))

	# get N and M for go list in database
	go_with_background, total_set = Gene_in_GO(args.go,go_alias_dict) #return dict
	N = len(total_set)
	print('The background gene number N = {0}'.format(N))
	
	# get n for candidate genes which is from input list and exits in genes with go annotation
	candidate_set = Read_gene(args.input)
	n = len(candidate_set & total_set)
	print('The candidate gene number n = {0}'.format(n))
	
	go_path = pd.read_table(args.path, low_memory=False, header=None, index_col=0) # col0 = ontology,col1= ontology accession
	go_path = go_path.drop(go_path.loc[:,[1]], axis=1) #drop col1= ontology accession
	
	all_go_with_gene = pd.DataFrame()
	all_go_with_count = pd.DataFrame()
	
	group_path = go_path.groupby(go_path.index,as_index=False) #group by ontology to lower resource used
	for name,group in group_path:
		if not name in Ontology :
			print('Illegal Ontology in class file !')
			exit()
		extend_go_with_background = Extend_Go(go_with_background,group)
		print('The extend_go_with_background number in {1}: {0}'.format(len(extend_go_with_background),name))
		print(time.ctime())

		#get M for go term in database
		extend_go_with_candidate,part_set = Candidate_in_Go(extend_go_with_background,candidate_set) #return dict
		print('The extend_go_with_candidate number in {1}: {0}'.format(len(extend_go_with_candidate),name))

		# get the union go and gene number m in each go term
		union_dict, count_dict = Merge_Go(extend_go_with_background,extend_go_with_candidate)
		union = pd.DataFrame.from_dict(union_dict, orient='index')
		union.columns =['Background_Gene','Background_Count','Candidate_Gene','Candidate_Count']
		gene_result = Extract_Class(args.classify,union,name)
		all_go_with_gene = pd.concat([all_go_with_gene,gene_result])
		
		count = Statistic_in_Go(count_dict,n,N)
		count_result = Extract_Class(args.classify,count,name)
		all_go_with_count = pd.concat([all_go_with_count,count_result])
	
	go_with_gene = '{0}.go_with_gene'.format(args.output)
	all_go_with_gene.to_csv(go_with_gene, sep='\t',header=True, index=False)
	go_with_count = '{0}.go_with_count'.format(args.output)
	all_go_with_count.to_csv(go_with_count, sep='\t',header=True, index=False)

if __name__=='__main__':
	print('Start : ',time.ctime())
	main()
	print('End : ',time.ctime())
