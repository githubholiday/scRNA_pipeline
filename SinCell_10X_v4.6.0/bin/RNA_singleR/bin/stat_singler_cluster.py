#! /usr/bin/env python3
import argparse
import time
import sys
import re
import os
import csv
import pandas as pd
bindir = os.path.abspath(os.path.dirname(__file__))

__author__ = ''
__mail__ = ''
__doc__ = 'the description of program'

pat1=re.compile('^s+$')
'''
统计每个cluster中细胞类型
'''
def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-c','--cell',help='cell file',dest='cell', required=True)
	parser.add_argument('-o','--output',help='output dir',dest='output',required=True)
	args=parser.parse_args()

	#info.columns=['Sample','Number of Cells','Mean Reads per Cell','Median Genes per Cell']
	df = pd.read_csv(args.cell,index_col = None,sep='\t')
	df = df[["Barcode","labels","pruned.labels","clusters","stim"]]
	df_tmp=df.groupby(["clusters","labels"])["Barcode"].count().sort_values(ascending=False).reset_index()
	df_tmp_new1=df_tmp.groupby("clusters")["Barcode"].sum().reset_index()
	df_tmp_new=df_tmp_new1.rename(columns={"Barcode":"cluster_totalcells"})
	df_merge=pd.merge(df_tmp,df_tmp_new,on="clusters",how='left')
	#print(df_merge)
	df_merge["precent(%)"]=df_merge["Barcode"]/df_merge["cluster_totalcells"]
	df_merge["precent(%)"]=df_merge["precent(%)"].apply(lambda x:round(x*100,2))
	df_merge=df_merge.sort_values(by=["clusters","Barcode"],ascending=[1,0]).reset_index()
	#df_merge=df_merge[df_merge.groupby("clusters")["Barcode"].head(3).reset_index()]
	#print(df_merge)
	df_merge=df_merge[["clusters","labels","Barcode","cluster_totalcells","precent(%)"]]
	df_merge.to_csv(args.output+"/singleR_celltype_stat.xls", sep='\t',index=0)
	#print(df_merge)
	#df_merge=df_merge[df_merge["precent(%)"]>30]
	grouped=df_merge.groupby(['clusters']).agg({'precent(%)':'max'}).reset_index()
	#print(grouped)
	df_merge_new=pd.merge(grouped,df_merge,on=["clusters",'precent(%)'],how='left')
	df_merge_new.to_csv(args.output+'/singleR_celltype_stat_top1.xls', sep='\t',index=0,columns=["clusters","labels","Barcode","cluster_totalcells","precent(%)"])

if __name__ == '__main__':
	main()
