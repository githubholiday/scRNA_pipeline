#!/usr/bin/evn python3
######################################################################### import ##########################################################
import argparse
import os
import sys
#import pymysql as mariadb
import mysql.connector as mariadb
import pandas as pd
import configparser
from itertools import islice
######################################################################### ___  ##########################################################
__author__ = 'Liu San Yang/Zhang yue'
__mail__ = 'sanyangliu@genome.cn'
__date__ = '2019年08月12日 星期一 11时30分40秒'
__version__ = '1.0'
######################################################################### main  ##########################################################

def insert_table(args,table_name,config,cursor,conn):
	with open (args.input) as f :
		for  i in islice(f,1,None):
			i=i.strip("\n")
			line=i.split("\t")
			column_use=[]
			try:
				column_use=list(map(int,config.get(args.key, "column").split("\n")))
			except configparser.NoOptionError as err:
				print (err)
			insert_line=[ line[i] for i in column_use ]
			if args.primary:
				insert_line.insert(0,args.primary)
			print(config[args.key]['insert'].format(table_name),insert_line)
			cursor.execute(config[args.key]['insert'].format(table_name),insert_line)

def create_table(table_name,config,cursor,conn):
	sql=config['TABLE_INFO']['create'].format(table_name)
	print(sql)
	cursor.execute(sql)
def clean_table(table_name,config,cursor,conn):
#	sql="DELETE FROM {0}".format(table_name)
	print(sql)
	cursor.execute(sql)
def alter_table(table_name,config,cursor,conn):
	sql=config['TABLE_INFO']['alter'].format(table_name,config['TABLE_INFO']['column_name'])
	print(sql)
	cursor.execute(sql)
def update_table(args,table_name,config,cursor,conn):
	with open (args.input) as f :
		for  i in islice(f,1,None):
			i=i.strip("\n")
			line=i.split("\t")
			column_use=[]
			try:
				column_use=list(map(int,config.get(args.key, "column").split("\n")))
			except configparser.NoOptionError as err:
				print (err)
			insert_line=[ line[i] for i in column_use ]
			insert_line.insert(0,table_name)
			sql=config[args.key]['update'].format(*insert_line)
			print(sql)
			cursor.execute(sql)

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter,
		epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
	parser.add_argument('-i',help='tsv/csv file',dest='input',type=str,required=True)
	parser.add_argument('-c',help='config file set to write db',dest='config',type=str,required=True)
	parser.add_argument('-p',help='primary key,if give will add to first column',dest='primary',type=str)
	parser.add_argument('-u',help='update table mode',action="store_true")
	parser.add_argument('-n',help='new table mode',action="store_true")
	parser.add_argument('-r',help='force recover table mode',action="store_true")
	parser.add_argument('-t',help='table name,not use yet',dest='table_name',type=str,required=True)
	parser.add_argument('-e',help='table exists mode',action="store_true")
	parser.add_argument('-k',help='tsv colum info section name',dest='key',type=str,required=True)
	args=parser.parse_args()
	config = configparser.ConfigParser()
	config.read(args.config)
	args=parser.parse_args()
	conn = mariadb.connect(user="sci_bioinfo",password="McIFBcxwU#9Y",host="192.168.169.37",db="sci_bioinfo",autocommit=True)
	cursor = conn.cursor()
	if not args.n and not args.e and not args.r and not args.u:
		exit("-n  -e -f -u 四个参数不能同时为空")
	if args.n or args.r or args.e:
		try:
			create_table(args.table_name,config,cursor,conn)
		except mariadb.OperationalError as err:
			if args.e:
				print ("表格已经存在,续写模式！")
			elif args.r:
				clean_table(args.table_name,config,cursor,conn)
				config['CSV_INFO']['insert']=config['CSV_INFO']['insert'].replace("insert","replace")
			else:
				exit(err)
		except KeyError as err:
			print("没有设置TABLE_INFO")
			print(err)
		insert_table(args,args.table_name,config,cursor,conn)
		cursor.close()
		conn.close()
	if args.u:
		try:
			alter_table(args.table_name,config,cursor,conn)
		except mariadb.OperationalError as err:
			print(err)
		except KeyError as err:
			print("config文中没定义key:{0},忽略".format(err))
		update_table(args,args.table_name,config,cursor,conn)
if __name__=="__main__":
	main()
