#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
使用说明
	该程序根据提供的质控等级表，对受质控的文件中样品的质控指标进行等级判断,并将判断结果存入数据库中.

	参数说明：
	
	  -i: 输入目录(项目分析目录)
	  -t: 质控标准文件
	  -p: 项目名称（作为数据库名称）
	  -n: 质控内容(质控内容以逗号分隔;仅包含质控指标名称的单列文件也可以)
	  -s: 样本名(质控样本名称,以逗号分隔;仅包含样本名称的单列文件也可以)
	  -e: 邮件发送(yes/no/B/C.'yes':不考虑结果均发邮件,'no':不考虑结果均不发邮件;B/C:判断发邮件,含有或低于此等级发邮件,如C:当结果存在C或者D时才发邮件,且标题为相应的等级(合格,让步,不合格,终止),如果结果只包含A,B则不发邮件.
	  -m: 邮件配置文件(设定发件人,收件人,格式和home目录下的~/.email/.email.txt一致)
	  -d: 数据库db文件输出路径
	  -f: 项目完成选项,如果添加该参数,将输出多次运行后数据库中所有指标和样品的汇总结果
	  -o: 质控结果输出目录+输出前缀
	  
	  示例1：python3 /annoroad/data1/bioinfo/PMO/zhoufulai/QC/QC/qc.py -i /my/dir2 -t QC_config -p HA-01 -n 'Clean_Reads_Rate,Clean_Bases_Number' -s A1,A2 -e no -m /my/dir2/email.txt -d /my/dir2 -o /my/dir2/result
	  示例2: python3 /annoroad/data1/bioinfo/PMO/zhoufulai/QC/QC/qc.py -i /my/dir2 -t QC_config -p HA-01 -n qc.list -s sample.list -e yes -m /my/dir2/email.txt -d /my/dir2 -o /my/dir2/result
	  示例3: python3 /annoroad/data1/bioinfo/PMO/zhoufulai/QC/QC/qc.py -i /my/dir2 -t QC_config -p HA-01 -n qc.list -s sample.list -e B -m /my/dir2/email.txt -d /my/dir2 -o /my/dir2/result
	
	质控标准文件：
		
	  1. 第一列为质控指标； 
	  2. 质控指标如果包含空格和'-'会被替换成'_',如果包含'%','()','>','<','='会被删除；
	  3. 质控指标名称需与文件中名称对应，包括大小写及空格，否则会出现错误提示;
	  4. 第二、三、四、五列为A,B,C,D档，分别代表合格,让步,终止,终止;
	  5. 数值中不包含百分号(%)和千分位(,);
	  6. 判断区间使用小括号，判断区间中必须包含上、下限;
	  7. 同一档中包含多个判断区间，判断区间之间使用分号分割;
	  8. 质控内容可以忽略大小写.如某一档内容为YES时,文件内为yes，认为相同;
	  9. 第六列表头必须为'Indir',该列填写质控指标对应的结果文件所在路径,程序将填写路径中的'Indir'关键字替换成-i参数中输入的路径;用于文件捕获；
	  10. 填写'Indir'那一列应注意，表头必须为'Header'如果填写路径包含样本名称目录，样本名称所在目录需用关键字'SAMPLE'代替；
	  11．最后一列为报告内文件样品放置形式，如果文件第一行为样品名则填写R，如果文件第一列为样品名则填写C;
	  
	更新说明：
	  2017/03/24
		1.将输出outdir改为output(outdir+prefix)
		2.数据库中添加表LOG记录多次运行结果，表QUALITY,SCORE和DIR只保留最新结果
	  2017/08/01
		1.替换发邮件模块，之前邮件模块开发使用perl语言，现更换为python版本
		2.增加html模块，邮件中不合格内容标红
		3.增加发邮件函数，实现从跳转机上发送邮件
		4.修改数据库，删除数据库中的子表，只保留总表；避免相同模块单个样本同时质控，子表写入冲突
	  2017/08/10
		1.修复html中samples和quals未进行sorted，导致邮件中表头显示异常
	 
'''
import os
import sys
import re
import glob
import argparse
import sqlite3
import logging
import datetime
import shutil
import getpass
import configparser
import subprocess
import zipfile

bindir=os.path.abspath(os.path.dirname(__file__))

__author__='zhoufulai and renxue'
__mail__='zhoufulai@genome.cn and renxue@genome.cn'

pat1=re.compile('^\s*$')
pat2=re.compile('(\d+,\d+)')
pat3=re.compile(';')
	
class Database():
	def __init__(self, path ):
		self.database = path
		self.conn = sqlite3.connect( self.database )
	def check_table_exists(self ,table_list, index_list):
		cmd = 'select name from sqlite_master where type = "table";'
		cursor = self.conn.execute(cmd)
		tables = [ j for i in cursor.fetchall() for j in i]
		r_value = True
		if "LOG" not in  tables:
			self.create_LOG()
		for i in table_list :
			if len(tables) ==  0 or  not i in tables:
				self.create_database( i , index_list)
				r_value = False
		return r_value 	
		
	def create_database(self , table , index_list):
		cmd = ''
		header=["Sample	TEXT	PARAMETER KEY"]
		for i in index_list:
			header.append('{0}\tTEXT'.format(i))
			header_new="," .join([str(i) for i in header])
		cmd = 'CREATE TABLE {1} ({0})'.format(header_new, table)
		self.conn.execute(cmd)
		self.commit()
	def create_LOG(self):
		cmd = ''
		cmd = '''CREATE TABLE LOG 
			(BLOCK	TEXT	NOT NULL , 
			QCINDEX	TEXT	NOT NULL , 
			SAMPLE	TEXT	NOT NULL , 
			SCORE	TEXT	NOT NULL , 
			QUALITY	TEXT	NOT NULL , 
			TIME	TEXT	NOT NULL , 
			DIR	TEXT	NOT NULL);'''
		self.conn.execute(cmd)
		self.commit()

	def check_max_table(self):
		cmd = 'select name from sqlite_master where type = "table";'
		cursor = self.conn.execute(cmd)
		tables = [ j for i in cursor.fetchall() for j in i]
		p = re.compile('QUALITY_(\d+)$')
		max=0
		for i in tables :
			m = re.search(p,i)
			if m:
				num=m.group(1)
				if int(num)>max:
					max=int(num)
		return max 	
	def check_field(self , table_name):
		cmd ='PRAGMA table_info({0});'.format(table_name)
		cursor = self.conn.execute(cmd)
		tt = cursor.fetchall()
		field=[]
		for i in tt:
			field.append(i[1])
		return field	
	def alter_field(self , table_name , field):
		cmd='ALTER TABLE {0} ADD {1} TEXT '.format(table_name,field)
		self.conn.execute(cmd)
		self.commit()
	def query(self , table_name , key , a_value , cols):
		cmd = 'SELECT {3} FROM {0} WHERE {1} = "{2}"; '.format(table_name , key , a_value , ",".join(cols))
		cursor = self.conn.execute(cmd)
		tt = cursor.fetchall()
		return tt
	def query_value(self , table_name, sample, index ):
		cmd = 'SELECT {3} FROM {0} WHERE {1} = "{2}"; '.format(table_name , "Sample" ,sample, index)
		cursor = self.conn.execute(cmd)
		tt = cursor.fetchall()
		value=''
		for i in tt:
			value=i[0]
			return value
	def query_sample(self , table_name ):
		cmd = 'SELECT Sample FROM {0}'.format(table_name)
		cursor = self.conn.execute(cmd)
		tt = cursor.fetchall()
		sample_list=[]
		for i in tt:
			sample_list.append(i[0])
		return sample_list
	def query_part(self,table,key,index,keys):
		#keys_new=[str(key)+"="+str(i) for i in keys]
		keys_new=['{0}="{1}"'.format(key,i) for i in keys]
		cmd = 'SELECT {0} from {1} where {2}'.format(",".join(index),table," or ".join(keys_new))
		cursor = self.conn.execute(cmd)
		tt = cursor.fetchall()
		return tt
	def query_all(self , table_name):
		cmd = 'SELECT * FROM {0};'.format(table_name)
		cursor = self.conn.execute(cmd)
		tt = cursor.fetchall()
		return tt
	def insert(self , table_name , name_list , value_list):
		names =",".join(['"' + str(i) + '"' for i in  name_list])
		values = ",".join(['"' + str(i) + '"' for i in  value_list])
		cmd = 'INSERT INTO {0} ({1}) VALUES ({2});'.format(table_name , names , values)
		self.conn.execute(cmd)
		self.commit()
	def delete(self , table_name , key , a_value):
		cmd = 'DELETE FROM {0} WHERE {1} = {2} ; '.format(table_name , key , a_value)
		self.conn.execute(cmd)
		self.commit()
	def delete_tables(self , table_name_list):
		for i in table_name_list:
			cmd = 'DROP TABLE {0}; '.format(i)
			self.conn.execute(cmd)
			self.commit()
	def rename(self , table_name , table_name_new):
		cmd = 'ALTER TABLE {0} RENAME TO {1} ;  '.format(table_name, table_name_new)
		self.conn.execute(cmd)
		self.commit()	
	def cptable(self , table_name , table_name_new):
		cmd = 'CREATE TABLE {0} as select * from {1} ;  '.format(table_name, table_name_new)
		self.conn.execute(cmd)
		self.commit()	
	def delete_all_content(self , table_name ):
		cmd = 'DELETE  FROM {0} ;  '.format(table_name )
		self.conn.execute(cmd)
		self.commit()
	def update(self , table_name , record , key , a_value):
		cmd = ''
		for i , j in record.items():
			cmd = 'UPDATE {0} set {1} = "{2}" WHERE {3} = "{4}";'.format(table_name , i, j , key , a_value )
			self.conn.execute(cmd)
			self.commit()
	def export_all(self , outfile , table_name):
		all_values=self.query_all(table_name)
		for i in all_values:
			for j in range(len(i)):
				outfile.write('{0}\t'.format(i[j]))
			outfile.write('\n')
	def export(self , outfile , table_name, key , keys , indexs):
		all_values=self.query_part(table_name,key,indexs,keys)
		for i in all_values:
			for j in range(len(i)):
				outfile.write('{0}\t'.format(i[j]))
			outfile.write('\n')
	def commit(self):
		self.conn.commit()
	def close(self):
		self.conn.close()

		
def string_convert(string):
	#replace the whitespace and - with _
	if re.search('[\s+|-]',string):
		string=re.sub('[\s+|-]','_',string)
		#replace () % （） > < = with ''
		if re.search('[\(|\)|%|\（|\）|>|<|=|\~|\!|\@|\#|\$|\^|\&|\*|\+|\:|\;|\,|\.|\'|\"]',string):
			string=re.sub('[\(|\)|%|\（|\）|>|<|=|\~|\!|\@|\#|\$|\^|\&|\*|\+|\:|\;|\,|\.|\'|\"]','',string)
	elif re.search('[\(|\)|%|\（|\）|>|<|=|\~|\!|\@|\#|\$|\^|\&|\*|\+|\:|\;|\,|\.|\'|\"]',string):
		string=re.sub('[\(|\)|%|\（|\）|>|<|=|\~|\!|\@|\#|\$|\^|\&|\*|\+|\:|\;|\,|\.|\'|\"]','',string)
	return string
	
def check_input(input_para):
	in_list=[]
	if os.path.isfile(input_para):
		in_file=open(input_para,'r')
		for line in in_file:
			if line.startswith('#') or pat1.search(line):continue
			line_content=line.strip()
			in_list.append(line_content)
		in_file.close()
	elif isinstance(input_para,str):
		if re.search('[,|;]',input_para):
			in_list=re.split('[,|;]',input_para)
		else:
			in_list.append(input_para)
	else:
		logging.error ("Please check your input parameter {0}".format(input_para))
		sys.exit()
	return in_list
	
def check_exists(input):
	if not os.path.exists(input):
		logging.error("{0} does not exist,please check it.".format(input))
		sys.exit()
		
def read_template(t_file):
	index_value={} 
	index_type={}
	num=1
	#debug_log.info('Read template')
	template=open(t_file,'r')
	for line in template:
		if num==1:
			temp1=line.rstrip().split('\t') #get the header of the table
			num +=1
			continue
		else:
			temp2=line.rstrip().split('\t')
			#if len(temp2)!=6:
			#	logging.error("{0} does not have 6 columns,please check your template.".format(temp2[0]))
			#	sys.exit()
			for j in range(1,len(temp2)): 		#
				if temp2[j]=='':
					logging.error("{0}_{1} is empty,please check you template.".format(temp2[0],temp1[j]))
					sys.exit()
				index_type[temp1[j]]=temp2[j]
			# replace the '-' and space of types into '_' 
			type=string_convert(temp2[0])
			index_value[type]=index_type
			index_type={}
			num +=1
			continue
	template.close()
	for i in index_value:
		if "Standard" in index_value[i]:
			index_value[i]["Standard"]=index_value[i]["Standard"].replace("%","%%")

	#print (index_value)
	return index_value
def read_file(r_file):
	sample_index_value={}
	sample_index_type={}
	num=1
	rfile=open(r_file,'r')
	for line in rfile:
		if num==1:
			samples=line.strip().split('\t') #get the header of the table
			num +=1
			continue
		else:
			temp=line.strip().split('\t')
			for i in range(1,len(temp)): #col
				if re.search(',',temp[i]):
					temp[i]=re.sub(',','',temp[i])
				sample_index_type[samples[i]]=temp[i]
			# replace the '-' and space of types into '_' 
			type=string_convert(temp[0])
			sample_index_value[type]=sample_index_type	
			sample_index_type={}
			num +=1
			continue
	rfile.close()
	return sample_index_value
def read_trans_file(r_file):
	sample_value={}
	num=1
	rfile=open(r_file,'r')
	for line in rfile:
		if num==1:
			types=line.strip().split('\t') #get the header of the table
			types=[string_convert(type) for type in types]
			num +=1
			continue
		else:
			temp=line.strip().split('\t')
			sample=temp[0]
			for i in range(1,len(temp)): #col
				if types[i] not in sample_value:
					sample_value[types[i]]={}
				# replace the '-' and space of types into '_'
				if re.search(',',temp[i]):
					temp[i]=re.sub(',','',temp[i])
				sample_value[types[i]][sample]=temp[i]
			num +=1
			continue
	rfile.close()
	#print (sample_value)
	return sample_value

def check_result(indir,qc_type,d_type,samples,result):
	check_table={}
	check_image={}
	check_upload={}
	check_flag="质控通过"
	image_flag=",图片请人工核查"
	upload_flag=''
	qual_noexist=[]
	for qual in qc_type:
		quality_failed=[]
		if not "png" in qual and not "upload" in qual:
			sample_flag=0
			if len(samples)==1 and samples[0]=="0":
				for s in result:
					if qual in result[s]:
						sample_flag=1
						if result[s][qual]["QUALITY"] in ["C","D"]:
							f_info=' {0}:{1} '.format(s,result[s][qual]["SCORE"])
							quality_failed.append(f_info)
					else:
						continue
			else:
				for s in samples:
					if qual in result[s]:
						sample_flag=1
						if result[s][qual]["QUALITY"] in ["C","D"]:
							f_info=' {0}:{1} '.format(s,result[s][qual]["SCORE"])
							quality_failed.append(f_info)
					else:
						continue
			if sample_flag==0:
				quality_failed.append("质控指标在数据库中不存在")
				qual_noexist.append(qual)
			if len(quality_failed)>=1:
				info=quality_failed[0]
				if len(quality_failed)>1:info=';'.join(quality_failed)
				check_table[qual]=["不通过",info]
				check_flag="质控不通过"	
			else:
				check_table[qual]=["通过",""]
		elif "png" in qual:
			png_file=d_type[qual]['Indir']
			file=png_file.replace('Indir',indir)
			if len(samples)==1 and samples[0]=="0":
				files=glob.glob(file)
				if not files:
					logging.error ("{0},no picture was found,please check the Indir cloumn in your template or the result is wrong\n.".format(qual))
					image_flag=", 图片质控不通过"
					check_image[qual]=["1",""]
				else:
					check_image[qual]=["0",files]
		elif "upload" in qual:
			upload_file=d_type[qual]['Indir']
			file=upload_file.replace('Indir',indir)
			if len(samples)==1 and samples[0]=="0":
				files=glob.glob(file)
				if not files:
					logging.error ("{0},no attachment was found,please check the Indir cloumn in your template or the result is wrong\n.".format(qual))
					upload_flag=", 附件存在问题"
					check_upload[qual]=["1",""]
				else:
					check_upload[qual]=["0",files]
	if check_image:
		if check_flag=="质控通过":
			check_flag="数值质控通过{0}".format(image_flag)
		else:
			check_flag="数值质控不通过{0}".format(image_flag)
	if check_upload:
		check_flag="{0}{1}".format(check_flag,upload_flag)

	return check_table,check_image,check_flag,check_upload,qual_noexist

def check_email(table_result,image_result,upload_result,check_flag,project_name,nowtime,indir,qc_type,template,d_type):
	email_subject ='【分析项目质控反馈】【{0}】{1}'.format(check_flag,project_name)

	user = getpass.getuser()
	text1 = '各位:<br> ' 
	text1 += '\t\t该项目质控信息汇总如下：</br> ' 
	text1 += '项目名称：{0}<br>'.format(project_name) 
	text1 += '分析人员：{0}<br>'.format(user) 
	text1 += '分析日期：{0}<br>'.format(nowtime) 
	text1 += '质控路径：{0}<br>'.format(indir) 
	text1 += '<font color=\'red\'>' + '质控结果：{0}'.format(check_flag) + '</font></br>'  
								 
	html = '<html><head></head><body> '		
	#table 
	table1 = '<table border="1"><tr> '
	table1 += '<caption><b>具体指标质控情况表</b></caption>'
	table1 += '<tr><th>质控指标</th><th>参考标准</th><th>审核情况</th><th>备注信息</th></tr>'
	
	for qual in qc_type:
		if not "png" in qual and not "upload" in qual:
			if table_result[qual][0]=="不通过":
				table1 += '<tr><td align=\"center\"><font color=\'red\'> {0} </font></td><td align=\"center\"><font color=\'red\'> {1} </font></td><td align=\"center\"><font color=\'red\'> {2} </font></td><td align=\"center\"><font color=\'red\'> {3} </font></td></tr>'.format(d_type[qual]["Description"],d_type[qual]["Standard"],table_result[qual][0],table_result[qual][1])			
			else:
				table1 += '<tr><td align=\"center\">{0}</td><td align=\"center\">{1}</td><td align=\"center\">{2}</td><td align=\"center\">{3}</td></tr>'.format(d_type[qual]["Description"],d_type[qual]["Standard"],table_result[qual][0],table_result[qual][1])
	
	table1 += '</table>'
		
	html += '<p>{0}</p>{1}</br> '.format(text1,table1)
	if image_result:
		for qual in image_result:
			if image_result[qual][0]=="0":
				tag1=qual.lower()
				html += '<p>{0}如下:<br> <img src="cid:{1}" width="400" height="400" border=0.3></p><br>'.format(d_type[qual]["Description"],tag1)
			else:
				html += '<p><font color=\'red\'>{0}如下:<br> 质控失败，图片可能不存在，或路径不对，请人工再次核查 </font></p><br>'.format(d_type[qual]["Description"])
	if upload_result:
		for qual in upload_result:
			if upload_result[qual][0]=="1":
				html+='<p><font color=\'red\'>{0}附件不存在，请确认</font></p><br>'.format(qual)
	html += '</br>'+'+++++++++++++++++++++++++++++++++++ 说明 +++++++++++++++++++++++++++++++++'+'</br>' 
	html += '(1) 附件中zip文件是具体质控数据，其他为重要图表文件；<br>'
	html += '(2) 该产品质控标准文件：{0}；<br>'.format(template)
	html += '(3) 此次质控分为4个等级：ABCD，分别对应合格，让步通过，不通过，终止流程；<br>'
	html += '(4) 质控通过包括A和B，质控失败包括C和D。<br>'
	html += '</body></html>'


	return html,email_subject

def check_log(table_result,image_result,upload_result,check_flag,project_name,nowtime,indir,qc_type,template,d_type,outdir,prefix):
	user = getpass.getuser()
	file='{0}/{1}_check.log'.format(outdir,prefix)
	log_file=open(file,"w")
	log_file.write('''
[info]
项目名称={0}
分析人员={1}
分析日期={2}
质控路径={3}
指控结果={4}
指控指标={5}

[table]
质控指标\t指标描述\t参考标准\t审核情况\t备注信息
'''.format(project_name,user,nowtime,outdir,check_flag,template))
	for qual in qc_type:
		if not "png" in qual and not "upload" in qual:
			log_file.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(qual,d_type[qual]["Description"],d_type[qual]["Standard"],table_result[qual][0],table_result[qual][1]))
	if image_result:
		log_file.write("\n[image]\n质控指标\t指标描述\t图片路径\n")
		for qual in image_result:
			if image_result[qual][0]=="0":
				log_file.write("{0}\t{1}\t{2}\n".format(qual,d_type[qual]["Description"],",".join(image_result[qual][1])))
			else:
				log_file.write("{0}\t{1}\t质控不通过，图片可能不存在，或路径不对\n".format(qual,d_type[qual]["Description"]))
	if upload_result:
		log_file.write("\n[attachments]\n附件指标\t指标描述\t文件路径\n")
		for qual in upload_result:
			if upload_result[qual][0]=="0":
				log_file.write("{0}\t{1}\t{2}\n".format(qual,d_type[qual]["Description"],",".join(upload_result[qual][1])))
			else:
				log_file.write("{0}\t{1}\t{2}\t不存在，或路径不对\n".format(qual,d_type[qual]["Description"],d_type[qual]["Indir"]))
	log_file.close()
				
def judge(indir,qc_type,d_type,samples):
	result={}
	#quals=list(d_type.keys()) #get the quality control types
	quals=qc_type
	for sample in samples:
		if sample not in result:
			result[sample]={}
		for qual in quals:
			if qual not in result[sample]:
				result[sample][qual]={}
			if qual not in d_type:
				logging.error("{0} does not in QC template file,please check it.".format(qual))
				sys.exit()
			if 'Indir' not in d_type[qual]:
				logging.error("{0} {1} column does not exist,please check you template.".format(qual,'Indir'))
				sys.exit()
			dir=d_type[qual]['Indir']
			if re.search('SAMPLE',dir):
				dir=dir.replace('Indir',indir)
				file=dir.replace('SAMPLE',sample)
			else:
				file=dir.replace('Indir',indir)
			files=glob.glob(file)
			if len(files)>1:
				logging.error ("{0},more than one file were found,please check the Indir cloumn in your template.".format(files))
				sys.exit()
			elif len(files)==1:
				i_file=files[0]
				if not os.path.exists(i_file):
					logging.error ("{0} does not exist.".format(i_file))
					sys.exit()
			else:
				logging.error("{0} does not exist,please check your template.".format(file))
				sys.exit()
			result[sample][qual]['DIR']=i_file
			result[sample][qual]['SCORE']=''
			result[sample][qual]['QUALITY']=''
			if 'Header' not in d_type[qual]:
				logging.error("{0} {1} column does not exist,please check you template.".format(qual,'Header'))
				sys.exit()
			if d_type[qual]['Header']=='R':
				d_sample=read_file(i_file)
			elif d_type[qual]['Header']=='C':
				d_sample=read_trans_file(i_file)
			else:
				logging.error("{0}_{1} should be \'R\' or \'C\',please check you template".format(qual,'Header'))
				sys.exit()
			if qual not in d_sample:
				logging.error ("QC type {0} does not exist in your QC file,please check it.".format(qual))
				sys.exit()
			if sample not in d_sample[qual]: 
				logging.error ("Sample name {0} does not exist in your QC file,please check your sample name.".format(sample))
				sys.exit()
			for rank in d_type[qual]:  
				value=d_type[qual][rank]
				if pat2.search(value): #Number judge
					if pat3.search(value): #Multi criteria judge
						judge1=value.split(';')
						for m in judge1:
							judge1_min,judge1_max=m.lstrip('(').rstrip(')').split(',')
							if float(d_sample[qual][sample])>=float(judge1_min) and float(d_sample[qual][sample])<float(judge1_max):
								result[sample][qual]['SCORE']=d_sample[qual][sample]
								result[sample][qual]['QUALITY']=rank
							else:
								continue
					else: #Single criteria judge
						#print value
						judge2_min,judge2_max=value.lstrip('(').rstrip(')').split(',')
						if float(d_sample[qual][sample])>=float(judge2_min) and float(d_sample[qual][sample])<float(judge2_max):
							result[sample][qual]['SCORE']=d_sample[qual][sample]
							result[sample][qual]['QUALITY']=rank
						else:
							continue
				else: #String judge
					judge3=value
					if d_sample[qual][sample].upper()==judge3.upper():
						result[sample][qual]['SCORE']=d_sample[qual][sample]
						result[sample][qual]['QUALITY']=rank
					else:
						continue
			if result[sample][qual]['SCORE']=='':
				logging.error ("The {0}_{1} value {2} not in any rank of {3},please check you template.".format(sample,qual,d_sample[qual][sample],qual))
	#print (result)
	return (result)	

def push_all_db(all_dict,db_qc,prefix,tables,nowtime):
	heads=sorted( list(list(all_dict.values())[0].keys()) )
	db_qc.check_table_exists(tables,heads)
	log_header=[]
	log_header=db_qc.check_field("LOG")
	for sample in sorted(all_dict):
		for index in sorted(all_dict[sample]):
			quality_value={}
			score_value={}
			dir_value={}
			quality_value[index]=all_dict[sample][index]["QUALITY"]
			score_value[index]=all_dict[sample][index]["SCORE"]
			dir_value[index]=all_dict[sample][index]["DIR"]
			for table in tables:  #check the summary table
				column=db_qc.check_field(table)  #get the column names
				if not index in column: #if the index not in columns,creat it 
					db_qc.alter_field(table,index)
				t_sample_list=db_qc.query_sample(table)
				if not sample in t_sample_list:
					if "QUALITY" in table:
						db_qc.insert(table, ["Sample",index],[sample,all_dict[sample][index]["QUALITY"]])
					elif "SCORE" in table:
						db_qc.insert(table, ["Sample",index],[sample,all_dict[sample][index]["SCORE"]])
					elif "DIR" in table:
						db_qc.insert(table, ["Sample",index],[sample,all_dict[sample][index]["DIR"]])
				else:
					old_value=db_qc.query_value(table , sample, index)
					if "QUALITY" in table:
						if old_value != quality_value[index]:
							db_qc.update(table,quality_value,"Sample",sample)
							if table=="QUALITY" and old_value !=None:
								logging.warning("{0} {1} quality is different from the last time, last time is {2}, this time is {3} ".format(sample,index,old_value,quality_value[index]))
							db_qc.update(table,quality_value,"Sample",sample)
					elif "SCORE" in table:
						if old_value != score_value[index]:
							db_qc.update(table,score_value,"Sample",sample)
							if table=="SCORE" and old_value !=None:
								logging.warning("{0} {1} score is different from the last time, last time is {2}, this time is {3} ".format(sample,index,old_value,score_value[index]))
					elif "DIR" in table:
						if old_value != dir_value[index]:
							db_qc.update(table,dir_value,"Sample",sample)
							if table=="DIR" and old_value !=None:
								logging.warning("{0} {1} dir is different from the last time, last time is {2}, this time is {3} ".format(sample,index,old_value,dir_value[index]))
			db_qc.insert("LOG",log_header,[prefix,index,sample,all_dict[sample][index]["SCORE"],all_dict[sample][index]["QUALITY"],nowtime,all_dict[sample][index]["DIR"]])		

def read_all_db(db_qc,tag):
	result={}
	new=db_qc.check_max_table()
	for i in [ 'QUALITY','SCORE',"DIR"]:
		if tag=="all":
			name="{0}_{1}".format(i,new)
		else:
			name=i
		all_values=db_qc.query_all(name)
		column=db_qc.check_field(name)
		for a in all_values:
			ind = a[0]
			if not ind in result:
				result[ind]={}
			for j in range(len(a)):
				if j !=0:
					if not column[j] in result[ind]:
						result[ind][column[j]]={}
					if a[j]:
						result[ind][column[j]][i]=a[j]
					else:
						result[ind][column[j]][i]="-"
	return result

def html(result, outdir):
	rank = {}
	samples = sorted( list(result.keys()) )
	quals=sorted( list( ( list(result.values())[0] ).keys() ) )
	
	text = '各位:' + '</br>' 
	text += '<b>质控结果输出目录：</b>' + outdir + '<br>' 
	html = '<html><head></head><body> '
	
	#tables
	table1 = '<table border="1"><tr> '
	table1 += '<caption><b>Table1: QC rank</b> </caption> '
	table1 += '<th>Sample</th> '

	table2 = '<table border="1"><tr> '
	table2 += '<caption><b>Table2: QC value</b> </caption> '
	table2 += '<th>Sample</th> '
	
	#table head
	for qual in quals:
		table1 += '<th>{0}</th> '.format(qual)
		table2 += '<th>{0}</th> '.format(qual)
	table1 += '</tr> '
	table2 += '</tr> '
	
	#table row
	ref ={'A':'合格','B':'让步','C':'不合格','D':'终止','-':'-','Standard':'合格'}
	for sample in sorted(result):
		outcome = ''
		text += '<b>' + sample + '</b>' + ' '
		table1 += '<tr> '
		table2 += '<tr> '
		table1 += '<td align=\"center\">{0}</td> '.format(sample)
		table2 += '<td align=\"center\">{0}</td> '.format(sample)
		for qual in sorted(result[sample]):
			quality = result[sample][qual]['QUALITY']
			score = result[sample][qual]['SCORE']
			if re.search('[B|C|D]',quality):
				outcome = 'Failed'
				text += '<font color=\'red\'>' + quality + '</font>' + ' '
				table1 += '<td align=\"center\" bgcolor=\'red\' >{0}</td> '.format(ref[quality])
				table2 += '<td align=\'center\' bgcolor=\'red\'>{0}</td> '.format(score)
			else:
				text += quality + ' '
				table1 += '<td align=\"center\">{0}</td> '.format(ref[quality])
				table2 += '<td align=\"center\">{0}</td> '.format(score)
			
			#stat the total rank
			if quality not in rank:
				rank[quality]=1
			else:
				rank[quality]+=1
			
		
		text += '<font color=\'red\'>' + outcome + '</font>' + '\n'	
		table1 += '</tr> '
		table2 += '</tr> '
	table1 += '</table>'
	table2 += '</table>'
	
	html += '<p>{0}</p>'.format(text)
	html += '{0} </br></br> {1} </body></html>'.format(table1,table2)
	return (html,rank)

def email(result, projectname, outdir, prefix, flag):
	qc_email = 0
	email_subject = ''
	
	#get the email_body and rank 
	(email_body, rank) = html(result,outdir)
	
	#get the email subject	
	if flag == 'yes':
		if 'D' in rank:
			email_subject = "【分析项目质控点反馈】{0}+{1}+{2}".format(projectname,prefix,'终止')
		elif 'C' in rank:
			email_subject = "【分析项目质控点反馈】{0}+{1}+{2}".format(projectname,prefix,'不合格')
		elif 'B' in rank:
			email_subject = "【分析项目质控点反馈】{0}+{1}+{2}".format(projectname,prefix,'让步')
		elif 'A' in rank:
			email_subject = "【分析项目质控点反馈】{0}+{1}+{2}".format(projectname,prefix,'合格')
		qc_email = 1
		
	elif re.search('[B|C]',flag):
		if flag == 'B':
			if 'D' in rank:
				email_subject = "【分析项目质控点反馈】{0}+{1}+{2}".format(projectname, prefix, '终止')
				qc_email = 1
			elif 'C' in rank:
				email_subject = "【分析项目质控点反馈】{0}+{1}+{2}".format(projectname, prefix, '不合格')
				qc_email = 1
			elif 'B' in rank:
				email_subject = "【分析项目质控点反馈】{0}+{1}+{2}".format(projectname, prefix, '让步')
				qc_email = 1
			else:
				logging.info("The qc type rank of all samples are A,don't send email.")
		elif flag == 'C':
			if 'D' in rank:
				email_subject = "【分析项目质控点反馈】{0}+{1}+{2}".format(projectname, prefix, '终止')
				qc_email = 1
			elif 'C' in rank:
				email_subject = "【分析项目质控点反馈】{0}+{1}+{2}".format(projectname, prefix, '不合格')
				qc_email = 1
			else:
				logging.info("The qc type rank of all samples are A or B,don't send email.")
	elif flag=='no':
		logging.info("don't send email.")
		sys.exit()
	else:
		logging.error("The parameter -e should be yes/no/[B/C],please check it.")
		sys.exit()
	return qc_email,email_subject, email_body
	
def get_email_config(email_cfg, email_subject, email_attachment, email_body, outdir):		
	
	#readin email config file
	config = configparser.RawConfigParser()
	config.read(email_cfg)
	user = getpass.getuser()
	
	#revise the email_config.ini
	sections=config.sections()
	if "BODY" not in sections:
		config.add_section('BODY')
		config.set('BODY','Subject',email_subject)
		config.set('BODY','Attachment',email_attachment)
		config.set('BODY','Body',email_body)
	else:
		config.set('BODY','Subject',email_subject)
		config.set('BODY','Attachment',email_attachment)
		config.set('BODY','Body',email_body)
			
	#Save the revised config
	config_file = open("{0}/email_config.ini".format(outdir),'w')
	config.write(config_file)
	config_file.close()

def get_email_config_image(email_cfg, image):	
	#readin email config file
	config = configparser.ConfigParser()
	config.read(email_cfg)
	sections=config.sections()
	if "IMAGE" not in sections:
		config.add_section('IMAGE')
		for i in image:
			config.set('IMAGE',i,image[i])
	
	config_file = open(email_cfg,'w')
	config.write(config_file)
	config_file.close()	
		

def jump_email(bindir,email_dir,outdir):
	# upload the email content and send email
	print("/annoroad/data1/software/bin/miniconda/envs/python3_base/bin/python3 {0}/jump_send_email.py -i {1} && cd {2}".format(bindir,email_dir,outdir))
	qc_mail_return = subprocess.call("/annoroad/data1/software/bin/miniconda/envs/python3_base/bin/python3  {0}/jump_send_email.py -i {1} && cd {2}".format(bindir,email_dir,outdir), shell=True)
		
	if qc_mail_return == 0:
		logging.info("如果上一行面没有Reason：比如类似list index out of range 提醒，则send email successfully。如果有报错提醒则是跳转机内发送邮件失败")
	else:
		logging.warning("send email failed")
	#Remove email_dir
	#shutil.rmtree(email_dir)
	if not os.path.exists(email_dir):
		logging.info("Remove email_dir sucessfully")

def send_email(result, projectname, outdir, prefix, flag, email_cfg):

	#get the qc_mail email_subject, email_body
	qc_mail, email_subject, email_body = email(result, projectname, outdir, prefix, flag)
	nowtime = datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
	
	#get the username
	user = getpass.getuser()
	
	if qc_mail:
		#creat email dir
		email_dir = "{0}/{1}_{2}".format(outdir, projectname, nowtime)
		if not os.path.exists(email_dir):
			os.mkdir(email_dir)
		
		#copy the *_Result_*.xls and script to email dir 
		quality_file = "{0}/{1}_Result_QUALITY.xls".format(outdir,prefix)
		score_file = "{0}/{1}_Result_SCORE.xls".format(outdir,prefix)
		shutil.copy(quality_file, email_dir)
		shutil.copy(score_file , email_dir)
		shutil.copy("{0}/myemail.py".format(bindir), email_dir)
		shutil.copy("{0}/send_email.py".format(bindir), email_dir)
		
		#zip files 
		filenames = [quality_file, score_file]
		zipfiles_name = "/{0}/{1}_QC_Result.zip".format(email_dir,prefix)
		zipfiles = zipfile.ZipFile(zipfiles_name,'w',compression=zipfile.ZIP_DEFLATED)
		for file in filenames:
			zipfiles.write(file,os.path.basename(file))
		zipfiles.close()
		
		#get the email attach
		email_attachment = "/home/{0}/{1}/{2}_QC_Result.zip".format(user, os.path.basename(email_dir), prefix)
		
		#generate the email_config.ini
		get_email_config(email_cfg, email_subject, email_attachment, email_body, email_dir)
		
		# upload the email content and send email
		jump_email(bindir,email_dir,outdir)

def send_check_email(table_result,image_result,upload_result,check_flag,projectname,title,outdir,prefix,indir,qc_type,template,email_cfg,d_type):

	nowtime = datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
	nowtime1 = datetime.datetime.now().strftime('%Y-%m-%d:%H:%M:%S')
	email_dir = "{0}/{1}_{2}_check".format(outdir, projectname, nowtime)
	user = getpass.getuser()

	if not os.path.exists(email_dir):
		os.mkdir(email_dir)	
	#copy the *_Result_*.xls and script to email dir
	quality_file = "{0}/{1}_Result_QUALITY.xls".format(outdir,prefix)
	score_file = "{0}/{1}_Result_SCORE.xls".format(outdir,prefix)

	shutil.copy(quality_file, email_dir)
	shutil.copy(score_file , email_dir)
	shutil.copy("{0}/myemail.py".format(bindir), email_dir)
	shutil.copy("{0}/send_email.py".format(bindir), email_dir)

	filenames = [quality_file, score_file]
	zip_attach = "/home/{0}/{1}/{2}_QC_Result.zip".format(user, os.path.basename(email_dir), prefix)
	attach_files =[]
	# copy png to email_dir 
	image_html={}
	for qual in image_result:
		if image_result[qual][0]=="0":
			num=0
			for f in image_result[qual][1]:
				num+=1
				name=os.path.basename(f)
				new_name='/home/{0}/{1}/{2}'.format(user,os.path.basename(email_dir),name)
				if num==1:
					image_html[qual]=new_name
					attach_files.append(new_name)
				else:
					filenames.append(f)
				shutil.copy(f, email_dir)
	#zip files 
	zipfiles_name = "/{0}/{1}_QC_Result.zip".format(email_dir,prefix)
	zipfiles = zipfile.ZipFile(zipfiles_name,'w',compression=zipfile.ZIP_DEFLATED)
	for file in filenames:
		zipfiles.write(file,os.path.basename(file))
	zipfiles.close()

	#upload files
	for qual in upload_result:
		if upload_result[qual][0]=="0":
			num=0
			for f in upload_result[qual][1]:
				num+=1
				name=os.path.basename(f)
				new_name='/home/{0}/{1}/{2}'.format(user,os.path.basename(email_dir),name)
				if num==1:
					attach_files.append(new_name)
				else:
					filenames.append(f)
				shutil.copy(f, email_dir)

	zipfiles_name = "/{0}/{1}_QC_Result.zip".format(email_dir,prefix)
	zipfiles = zipfile.ZipFile(zipfiles_name,'w',compression=zipfile.ZIP_DEFLATED)
	for file in filenames:
		zipfiles.write(file,os.path.basename(file))
	zipfiles.close()
	attach_files.append(zip_attach)
	#get the email attach
	email_attachment= ";".join(attach_files)

	#get the email body
	email_body,email_subject = check_email(table_result,image_result,upload_result,check_flag,title,nowtime1,indir,qc_type,template,d_type)
	check_log(table_result,image_result,upload_result,check_flag,projectname,nowtime1,indir,qc_type,template,d_type,outdir,prefix)
	#generate the email_config.ini
	get_email_config(email_cfg, email_subject, email_attachment, email_body, email_dir)
	new_email_cfg = '{0}/email_config.ini'.format(email_dir)
	get_email_config_image(new_email_cfg, image_html)	
	# upload the email content and send email
	jump_email(bindir,email_dir,outdir)

def db_finish(db_qc,output,prefix,time):
	old=db_qc.check_max_table()
	new=old+1
	new_tables=[]
	for i in [ 'QUALITY' ,'SCORE', 'DIR' ]:
		name = i
		if new >1:
			name = "{0}_{1}".format(i,old)
		new_name = "{0}_{1}".format(i,new)
		new_tables.append(new_name)
		db_qc.cptable(new_name,name)

	all_values=read_all_db(db_qc,"qual")
	push_all_db(all_values,db_qc,prefix,new_tables,time)
	for i in [ 'QUALITY' ,'SCORE', 'DIR' ]:
		name = i
		if new >=1:
			new_name = "{0}_{1}".format(i,new)
			column=db_qc.check_field(new_name)
			outfile= open("{0}_Result_{1}.xls".format(output,i),'w')
			outfile.write("{0}\n".format("\t".join(column)))
			db_qc.export_all(outfile,new_name)
			outfile.close()


def main():
	parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-i','--indir',help='Indir',dest='indir',required=True)
	parser.add_argument('-t','--template',help='QC index table',dest='template',required=True)
	parser.add_argument('-p','--project',help='Project id',dest='project',required=True)
	parser.add_argument('-n','--name',help='QC module name',dest='name',required=True)
	parser.add_argument('-s','--sample',help='Sample name',dest='sample',required=True)
	parser.add_argument('-e','--email',help='Email',dest='email',required=True)
	parser.add_argument('-m','--emailconfig',help='Email config',dest='emailconfig',required=True)
	parser.add_argument('-d','--dbdir',help='Database outdir',dest='dbdir',required=True)
	parser.add_argument('-f','--finish',help='Project finish,export all results',dest='finish',action='store_true')
	parser.add_argument('-c','--check',help='check',dest='check',action='store_true')
	parser.add_argument('-o','--output',help='outdir+prefix',dest='output',required=True)
	parser.add_argument('-pn','--projectname',help='Project name',dest='title',default=' ')
	args=parser.parse_args()
	
	result={}
	samples=[]
	qc_type=[]
	
	#set the logging
	logging.basicConfig(level=logging.DEBUG,format="%(asctime)s - %(filename)s - %(levelname)s - %(message)s")
	
	# Get the local time
	time=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
	
	#get the outdir and prefix
	output=os.path.abspath(args.output)
	outdir=os.path.dirname(output)
	prefix=os.path.basename(output)
	
	projectname=args.project
	title=args.title
	title=title.replace("结题报告","")
	title='{0}{1}'.format(projectname,title)
	indir=os.path.abspath(args.indir)
	dbdir=os.path.abspath(args.dbdir)
	#check whether the dir exists
	check_exists(indir)
	check_exists(dbdir)
	#check_exists(outdir)
	os.system('mkdir -p {0}'.format(outdir))
	
	#export all info in database
	if args.finish:
		db_qc=Database( '{0}/{1}.db'.format(dbdir,projectname))
		db_finish(db_qc,output,prefix,time)
		result=read_all_db(db_qc,"all")
		if not args.check:
			email_cfg = os.path.abspath(args.emailconfig)
			send_email(result, projectname, outdir, prefix, args.email, email_cfg)
			logging.info('Total QC data print out; -f exist, so print all  result in db and exit 0\n')
			sys.exit(0)
			
	#get the sample names
	samples=check_input(args.sample)
	logging.info('get the QC sample name.')
	
	#get the Quality types
	types=check_input(args.name)
	qc_type=[string_convert(type) for type in types]
	logging.info('get the QC index.')
	
	#readin template
	template=os.path.abspath(args.template)
	d_type=read_template(template) 
	logging.info('readin the QC template.')
	
	db_qc=Database( '{0}/{1}.db'.format(dbdir,projectname))

	##用于一些指控指标不存在，其他需要输出判断
	qual_noexist=[]
	if not args.check:
		# Quality control
		logging.info('QC indexs judge start.')
		result=judge(indir,qc_type,d_type,samples)
		logging.info('QC indexs judge end.')	

		#push the Quality control into database
		logging.info('save the QC result into database.')
		#print(result)
		tables = [ 'QUALITY' ,'SCORE', 'DIR' ]
		push_all_db(result,db_qc,prefix,tables,time)

	else:
		logging.info('QC indexs check start.')
		result = read_all_db(db_qc,"all")
		table_result,image_result,check_flag,upload_result,qual_noexist = check_result(indir,qc_type,d_type,samples,result)
	 
	#export the results
	logging.info('export the QC result.')
	tn=db_qc.check_max_table()
	column=["Sample"]
	for i in qc_type:
		if not "png" in i and not "upload" in i and not i in qual_noexist: 
			  column.append(i)
	if len(samples)==1 and samples[0]=="0":
		samples=[]
		for s in result:
			samples.append(s)
	for i in [ 'QUALITY' ,'SCORE', 'DIR' ]:
		tab=i
		name = '{0}_Result_{1}'.format(prefix, i)
		if args.check and tn >=1:
			tab = "{0}_{1}".format(tab,tn)			
		outfile1=open("{0}/{1}.xls".format(outdir,name),'w')
		outfile1.write("{0}\n".format("\t".join(column)))
		db_qc.export(outfile1,tab,"Sample",samples,column)
		outfile1.close()
		
	#send email
	email_cfg = os.path.abspath(args.emailconfig)
	if not args.check:
		flag=args.email		
		send_email(result, projectname, outdir, prefix, flag, email_cfg)
	else:
		send_check_email(table_result,image_result,upload_result,check_flag,projectname,title,outdir,prefix,indir,qc_type,template,email_cfg,d_type)
	
	logging.info('All the analysis finished.\n')
	
				
if __name__=="__main__":
	main()

