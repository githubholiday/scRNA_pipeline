'''
连接lims数据库，从流程的config文件中获取到Lims数据库的连接信息，进行连接
执行单个功能的时候，按照table_name进行即可
'''
import os
import sys
import argparse
import mysql.connector
import re
pat1=re.compile('^\s+$')

__author__ = 'suyanxun'
__mail__ = 'yanxunsu@genome.cn'
pat1=re.compile('^\s+$')

class LIMS():
	def __init__(self, config_file, port = 3306, charset='utf8' ):
		self.config_file = config_file
		self.config_dic = self.read_config( )
		usr  = self.config_dic['sql_usr']
		pwd = self.config_dic['sql_pwd']
		port = self.config_dic['sql_port']
		host = self.config_dic['sql_host']
		database = self.config_dic['sql_db']
		try:
			self.cnx = mysql.connector.connect(user=usr, password=pwd, host=host, database=database, port = port, charset = charset)
			print( 'connect db successed' )
			self.cursor = self.cnx.cursor()
		except mysql.connector.Error as err:
			print ( 'connect db failed' )
			print ( 'Error: {0}'.format( err ) )
			sys.exit()
		#self.table = table
		self.charset = charset
	def read_config(self):
		config_dic = {}
		with open(self.config_file, 'rt', encoding='utf-8') as infile:
			for line in infile.readlines():
				if re.search(pat1, line ) or line.startswith('#') : continue
				tmp = line.rstrip().split('=',1)
				target = tmp[0].rstrip(' ')
				value = tmp[1].lstrip(' ')
				if target not in config_dic :
					config_dic[ target ] = value
				else :
					print("{0} is repeat in {1}".format(target, self.config_file))
		return config_dic
	
	def colname( self, table ):
		cmd = 'describe {0}'.format( table )
		self.execute( cmd )
		name_list = [i[0] for i in self.cursor.fetchall() ]
		return name_list
	
	def select( self, table_name, col_list = '*', conditions = None ):
		'''
		默认返回所有值
		col_list: 查询的列名,list
		conditions: 条件,形如[ (colname1,value1), (colname2,value2), ... ]
		'''
		table = self.config_dic[table_name]
		if col_list == '*':
			col = col_list
		else:
			col = ','.join( col_list )
		cmd = 'SELECT {0} from {1}'.format( col, table )
		if conditions:
			cmd = cmd + ' where {0} = "{1}"'.format( conditions[0][0], conditions[0][1] )
			if len( conditions ) > 1:
				for n,v in conditions[1:]:
					cmd = cmd + ' and {0} = "{1}"'.format( n,v )
		self.execute( cmd )

		return self.cursor.fetchall()
	
	def insert( self, table_name, col_list, value_list ):
		'''
		col_list: 要插入的列名list
		value_list: 要插入的值,形如[(value1.1,value1.2,...), (value2.1,value2.2,...)]
		'''
		table = self.config_dic[table_name]
		cmd = 'INSERT INTO {0}({1}) VALUES ({2})'.format( table, ','.join( col_list ), ','.join( [str(i) for i in value_list[0]] ) )
		if len( value_list ) > 1:
			for i in value_list[1:]:
				cmd = cmd + ',({0})'.format( ','.join([str(j) for j in i]) )
		self.execute( cmd )
	
	def update( self, table_name, value_list, conditions = None ):
		'''
		value_list: 要更新的值,形如[ (colname1,value1), (colname2,value2), ... ]
		conditions: 更新条件: 形如[ (colname1,value1), (colname2,value2), ... ]
		'''
		table = self.config_dic[table_name]
		cmd = 'UPDATE {0} SET '.format( table )
		for n,v in value_list:
			cmd = cmd + '{0}="{1}", '.format( n,v )
		cmd = cmd.rstrip( ', ' )
		if conditions:
			cmd = cmd + ' where {0} = "{1}"'.format( conditions[0][0], conditions[0][1] )
			if len( conditions ) > 1:
				for n,v in conditions[1:]:
					cmd = cmd + ' and {0} = "{1}"'.format( n,v )
		self.execute( cmd )
	
	def delete( self, table_name, conditions ):
		'''
		conditions: 删除条件: 形如[ (colname1,value1), (colname2,value2), ... ]
		'''
		table = self.config_dic[table_name]
		cmd = 'delete from {0} where {1}="{2}"'.format( table, conditions[0][0], conditions[0][1] )
		if len( conditions ) > 1:
			for n,v in conditions[1:]:
				cmd = cmd + ' and {0} = "{1}"'.format( n,v )
		self.execute( cmd )
	
	def execute( self, cmd, times = 3 ):
		if times > 0:
			try:
				self.cursor.execute(cmd)
				if cmd.startswith( ('INSERT' ,'UPDATE' , 'DELETE')):
					self.cnx.commit()
				#print ( '{0} 运行成功'.format(cmd) )
			except mysql.connector.Error as err:
				print ('{0} 尝试倒数第{1}次失败'.format( cmd, times ) )
				print ( 'Error: {0}'.format( err ) )
				self.execute( cmd, times-1 )
		else:
			self.close()
			sys.exit()
	
	def close( self ):
		self.cnx.commit()
		self.cursor.close()
		self.cnx.close()
