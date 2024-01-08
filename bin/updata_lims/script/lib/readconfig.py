#! /usr/bin/env python3
import configparser

__author__='Tu Chengfang'
__mail__= 'chengfangtu@genome.cn'

class myconf(configparser.ConfigParser):
	def __init__(self,defaults=None):
		configparser.ConfigParser.__init__(self,defaults=None,allow_no_value=True)
	def optionxform(self, optionstr):
		return optionstr


class Config_Parser( ) :
	def __init__( self, config_file ) :
		self.config_file = config_file 
		self.config = self.read_config()

	def read_config( self ) :
		config = myconf()
		config.read( self.config_file )
		return config

	def data_block( self ) :
		data = self.config['DATA']
		for rn in data :
			info = data[rn].rstrip().split('\t')
			yield info

	def all_block( self, title, head ) :
		return self.config[title][head]

	def return_block_list( self, title ) :
		try :
			data = self.config[ title]
		except :
			return []
		for rn in data :
			info = data[rn].rstrip().split('\t')
			yield info

	def return_head_value( self, title, head ) :
		return self.config[title][head]

	def return_title_block(self, title ) :
		try :
			data = self.config[ title]
		except :
			data = []
		return data
