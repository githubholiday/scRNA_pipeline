#! /usr/bin/env python3
'''
Created on 20190911

程序功能:
绘制差异基因火山图动态图，以log2FoldChange值为横坐标，-lg(padj)为纵坐标，up/down信息区分颜色。

程序参数:
-i --input：  【必选】输入文件，包括4列，分别为gene，log2FoldChange，-lg(padj)，up/down信息
-o --outdir： 【必选】输出文件，以html结尾，路径不存在时进行构建目录
-x --datax:   【可选】输入参数，点图的横坐标，对应输入文件的第二列表头，用于自定义绘制点图
-y --datay:   【可选】输入参数，点图的纵坐标，对应输入文件的第三列表头，用于自定义绘制点图
-t --title:   【可选】输入参数，点图的标题，用于自定义绘制点图

使用方法：
/annoroad/share/software/install/miniconda3/bin/python3 volcano_html.py -i A_B.volcano.tmp -o outdir/out.html
'''
import argparse
import sys
import re
import pandas as pd
import os
import logging
bindir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(bindir)
import plotly.graph_objs as go
import plotly.express as px
import plotly.offline as py

__author__='chen pengyan'
__mail__= 'pengyanchen@genome.cn'
__date__= '2019/9/11'

pat1=re.compile('^\s+$')

LOG = os.path.basename(__file__)

def my_log( level, message ) :
	logging.basicConfig(level = logging.INFO,format = '%(asctime)s - %(filename)s - %(levelname)s - %(message)s')
	logger = logging.getLogger(__name__)
	if level == 'info' :
		return logger.info( message )
	if level == 'warning' :
		return logger.warning( message )
	if level == 'debug' :
		return logger.debug( message )
	if level == 'error' :
		return logger.error( message )

def check_file_exists( *file_list ) :
	for file in file_list :
		if os.path.exists( file ) :
			my_log( 'info', 'file : {0}'.format( file ) )
		else :
			my_log( 'error', 'file is not exists : {0}'.format( file ) )
			sys.exit(1)

def make_dir( dir ) :
	try :
		os.makedirs( dir )
		time.sleep(1)
		my_log( 'info', 'mkdir {0} sucessful!'.format( dir) )
	except :
		my_log( 'error', 'mkdir {0} failed!'.format( dir) )

def scatter(data_frame, x=None, y=None, color=None, symbol=None, size=None,
			hover_name=None, hover_data=None, text=None, facet_row=None,
			facet_col=None, error_x=None, error_x_minus=None, error_y=None,
			error_y_minus=None, animation_frame=None, animation_group=None, 
			category_orders={}, labels={}, color_discrete_sequence=None, 
			color_discrete_map={}, color_continuous_scale=None, 
			range_color=None, color_continuous_midpoint=None, 
			symbol_sequence=None, symbol_map={}, opacity=None,
			size_max=None, marginal_x=None, marginal_y=None, trendline=None,
			trendline_color_override=None, log_x=False, log_y=False,
			range_x=None, range_y=None, render_mode='auto', title=None,
			template=None, width=None, height=None):
	return

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--input',help='input file',dest='input',required=True)
	parser.add_argument('-o','--outfile',help='output file',dest='outfile',required=True)
	parser.add_argument('-x','--datax',help='x data',dest='datax')
	parser.add_argument('-y','--datay',help='y data, if more than 1, separated by ","',dest='datay')
	parser.add_argument('-t','--title',help='title',dest='title',default='Line graph')
	args=parser.parse_args()

	check_file_exists( args.input )
	outdir = os.path.dirname( args.outfile )
	if not os.path.exists( outdir ) :
		my_log( 'info', '{0}  directory is not exists, i creat it'.format(outdir) )
		make_dir( outdir )
	matrix = pd.read_table(args.input, header=0, index_col=None, encoding='utf-8')
	fig=px.scatter(matrix, x="log2FoldChange", y="-lg(padj)", color="up/down",
				color_discrete_map={"up":"yellow","down":"blue","none":"grey"},
				render_mode='auto',width=600, height=600)
	#layout = go.Layout(margin=go.layout.Margin(l=400,r=400,b=200,t=200,pad=0))
	filename=args.outfile
	if filename :
		py.plot(fig, filename=filename, auto_open= False)
	else:
		py.init_notebook_mode(connected=True)
		py.iplot(fig )
	my_log( 'info', '{0} is generated'.format(filename) )
if __name__ == '__main__':
	main()
