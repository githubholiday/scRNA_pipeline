import numpy as np 
import pylab 
import math
import scipy.stats as stats
import pandas as pd
import sys
from pandas import Series,DataFrame #此两个为pandas中的数据结构，特别是DataFrame 非常实用
import glob,os
import plotly.offline as py
from plotly.graph_objs import Scatter, Layout
import pprint as pp
import calendar
import plotly.graph_objs as go 
import colorlover as cl


def qqplot(series , name):
	tt = series[pd.notna(series )]
	res = stats.probplot(tt, dist="norm", plot=pylab)
	pylab.title('qqplot of {0}'.format(name))
	pylab.show()
	slope, intercept, r = res[1]
	r_square = r * r
	return(r_square)
	##https://www.programcreek.com/python/example/1270/pylab.show pylab example

class myPlot():
	def __init__(self, df , xlabel='' , ylabel='' , title='' , color_by = '' ) :
		''' 以index为x轴，column name为y轴进行绘图，数据转换有fromat完成'''
		self.df = df
		self.xlabel = xlabel
		self.ylabel = ylabel
		self.title = title
		self.color_by = pd.DataFrame(df[color_by]) if color_by else ''
		self.color_name = color_by 

	def redef_x(self , x , xlabel ):
		if x : 
			self.xlabel = x
			self.x = self.df[x]
		else:
			self.xlabel = 'all_index'
			self.x = self.df.index
		if xlabel: self.xlabel = xlabel

	def redef_y(self , y , ylabel):
		if y and not isinstance(y , list) : sys.exit('y must be list')
		if y : 
			self.ylabel = "&".join(y)
			self.y = y
		else:
			self.y = set(self.df.columns.get_level_values(0))
			self.ylabel = 'all_column'
		if ylabel: self.ylabel = ylabel

	def choose_col(self):
		category = self.color_by.iloc[: , 0].unique()
		#print(category)
		#print(len(category))
		if len(category) > 11:
			col_list = cl.interp(cl.scales['11']['div']['RdYlBu'], 500)
		elif len(category)>2:
			col_list = cl.scales[len(category)]['div']['RdYlBu']
		elif len(category) == 2:
			col_list = ['rgb(252,141,89)', 'rgb(145,191,219)']
		else:
			col_list = ['rgb(252,141,89)']
		category_dict = { j:col_list[i] for i,j in enumerate(category)}
		category_flag = { j : 0 for j in category}  ## use to store whether used 
		return category_dict,category_flag

	def set_layout(self , xtype):
		x_axis_template=dict(
			showgrid=False,  #网格
			zeroline=False,  #是否显示基线,即沿着(0,0)画出x轴和y轴
			nticks=20,
			showline=True,
			title=self.xlabel,
			mirror='all',
		)
		if xtype: x_axis_template['type'] = xtype
		
		y_axis_template=dict(
			showgrid=False,  #网格
			zeroline=False,  #是否显示基线,即沿着(0,0)画出x轴和y轴
			nticks=20,
			showline=True,
			title=self.ylabel, #y轴标题
			mirror='all'
		)
		

		layout=go.Layout(
			title=self.title,
			#margin=go.layout.Margin(l=400,r=400,b=100,t=150,pad=0),
			xaxis=x_axis_template,
			yaxis=y_axis_template
			)
		return (layout)

	def line_plot(self , xlabel='' , ylabel='' , title='' , x = '' , xtype='' , y='', filename=''):
		myPlot.redef_x(self, x, xlabel)
		myPlot.redef_y(self , y ,ylabel)
		if title : self.title  = title

		traces = []
		if isinstance( self.color_by , pd.core.frame.DataFrame): 
			category_dict , category_flag = myPlot.choose_col(self)

		for j in self.y: 
			#print(df.index)
			if isinstance( self.color_by ,pd.core.frame.DataFrame):
				for group_name in category_dict :
					#group_name  = self.color_by.iloc[:,0][j]
					group_color = category_dict[ group_name ]
					index_list = self.color_by[self.color_name] == group_name
					if category_flag[group_name] : 
						trace_tmp=go.Scatter( 
							name= group_name ,
							x=self.x.loc[index_list] ,
							y=self.df[j].loc[index_list],
							legendgroup=group_name,
							showlegend=False,
							line = dict(color= group_color),
						)
					else:
						trace_tmp=go.Scatter( 
							name= group_name ,
							x=self.x.loc[index_list],
							y=self.df[j].loc[index_list],
							legendgroup=group_name,
							showlegend=True,
							line = dict(color= group_color),
						)
						category_flag[group_name] += 1 
					traces+=[trace_tmp]
			else:
				trace_tmp=go.Scatter( 
					name=j, 
					x=self.x,
					y=self.df[j],
				)
				traces+=[trace_tmp]

		layout = myPlot.set_layout(self, xtype)
		fig=go.Figure( data=traces, layout=layout)
		if filename : 
			py.plot(fig, filename=filename, auto_open= False)
		else:
			py.init_notebook_mode(connected=True)
			py.iplot(fig )

	def scatter_plot(self , xlabel='' , ylabel='' , title='' , x = '' , xtype='' , y='' , alpha=1 , filename='' ):
		myPlot.redef_x(self, x, xlabel)
		myPlot.redef_y(self , y ,ylabel)
		if title : self.title  = title

		traces = []
		if isinstance( self.color_by , pd.core.frame.DataFrame): 
			category_dict , category_flag = myPlot.choose_col(self)

		for j in self.y: 
			#print(df.index)
			if isinstance( self.color_by ,pd.core.frame.DataFrame):
				for group_name in category_dict :
					#group_name  = self.color_by.iloc[:,0][j]
					group_color = category_dict[ group_name ] 
					index_list = self.color_by[self.color_name] == group_name
					if category_flag[group_name] : 
						trace_tmp=go.Scatter( 
							name= group_name ,
							x=self.x.loc[index_list],
							y=self.df[j].loc[index_list],
							legendgroup=group_name,
							showlegend=False,
							mode = 'markers',
							marker = dict(color= group_color,opacity= alpha),
						)
					else:
						trace_tmp=go.Scatter( 
							name= group_name ,
							x=self.x.loc[index_list],
							y=self.df[j].loc[index_list],
							legendgroup=group_name,
							showlegend=True,
							mode = 'markers',
							marker = dict(color= group_color ,opacity= alpha),
						)
						category_flag[group_name] += 1 
					traces+=[trace_tmp]
			else:
				trace_tmp=go.Scatter( 
					name=j, 
					x=self.x,
					y=self.df[j],
					mode = 'markers',
					marker = dict(opacity= alpha),
				)
				traces+=[trace_tmp]
		layout = myPlot.set_layout(self, xtype)
		fig=go.Figure( data=traces, layout=layout)
		if filename : 
			py.plot(fig, filename=filename, auto_open= False)
		else:
			py.init_notebook_mode(connected=True)
			py.iplot(fig )
	def heat_plot(self , xlabel='' , ylabel='' , title='' , x = '' , xtype='' , y='' , alpha=1 ,color='Hot' ,filename=''):
		''' Greys,YlGnBu,Greens,YlOrRd,Bluered,RdBu,Reds,Blues,Picnic,Rainbow,Portland,Jet,Hot,Blackbody,Earth,Electric,Viridis,Cividis. '''
		
		fig = [go.Heatmap( z=self.df.values.tolist(), 
		                    x = self.df.index,
							y = self.df.columns,
		                    colorscale=color)]
		if filename : 
			py.plot(fig, filename=filename, auto_open= False)
		else:
			py.init_notebook_mode(connected=True)
			py.iplot(fig )
	def hist2D_plot(self , xlabel='' , ylabel='' , title='' , x = '' , xtype='' , y='' , alpha=1 ,color='Hot' ,filename=''):
		myPlot.redef_x(self, x, xlabel)
		myPlot.redef_y(self , y ,ylabel)
		if title : self.title  = title
		fig = [go.Histogram2d( x = self.x,
							  y = self.df[ self.y[0] ],
							  #autobinx=False,
							  #xbins=dict(start=0, end=1000, size=10),
							  #autobiny=False,
							  #ybins=dict(start=0, end=1000, size=10),
							  #colorscale=[[0, 'rgb(12,51,131)'], [0.25, 'rgb(10,136,186)'], [0.5, 'rgb(242,211,56)'], [0.75, 'rgb(242,143,56)'], [1, 'rgb(217,30,30)']]
							  )]
		if filename : 
			py.plot(fig, filename=filename, auto_open= False)
		else:
			py.init_notebook_mode(connected=True)
			py.iplot(fig )

	def box_plot(self , xlabel='' , ylabel='' , title='' , x = '' , xtype='' , y='' , alpha=1,filename=''):
		myPlot.redef_x(self, x, xlabel)
		myPlot.redef_y(self , y ,ylabel)
		if title : self.title  = title
		traces = []
		if isinstance( self.color_by , pd.core.frame.DataFrame): 
			category_dict , category_flag = myPlot.choose_col(self)
		for j in self.y: 
			#print(df.index)
			if isinstance( self.color_by ,pd.core.frame.DataFrame):
				print('Error: donot support group')
				pass
			else:
				trace_tmp=go.Box( 
					name=j, 
					#x=x,
					y=self.df[j],
					#mode = 'markers',
					#marker = dict(opacity= alpha),
				)
			traces+=[trace_tmp]
		layout = myPlot.set_layout(self, xtype)
		fig=go.Figure( data=traces, layout=layout)
		if filename : 
			py.plot(fig, filename=filename, auto_open= False)
		else:
			py.init_notebook_mode(connected=True)
			py.iplot(fig )
	def pie_plot(self, xlabel='' , ylabel='' , title='' , x = '' , xtype='' , y='' , alpha=1,filename=''):
		myPlot.redef_x(self, x, xlabel)
		myPlot.redef_y(self , y ,ylabel)
		if title : self.title  = title
		traces = []
		if isinstance( self.color_by , pd.core.frame.DataFrame): 
			category_dict , category_flag = myPlot.choose_col(self)
		if len(y) != 1 : 
			sys.exit('y should be in length 1')
		else : 
			#print(df.index)
			if isinstance( self.color_by ,pd.core.frame.DataFrame):
				print('Error: donot support group')
				pass
			else:
				trace_tmp=go.Pie( 
					labels=self.df.index, 
					values=self.df[y[0]],
				)
			traces+=[trace_tmp]
		layout = myPlot.set_layout(self, xtype)
		fig=go.Figure( data=traces, layout=layout)
		if filename : 
			py.plot(fig, filename=filename, auto_open= False)
		else:
			py.init_notebook_mode(connected=True)
			py.iplot(fig )



class myDF():
	def __init__(self , df):
		self.df = df
	
	def normalization(self):
		df_norm = (self.df - self.df.mean()) / (self.df.max() - self.df.min())
		return (df_norm)
