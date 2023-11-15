#! /usr/bin/env python3
import argparse
import os
import sys
import pandas as pd
import plotly as py
import plotly.graph_objs as go
#import plotly.express as px
pyplt = py.offline.plot


def main():
	parser=argparse.ArgumentParser(description='Bar plot using plotly')
	parser.add_argument('-i','--input',help='input file or stdin',dest='input',required=True)
	parser.add_argument('-o','--outfile',help='output file',dest='outfile',required=True)
	parser.add_argument('-t','--title',help='title',dest='title',default='DE Gene Count')
	parser.add_argument('-x','--xaxis',help='xaxis',dest='xaxis',default='')
	parser.add_argument('-y','--yaxis',help='yaxis',dest='yaxis',default='Count')
	parser.add_argument('-m','--barmode',help='barmode',dest='barmode',choices=['group','stack'])
	parser.add_argument('-W','--width',help='width',dest='width',default='600')
	parser.add_argument('-H','--height',help='height',dest='height',default='600')
	args=parser.parse_args()
	
	# handle
	if os.path.isfile(args.input):
		handle = args.input
	else:
		handle = sys.stdin
	
	# read table
	df = pd.read_csv(handle, header=0, index_col=None, sep='\t', encoding='utf-8')
	
	# bar traces
	sample_head = list(df.columns[1:])
	item = list( df.iloc[:,0] )
	traces = []
	for i in range( len(item) ):
		print(i)
		if i <= 1:
			col_list = list(df.iloc[i,:][1:])
			up_color = []
			down_color = []
			if i == 0:
				for j in col_list:
					up_color.append("rgb(255, 106, 106)")
				upcount_list = list(df.iloc[i,:][1:])
				trace = go.Bar( x = sample_head, y = list(df.iloc[i,:][1:]), name = item[i] ,opacity=0.8, marker=dict(color=up_color) )
				traces.append( trace )
			elif i == 1:
				for j in col_list:
					down_color.append("rgb(0, 238, 238)")
				downcount_list = list(df.iloc[i,:][1:])
				trace = go.Bar( x = sample_head, y = list(df.iloc[i,:][1:]), name = item[i] ,opacity=0.8, marker=dict(color=down_color) )
				traces.append( trace )
		else:
			xxx = 1
			#print(list(df.iloc[i,:][1:]))
			Onto_list = list(df.iloc[i,:][1:])
			print(upcount_list)
			max_list = []
			max_list.append(max(upcount_list))
			max_list.append(max(downcount_list))
			y_max = max(max_list)
			term_y = 0 - (int(y_max)/4)
			color_list = []
			cc_x_list = []
			bp_x_list = []
			mf_x_list = []
			cc_y_list = []
			bp_y_list = []
			mf_y_list = []
			bp_color_list = []
			cc_color_list = []
			mf_color_list = []
			for i,e in enumerate(Onto_list):
				if e == "cellular_component":
					cc_color_list.append("rgb(0, 191, 255)")
					cc_x_list.append(sample_head[i])
					cc_y_list.append(term_y)
				elif e == "biological_process":
					bp_color_list.append("rgb(255, 62, 150)")
					bp_x_list.append(sample_head[i])
					bp_y_list.append(term_y)
				elif e == "molecular_function":
					mf_color_list.append("rgb(255, 165, 0)")
					mf_x_list.append(sample_head[i])
					mf_y_list.append(term_y)
				
			trace =  go.Scatter(x = cc_x_list, y = cc_y_list, mode='lines+markers', name="CC", marker=dict(size=10, color="white", line=dict(width=7 , color=cc_color_list)))
			traces.append( trace )
			trace =  go.Scatter(x = bp_x_list, y = bp_y_list, mode='lines+markers', name="BP", marker=dict(size=10, color="white", line=dict(width=7 , color=bp_color_list)))
			traces.append( trace )
			trace =  go.Scatter(x = mf_x_list, y = mf_y_list, mode='lines+markers', name="MF", marker=dict(size=10, color="white", line=dict(width=7 , color=mf_color_list)))
			traces.append( trace )
			
	print ( traces )
	
	# layout
	layout = go.Layout( dict(title=args.title,
			xaxis=dict(title=args.xaxis,
						tickangle = -45,
						tickfont = dict(size = 10, color = 'rgb(55, 83, 109)'),
						#titlefont = dict(size = 18, color = 'rgb(26, 118, 255)')),
						titlefont = dict(size = 18, color = 'rgb(55, 83, 109)')),
			yaxis=dict(title=args.yaxis, showticklabels=True,
						tickfont = dict(size = 15, color = 'rgb(26, 118, 255)'),
						titlefont = dict(size = 18, color = 'rgb(55, 83, 109)')),
			legend=dict(#x=0.03,
						#y=1.03,
						bgcolor='rgba(255, 255, 255, 0)',
						bordercolor='rgba(255, 255, 255, 0)',
						borderwidth = 1
						), # showlegend=False
			barmode = args.barmode,
			bargap = 0.01,
			bargroupgap = 0.1,
			autosize = False,
			width = int(args.width),
			height = int(args.height),
			margin = dict(autoexpand=False, l=100, r=50 , t=60 , b=250)
			))
	# plot
	
	fig = go.Figure(data = traces, layout = layout)
	pyplt( fig, filename = args.outfile ,auto_open = False)


if __name__ == '__main__':
	main()
