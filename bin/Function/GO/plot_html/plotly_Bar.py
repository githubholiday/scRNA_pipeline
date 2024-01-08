#! /usr/bin/env python3
import argparse
import os
import sys
import pandas as pd
import plotly as py
import plotly.graph_objs as go
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
		trace = go.Bar( x = sample_head, y = list(df.iloc[i,:][1:]), name = item[i] ,opacity=0.8 )
		traces.append( trace )
	print ( traces )
	
	# layout
	layout = go.Layout( dict(title=args.title,
			xaxis=dict(title=args.xaxis,
						tickfont = dict(size = 15, color = 'rgb(55, 83, 109)'),
						titlefont = dict(size = 18, color = 'rgb(26, 118, 255)')),
			yaxis=dict(title=args.yaxis, showticklabels=True,
						tickfont = dict(size = 15, color = 'rgb(26, 118, 255)'),
						titlefont = dict(size = 18, color = 'rgb(55, 83, 109)')),
			legend=dict(#x=0.03,
						#y=1.03,
						bgcolor='rgba(255, 255, 255, 0)',
						bordercolor='rgba(255, 255, 255, 0)',
						borderwidth = 1), # showlegend=False
			barmode = args.barmode,
			bargap = 0.05,
			bargroupgap = 0.1,
			#autosize = False,
			#width = int(args.width),
			#height = int(args.height),
			#margin = dict(autoexpand=False, l=300, r=300, b=100, t=100, pad=0)
			))
	# plot
	fig = go.Figure(data = traces, layout = layout)
	pyplt( fig, filename=args.outfile, auto_open= False)


if __name__ == '__main__':
	main()
