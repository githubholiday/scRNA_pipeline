#! /usr/bin/env python3
import argparse
import os
import sys
import pandas as pd
import numpy as np
import plotly as py
import string
import plotly.graph_objs as go
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as cols
pyplt = py.offline.plot

#offline.init_notebool_mode(connect=False)
def main():
    parser=argparse.ArgumentParser(description='Bar plot using plotly')
    parser.add_argument('-i','--input',help='input file or stdin',dest='input',required=True)
    parser.add_argument('-o','--outfile',help='output file',dest='outfile',required=True)
    parser.add_argument('-t','--title',help='title',dest='title',default='DE Gene Count')
    parser.add_argument('-x','--xaxis',help='xaxis',dest='xaxis',default='')
    parser.add_argument('-y','--yaxis',help='yaxis',dest='yaxis',default='Count')
    parser.add_argument('-m','--barmode',help='barmode',dest='barmode',choices=['group','stack'])
    parser.add_argument('-T','--trans',help='trans or not',dest='trans',default='v')
    parser.add_argument('-C','--colorscale',help='colorscale',dest='colorscale')
    parser.add_argument('-a','--angle',help='angle of xaxis words',dest='angle',default='0')
    parser.add_argument('-s','--size',help='words size of xaxis',dest='size',default='16')
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
    colors = np.linspace(0,1,len(sample_head))

    traces = []
    if args.trans=='h' and args.colorscale: 
        for i in range( len(item) ):
            trace = go.Bar( x = list(df.iloc[i,:][1:]), y = sample_head, name = item[i], 
                    orientation= 'h', ##调整图形转置
                    marker= dict(color=colors,
                        colorscale = args.colorscale), ##添加颜色
                    opacity=0.8 )
            traces.append( trace )
        print ( traces )
    elif args.trans=='v' and  args.colorscale:
        for i in range( len(item) ):
            trace = go.Bar( y = list(df.iloc[i,:][1:]), x = sample_head, name = item[i],
                    marker= dict(color=colors,
                        colorscale = args.colorscale),
                    opacity=0.8 )
            traces.append( trace )
        print ( traces )
    else:
        for i in range( len(item) ):
            trace = go.Bar( y = list(df.iloc[i,:][1:]), x = sample_head, name = item[i],opacity=0.8 )
            traces.append( trace )
        print ( traces )
    
    # layout
    layout = go.Layout( dict(title=args.title,
            xaxis=dict(title=args.xaxis, 
                        tickangle = int(args.angle),
                        tickfont = dict(size = 16, color = 'rgb(55, 83, 109)'),
                        titlefont = dict(size = 18, color = 'rgb(26, 118, 255)')),
            yaxis=dict(title=args.yaxis, showticklabels=True,
                        tickfont = dict(size = int(args.size), color = 'rgb(26, 118, 255)'),
                        titlefont = dict(size = 18, color = 'rgb(55, 83, 109)')),
            legend=dict(#x=0.03,
                        #y=1.03,
                        bgcolor='rgba(255, 255, 255, 0)',
                        bordercolor='rgba(255, 255, 255, 0)',
                        borderwidth = 1), # showlegend=False
            barmode = args.barmode,
            bargap = 0.05,
            bargroupgap = 0.1,
            #margin =go.layout.Margin(l=500,r=500,b=100,t=100,pad=0),
            ),
            )
    fig = go.Figure(data = traces, layout = layout)
    pyplt( fig, filename=args.outfile ,auto_open=False)


if __name__ == '__main__':
    main()
