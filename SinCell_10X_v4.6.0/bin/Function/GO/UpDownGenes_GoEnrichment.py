"""
For plot up down genes go enrichment bar
"""
import argparse
import os
import sys
import re
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')

#import warnings
#warnings.filterwarnings("ignore")

Bin = os.path.abspath(os.path.dirname(__file__))

__author__='Yuan Zan'
__mail__= 'seqyuan@gmail.com'
__date__= '20171219'
__modifier__='Huiling Liu'

class up_down_go_enrich:
    df = None
    legend_df = None
    ax = None
    Up_color = '#951F2B'
    Down_color = '#9ABBA6'

    def Init(self,gene_up_down_enrich_go_file,ax,sample,outdir):
        df = pd.read_table(gene_up_down_enrich_go_file, header = 0, index_col=False, encoding='utf-8')
        self.df = df
        df['termcolor'] = '#DAD9DE'        
        df.loc[df['Ontology']=='cellular_component','termcolor'] = '#FEF3C6'
        df.loc[df['Ontology']=='molecular_function','termcolor'] = '#E0BCC8'
        self.df = df
        self.ax = ax

    def plot_bg(self):
        for i in ['biological_process','cellular_component','molecular_function']:
            df = self.df[self.df['Ontology']==i]
            self.ax.bar([df.index[0]], [100], df.index[-1]-df.index[0]+1, alpha=1, color=df.iloc[-1,6],linewidth=0,edgecolor=df.iloc[-1,6],align='edge')
            self.ax.text((df.index[-1]+df.index[0]+1)/2, 92, i, ha='center', va= 'bottom',color='k')      

    def plot_bar(self):
        bar_width = 0.3
        rects1 = self.ax.bar(self.df.index+0.2, self.df['Up_Percent'], bar_width, alpha=1, color=self.Up_color,linewidth = 0,edgecolor=self.Up_color, label='Up',align='edge')
        rects2 = self.ax.bar(self.df.index+0.2+bar_width, self.df['Down_Percent'], bar_width, alpha=1, color=self.Down_color,linewidth = 0,edgecolor=self.Down_color, label='Down',align='edge')
           
        for i in ['bottom','left','top','right']:
            self.ax.spines[i].set_color('black')
            self.ax.spines[i].set_linewidth(0.5)

        self.ax.set_xticks(list(self.df.index+0.5))
        self.ax.set_xticklabels(list(self.df['Term_name']),rotation=70,fontsize='smaller',ha='right')
        self.ax.set_ylabel('Percent of Genes')
        self.ax.tick_params(bottom ='on',top='off',left='on',right='off')
        self.ax.xaxis.set_ticks_position('bottom')
        self.ax.set_xlim([0,self.df.shape[0]])
        self.ax.set_ylim([0,100])
        self.ax.grid(False)

        legend = self.ax.legend(bbox_to_anchor=(1, 1.07), loc=1, borderaxespad=0.,prop={'size':8},ncol=2,frameon=False)        
        #xlabels = self.ax.get_xticklabels()
        #for xl in xlabels:
        #    xl.set_rotation(15)

    def plot_right_txt(self):
        df = self.df
        #[(self.df['Up_Percent'] != 0) & (self.df['Down_Percent'] != 0)]
        all_Up=0;all_Down=0
        #if df['Up_Percent'].sum() != 0 : all_Up = int(list(df['Up_Count']/df['Up_Percent'])[0]*100)
        if df['Up_Percent'].sum() != 0 : all_Up = int([x for x in list(df['Up_Count']/df['Up_Percent']) if str(x) != 'nan'][0]*100)
        #if df['Down_Percent'].sum() != 0 : all_Down = int(list(df['Down_Count']/df['Down_Percent'])[0]*100)
        if df['Down_Percent'].sum() != 0 : all_Down = int([x for x in list(df['Down_Count']/df['Down_Percent']) if str(x) != 'nan'][0]*100)
        self.ax.text(self.df.index[-1]+1.5, 100, str(all_Up), ha='left', va= 'bottom',color=self.Up_color)      
        self.ax.text(self.df.index[-1]+1.5, 100, str(all_Down), ha='left', va= 'top',color=self.Down_color) 
        self.ax.text(self.df.index[-1]+1.5, 0, '0', ha='left', va= 'bottom',color=self.Up_color)      
        self.ax.text(self.df.index[-1]+1.5, 0, '0', ha='left', va= 'top',color=self.Down_color) 
        self.ax.text(self.df.index[-1]+1.5, 50, str('%.1f' %(all_Up/2)), ha='left', va= 'bottom',color=self.Up_color)      
        self.ax.text(self.df.index[-1]+1.5, 50, str('%.1f' %(all_Down/2)), ha='left', va= 'top',color=self.Down_color)

        self.ax.text(self.df.index[-1]+7, 50, 'Number of Genes', ha='center', va= 'center',color='k',rotation='vertical')
        line1, = self.ax.plot([self.df.index[-1]+0.5,self.df.index[-1]+1], [50,50], '-', linewidth=0.2,color='k')


def main():
    parser=argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\n'.format(__author__,__mail__,__date__))
    parser.add_argument('-f','--updown',help='D_G_Up_Down.txt file',dest='updown',type=str,required=True)
    parser.add_argument('-s','--sample',help='sample',type=str,default='sample')
    parser.add_argument('-o','--outDir',help='outDir',dest='outDir',type=str,default=Bin)

    args=parser.parse_args()

    fig = plt.figure(figsize=(9,7),facecolor='white')
    [ax_x, ax_y, ax_w, ax_h] = [0.1,0.57,0.8,0.35]   #[0.05,0.07,0.07,0.66] 
    ax = fig.add_axes([ax_x, ax_y, ax_w, ax_h], frame_on=True,axisbg = 'white')
    ax.set_title(args.sample)
    
    udge = up_down_go_enrich()
    udge.Init(args.updown,ax,args.sample,args.outDir)    
    udge.plot_bg()
    udge.plot_bar()
    udge.plot_right_txt()

    fig.savefig(os.path.join(args.outDir,"{0}_Up_Down.pdf".format(args.sample)),facecolor = fig.get_facecolor(),edgecolor ='none',bbox_inches='tight', dpi=100)
    fig.savefig(os.path.join(args.outDir,"{0}_Up_Down.png".format(args.sample)),facecolor = fig.get_facecolor(),edgecolor ='none',bbox_inches='tight', dpi=100)

if __name__=="__main__": 
    main()
