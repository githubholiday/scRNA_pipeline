#!/usr/bin/env python3
import os
import sys
import re
import argparse
bindir = os.path.abspath(os.path.dirname(__file__))

__author__='Yang Zhang'
__mail__= 'yangzhang@genome.cn'
__doc__='used for report to add ko' 

def main():
        parser=argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter,
                        epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
        parser.add_argument('-i','--report',help='input file',dest='report',required=True)
        parser.add_argument('-k','--kolist',help='input file',dest='kolist',required=True)
        parser.add_argument('-o','--output',help='output file',dest='output',required=True)
        args=parser.parse_args()
        ko_dict={}
        with open (args.kolist , 'r') as ko:
            for line in ko:
                newline = line.strip().split('\t')
                ko_dict[newline[0]] = newline[1]
        OUT = open( args.output , 'w')
        with open ( args.report , 'r') as report:
            for line in report:
                if line.startswith("Map"):
                    OUT.write(line)
                else:
                    tmp = line.strip().split('\t')
                    Up_Gene , Down_Gene , url = tmp[11] , tmp[13] , tmp[15]
                    ko_list = [url]
#                    print("Up_Gene:\n",Up_Gene,"\nDown_Gene:\n",Down_Gene)
                    up_gene_list = Up_Gene.split('|')
                    new_up_gene_list = []
                    for i in up_gene_list:
                        if i in ko_dict:
                            gene_ko = i+"|"+ko_dict[i]
                            new_up_gene_list.append( gene_ko )
                            ko_list.append( ko_dict[i] )
                    new_up_gene = ";".join( new_up_gene_list )
                    down_gene_list = Down_Gene.split('|')
                    new_down_gene_list = []
                    for i in down_gene_list:
                        if i in ko_dict:
                            gene_ko = i+"|"+ko_dict[i]
                            new_down_gene_list.append( gene_ko )
                            ko_list.append( ko_dict[i] )
                    new_down_gene = ";".join( new_down_gene_list )
                    new_url = "+".join(ko_list)
                    out_line ="\t".join( [tmp[0],tmp[1],tmp[2],tmp[3],tmp[4],tmp[5],tmp[6],tmp[7],tmp[8],tmp[9],tmp[10],new_up_gene,tmp[12],new_down_gene,tmp[14],new_url,tmp[16]] )+"\n"
                    OUT.write( out_line )
        OUT.close()


if __name__=="__main__":
        main()
