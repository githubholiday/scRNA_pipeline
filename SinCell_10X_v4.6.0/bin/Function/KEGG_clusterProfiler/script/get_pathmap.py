#!/usr/bin/env python3
"""
do something
"""
import argparse
import os
import sys
import re
import logging
bin = os.path.abspath(os.path.dirname(__file__))
__author__='Ren Xue'
__mail__= 'xueren@genome.cn'
__date__= '2023年02月02日 星期四 13时42分14秒'


def read_map(map):
	maplist={}
	n=0
	with open(map,"r") as file:
		for line in file:
			n+=1
			if n==1:continue
			else:
				tmp =line.rstrip("\n").split("\t")
				mapid=tmp[0].replace("map","")
				maplist[tmp[0]]=mapid
	return maplist

def read_db(db):
	ordb={}
	with open(db,"r") as file:
		for line in file:
			tmp =line.rstrip("\n").split("\t")
			ordb[tmp[2]] = tmp[1]
	return ordb
def main():
	parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter,epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\n'.format(__author__,__mail__,__date__))
	parser.add_argument('-o','--outdir',help='outfile name',dest='outdir',type=str,required=True)
	parser.add_argument('-k','--kegg_map',help='kegg_map',dest='kegg_map',type=str,required=True)
	parser.add_argument('-d','--orgdb',help='orgdb',dest='orgdb',type=str,required=True)
	parser.add_argument('-s','--species',help='species',dest='species',type=str,required=True)
	parser.add_argument('-m','--maplist',help='maplist',dest='maplist',type=str,required=True)
	args=parser.parse_args()
	maplist = read_map(args.maplist)
	ordb = read_db(args.orgdb)
	if not args.species in ordb:
		print("{0} is not in Orgdb list {1}".format(args.species,args.orgdb))
		sys.exit(1)
	else:
		os.system("mkdir -p {0}".format(args.outdir))
		spe = ordb[args.species]
		for i in maplist:
			id ='{0}{1}'.format(spe,maplist[i])
			old_png = '{0}/{1}.png'.format(args.kegg_map,i)
			new_png = '{0}/{1}.png'.format(args.outdir,id)
			old_xml = '{0}/{1}.xml'.format(args.kegg_map,i)
			new_xml = '{0}/{1}.xml'.format(args.outdir,id)
			if os.path.exists(old_png) and os.path.exists(old_xml):
				os.system("cp -f {0} {1}".format(old_png,new_png))
				os.system("cp -f {0} {1}".format(old_xml,new_xml))
		print("kegg map 拷贝处理完成")
if __name__=="__main__":
	main()
