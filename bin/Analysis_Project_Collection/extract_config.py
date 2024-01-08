#! /usr/bin/env python3

import argparse
import os
import sys
import re

__author__='Simon lee'
__mail__='huayunli@genome.cn'
__doc__=''


def read_file(file):
	if os.path.exists(file):
		with open(file,'r')as fl:
			for line in fl:
				if not line.strip(''):continue
				yield line.strip('\n')
	else:
		print ("file {0} didn't exist!".format(file))


def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-c','--config',help='input pepline parameter config',dest='cfg',required=True)
	parser.add_argument('-i','--indir',help='pepline result dir',dest='indir',required=True)
	parser.add_argument('-o','--outdir',help='human single cell database dir',dest='outdir',required=True)

	args=parser.parse_args()
	###############参数梳理###############
	cfg, indir, outdir = read_file(args.cfg), args.indir, args.outdir
	###############获取ms文件列表##############
	cfg_dict = {}
	for line in cfg:
		if not line.strip(""):continue
		if line.startswith("[sample]"):
			key = "sample"
			continue
		if line.startswith("[group]"):
			key = "group"
			continue
		if line.startswith("[cmp]"):
			key = "cmp"
			continue
		if line.startswith("[Para]"):
			key = "Para"
			continue
		if key in cfg_dict:
			cfg_dict[key] = cfg_dict[key] + "\n" + line
		else:
			cfg_dict[key] = line
	para_list = cfg_dict["Para"].split("\n")
	num = 0
	for i,e in enumerate(para_list):
		if "Para_ref" in e:
			human_ref = ["GRCh38","human","hg19"]
			for ref in human_ref:
				if ref in e:
					num = num + 1
		if "Para_project" in e:
			pro_id = e.split(" = ")[1]
	if num == 0:
		print("此项目没有人类样品，无需收集到数据库！\n")
		exit(0)
	else:
		print("\n")
		print("此项目含有人类样品，cellranger和seurat结果将收集到：{0}\n".format(outdir))
		print("\n")
		sfile = outdir + "/3_Sample_Info/" + pro_id + ".xls"
		with open(sfile, "w")as sf:
			sf.write(cfg_dict["sample"])
	cellranger_dir = indir + "/CellRanger_Count"
	out_cellranger_dir = outdir + "/1_Cellranger"
	seurat_dir = indir + "/DiffExpr"
	out_seurat_dir = outdir + "/2_Single_Sample_Seurat" ###20211104
	cmd1 = "cp -rf {0}/* {1}/.".format(cellranger_dir, out_cellranger_dir)
	cmd2 = "cp -rf {0}/* {1}/.".format(seurat_dir, out_seurat_dir)
	#cmd1 = "ln -sf {0}/* {1}/.".format(cellranger_dir, out_cellranger_dir)
	#cmd2 = "ln -sf {0}/* {1}/.".format(seurat_dir, out_seurat_dir)
	os.system(cmd1)
	os.system(cmd2)


if __name__ == '__main__':
	main()
