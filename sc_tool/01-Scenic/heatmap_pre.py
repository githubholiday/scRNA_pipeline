#!/usr/bin/env python3
import os
import sys
import re
import argparse
import numpy as np

def main():
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-a','--auc',help='auc_mtx.csv file',dest='auc',required=True,type=str)
    parser.add_argument('-g','--group',help='group file',dest='group',required=True,type=str)
    parser.add_argument('-o','--heatmap',help='heatmap file',dest='heatmap',required=True,type=str)
    args=parser.parse_args()

    with open(args.auc, 'r') as fa, open(args.group, 'r') as fg:
        file1_data = [line.strip().replace("(+)", "")  for line in fa]
        file2_data = [line.strip().split("\t") for line in fg]

## 根据文件2中每行的第一列信息找到文件1中对应行的索引
        file1_indices = []
        for row in file2_data:
            index = [i for i, x in enumerate(file1_data) if x.startswith(row[0])]
            file1_indices.extend(index)

# 将文件1中对应的行按照文件2中的顺序进行排序

        sorted_file1 = [file1_data[i].split(",") for i in file1_indices]
        sorted_file1 = sorted(sorted_file1, key=lambda x: file2_data[file1_indices.index([i for i, x in enumerate(file1_data) if x.startswith(x[0])][0])][1:])

# 进行行列转换
        transposed_data = np.transpose(sorted_file1)

# 将处理后的数据写入输出文件中
    with open(args.heatmap, 'w') as fh:
        for row in transposed_data:
            fh.write('\t'.join(row) + '\n')

if __name__=="__main__":
    main()
