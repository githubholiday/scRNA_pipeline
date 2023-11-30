#!/usr/bin/env python3
import os
import sys
import re
import argparse

def main():
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-m','--meta',help='meta.csv file',dest='meta',required=True,type=str)
    parser.add_argument('-o','--group',help='group file',dest='group',required=True,type=str)
    args=parser.parse_args()

    with open(args.meta, 'r') as fm, open (args.group, 'w') as fg:
        lines = fm.readlines()[1:]
        lines.sort(key=lambda line: (line.split(',')[2], line.split(',')[1]))
        result = []
        for line in lines:
            fields = line.strip().split(',')[:3]
            new_fields = '\t'.join(fields)
            result.append(new_fields)
        output = "Cell\tcell_type\tsample\n" + "\n".join(result) + "\n"
        fg.write(output)

if __name__=="__main__":
    main()
