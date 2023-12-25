#!/usr/bin/env python3
import os
import sys
import re
import argparse
import numpy as np
import pandas as pd

def main():
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i','--infile',help='infile',dest='infile',required=True)
    parser.add_argument('-o','--outdir',help='outdir',dest='outdir',required=True)
    parser.add_argument('-p','--prefix',help='prefix',dest='prefix',required=True)
    args=parser.parse_args()
    cell_type_out = '{0}/{1}.celltype.xls'.format( args.outdir, args.prefix)
    cluster_out = '{0}/{1}.cluster.xls'.format( args.outdir, args.prefix)
    
    data = pd.read_csv( args.infile,sep='\t' )
    cell_type_mean = data.groupby('cell_type').mean()
    cell_type_mean_drop = cell_type_mean.drop(labels=['cell_cluster'], axis=1)
    cell_type_mean_drop.to_csv( cell_type_out,sep='\t' )

    cluster_mean = data.groupby('cell_cluster').mean()
    cluster_mean.to_csv( cluster_out,sep='\t' )
        

if __name__=="__main__":
    main()
