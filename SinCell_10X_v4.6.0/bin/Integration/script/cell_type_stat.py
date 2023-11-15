#!/usr/bin/evn python3
######################################################################### import ##########################################################
import argparse
import os
import sys
import time
import re
import pandas
######################################################################### ___  ##########################################################
__author__ = 'qianzhao'
__mail__ = 'qianzhao@genome.cn'
__date__ = '2021年6月8日 星期二'
__version__ = '3.0'
###################################### main  #######################################
def Read_scibet(csv):
    sci_dict={}
    with open(csv) as F:
        for i,line in enumerate(F):
            if i ==0:continue
            if line.startswith(','):continue
            arrays = line.rstrip().split(',')
            cluster,celltype=arrays[2],arrays[3]
            #sci_dict[cluster]['total'] = 0
            if cluster not in sci_dict :
                sci_dict[cluster] = {}
                sci_dict[cluster]['total'] = 1
                if celltype not in sci_dict[cluster]:
                    sci_dict[cluster][celltype]= 1
            else:
                sci_dict[cluster]['total'] += 1
                if celltype not in sci_dict[cluster]:
                    sci_dict[cluster][celltype]= 1
                else:
                    sci_dict[cluster][celltype] += 1
    
    #print(sci_dict)
    #for clu in sci_dict.keys():
    #    print(clu)
    #    print(sorted(sci_dict[clu].items(), key=lambda item:item[1],reverse=True))
    #    sci_dict[clu]=(sorted(sci_dict[clu].items(), key=lambda item:item[1],reverse=True))[0]
    return sci_dict

def Read_marker(csv):
    mark_dict={}
    with open(csv) as F:
        for i,line in enumerate(F):
            if i ==0:
                title=line.rstrip().split(',')
                gene_index=title.index('gene_name')
                cluster_index=title.index('cluster')
            else:
                arrays=line.rstrip().split(',')
                Gene_name,cluster=arrays[gene_index],arrays[cluster_index]
                mark_dict.setdefault(cluster,[]).append(Gene_name)
    return mark_dict

def Read_list(annofile):
    delltype_dict={}
    with open(annofile) as F:
        for line in F:
            arrays=line.rstrip().split('\t')
            cell_type,marker=arrays[0],arrays[1]
            markers=marker.split(',')
            marker2=[]
            marker2=[mar.upper() for mar in markers if isinstance(mar,str)]
            delltype_dict[cell_type]=marker2
    return delltype_dict

def main():
    function="this progrem is used to merge scibet result 、 annolist and marker gene"
    parser=argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}\nfunction:\t{4}'.format(__author__,__mail__,__date__,__version__,function))
    parser.add_argument('-m',help='cluster marker gene file',type=str,required=True)
    parser.add_argument('-s',help='scibet file',type=str,required=True)
    parser.add_argument('-o',help='output file',type=str,required=True)
    parser.add_argument('-l',help='anno maker list',type=str,required=False)
    parser.add_argument('-r',help='scibet stat',type=str,required=True)
    args=parser.parse_args()

    sci_dict = Read_scibet(args.s)
    clust_list=[int(x) for x in sci_dict.keys()]
    clust_list.sort()
    dic_scibet={}
    #####对scibet的结果进行统计，统计前3种细胞类型的数量
    with open(args.r,'w')as OUT:
        OUT.write("CLUSTER\tCellType1\tCellType1_Count\tCellType2\tCellType2_Count\tCellType3\tCellType3_Count\n")
        for clu in clust_list:
            new_arr=[]
            new_arr=sorted(sci_dict[str(clu)].items(), key=lambda item:item[1],reverse=True)
            if new_arr[0][0] != 'total': ###当细胞均分为1类时，total和单类细胞的数量一致
                new_arr[0],new_arr[1]=new_arr[1],new_arr[0]
            dic_scibet[str(clu)]=new_arr[1]
            celltyperatio={}
            for i,line in enumerate(new_arr):
                celltyperatio[i]=str(round(new_arr[int(i)][1]/new_arr[0][1]*100,2))
            if len(new_arr)==2:
                celltypeline=str(clu)+'\t'+new_arr[1][0]+'\t'+str(new_arr[1][1])+'('+ celltyperatio[1] +')'+'\t-'*4 +'\n'
            elif len(new_arr)==3:
                celltypeline=str(clu)+'\t'+new_arr[1][0]+'\t'+str(new_arr[1][1])+'('+ celltyperatio[1] +')'+'\t'+new_arr[2][0]+'\t'+str(new_arr[2][1])+'('+ celltyperatio[2] +')'+'\t-'*2 +'\n'
            elif len(new_arr) >=4:
                celltypeline=str(clu)+'\t'+new_arr[1][0]+'\t'+str(new_arr[1][1])+'('+ celltyperatio[1] +')'+'\t'+new_arr[2][0]+'\t'+str(new_arr[2][1])+'('+ celltyperatio[2] +')'+'\t'+new_arr[3][0]+'\t'+str(new_arr[3][1]) +'('+ celltyperatio[3] +')'+'\n'
            OUT.write(celltypeline)
    #### 对cluster marker进行验证
    dic_mark = Read_marker(args.m)
    with open(args.o,'w')as OUT:
        if args.l:
            anno_mark_dict=Read_list(args.l)
            OUT.write("CLUSTER\tANNO\tExpected_Count\tCovered_Count\tRatio(%)\tExpected_Marker\tCovered_Marker\n")
            for i in clust_list:
                ANNO=dic_scibet[str(i)][0]
                All_maker=dic_mark[str(i)]
                Covered_Count=0;covered_Marker=[]
                if ANNO in anno_mark_dict:
                    Expected_Count=len(anno_mark_dict[ANNO])
                    Expected_Marker=','.join(anno_mark_dict[ANNO])
                    for mark in anno_mark_dict[ANNO]:
                        if mark in All_maker:
                            Covered_Count +=1
                            covered_Marker.append(mark)
                        else:continue
                    if Covered_Count!=0:
                        Ratio='%.2f' %((Covered_Count/Expected_Count)*100)
                        Covered_Marker=','.join(covered_Marker)
                    else:
                        Ratio=0
                        Covered_Marker=','.join(All_maker[0:10])

                else:
                    Expected_Count=0
                    Expected_Marker='-'
                    Ratio=0
                    Covered_Marker=','.join(All_maker[0:10])
                OUT.write(str(i) +"\t"+ANNO+"\t"+str(Expected_Count)+"\t"+str(Covered_Count)+"\t"+str(Ratio)+"\t"+Expected_Marker+"\t"+Covered_Marker+'\n')
        else:
            OUT.write("CLUSTER\tANNO\tCovered_Marker\n")
            for i in clust_list:
                #print(clust_list)
                print(dic_scibet[str(i)])
                ANNO=dic_scibet[str(i)][0]
                if str(i) in dic_mark:  ##可能存在没有显著差异基因的情况
                    All_maker=dic_mark[str(i)]
                    Covered_Marker=','.join(All_maker[0:10])
                else:
                    Covered_Marker='-'
                OUT.write(str(i) +"\t"+ANNO+"\t"+Covered_Marker+'\n')

if __name__=="__main__":
    main()
