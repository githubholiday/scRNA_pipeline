#! /bin/bash
basepath=$(cd `dirname $0`;pwd)

help(){
	echo "生成蛋白互作网络动态图"
        printf "%-5s %-30s\n" "-n" "蛋白节点描述文件,如/nas/sentieon/wanghong02/networkD3/Si1_NC.ppi.node.txt"
        printf "%-5s %-30s\n" "-l" "蛋白关系描述文件,如/nas/sentieon/wanghong02/networkD3/Si1_NC.ppi.edge.txt"
	printf "%-5s %-30s\n" "-o" "输出文件，如./test.html"
        echo "eg: $0 -n /nas/sentieon/wanghong02/networkD3/Si1_NC.ppi.node.txt -l /nas/sentieon/wanghong02/networkD3/Si1_NC.ppi.edge.txt -o ./test.html"
        exit 1	
}
while getopts "n:l:o:" arg
do
        case $arg in
        n)
                node=$OPTARG
                ;;
        l)
                link=$OPTARG
                ;;
        o)
                output=$OPTARG
                ;;
        h|*)
                help
                ;;
        esac
done
if [ $# -ne 6 ];then
	help
fi
##处理node文件和link文件
outpath=`dirname $output`
node_name=`basename $node`.tmp
link_name=`basename $link`.tmp
node_tmp=$outpath/$node_name
link_tmp=$outpath/$link_name
##处理node文件
echo -e "\"name\"\t\"group\"\t\"size\"" > $node_tmp
cat $node | sed "s/^/\"/g" | sed "s/\t/\"\\t\"/" | sed "s/$/\"\t20/" >> $node_tmp
##处理link文件
echo -e "\"source\"\t\"target\"\t\"value\"" > $link_tmp
cat $link | while read line;do
	node1=`echo $line | awk '{print $1}'`
	node2=`echo $line | awk  '{print $2}'`
	node1_id=`cat $node| grep -n $node1| awk -F ':' '{print $1}'`
	node1_id=`expr $node1_id - 1`
	node2_id=`cat $node| grep -n $node2| awk -F ':' '{print $1}'`
	node2_id=`expr $node2_id - 1`
	echo -e "$node1_id\t$node2_id\t1" >> $link_tmp
done
#调用R脚本
$basepath/networkD3.R $node_tmp $link_tmp $output

