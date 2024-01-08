#! /usr/bin/env python3
# -*- coding=utf-8 -*-
import argparse
import sys
import os
import re
bindir = os.path.abspath(os.path.dirname(__file__))

__author__='Su Lin'
__mail__= 'linsu@annoroad.com'
__doc__='the program is to create HTML table '

pat1=re.compile('^\s+$')

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--input',help='input kegg table file',dest='input',type=open,required=True)
	parser.add_argument('-o','--output',help='output html file',dest='output',required=True)
	parser.add_argument('-c','--col',help='the index ko colnm',type=int,required=True)
	args=parser.parse_args()
	cmd = 'cp -r {0}/kegg_html/src {1}'.format(bindir,os.path.dirname(args.output))
	os.system(cmd)
	name = os.path.splitext(os.path.basename(args.output))[0]
	output = open(args.output,'w',encoding='utf-8')
	out_content = '''
<!DOCTYPE html>
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=gb2312" />
<title>Pathway Enrichment</title>
<SCRIPT type="text/javascript" src="src/js/jquery-1.4.2.min.js"></SCRIPT>
<link href="src/base.css" type="text/css" rel="stylesheet">
</head>

<body>
<div>
<img style="height:50px;padding-top:10px" src="src/logo2.jpg"  />
<h1>The most enriched pathway terms</h1>
</div>

<p>Statistic method: hypergeometric test</p>
<p>FDR correction method: Benjamini and Hochberg</p>
'''
	for index,line in enumerate(args.input):
		if index == 0:
			lines = line.rstrip().split('\t')
			out_content += '<table class="gy" id="tblMain" border="1">\n<tr>'
			for item in lines:
				out_content += '<th>{0}</th>\n'.format(item)
			out_content += '</tr>'
			continue
		lines = line.rstrip().split('\t')
		if lines[9] == '-': continue
		map = lines[args.col]#.replace('map','ko')
		out_content += '<tr><td>\n<a href="{0}/{1}.html" target="_blank">{2}</a></td>'.format(name,map,lines[0])
		for item in lines[1:]:
			out_content += '<td>{0}</td>'.format(item[0:40])
		out_content += '</tr>\n'
	out_content += '''
</table>
<div id="shuoming" style="text-align:left;">
由于KEGG注释结果内容较多，所以KEGG结果展示了一部分内容，全部注释结果请见相应的文件夹中的*kegg_report.xls文件<br>
结果说明如下：<br>
（1）Name：通路的名称；<br>
（2）Map：进行富集的map条目；<br>
（3）Count1，Count2，Count3，Count4：进行Fisher检验的四个数据，分别为上面公式里的m，M-m，n-m，N-n-M+n；<br>
（4）p：检验后的p值；<br>
（5）q：多重检验校正的p值；<br>
（6）Gene_in_background：在背景中的基因的名称；<br>
（7）Gene_in_DE：差异表达基因中注释到该通路的差异表达基因名称；<br>
（8）Up_Gene：差异表达基因中注释到该通路的上调基因；<br>
（9）Up_Count：差异表达基因中注释到该通路的上调基因的个数；<br>
（10）Down_Gene：差异表达基因中注释到该通路的下调基因；<br>
（11）Down_Count:差异表达基因中注释到该通路的下调基因的个数；<br>
（12）Links：该通路的KEGG数据库链接；<br>
（13）Result：该通路是否显著。<br>
</div>

<div id="bot">
      电话 : 4008-986-980　　　网址 : www.genome.cn　　　邮箱 : service@genome.cn       <font style="margin-left:50px"> 微信 :</font>  <IMG id="weixin" style="margin-bottom:-5px" src="src/8.png"> 
</div>
<div id="wx" style="display:none"><IMG id="code" src="src/wx2.png"></div>
<!--<IMG id="backtop" src="src/backtop.gif">-->
<script>
$(function(){
//二维码
$("#weixin").mouseover(function () {
$("#wx").css("display","block");
});
$("#weixin").mouseleave(function () {
$("#wx").css("display","none");
});
$('#closecode').click(function(){//关闭二维码
$('#code').fadeOut('slow');
$(this).fadeOut('slow');
});
})
</script>
<!--隔行变色-->
<script type='text/javascript'>
var trs=document.getElementsByTagName('tr');
for(var i=0; i<trs.length;i++){
if(i%2==0){
trs[i].style.background="#DBEAF9";
}else{
trs[i].style.background="#F2F4FF";
}
}
</script>
</body></html>
'''
	output.write(out_content.encode('utf-8').decode('utf-8'))
	output.close()
	#cmd = '/usr/bin/iconv -f UTF-8 -t GBK -o {0} {1} && rm {1}'.format(args.output,args.output+'_tmp')
	#os.system(cmd)


if __name__ == '__main__':
	main()
