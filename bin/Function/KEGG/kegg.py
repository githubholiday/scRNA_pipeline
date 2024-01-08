#! /usr/bin/env python3
import re,os,sys
import argparse
import glob
import time
import numpy as np
bindir = os.path.abspath(os.path.dirname(__file__))

__author__='Su Lin'
__modifier__= 'jiangdezhi'
__modified_time = '2016.12.5'
__mail__= 'linsu@annoroad.com'
__doc__='the decription of program'

def LoadHtml(ko,mapcolor_dir,outdir,ko_dic,gene_des):
	patt = re.compile('\<img src\=\"(.*.png)\"')
	html_file = '{0}/{1}.html'.format(mapcolor_dir,ko)
	html = open(html_file,'r') 
	print(html_file)
	new_html = open('{0}/{1}.html'.format(outdir,ko),'w') 
	
	content = html.read()	
	content = content.replace('</head>','''
<script type="text/javascript">
<!--
function showInfo(info) {
obj = document.getElementById("detail");
obj.innerHTML = "<div style='cursor: pointer; position: absolute; right: 5px; color: #7a9c14;' onclick='javascript: document.getElementById(\\"detail\\").style.display = \\"none\\";' title='close'>Close</div>" + info;
obj.style.display = "";
}
//-->
</script>
</head>''')

	content = content.replace('</body>','''
<div id='detail' style='position: fixed; left:25%; top:0; width: 50%; border: 3px solid #A7C942; background-color: #EFFEC1; filter: alpha(opacity=95); opacity: 0.95; font-size: 11px; padding-right: 20px; display: none;' onmouseover="javascript: this.style.filter = 'alpha(opacity=100)'; this.style.opacity = 1;" onmouseout="javascript: this.style.filter = 'alpha(opacity=95)'; this.style.opacity = 0.95;"></div>
</body>''')
	content = content.replace('<table><tr><td>\n<form name="selmenu" method="get">','<!--\n<table><tr><td>\n<form name="selmenu" method="get">')
	content = content.replace('[\n<a','<!--\n[\n<a')
	content = content.replace('\n<table id="description"','-->\n<table id="description"')
	content = content.replace('style="display: none;"','')
	content = content.replace('/Fig/bget/kegg3.gif','kegg3.gif')
	content = content.replace('<td valign="bottom"','<!--\n<td valign="bottom"')
	content = content.replace('onmouseup="btn(this,\'Hb\')" /></a>\n  </td>','onmouseup="btn(this,\'Hb\')" /></a>\n  </td>\n-->\n')
	content = content.replace('<!--</form>-->','</form>')
	content = content.replace('href="','href="http://www.genome.jp')
	contents = content.split('\n')

	index = 0
	for line in contents:
		if 'pathwayimage' in line:
			line = line.replace('<img','-->\n<img')
			line = line.replace('/kegg/pathway/map/','')

		if line.startswith('<div id="poplay"'):
			line = ''

		if line.startswith('\t\t<area id=') and '/entry/K' in line:
			ko_patt = re.compile('title="(.*)" />')
			pp=ko_patt.search(line)
			ko_info = pp.group(1)
			ko_list = line.split('entry/')[1].split('"\ttitle=')[0].split('+')
			ko_list = [i for i in ko_list if i.split()[0] in ko_dic]
			if ko_list == []:mouse_info=''
			else:mouse_info = GetMouse(ko_info,ko_list,ko_dic,gene_des)
			line = line.replace('/>','onmouseover=\'javascript: showInfo("'+mouse_info+'");\'/>')

		index += 1
		new_html.write(line +"\n") 

	# cp mapid.png to outdir
	print ("cp " + mapcolor_dir + "/" + ko + ".png " + outdir )
	os.system("cp " + mapcolor_dir + "/" + ko + ".png " + outdir ) 	
	
	html.close()
	new_html.close()

def GetMouse(ko_info,ko_list,ko_dic,gene_des):
	info = '<ul><li>{0}</li><ul>'.format(ko_info)
	info += '<li><font color=\\"red\\">Up Regulated Genes</font><ul><font color=\\"red\\">'
	up = ''
	down = ''
	for ko in ko_list:
		ko  = ko.split(' ')[0]

		up_gene_list = ko_dic[ko]['up']
		if len(up_gene_list) == 0: up_gene_info=[]
		for up_gene in up_gene_list:
			up_gene_info = ['<li>'+i+':'+gene_des[up_gene]+'</li>' for i in up_gene_list]
		up += ' '.join(up_gene_info)

		down_gene_list = ko_dic[ko]['down']
		if len(down_gene_list) == 0: down_gene_info=[]
		for down_gene in down_gene_list:
			down_gene_info = ['<li>'+i+':'+gene_des[down_gene]+'</li>' for i in down_gene_list]
		down += ' '.join(down_gene_info)
	
	if up == '': up = "None"
	if down == '': down = "None"

	info += up+'</font></ul>'
	info += '<li><font color=\\"green\\">Down Regulated Genes</font></li><ul><font color=\\"green\\">'
	info +=  down+'</font></ul></ul></ul>'
	info = info.replace("\'","") #remove single quote of gene description  

	return(info)

def GetFeature(feature_file):
	ko_dic = {}
	gene_des = {}
	for line in feature_file:
		if line.startswith('#'):continue
		gene,kegg,map,up_down,descri = line.rstrip().split('\t')
		keggs = kegg.split(",")
		for keg in keggs:	
			if keg not in ko_dic:
				ko_dic[keg] = {}
				ko_dic[keg]['up'] = []
				ko_dic[keg]['down'] = []
			ko_dic[keg][up_down].append(gene)
			gene_des[gene] = descri

	for ko in ko_dic.keys():
		for key in ko_dic[ko].keys():
			if ko_dic[ko][key]:
				ko_dic[ko][key] = np.unique(ko_dic[ko][key])
	
	return(ko_dic,gene_des)

def MKDIR(dir):
	
	patt = re.compile("/")
	
	if dir.startswith("/"):
		dir = dir 
	
	elif not patt.search(dir):
		dir = '{0}/{1}'.format(os.getcwd(),dir)

	else:
		print ("Error! Please provided right directory format,such as /*/*/dir or dir!\n")
		os.exit()

	# create dir

	if not os.path.isdir(dir):
		 #os.system("mkdir -p " + str(dir))
		 os.makedirs(dir)
	else:
		 #print ( " rm -rf ###" + str(dir)) 
		 os.system( " rm -rf " + str(dir)) 
		 #os.system("mkdir -p " + str(dir))
		 os.makedirs(dir)


def Stat(mapcolor_dir,outdir,kegg_report,map_num):

	png_in_mapcolor = glob.glob(mapcolor_dir + "/*.png")
	png_in_outdir = glob.glob(outdir + "/*.png")

	png_count_in_mapcolor,png_count_in_outdir = len(png_in_mapcolor),len(png_in_outdir)

	if ( png_count_in_mapcolor == png_count_in_outdir == map_num):
		print ('{0} png count : {1}'.format(kegg_report.name, map_num)) 		
		print ('{0} png count : {1}'.format(mapcolor_dir, png_count_in_mapcolor)) 		
		print ('{0} png count : {1}'.format(outdir, png_count_in_outdir)) 		
	else:
		print ("Warning: Png count among mapcolor_dir、outdir and kegg_report file is not equal!Please make sure where is wrong!")
		print ('{0} png count : {1}'.format(kegg_report.name, map_num)) 		
		print ('{0} png count : {1}'.format(mapcolor_dir, png_count_in_mapcolor)) 		
		print ('{0} png count : {1}'.format(outdir, png_count_in_outdir)) 		

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-r','--report',help='input kegg report file',type=open,required=True)
	parser.add_argument('-f','--feature',help='input feature file',type=open,required=True)
	parser.add_argument('-o','--outdir',help='outdir',dest='outdir',required=True)
	parser.add_argument('-m','--mapcolor_dir',help='mapcolor_dir',dest='mapcolor_dir',required=True)
	args=parser.parse_args()
	ko_dic,gene_des = GetFeature(args.feature)
	filter_ko = ["map01100","map01110","map01120","map01130","map00312"] # map00312 had been deprecated by KEGG database

	# create outdir
	report,outdir,mapcolor_dir = args.report, args.outdir, args.mapcolor_dir
	MKDIR(outdir)
	
	map_num = 0
	
	for index,line in enumerate(report):
		if index == 0 or line.startswith('#'):continue
		ko = line.split("\t")[0]
		if ko in filter_ko: continue
	
		map_num += 1; 
		ko_file ,png_file = mapcolor_dir+'/{0}.html'.format(ko) ,mapcolor_dir+'/{0}.png'.format(ko)
		
		if not os.path.exists(ko_file): continue
		
		LoadHtml(ko,mapcolor_dir,outdir,ko_dic,gene_des)

	os.system("cp " + bindir + "/kegg_html/kegg3.gif " + outdir ) 	

	# end 
	print ("The onmouse of htmls have been completed successfully!\n")
	print ("The KEGG map colored and html onmouse have been completed successfully!\n")
	Stat(mapcolor_dir,outdir,args.report,map_num)

if __name__ == '__main__':
	main()
