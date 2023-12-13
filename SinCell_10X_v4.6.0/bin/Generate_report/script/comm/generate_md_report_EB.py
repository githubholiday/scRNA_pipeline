#! /usr/bin/env python3
'''
## 参数
* -d , --indir : 输入目录，可以是包含一个模块，也可以包含多个模块的结果，每个模块的目录下面应该包含"input.json" 和 "template.md"这两个文件，才会被检索
* -o , --output : 输出文件，主要是md文件，可以用于后续产生html或者pdf。 
                  同时在同级目录下，会生成一个mapping.json，用于动态表格的json检索；
                  还有一个tmp_template.tmp, 可以用来拍错
* -m , --module: 指定模块的顺序，默认以glob的。 每行一个模块名称。
* -pdf, --pdf ： 如果指定，生成的报告为pdf转化用的，去掉动态表，下载链接，动态图、结果目录等pdf不支持的格式
* -l, --local : 如果指定，则使用服务器的相对路径来表示文件，否则会加上云端的路径前缀来表示文件；添加方法受到-pipeline和-smalltool的区别影响；
* -u , --upload : 如果指定，则会重新整理一个新的upload目录。
* -pipline , --pipeline : 提供项目的项目号
* -smalltool , --smalltool: 提供小工具运行的路径和小工具id两个参数

## 功能
输入md模板文件，然后使用json作为变量的输入，导出最终的md的文件。
md模板中有以下类型：
* 字符串：{{string}}
* 图片： ![title,w,h]{{image_}},支持单图、多图以及html结尾的动图， 末尾带* ，则自动添加下载链接, 可以加一个image_*_zip 来表示zip文件路径
* 表格： ![title,column_number,row_number]{{table_}}，自动将第一个表格转为md，末尾带*，自动添加下载链接，是否需要压缩后下载，可以加一个table_*_zip 来表示zip文件路径
* 下载： ![title]{{download_}}
* 跳转链接: ![title]{{href_}} , 构建跳转链接。
* json: ![title]{{json_}} , 后面的注释以@开头
* 小工具：![参数]{{smalltool_}}
* note : !!! , 以三个感叹号开头

## 功能：
* 添加图片width的维持原比例缩放
* 整理文件到upload目录下，并进行压缩成文件
* 支持跳转  SmallTool:(输入文件:upload/*_differential/*venn/*.report.xls),(首行为表头:否),(边框类型:虚线),(边框颜色:粉红,紫罗兰,天蓝色,纯黄,马鞍
    色),(填充颜色:粉红,紫罗兰,天蓝色,纯黄,马鞍色),(透明度:5),(id:bc088020d5c64df9a712d5f5b0f94487),(version:2.0.1),0,,,>
* 自动增加标题级别
* 增加公司信息
* 图片增加自定义功能，也可以不写，默认进行缩放
* 表格增加自定义功能，也可以不写，默认4列，10行
* 根据文件的多少来判断选用单图还是多图
* html图片下面必须有一个 png图，即file.replace('html','png')

## json和md文件注意事项
* json需要为二级dict，第一级为模块名称，第二级为变量和值
* md中需要有@@@@模块名 作为起始，且模块名称和json中应当一致，否则会报错 。
示例见 ： /annoroad/data1/bioinfo/PMO/liutao/bin/report/bin/test/Analysis

## 注意
* 为了避免风险，所有文件的文件名最好带有样品名，或者保持唯一；
json里有*/a.pdf ，这种会导致copy 文件到upload的时候，相同文件名覆盖 ，因此建议最好文件名是不一样， 加上样品名
* input.json 和template.md 必须平级，同时json里用相对路径
* 涂涂吐槽我说明写的少，我就再写点

## 如何添加一个新功能
1.  134和135  type_of_variable 添加一个新类型
2.  279  replace_content 添加一个类型
3.  264  output 增加一个动作。


## reference
* https://c.solargenomics.com//final-report/test.html#1.6%20%20%E5%B0%8F%E5%B7%A5%E5%85%B7

'''
import argparse
import sys
import os
import re
from jinja2 import Environment, FileSystemLoader, StrictUndefined, meta, DebugUndefined
import json
from pytablewriter import MarkdownTableWriter
from pytablewriter.style import Style
import glob
from shutil import copyfile
import zipfile
from PIL import Image
import collections
from pprint import pprint

bindir = os.path.abspath(os.path.dirname(__file__))

default_picture = '{0}/element/no_picture.jpg'.format(bindir)
default_table = '{0}/element/no_table.xls'.format(bindir)
default_download = '{0}/element/download.png'.format(bindir)
default_file = '{0}/element/file.png'.format(bindir)
default_js1 = '{0}/element/jquery.min.js'.format(bindir)
default_js2 = '{0}/element/jquery.albumSlider.min.js'.format(bindir)
default_html = '{0}/element/albumSlider.html'.format(bindir)
default_html_all = '{0}/html'.format(bindir)

__author__ = 'Liu Tao'
__mail__ = 'taoliu@annoroad.com'

pat1 = re.compile('{{(\S*)}}')
pat2 = re.compile('\[(.*?)\]')


class Content():
    coder = {1: 0, 2: 0, 3: 0}

    def __init__(self, line, config, bool_pdf=False):
        self.line = line.lstrip()
        self.downloadable = False
        self.config = config
        self.description = []
        self.note_content = []
        self.out_line = ''
        self.bool_pdf = bool_pdf
        self.add_title_number()

    def add_title_number(self):
        '''
        add title coding number after ### symbols
        '''
        if self.line.startswith('#'):
            count = 0
            for i in self.line:
                if i == '#':
                    count += 1
            if count <= 4:
                header = self.code_header(count)
                self.line = '\n\n' + '#' * count + ' ' + header + self.line[count:]
            else:
                pass
        else:
            pass

    def code_header(self, count):
        '''
        coding title number
        '''
        coder = Content.coder
        if count == 1:
            coder[1] += 1
            coder[2] = 0
            coder[3] = 0
            return ('{0}'.format(coder[1]))
        elif count == 2:
            coder[2] += 1
            coder[3] = 0
            return ('{0}.{1}'.format(coder[1], coder[2]))
        elif count == 3:
            coder[3] += 1
            return ('{0}.{1}.{2}'.format(coder[1], coder[2], coder[3]))
        else:
            pass

    def is_contain_variable(self ):
        '''
        check a line whether contain a variable
        '''
        s1 = pat1.search(self.line)  ## for {{}}
        s2 = pat2.search(self.line)  ## for []
        if s1:
            self.variable = s1.group(1)
            if self.line.startswith('!') :
                if s2:
                    self.title, *self.parameters = s2.group(1).split(',')
                    assert len(self.parameters) in [0, 2], '参数设置的不对，要么是空，要么是2个'
                    if self.line.rstrip().endswith('*'):
                        self.downloadable = True
            return (True)
        elif self.line.startswith('!!!') :
            self.variable = "!!!"
            return (True)
        else:
            self.variable = None
            self.title = None
            return (False)

    def type_of_variable(self):
        if self.variable.startswith('table_'):
            self.type = 'table'
        elif self.variable.startswith('image_'):
            self.type = 'image'
        elif self.variable.startswith('download_'):
            self.type = 'download'
        elif self.variable.startswith('json_'):
            self.type = 'json'
        elif self.variable.startswith('smalltool_'):
            self.type = 'smalltool'
        elif self.variable.startswith('href_'):  # ![title]{{href_}}
            self.type = 'href'
        elif self.variable.startswith('!!!'):
            self.type = 'note'
        else:
            self.type = 'string'

    def output_table(self, num):
        r_str = '\n\n<center>表{0} {1}</center> \n\n{{{{{2}}}}}\n\n'.format(num, self.title, self.variable)
        self.out_line = r_str

    def output_html_image(self, num):
        self.out_line = '''<p align="center"><iframe height="400px" width="300px" \
            src="{{{{{0}}}}}" frameborder=0 allowfullscreen></iframe></p>\n\n<center> \
            图{1} {2}</center>\n\n'''.format(self.variable, num, self.title)

    def output_multi_html_image(self, a_list, num):
        if self.bool_pdf:
            self.out_line = '''<div><img  \
               src="{0}.png" style="display: block;margin: 0 auto;" width="300px" height="400px"  >\
               </div>\n\n<center>图{1} {2}</center>\n\n'''.format(os.path.splitext(a_list[0])[0], num, self.title)
        else:
            for i, j in enumerate(a_list):
                if i == 0:
                    self.out_line = '''<div class="albumSlider"><div class="album"><div class="fullview"> \
            <img src="{0}.png" style=""></div> <div class="slider"> <div class="button movebackward" title="向上滚动"> \
            &lt;&lt;</div><div class="imglistwrap"> <ul class="imglist" style="top: 0;"> \
            <li class="current" html="{0}.html"><img src="{0}.png" alt="{1}"> </li> '''.format(os.path.splitext(j)[0],
                                                                                               self.title)
                else:
                    self.out_line += '''<li class="" html="{0}.html"><img src="{0}.png" alt="{1}"></li>'''.format(
                        os.path.splitext(j)[0], self.title)
            self.out_line += '''</ul> </div> <div class="button moveforward" 
            title="向下滚动">&gt;&gt;</div></div> </div>\n\n<center>图{0} {1}</center>\n\n</div>'''.format(num, self.title)

    def output_multi_image(self, a_list, num, w, h):
        if self.bool_pdf:
            self.out_line = '''<div><img  \
               src="{0}" style="display: block;margin: 0 auto;" width="{4}px" height="{3}px"  >\
               </div>\n\n<center>图{1} {2}</center>\n\n'''.format(a_list[0], num, self.title, h, w)
        else:
            all_path = []
            all_name = []
            for i, j in enumerate(a_list):
                all_path.append(j)
                all_name.append(os.path.basename(j))
            
            self.out_line = '''<an-multiimg
imgs="{0}"
titles="{1}"
index="">
</an-multiimg>
                '''.format(",".join(all_path) ,",".join(all_name))
            self.out_line += '''\n<center>图{0} {1}</center>\n\n</div>'''.format(num, self.title)


    def output_single_image(self, num, w, h):
        self.out_line = '''<div><img  \
        src="{{{{{0}}}}}" style="display: block;margin: 0 auto;" width="{4}px" height="{3}px"  >\
        </div>\n\n<center>图{1} {2}</center>\n\n'''.format(self.variable, num, self.title, h, w)

    def output_download(self, name=False):
        if self.bool_pdf:
            pass
        else:
            if not name: name = self.variable
            # self.out_line += '''<div>&emsp;&emsp;{0}下载链接：</div>\n'''.format(self.title)
            
            self.out_line += '''<an-multidownload 
    titles="" download="{{{{{0}}}}}">
</an-multidownload>'''.format( name )

    def output_json(self, project, keys, group):
        if self.bool_pdf:
            pass
        else:
            self.out_line = '''<div class="fileListBox"><div class="downloadBtn"></div>\
            <span class="leftBtn"></span><div class="fileListContent">\
                <img src="{{_Download_file}}" alt="文件图标"><ul>'''
            for i in keys:
                self.out_line += '<li url="table.html?idCode={0}&sourceKey={1}&group=g{2}"><p>{3}</p></li> '.format(
                    project, i, group, keys[i])
            self.out_line += '</ul></div><span class="rightBtn"></span></div>'

    def output_smalltool(self, a_list):
        if self.bool_pdf:
            pass
        else:
            self.out_line = ''' SmallTool:(输入文件:{{{{{0}}}}}){1},\n '''.format(self.variable, self.title)
            tool_config = self.config['smalltools'][self.title]
            self.out_line = '''<div class='tools'> 
                                  <label class='title'>{1}</label> 
                                  <div class='toolCode' alias="{0[alias]}">
                                  <label>输入文件:</label>
                                  <select class="c-tools" name='{0[input_name]}'>'''.format(tool_config, self.title)
            for i in a_list:
                self.out_line += '<option value="{0}">{0}</option>'.format(i)
            self.out_line += '''</select> <button class="c-jump-btn">自助分析</button></div></div>\n'''.format(self.title)

    def output_href(self):
        self.out_line = '''&emsp;&emsp;<a href="{{{{{0}}}}}" target="_blank">{1}</a>\n '''.format(self.variable,
                                                                                                  self.title)

    def output_note(self):
        if self.bool_pdf:
            for i in self.note_content:
                self.out_line += i
        else:
            self.out_line = '''<div class="notes">'''
            for i in self.note_content:
                self.out_line += i
            self.out_line += '</div>'


    def replace_content(self, inputs_dict, count, bool_local, project):
        a_path = inputs_dict[self.variable]
        if self.type == 'table':
            count['table'] += 1
            inputs_dict[self.variable] = a_path.file_to_table(self.parameters)  ### 读取表格文件，直接替换成文件内容
            self.output_table(count['table'])  ## 添加个 标题
            if self.downloadable:
                a_path.zip_files( inputs_dict )  ### 打包文件，存储zip文件路径
                inputs_dict['{0}_zip'.format(self.variable)] = a_path.add_prefix('zip', project, bool_local)  ### 获得文件的路径，本地就是路径
                self.output_download('{0}_zip'.format(self.variable))  ## 添加一行html
        elif self.type == 'image':
            count['image'] += 1
            files_number = len(a_path.get_all_files())
            if a_path.path.endswith('html') or a_path.path.endswith('svg'):
                if files_number > 0:
                    # if '*' in a_path.path:
                    self.output_multi_html_image(a_path.add_prefix('all', project, bool_local), count['image'])
                else:
                    self.output_html_image(count['image'])
                    inputs_dict[self.variable] = a_path.add_prefix(0, project, bool_local)
            else:
                if files_number > 1:
                    # if '*' in a_path.path:
                    w, h = a_path.adjust(self.parameters)
                    self.output_multi_image(a_path.add_prefix('all', project, bool_local), count['image'], w, h)
                else:
                    inputs_dict[self.variable] = a_path.add_prefix(0, project, bool_local)
                    w, h = a_path.adjust(self.parameters)
                    self.output_single_image(count['image'], w, h)
            ## inputs_dict['{0}_zip'.format(self.variable)] = a_path.file1
            if self.downloadable:
                bool_image = True
                a_path.zip_files( inputs_dict , bool_image)
                inputs_dict['{0}_zip'.format(self.variable)] = a_path.add_prefix('zip', project, bool_local)
                self.output_download('{0}_zip'.format(self.variable))
        elif self.type == 'download':
            a_path.zip_files( inputs_dict )
            inputs_dict['{0}_zip'.format(self.variable)] = a_path.add_prefix('zip', project, bool_local)
            self.output_download('{0}_zip'.format(self.variable))
        elif self.type == 'json':
            keys, group = a_path.to_json(self.description)
            self.output_json(project, keys, group)
            if self.downloadable:
                a_path.zip_files( inputs_dict )
                inputs_dict['{0}_zip'.format(self.variable)] = a_path.add_prefix('zip', project, bool_local)
                self.output_download('{0}_zip'.format(self.variable))
        elif self.type == 'smalltool':
            self.output_smalltool(a_path.add_prefix('st', project, bool_local))
            # inputs_dict['{0}'.format(self.variable)] = a_path.add_prefix('st', project , bool_local)
            if self.downloadable:
                a_path.zip_files( inputs_dict )
                inputs_dict['{0}_zip'.format(self.variable)] = a_path.add_prefix('zip', project, bool_local)
                self.output_download('{0}_zip'.format(self.variable))
        elif self.type == 'href':
            inputs_dict[self.variable] = a_path.add_prefix(0, project, bool_local)
            self.output_href()
        elif self.type=="note":
            self.output_note()
        else:
            self.out_line = self.line


class Pathway():
    sourceKey = 0
    group = 0
    file_mapping = []

    def __init__(self, name, value, module_name, dirname, url, mode):
        self.module = module_name  ### module name 
        self.name = name  ## key
        self.value = value  ## value
        self.dirname = dirname  ### old file dirname
        self.is_file_miss = False
        self.judge_type()  ## set self.type , true is a file ,  then add old dirname into path to access file
        self.add_path()  ## 加上全路径
        self.choose_file()  ## glob all files to get a list with origin directory
        self.parameter = url
        self.mode = mode

    def judge_type(self):
        name = self.name
        if name.startswith('image') or name.startswith('table') or name.startswith('download') \
                                    or name.startswith('json') or name.startswith('smalltool') \
                                    or name.startswith('href'):
            self.type = True
        else: 
            ### 字符串类型
            self.type = False

    def add_path(self):
        '''
        如果是文件，添加全路径 
        '''
        if self.type:
            if self.value.startswith('http'):
                self.path = self.value
            else:
                self.path = '{0}/{1}'.format(self.dirname, self.value)
        else:
            self.path = self.value
        print("generate_md_report path:" + self.path)

    def choose_file(self):
        '''
        选择filelist的第一个作为 file1
        '''
        if '*' in self.path:
            self.filelist = self.get_all_files()
        else:
            self.filelist = [self.path]
        self.check_filelist()
        self.file1 = self.filelist[0]

    def get_all_files(self):
        '''get all file from old path
        根据path获得所有的文件列表，存储到filelist中

        '''
        all_files = sorted(glob.glob(self.path))
        if len(all_files):
            return (all_files)
        else:
            print("no file in {0}".format(self.path))
            return ([])

    def check_filelist(self):
        '''
        检查文件是否存在，如果不存在，就把默认空图空表copy过去，
        '''
        # print(self.filelist)
        if len(self.filelist) == 0:
            self.is_file_miss = True
            if self.name.endswith('zip'):
                pass
            elif self.name.startswith('image'):
                dest_file = "{0}/{1}".format(glob.glob(os.path.dirname(self.path))[0],
                                             os.path.basename(default_picture))
                copyfile(default_picture, dest_file)
                self.filelist.append(dest_file)
            elif self.name.startswith('table') or self.name.startswith('download') or self.name.startswith('json'):
                try:
                    dest_file = "{0}/{1}".format(glob.glob(os.path.dirname(self.path))[0], os.path.basename(default_table))
                    copyfile(default_table, dest_file)
                    self.filelist.append(dest_file)
                except:
                    print("no file in {0}".format(self.path))
                    sys.exit(1)
            else:
                pass
        else:
            for i, j in enumerate(self.filelist):
                if not os.path.isfile(j):
                    self.is_file_miss = True
                    if self.name.endswith('zip'):
                        pass
                    elif self.name.startswith('image'):
                        copyfile(default_picture, self.filelist[i])
                    elif self.name.startswith('table') or self.name.startswith('download') or self.name.startswith(
                            'json'):
                        copyfile(default_table, self.filelist[i])
                    # print(i)
                    else:
                        pass
                        # print('tttttttttt')
                else:
                    pass
                    # print('mmmmmmmmmmm')
        # print(self.filelist)

    def file_to_table(self, parameters, ts=False):
        '''
        控制markdown中表格格式
        output file1 into markdown table format
        '''
        header = []
        content = []
        number = 4  ## 列数
        line_cutoff = 15 ## 行数
        if len(parameters):
            number, line_cutoff = int(parameters[0]), int(parameters[1])
        with open(self.file1) as f_in:
            for count, line in enumerate(f_in):
                tmp = line.rstrip().split('\t')
                if count == 0:
                    header = tmp[:number]
                else:
                    content += [tmp[:number]]
                if count >= line_cutoff - 1: break
        
        #content += [["a"] * len(header)]
        writer = MarkdownTableWriter()
        # writer.table_name = "example_table"
        writer.headers = header
        writer.type_hints = ["str"] * len(header)

        writer.value_matrix = content
        # writer.margin = 1
        for i in range(len(header)):
            writer.set_style(i, Style(align="center"))
        #out_table_md = "\n".join(writer.dumps().split('\n')[:-2])
        print( writer.dumps())
        return ('\n' + writer.dumps() + '\n')

    def add_prefix(self, position, project, local):
        '''
        position: which file , 0- single file; all: modify all files ; st: modify samtools files ;other: return zip file
        project :
        local : used for local or url
        '''
        # url = 'https://source.solargenomics.com/user/report'
        url = ''
        if self.mode == 'pipeline':
            url = self.parameter[0]
        else:
            url = "/".join(self.parameter[:2])

        f_modify = lambda p: p[p.find('upload', ):] if 'upload' in p else p

        if position == 0:
            #print('************',self.file1)
            if local:
                return (self.file1)
            else:
                if self.file1.startswith('http'):
                    return (self.file1)
                else:
                    if self.mode == 'pipeline':
                        return ('{0}/{1}/{2}'.format(url, project, f_modify(self.file1)))
                    else:
                        return ('{0}/{1}'.format(url, f_modify(self.file1)))
                        # self.newfilelist.append('{0}/{1}'.format(url , f_modify(i)))
        elif position == 'all':
            self.newfilelist = []
            for i in self.filelist:
                if local:
                    self.newfilelist.append(i)
                else:
                    if self.mode == 'pipeline':
                        self.newfilelist.append('{0}/{1}/{2}'.format(url, project, f_modify(i)))
                    else:
                        self.newfilelist.append('{0}/{1}'.format(url, f_modify(i)))
            return (self.newfilelist)
        elif position == 'st':
            self.newfilelist = []
            # print(self.filelist)
            for i in self.filelist:
                if local:
                    self.newfilelist.append(i)
                else:
                    if self.mode == 'pipeline':
                        self.newfilelist.append("结题报告/{0}/{1}".format(project, f_modify(i).replace('./upload/', '')))
                        # self.newfilelist.append(  "{0}".format( f_modify(i).replace('./upload/' , '')))

                    else:
                        self.newfilelist.append("{0}/{1}".format(url, f_modify(i).replace('./upload/', '')))
            return (self.newfilelist)
        else:  ####for zip file
            if local:
                return (self.zip_path)
            else:
                if self.mode == 'pipeline':
                    return ('{0}/{1}/{2}'.format(url, project, f_modify(self.zip_path)))
                else:
                    return ('{0}/{1}'.format(url, f_modify(self.zip_path)))

    def adjust(self, parameters):
        '''
        adjust image size and keep w/h raito unchange
        '''
        with Image.open(self.file1) as img:
            if len(parameters):
                return (parameters[0], parameters[1])
            else:
                width, height = img.size
                if width / height > 4 / 3:
                    return (400, height * 400 / width)
                else:
                    return (width * 400 / height, 400)

    def prepare_directory(self, outdir, bool_upload):
        """
        copy files into upload file , and modify new file path
        """
        if bool_upload and self.type:
            new_file_path = '{0}/{1}'.format(outdir, os.path.basename(self.path))
            self.copy(new_file_path)
            self.path = new_file_path
        else:
            pass

    def zip_files(self, input_dict , bool_image=False ):
        '''
        zip files and return zip file path
        '''
        output_dir = os.path.dirname(self.path)
        if '*' in output_dir:
            output_dir_tmp = glob.glob(output_dir)
            if len(output_dir_tmp)==1:
                output_dir = output_dir_tmp[0]
            else:
                output_dir = os.path.dirname(output_dir)
            #output_dir = glob.glob(output_dir)[0]
        a_key = "{0}_zip".format(self.name)
        if a_key in input_dict:
            out_zip = input_dict[a_key].value
        else: 
            out_zip = '{0}/{1}.zip'.format(output_dir, self.name)

        with zipfile.ZipFile(out_zip, 'w', zipfile.ZIP_DEFLATED) as zf:
            for i in self.filelist:
                if bool_image:
                    for _atype in ['png', 'pdf']:
                        for j in glob.glob('{0}*{1}'.format(i.rsplit('.', 1)[0], _atype)):
                            zf.write(j, os.path.basename(j))
                else:
                    zf.write(i, os.path.basename(i))
        zf.close()
        self.zip_path = out_zip
        return (self.zip_path)

    def copy(self, dest):
        '''
        copy file into upload dirctory
        '''
        dirname = os.path.dirname(os.path.abspath(dest))
        r_list = []
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        for i in self.filelist:
            # print(i , dest)
            dest = '{0}/{1}'.format(os.path.dirname(dest), os.path.basename(i))
            copyfile(i, dest)
            r_list.append(dest)
        self.filelist = r_list
        self.file1 = r_list[0]

    def to_json(self, description):
        ''' dumps table into json file with title ,data , label ,descrition
        '''
        tt = {}
        Pathway.group += 1
        print(self.name, self.filelist)
        for i in self.filelist:
            pp = {"title": [], "data": [], "label": [], "description": description}
            header = collections.OrderedDict()
            with open(i) as fin:
                for m, n in enumerate(fin):
                    tmp = n.rstrip().split('\t')
                    if m == 0:
                        header = {a.replace('.', '!'): "int" for a in tmp}
                    else:
                        pp["data"].append(tmp)
                        header = self.modify_type(header, tmp)
            out_json = i.rsplit('.', 1)[0] + '.json'
            self.write_to_json(out_json, pp, header)
            Pathway.sourceKey += 1
            tt[Pathway.sourceKey] = os.path.basename(out_json)
            Pathway.file_mapping.append({"source": out_json.replace('./upload/', '', 1),
                                         "sourceKey": "{0}".format(Pathway.sourceKey),
                                         "group": "g{0}".format(Pathway.group)})
        return (tt, Pathway.group)

    def isdecimal(self, value):
        try:
            float(value)
            return (True)
        except ValueError:
            return (False)

    def modify_type(self, header, tmp):
        '''
        Determine the type of title
        '''
        header_list = list(header.keys())
        for i, j in enumerate(tmp):
            a_key = header_list[i]
            if header[a_key] == 'int':
                if j.isnumeric():
                    pass
                elif self.isdecimal(j):
                    header[a_key] = 'double'
                else:
                    header[a_key] = 'string'
            elif header[a_key] == 'double':
                if not (j.isnumeric() or self.isdecimal(j)):
                    header[a_key] = 'string'
            else:
                pass
        return (header)

    def write_to_json(self, out_json, pp, header):
        '''
        dumps json
        '''
        for i in header:
            if header[i] == 'string':
                pp["title"].append({"ts": "false", "name": i, "type": header[i]})
            elif header[i] == 'double':
                pp["title"].append({"ts": "false", "name": i, "type": header[i]})
                pp["label"].append(i.replace('!', '.'))
            else:
                pp["title"].append({"ts": "true", "name": i, "type": header[i]})
                pp["label"].append(i.replace('!', '.'))
        with open(out_json, 'w') as f_out:
            f_out.write(json.dumps(pp, indent=4))


def get_md_json(indir, bool_upload, output_dir, url, mode):
    '''
    目的：读取每个目录下的input.json 和 template.md文件，并且按section进行读取，存到对应的字典中
    
    indir: args.indir 
    bool_upload:
    output_dir  : upload directory
    r_dict->module_name -> json -> key -> pathway object
                        -> md -> md contents list
    input.json
    {
    "GO_clusterProfiler": {
        "table_go_report": "upload/GO/*/go.example.report.xls",
        "download_go_report": "upload/GO/*/*go.report.xls"
    },
    "clusterProfiler":{},
    "GO_enrich":{
        "table_go_enrich":"upload/GO/*/*enrichment.xls",
        "image_go_bar": "upload/GO/*/*barplot*png",
        "image_go_dot": "upload/GO/*/*dotplot*png"
    }
    }

    '''
    r_dict = {}
    list_of_files = []
    module_sequences = []
    miss_file_module = []
    for root, dirs, files in os.walk(indir):
        list_of_files += [os.path.join(root, afile) for afile in files]
    for afile in list_of_files:
        if afile.endswith('input.json'):
            with open(afile) as fin:
                json_load = json.load(fin)
                if len(json_load):
                    for name in json_load:  ### name:模块名
                        json_obj_dict, miss_file_name = json_modify(name, 
                                                                    json_load[name], 
                                                                    os.path.dirname(afile),
                                                                    output_dir, 
                                                                    bool_upload, 
                                                                    url, 
                                                                    mode)
                        if not name in r_dict:
                            r_dict[name] = {}
                            r_dict[name]['json'] = json_obj_dict
                        miss_file_module += miss_file_name
                else:
                    raise Exception(" no target in {0} ".format(afile))

            ### 处理md文件
            md_file = afile.replace('input.json', 'template.md')
            if md_file in list_of_files:
                with open(md_file) as infile:
                    for i in infile:
                        if i.startswith('@@@@'):  ## @@@@filter
                            name = i.replace('@', '').rstrip()
                            if not name in module_sequences:
                                module_sequences.append(name)
                            else:
                                raise Exception("{0} is used repeatly in {1}".format(name, md_file))
                            if not name in r_dict:
                                raise Exception(' section {0} in {1} is not existed '.format(name, md_file))
                            r_dict[name]['md'] = []
                        else:
                            r_dict[name]['md'].append(i)
            else:
                raise Exception('Error:{0} is not found'.format(md_file))
    print(miss_file_module)
    ### 
    return (r_dict, module_sequences, miss_file_module)


def json_modify(name, json_dict, dirname, outdir, bool_upload, url, mode):
    '''
    name: module name
    outdir: basename + upload
    dirname:  old file dirname, from work directory to json file 
    json_dict : json load dict
    bool_upload:  whether prepare file into upload
    '''
    # print(outdir , bool_upload)
    r_dict = {}
    miss_file_module = []
    for i in json_dict:
        ## 判断类型，添加全路径，glob所有文件
        a_path = Pathway(i, json_dict[i], name, dirname, url, mode)  ### dirname is old directory
        a_path.prepare_directory('{0}/{1}'.format(outdir, name), bool_upload)   ## 把文件拷贝到新目录，同时修改pathlist，path，file1
        if a_path.is_file_miss:
            miss_file_module.append(a_path.name)
        # print('****',a_path.path)
        r_dict[a_path.name] = a_path
    return (r_dict, miss_file_module)


def replace_template(md_json_dict, count, bool_local, project, bool_pdf, config, icon_dir, function_name):
    '''
    template_list: a list of template file content 
    inputs_dict : a list of pathway object
    bool_local : 
    project : project name 
    '''
    r_list = []
    bool_json = False  ### used to store line for json table
    bool_note = False
    a_line = ''
    pat3 = re.compile('^\s*$')
    template_list = md_json_dict['md']  ### md文件得list
    inputs_dict = md_json_dict['json']
    inputs_dict["!!!"] = "" 
    if bool_local:
        inputs_dict['_Download_png'] = '{0}/download.png'.format(icon_dir)
        inputs_dict['_Download_file'] = '{0}/file.png'.format(icon_dir)
    else:
        inputs_dict['_Download_png'] = "https://source.solargenomics.com/source/report/download.png"
        inputs_dict['_Download_file'] = "https://source.solargenomics.com/source/report/file.png"
    for i in template_list:
        if bool_json:
            if re.search(pat3, i): continue
            if i.lstrip().startswith('@'):  ## 对json表格进行注释
                a_line.description.append(i.lstrip()[1:])
                continue
            else:
                bool_json = False
                modify_content(a_line, inputs_dict, count, bool_local, project, function_name)
                r_list.append(a_line.out_line)
        if bool_note:
            if re.search(pat3, i): continue
            if i.lstrip().startswith('!!!'):  ## 对json表格进行注释
                a_line.note_content.append(i.lstrip()[3:])
                continue
            else:
                bool_note = False
                modify_content(a_line, inputs_dict, count, bool_local, project, function_name)
                r_list.append(a_line.out_line)
    
        a_line = Content(i, config, bool_pdf)  ## 内容  配置文件 是否pdf

        if a_line.is_contain_variable():  ### 获得变量的名称 以及title，和属性
            a_line.type_of_variable()  ### 添加该行的类型
            if a_line.type == 'json':  
                bool_json = True
                continue
            if a_line.type == 'note':  
                a_line.note_content.append(i.lstrip()[3:])
                bool_note = True
                continue
            modify_content(a_line, inputs_dict, count, bool_local, project, function_name )
            if a_line.out_line[0].isalnum():
                r_list.append('\n\n&emsp;&emsp;' + a_line.out_line)
            else:
                r_list.append(a_line.out_line)
        else:
            if a_line.line:
                if a_line.line[0].isalnum():
                    r_list.append('\n\n&emsp;&emsp;' + a_line.line)
                else:
                    r_list.append(a_line.line)
            else:
                r_list.append(a_line.line)
    if bool_json:
        modify_content(a_line, inputs_dict, count, bool_local, project, function_name)
        r_list.append(a_line.out_line)
    if bool_note:
        modify_content(a_line, inputs_dict, count, bool_local, project, function_name)
        r_list.append(a_line.out_line)
    return r_list


def modify_content(a_line, inputs_dict, count, bool_local, project, function_name):
    print(inputs_dict)
    if a_line.variable in inputs_dict:  
        a_line.replace_content(inputs_dict, count, bool_local, project)
    else:
        raise Exception('{0} is not defined in {1}'.format(a_line.variable,  function_name))


def find_defined(infile):
    THIS_DIR = os.path.dirname(os.path.abspath(infile))
    env = Environment(loader=FileSystemLoader(THIS_DIR), trim_blocks=True, undefined=StrictUndefined)
    template_source = env.loader.get_source(env, os.path.basename(infile))[0]
    parsed_content = env.parse(template_source)
    a_set = meta.find_undeclared_variables(parsed_content)
    return (a_set)


def dumps(a_dict):
    '''
    output json from Pathway instance and string 
    '''
    r_dict = {}
    for i in a_dict:
        if isinstance(a_dict[i], str):
            r_dict[i] = a_dict[i]
        else:
            r_dict[i] = a_dict[i].path
    return (r_dict)


def range_module(fin, alist):
    r_list = []
    pat9 = re.compile('^\s*$')
    if fin:
        for line in fin:
            if line.startswith('#') or re.search(pat9, line): continue
            tmp = line.rstrip()
            r_list.append(tmp)
    else:
        r_list = alist
    return (r_list)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog='author:\t{0}\nmail:\t{1}'.format(__author__, __mail__))
    parser.add_argument('-d', '--indir', help='input directory ', dest='indir', required=True)
    parser.add_argument('-o', '--output', help='output file', dest='output', required=True)
    parser.add_argument('-m', '--module', help='module list', dest='module', type=open)
    parser.add_argument('-pdf', '--pdf', help='output pdf file', dest='pdf', action='store_true')
    parser.add_argument('-l', '--local', help='output local md file , or using cloud url insted of local path', dest='local', action='store_true')
    parser.add_argument('-u', '--upload', help='whether prepare upload directory', dest='upload', action='store_true')
    parser.add_argument('-pipeline', '--pipeline', help='parameters for pipeline , [project_id] ', dest='pipeline', nargs=1)
    parser.add_argument('-smalltool', '--smalltool', help='parameters for smalltool , [path , st_id]', dest='smalltool', nargs=2)
    args = parser.parse_args()

    mode = ''
    parameter = []
    if args.pipeline:
        mode = 'pipeline'
        parameter = ['https://source.solargenomics.com/user/report'] + args.pipeline  ### url , project
    elif args.smalltool:
        mode = 'smalltool'
        parameter = ['https://annoroad-cloud-product.oss-cn-beijing.aliyuncs.com'] + args.smalltool  ### url , path , smalltool id
    else:
        sys.exit()

    if args.output.startswith(r'/'): 
        sys.exit('Error: please -o use relative pathway')
    final_output = args.output
    with open('{0}/config.json'.format(bindir)) as fin:
        config = json.load(fin)

    if not args.output.startswith(r'.'): final_output = './' + args.output
    output_dir = os.path.dirname(final_output)
    if not os.path.exists(output_dir + '/upload'): os.makedirs(output_dir + '/upload')
    if not os.path.exists(output_dir + '/upload/_icon'):
        os.makedirs(output_dir + '/upload/_icon')
    if os.path.exists(output_dir + '/upload/_icon'):
        copyfile(default_download, output_dir + '/upload/_icon/download.png')
        copyfile(default_picture, output_dir + '/upload/_icon/file.png')
        copyfile(default_js1, output_dir + '/upload/_icon/jquery.min.js')
        copyfile(default_js2, output_dir + '/upload/_icon/jquery.albumSlider.min.js')
        copyfile(default_html, output_dir + '/upload/_icon/albumSlider.html')
    #print("cp -rf {} {}".format(default_html_all, output_dir))
    os.system("cp -rf {} {}".format(default_html_all, output_dir))

    ## iteration on all directory
    all_md_json, module_sequences, miss_file_module = get_md_json(args.indir,    ## 输入目录
                                                                  args.upload,   ## 是否整理新目录
                                                                  output_dir + '/upload',  ## 输出路径
                                                                  parameter,   ## 相关参数
                                                                  mode)  ## 流程 还是 小工具
    # sys.exit()
    ### 如果提供了模块顺序文件，那么报告按照提供的文件来
    function_list = range_module(args.module, module_sequences)  

    count = {'table': 0, 'image': 0}
    ### 输出报告文件
    f_report_out = open(final_output, 'w')

    #r_dict->module_name -> json -> key -> pathway object
    #                    -> md -> md contents list

    for name in function_list:
        new_template_list = replace_template(all_md_json[name], 
                                              count, 
                                              args.local,  ## 是否是本地版
                                              parameter[-1], ### 流程还是小工具
                                              args.pdf,  ##  是否是pdf
                                              config,  ## 一些配置
                                              output_dir + '/upload/_icon', name)

        with open('{0}/tmp_template.tmp'.format(output_dir), 'w') as f_out:
            for i in new_template_list:
                f_out.write(i)

        inputs_dict = dumps(all_md_json[name]['json'])
        #pprint(inputs_dict)

        env = Environment(loader=FileSystemLoader(output_dir), trim_blocks=True, undefined=DebugUndefined)
        template = env.get_template('tmp_template.tmp')
        f_report_out.write(template.render(**inputs_dict).replace("\n![", "\n\n!["))

    Content.coder[2] += 1
    if not args.pdf:
        #f_report_out.write('\n\n\n### {0}.{1} 结果目录\n'.format(Content.coder[1], Content.coder[2]))
        #f_report_out.write(config['dir_tree'].format(parameter[-1]))
        f_report_out.write('\n\n{0}'.format(config['company']))
    f_report_out.close()


    ## 给json table用的
    with open('{0}/mapping.json'.format(output_dir), 'w') as f_out:
        f_out.write(json.dumps({'file_mapping': Pathway.file_mapping}, indent=4))

    print("\n\nstep2 : check result ")
    a_set = find_defined(final_output)
    if len(a_set):
        print("[Error]: undefined variable is \n{0}\n".format("\n".join(a_set)))


if __name__ == '__main__':
    main()
