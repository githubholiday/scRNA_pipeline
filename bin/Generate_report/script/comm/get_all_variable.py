#! /usr/bin/env python3
import argparse
import sys
import os
import re
from jinja2 import Environment, FileSystemLoader,StrictUndefined,meta
import json
  
bindir = os.path.abspath(os.path.dirname(__file__))

__author__='Liu Tao'
__mail__= 'taoliu@annoroad.com'

pat1=re.compile('^\s*$')

def find_defined(infile):
    THIS_DIR = os.path.dirname(os.path.abspath(infile))
    env = Environment(loader=FileSystemLoader(THIS_DIR),trim_blocks=True , undefined=StrictUndefined)
    template_source = env.loader.get_source(env,   infile)[0]
    parsed_content = env.parse(template_source)
    a_set = meta.find_undeclared_variables(parsed_content) 
    return(a_set )

def merge_file(file_list , outdir ,name ):
    with open('{0}/{1}'.format(outdir , name), 'w') as outfile:
        for fname in file_list:
            with open(fname) as infile:
                outfile.write(infile.read())
    return('{0}/{1}'.format(outdir , name))

def main():
    parser=argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
    parser.add_argument('-t','--template',help='template file',dest='template',required=True)
    args=parser.parse_args()

    a_set = find_defined( args.template )
    tt_dict = dict.fromkeys(a_set ,'')
    print(json.dumps(tt_dict,sort_keys=True, indent=4, separators=(',', ':')))

if __name__ == '__main__':
    main()
