#!/usr/bin/env python3
'''
Information:
        this script is method used for Pipline.
Function:
        1) generateShell
        2) decorateShell
        3) myconf (read ini)
        4) mkdir
usage:
        sys.path.append(os.path.dirname(sys.argv[0]) + '/../lib')
        from PipMethod import mkdir
        from PipMethod import generateShell
        from PipMethod import myconf 
Modify Date:
        2015-09-26
'''
import argparse
import sys
import os
import re
import glob
import configparser
import time

#from multiprocessing import Pool

__author__='Yuan Zan'
__mail__= 'zanyuan@annoroad.com'


def generateShell(shell, content, finish_string="Live_long_and_prosper"):
    shell = str(shell)
    for file in glob.glob(shell + '.*'):
        os.remove(file)
    f=open(shell,'w')
    f.write('#!/bin/bash\n')
    #f.write('echo ==========start at : `date` ==========\n')
    #f.write(content + ' && ' + '\\\n')
    f.write(content + '\n')
    #f.write('echo ==========end at : `date` ========== && \\\n')
    #f.write('echo ' + finish_string + ' 1>&2 && \\\n')
    #f.write('echo ' + finish_string + ' > ' + shell + '.sign\n')
    f.close()

def decorateShell(shell, finish_string="Live_long_and_prosper"):
    shell = str(shell)
    for file in glob.glob(shell + '.*'):
        os.remove(file)
    cmd = 'cat ' + shell
    content = os.popen(cmd).read().rstrip()
    f=open(shell,'w')
    f.write('#!/bin/bash\n')
    f.write('echo ==========start at : `date` ==========\n')
    f.write(content + ' && ' + '\\\n')
    f.write('echo ==========end at : `date` ========== && \\\n')
    f.write('echo ' + finish_string + ' 1>&2 && \\\n')
    f.write('echo ' + finish_string + ' > ' + shell + '.sign\n')
    f.close()

class myconf(configparser.ConfigParser):
    def __init__(self,defaults=None):
        configparser.ConfigParser.__init__(self,defaults=None,allow_no_value=True)
    def optionxform(self, optionstr):
        return optionstr
def mkdir(inDirs):
    for i in inDirs:
        if os.path.exists(i) == False:
            os.makedirs(i)

def main():
    pass

if __name__ == '__main__':
        main()
