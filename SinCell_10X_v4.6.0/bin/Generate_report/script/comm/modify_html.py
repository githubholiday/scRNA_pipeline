#! /usr/bin/env python3
import argparse
import sys
import os
import shutil
bindir = os.path.abspath(os.path.dirname(__file__))
__author__ = 'Liu Tao'
__mail__ = 'taoliu@annoroad.com'

def copy(source, dest):
    os.system(f"cp -r {source} {dest}")

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     epilog='author:\t{0}\nmail:\t{1}'.format(__author__, __mail__))
    parser.add_argument('-i', '--input', help='input html ', dest='input', required=True)
    parser.add_argument('-o', '--output', help='output file', dest='output', required=True)
    args = parser.parse_args()

    copy(f"{bindir}/html" , "{0}".format(os.path.abspath(os.path.dirname(args.output))))
    with open(args.input) as fin:
        with open(args.output, 'w') as fout:
            for line in fin:
                content = line.rstrip()
                if "</title>" in content:
                    tt = '''{0} <link href="./html/css/jquery.tocify.css" rel="stylesheet" type="text/css"/>
  <link href="./html/css/datatables.min.css" rel="stylesheet" type="text/css"/>
  <link href="./html/css/markdown.css"type="text/css" rel="stylesheet">\n'''.format(line)
                elif "</body>" in content:
                    tt = '''{0} <script src="./html/js/jquery.min.js"></script>
<script src="./html/js/jquery-ui.min.js"></script>
<script src="./html/js/jquery.tocify.min.js"></script>
<script src="./html/js/datatables.min.js"></script>
<script src="./html/js/render.js"></script>
'''.format(line)

                else:
                    tt = line

                fout.write(tt)








    




if __name__ == '__main__':
    main()
