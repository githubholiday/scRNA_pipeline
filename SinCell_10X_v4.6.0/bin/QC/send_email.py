#!/usr/bin/env python3
'''

'''
import os
import sys
import getopt
import myemail

def usage():
	print ("Usage: python3 " + sys.argv[0] + " <email_config.ini>" )
	sys.exit()
def main():
	if len(sys.argv[1:]) !=1 :
		usage()
	config = os.path.abspath(sys.argv[1])
	email = myemail.Email([config])
	email.send_email()
if __name__ == "__main__":
	main()
