#! /usr/bin/env python3
# -*- coding: utf-8 -*-  
import argparse
import sys
import os
import re
import paramiko
import sqlite3
import getpass
bindir = os.path.abspath(os.path.dirname(__file__))
import myemail

__author__='Liu Tao'
__mail__= 'taoliu@annoroad.com'

pat1=re.compile('^\s+$')
IP='192.168.60.188'
IP2 = '192.168.13.20'
PORT=50732
PYTHON='/usr/local/bin/python3'

def check_jump_success(user):
	try:
		ssh = paramiko.SSHClient()
		ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
		ssh.connect(IP , PORT , user )
		ssh.close()
		return(True)
	except Exception  as e :
		ssh.close()
		print('Reason:', e)
		return (False)

def exists_email_config(user):
	try:
		ssh = paramiko.SSHClient()
		ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
		ssh.connect(IP , PORT , user )
		cmd = "test -f ~/.email/.email.txt && echo ok"
		flag = myrun(ssh , cmd)
		ssh.close()
		return flag
	except Exception  as e :
		ssh.close()
		print('Reason:', e)
		return False

def myrun(ssh , cmd):
	print(cmd)
	stdin , stdout , stderr = ssh.exec_command(cmd , timeout = 120)
	if stdout.readlines()[-1].rstrip().endswith('ok'):
		return True
	else:
		stdin , stdout , stderr = ssh.exec_command(cmd)
		print(stdout.readlines())
		return False
	
def scp_dir( user , upload_dir ):
	try:
		ssh = paramiko.SSHClient()
		ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
		ssh.connect(IP , PORT , user )
		stdin , stdout , stderr = ssh.exec_command('echo $$')
		# pid = stdout.readline().rstrip()
		# if not myrun( ssh , "mkdir -p ~/{0} && echo ok".format(pid) ):
			# sys.exit("目录创建失败")
		
		if not myrun(ssh , "scp -r {0}@{1}:{2} /home/{0}  && echo ok ".format(user , IP2, upload_dir)):
			sys.exit("拷贝失败")
		if not myrun(ssh, "{0} /home/{1}/{2}/{3} /home/{1}/{2}/{4} && echo ok".format(PYTHON, user, os.path.basename(upload_dir), "send_email.py", "email_config.ini" )):
			sys.exit("发送邮件失败")
		if not myrun( ssh , "rm -rf /home/{0}/{1} && echo ok".format(user,os.path.basename(upload_dir) ) ):
			sys.exit("目录删除失败")
			
		# if full_flag:
			# if not myrun( ssh  , "{0}  {1} -i ~/{3}/{2} -f   && echo ok ".format(PYTHON , SCRIPT , os.path.basename(old_db) ,pid)):
				# sys.exit("发送邮件失败")
		# else:
			# if not myrun( ssh ,  "{0}  {1} -i ~/{3}/{2}  && echo ok ".format(PYTHON , SCRIPT , os.path.basename(old_db) ,pid)):
				# sys.exit("发送邮件失败")
		#if not myrun( ssh , "rm -rf ~/{0} && echo ok".format(pid) ):
		#	sys.exit("目录删除失败")
		ssh.close()
	except Exception  as e :
		ssh.close()
		print('Reason:', e)
		return False
def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--input',help='input file',dest='input',required=True)
	#parser.add_argument('-full','--full',help='directory full ',dest='full', action='store_true')
	#parser.add_argument('-o','--output',help='output file',dest='output',type=argparse.FileType('w'),required=True)
	#parser.add_argument('-m','--mm',help='output file',dest='mm',action='store_false')
	args=parser.parse_args()
	user = getpass.getuser()
	#user = 'sci-qc'
	if not check_jump_success(user):
		sys.exit("免密跳转失败，请配置")
	else:
		print('免密跳转成功')
	# if not exists_email_config(user):
		# print('''please config your email address first :
# [DEFAULT]
# Max_count = 5 #发送失败时重复发送次数
# Sleep_time = 60 #发送失败时睡眠时间

# [HEADER]
# Addressor = **@genome.cm #邮件发送者
# Password = **** #邮箱密码
# Receiver = chengfangtu@annoroad.com; #接收者
# Copy = yanxunsu@annoroad.com #抄送者
# Server = smtp.exmail.qq.com #登陆服务器
# Receive_server = pop.exmail.qq.com
# ''')
		# sys.exit()
	# else:
		# print("发信配置成功")
	scp_dir(user , os.path.abspath(args.input))

if __name__ == '__main__':
	main()
