import sys
import datetime

LOG_FILE = sys.stdout

def set_log_file(file_name):
  global LOG_FILE
  if file_name is not None:
    LOG_FILE = open(file_name, 'w')
  else:
    LOG_FILE = None

def log(s=''):
  if LOG_FILE is not None:
    LOG_FILE.write(datetime.datetime.now().strftime("%H:%M:%S: ") + str(s) + '\n')

def new_line():
  LOG_FILE.write('\n')