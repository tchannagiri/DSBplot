import sys
import datetime

  
def log(s):
  print(datetime.datetime.now().strftime("%H:%M:%S: ") + str(s))

def log_input(s):
  log('(input) ' + str(s))

def log_output(s):
  log('(output) ' + str(s))

def new_line():
  print()