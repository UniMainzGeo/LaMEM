import sys

def py23input(prompt) :
  if sys.version_info[0] == 2 :
    v = raw_input(prompt)
  else :
    v = input(prompt)
  return(v)
