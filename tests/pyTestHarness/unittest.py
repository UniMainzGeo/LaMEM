
import os
import numpy as np
import math as math
import re

from pyTestHarness.colors import pthNamedColors as bcolors

def compareLiteral(input,expected):
  status = True
  err = ''
  
  if len(input) != len(expected):
    status = False
    err = err + "compareLiteral [failed]\nReason: input and expected are of different length\n"
    err = err + ("  expected: %s\n" % expected)
    err = err + ("  input:    %s\n" % input)
    return status,err

  for index in range(0,len(expected)):
    if input[index] != expected[index]:
      status = False
      err = err + "compareLiteral [failed]\nReason: strings are different\n"
      err = err + ("  expected: %s\n" % expected)
      err = err + ("  input:    %s\n" % input)
      err = err + "  index[" + str(index) +  "]" + " input \"" +  input[index] + "\" != expected \"" + expected[index] + "\"\n"

  return status,err

def compareFloatingPoint(input,tolerance,expected):
  status = True
  err = ''

  tmp = np.array(input)
  i_f = tmp.astype(np.float)
  tmp = np.array(expected)
  e_f = tmp.astype(np.float)
  tol_f = float(tolerance)

  if len(input) != len(expected):
    status = False
    err = err + "compareFloatingPoint [failed]\nReason: input and expected are of different length\n"
    err = err + ("  expected: %s\n" % e_f)
    err = err + ("  input:    %s\n" % i_f)
    return status,err

  for index in range(0,len(e_f)):
    absdiff = np.abs(i_f[index] - e_f[index]);
    if absdiff > tol_f:
      status = False
      err = err + "compareFloatingPoint [failed]\nReason: tolerance " + ("%1.4e" % tol_f) + " not satisifed\n"
      err = err + ("  expected: %s\n" % e_f)
      err = err + ("  input:    %s\n" % i_f)
      err = err + "  index[" + str(index) + "]" + (" input \"%1.6e\"" %  i_f[index])  + (" != expected \"%1.6e\"" % e_f[index]) + " (+/-" + ("%1.4e" % tol_f)+")\n"

  return status,err

def compareInteger(input,tolerance,expected):
  status = True
  err = ''

  tmp = np.array(input)
  i_i = tmp.astype(np.int)
  tmp = np.array(expected)
  e_i = tmp.astype(np.int)
  tol_i = int(tolerance)

  if len(input) != len(expected):
    status = False
    err = err + "compareInteger [failed]\nReason: input and expected are of different length\n"
    err = err + ("  expected: %s\n" % e_i)
    err = err + ("  input:    %s\n" % i_i)
    return status,err
  
  for index in range(0,len(e_i)):
    absdiff = np.abs(i_i[index] - e_i[index]);
    if absdiff > tol_i:
      status = False
      err = err + "compareInteger [failed]\nReason: tolerance " + str(tol_i) + " not satisifed\n"
      err = err + ("  expected: %s\n" % e_i)
      err = err + ("  input:    %s\n" % i_i)
      err = err + "  index[" + str(index) + "]" + (" input \"%s\"" %  i_i[index]) + (" != expected \"%s\"" % e_i[index]) + " (+/-"+str(tol_i)+")\n"

  return status,err


def parseFile(filename,keywords):

  flat = ''
  contents = []
  if filename:
    file = open(filename,"r")
    for line in file:
      if line.rstrip():
        if not any(keywords in line for keywords in keywords):
          rm_lb = line.lstrip()
          rm_lb = rm_lb.rstrip()
          flat = flat + (rm_lb + ' ')

          stripped_line = line.rstrip()
          stripped_line = stripped_line.lstrip()
          contents.append(stripped_line)

    file.close()
  return(contents,flat)


def getKeyValues(contents,keyword):

  c1 = contents
  f_c = []
  f_r = []
  f_s = []

  # Look for arrays, then strip them out
  #f1 = re.findall("\s*" + keyword + "\s*=\s*[\[\(\{].*?[\}\)\]]", c1)
  find_curly  = re.findall("\s*" + keyword + "\s*=?\s*[\{].*?[\}]", c1)
  if find_curly:
    #print('find_curly',find_curly)
    f_c  = re.findall("\s*" + keyword + "\s*=?\s*[\{](.*?)[\}]", c1)
    #print(f_c)
    for item in find_curly:
      c1 = c1.replace(item,' ')
  #print('c1 - curly:',c1)

  find_round  = re.findall("\s*" + keyword + "\s*=?\s*[\(].*?[\)]", c1)
  if find_round:
    #print('find_round',find_round)
    f_r  = re.findall("\s*" + keyword + "\s*=?\s*[\(](.*?)[\)]", c1)
    #print(f_r)
    for item in find_round:
      c1 = c1.replace(item,' ')
  #print('c1 - round:',c1)

  find_square  = re.findall("\s*" + keyword + "\s*=?\s*[\[].*?[\]]", c1)
  if find_square:
    #print('find_square',find_square)
    f_s  = re.findall("\s*" + keyword + "\s*=?\s*[\[](.*?)[\]]", c1)
    #print(f_s)
    for item in find_square:
      c1 = c1.replace(item,' ')
  #print('c1 - square:',c1)


  # Look for singletons, then strip them out
  #f3 = re.findall("\s*" + keyword + "\s*=\s*([0-9a-zA-Z-\_\.]*)", c1)
  #f4 = re.findall("\s*" + keyword + "\s+([0-9a-zA-Z-\_\.]*)", c1)

  find_single = []
  #find_single = re.findall("\s*" + keyword + "\s*(.*\s)", c1)
  f_single = re.findall("\s*" + keyword + "\s*=?\s*(.*?)\s", c1)
  #print('f_single',f_single)


  list = f_c + f_r + f_s + f_single

  filtered = []
  if list:
    for item in list:
      if item != '':
        trimmed = item.lstrip()
        trimmed = trimmed.rstrip()
        filtered.append(trimmed)

  return filtered


def getKeyValuesAsInt(contents,keyword):
  result = getKeyValues(contents,keyword)

  #print('result = ',result)
  flattened = []
  for sublist in result:

    if type(sublist) is list:
        pass
    else:
      if sublist != '':
        numlist = sublist.replace(',',' ')
        numlist = numlist.replace(';',' ')
        nums = numlist.split(' ')
        #print(nums)
        for num in nums:
          if num != '':
            flattened.append(num)
    

  #print('flat',flattened)
  tmp = np.array(flattened)
  values = tmp.astype(np.int)
  return values


def getKeyValuesAsFloat(contents,keyword):
  #print('DEBUG: contents = ',contents)

  result = getKeyValues(contents,keyword)

  #print('DEBUG: result = ',result)

  flattened = []
  for r in result:
    r = r.replace(',',' ')
    r = r.replace(';',' ')
    val = r.split(' ')
    for v in val:
      if v != '':
        flattened.append(v)

  #print('DEBUG: flattened = ',flattened)

  tmp = np.array(flattened)
  values = tmp.astype(np.float)
  return values

def getKeyValuesNLinesInclusive(contents,keyword,numlines):
  result = []
  counter = 0
  start = 0
  for line in contents:
    x = line.find(keyword)
    if x != -1:
      start = counter
      break
    counter = counter + 1

  for index in range(start,start+numlines):
    result.append(contents[index])

  return result

def getKeyValuesNLinesExclusive(contents,keyword,numlines):
  result = []
  counter = 0
  start = 0
  for line in contents:
    x = line.find(keyword)
    if x != -1:
      start = counter + 1
      break
    counter = counter + 1
  
  for index in range(start,start+numlines):
    result.append(contents[index])
  
  return result

def test1():
  (contents,flatcontents) = parseFile('test.std',['!'])
  print('File')
  for row in contents:
    print(row)

  #print(flatcontents)
  key = '\$verify'

  values = getKeyValues(flatcontents,key)
  print('values')
  print(values)

  values = getKeyValuesAsFloat(flatcontents,key)
  print(values)

def test2():
  (contents,flatcontents) = parseFile('test03c.expected',['!'])
  
  #print(flatcontents)
  key = 'KSP Residual norm'
  
  values = getKeyValues(flatcontents,key)
  print('values')
  print(values)
  
  values = getKeyValuesAsFloat(flatcontents,key)
  print(values)

  key = 'dt_courant ='
  values = getKeyValues(flatcontents,key)
  print('values')
  print(values)

  values = getKeyValuesNLinesInclusive(contents,'_DataExCompleteCommunicationMap',8)
  print(values)
  values = getKeyValuesNLinesExclusive(contents,'_DataExCompleteCommunicationMap',7)
  print(values)

def test3():
  (contents,flatcontents) = parseFile('test03c.expected',['!'])
  
  expected = getKeyValuesNLinesExclusive(contents,'_DataExCompleteCommunicationMap',7)
  output = getKeyValuesNLinesExclusive(contents,'_DataExCompleteCommunicationMap',7)
  output[2] = 'aaa'
  status,err = compareLiteral(output,expected)
  print(err)

  status,err = compareFloatingPoint([4.4,3.3],0.01,[4.4,3.4])
  print(err)

  status,err = compareInteger(['4','3'],'2',['3'])
  print(err)


class pthUnitTest:
  
  def __init__(self, name,ranks,execute,expected_file):
    self.passed = -1
    self.walltime = 2.0 # minutes
    self.errormessage = ''
    self.errno = -1
    self.name = name
    self.ranks = ranks
    self.execute = execute
    self.expected_file = expected_file
    self.keywords = [ '#', '!', '//' ]
    self.output_file = name + '-p' + str(ranks) + '.output'
    self.comparison_file = ''
    self.output_path = ''
    self.ignore = False
    self.verbosity_level = 1

  def setVerbosityLevel(self,value):
    self.verbosity_level = value
  
  def verify(self,junk):
    raise RuntimeError('[pth] A valid verification method for unit-test \"' + self.name + '\" was not found.\n\
              [pth] You must provide each unit-test with a method to verify the output.\n\
              [pth] The method is set via calling test.setVerifyMethod()')

  def setOutputPath(self,opath):
    self.output_path = opath

  def setComparisonFile(self,fname):
    self.comparison_file = fname

  def setVerifyMethod(self,verify):
    self.verify = verify

  def setWalltime(self,mins):
    self.walltime = mins

  def appendKeywords(self,keywords):
    self.keywords.append(keywords)
  #print(self.keywords)

  def getExitCode(self):
    return self.errno

  def getErrorReport(self):
    return self.errormessage

  def getErrorStatus(self):
    return self.passed

  def verifyOutput(self):
    
    if self.comparison_file == '':
      self.comparison_file = os.path.join(self.output_path,self.output_file)
    
    if self.verbosity_level > 0:
      print('[Parsing file]',self.expected_file)
    (self.expected_contents,self.expected_flatcontents) = parseFile(self.expected_file,self.keywords)
    
    if self.verbosity_level > 0:
      print('[Parsing file]',self.comparison_file)
    (self.output_contents,self.output_flatcontents) = parseFile(self.comparison_file,self.keywords)
    self.verify(self)


  def getOutput(self):
    return self.output_contents,self.output_flatcontents


  def getExpected(self):
    return self.expected_contents,self.expected_flatcontents


  def updateStatus(self,status,err):
    if self.passed == -1:
      self.passed = status
      if err != '':
        self.errormessage = self.errormessage + err
      return

    if status == False:
      self.passed = False
      if err != '':
        self.errormessage = self.errormessage + err


  def report(self,type):
    if type == 'summary':
      
      if self.ignore == True:
        print(bcolors.WARNING + ' [' + self.name + ']   skipped' + bcolors.ENDC)
      else:
        if self.passed == False:
          print(bcolors.FAIL + ' [' + self.name + ']   *** FAILED ***' + bcolors.ENDC)
        else:
          print(bcolors.OKGREEN + ' [' + self.name + ']   passed' + bcolors.ENDC)

    if type == 'log':
      if self.ignore == False:
        if self.passed == False:
          #print(bcolors.FAIL +  '[' + self.name + '] reason for failure\n' + '--------------------------------------------------------------\n' + self.errormessage + bcolors.ENDC)
          print('[' + self.name + '] reason for failure\n' + '--------------------------------------------------------------\n' + self.errormessage)


  def compareFloatingPoint(self,key,tolerance):
    expected,expected_flat = self.getExpected()
    output,output_flat = self.getOutput()
    values_e = getKeyValuesAsFloat(expected_flat,key)
    if len(values_e) == 0:
      raise RuntimeError('[pth][VerificationError] Test \"' + self.name + '\" queried the expected file \"' + self.comparison_file + '\" for key \"' + key + '\" which was not found. \n\t\t    Users verification code is likely incorrect (contains a typo in the key name)' )
    
    values   = getKeyValuesAsFloat(output_flat,key)
    status,err = compareFloatingPoint(values,tolerance,values_e)
    kerr = ''
    if status == False:
      kerr = 'Key = \"' + key + '\" --> ' + err
    self.updateStatus(status,kerr)


  def compareInteger(self,key,tolerance):
    expected,expected_flat = self.getExpected()
    output,output_flat = self.getOutput()
    values_e = getKeyValuesAsInt(expected_flat,key)
    if len(values_e) == 0:
      raise RuntimeError('[pth][VerificationError] Test \"' + self.name + '\" queried the expected file \"' + self.comparison_file + '\" for key \"' + key + '\" which was not found. \n\t\t    Users verification code is likely incorrect (contains a typo in the key name)' )

    values   = getKeyValuesAsInt(output_flat,key)
    status,err = compareInteger(values,tolerance,values_e)
    kerr = ''
    if status == False:
      kerr = 'Key = \"' + key + '\" --> ' + err
    self.updateStatus(status,kerr)

  def compareLiteral(self,key):
    expected,expected_flat = self.getExpected()
    output,output_flat = self.getOutput()
    values_e = getKeyValues(expected_flat,key)
    if len(values_e) == 0:
      raise RuntimeError('[pth][VerificationError] Test \"' + self.name + '\" queried the expected file \"' + self.comparison_file + '\" for key \"' + key + '\" which was not found. \n\t\t    Users verification code is likely incorrect (contains a typo in the key name)' )
    
    values   = getKeyValues(output_flat,key)
    status,err = compareLiteral(values,values_e)
    kerr = ''
    if status == False:
      kerr = 'Key = \"' + key + '\" --> ' + err
    self.updateStatus(status,kerr)

  def compareUnixDiff(self):
    expected,expected_flat = self.getExpected()
    output,output_flat = self.getOutput()
    
    status,err = compareLiteral(expected,output)
    kerr = ''
    if status == False:
      kerr = 'Key = \"' + '<entire file>' + '\" --> ' + err
    self.updateStatus(status,kerr)
  

  def clean(self):
    outfile = os.path.join(self.output_path,self.output_file)
    cmpfile = self.comparison_file
    print('<test> rm -f ' + outfile)
    if outfile != cmpfile and cmpfile:
      print('<test> rm -f ' + cmpfile)


