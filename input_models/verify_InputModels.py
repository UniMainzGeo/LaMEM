#!/usr/bin/env python
import os,sys
import subprocess
sys.path.insert(0, "../tests/pythontestharness/lib")

# This is a testing framework to ensure that the files in the LaMEM repository
# in /input_models remain giving consistent results

import argparse

#import pyTestHarness.unittest as pth
import pyTestHarness.harness as pthharness
import pyTestHarness.launcher as launch

# Build optimized and debug versions of LaMEM
os.system('cd ../src/;  make mode=opt all; cd ../input_models')
os.system('cd ../src/;  make mode=deb all; cd ../input_models')

# Import the tests that are in each of the directories in /input_models
sys.path.append(os.path.join(os.environ['PWD'], 'BuildInSetups/Tests'))

# add matlab-tests if matlab is available as ENVIRONMENTAL variable MATLAB
if os.environ.get('MATLAB') != None:
    sys.path.append(os.path.join(os.environ['PWD'], 'SubductionWithMATLABParticles/Tests'))
else:
   print('MATLAB tests cannot be executed, as the environmental variable $MATLAB is not set')

import test_BuildInSetup as BuildIn

if os.environ.get('MATLAB') != None:
  import test_SubductionWithMATLABParticles as SubdMAT  # requires MATLAB to run


#if os.environ.get('MATLAB') != None:
#  import test_3_Subduction1               as Sub1 # import test that requires MATLAB
 
def makeLocalPathAbsolute(localRelPath) :
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  return(os.path.join(thisDir,localRelPath))

def run_tests():
  os.environ['PYTHONUNBUFFERED'] = str('1')

  registeredTests = [ BuildIn.test_1(), BuildIn.test_2(),  BuildIn.test_3(), BuildIn.test_4(),
                      BuildIn.test_5(), BuildIn.test_6(),  BuildIn.test_7(), BuildIn.test_8(),
                      BuildIn.test_9(), BuildIn.test_10(), BuildIn.test_11(), BuildIn.test_12()];
 

  # Add matlab tests (There should be a better way to do this for a range of files at the same time)
  #if os.environ.get('MATLAB') != None:
    #registeredTests.append(SubdMAT.test_1());


  if os.path.isdir('TestOutput') == False:
    os.mkdir('TestOutput')

  for test in registeredTests:
    test.setOutputPath(makeLocalPathAbsolute('TestOutput'))

  # Run the tests:
  h = pthharness.Harness(registeredTests)
  h.execute()
  h.verify()

  os.system('cp pthErrorReport.log pthErrorReport_Temp.log')

  h.clean()   # clean output (remove sandbox test directories); as this also removes the logfile, we have the commands above/below to circumvent this)

  os.system('mv pthErrorReport_Temp.log pthErrorReport.log')


if __name__ == "__main__":
  run_tests()

