#!/usr/bin/env python
import os,sys
import subprocess
sys.path.append(os.path.join(os.environ['PWD'], 'pyTestHarness'))

import argparse

#import pyTestHarness.unittest as pth
import pyTestHarness.launch as launch
import pyTestHarness.harness as harness

# Import separate tests
sys.path.append(os.path.join(os.environ['PWD'], 't1_FB1_Direct'))
sys.path.append(os.path.join(os.environ['PWD'], 't2_FB2_MG'))
sys.path.append(os.path.join(os.environ['PWD'], 't4_Loc'))

# add matlab-tests if matlab is available as ENVIRONMENTAL variable MATLAB
if os.environ.get('MATLAB') != None:
    sys.path.append(os.path.join(os.environ['PWD'], 't3_SubductionMATLABinput'))
else:
    print('MATLAB tests cannot be executed, as the environmental variable $MATLAB is not set')

import test_1_FB1 as FB1
import test_2_FB2 as FB2
import test_4_localization as Loc1

if os.environ.get('MATLAB') != None:
  import test_3_Subduction1 as Sub1 # import test that requires MATLAB

def run_unittests_example1():
  os.environ['PYTHONUNBUFFERED'] = str('1')

  if os.path.isdir('output') == False:
    os.mkdir('output')

  # Register all non-MATLAB tests
  registeredTests = [ FB1.test_a(),  FB1.test_b(),  FB1.test_c(),  FB1.test_d(),
                      FB2.test_a(), Loc1.test_a(), Loc1.test_b(), Loc1.test_c()];

  # Add matlab tests (There should be a better way to do this for a range of files at the same time)
  if os.environ.get('MATLAB') != None:
    registeredTests.append(Sub1.test_a());
    registeredTests.append(Sub1.test_b());
    registeredTests.append(Sub1.test_c());
	registeredTests.append(Sub1.test_d());

  # Force output to written somewhere else, can be invoked using -o <path>
  for test in registeredTests:
    test.setOutputPath('output')

  # Build optimized and debug versions of LaMEM
  os.system('cd ../src/;  make mode=opt all; cd ../tests')
  os.system('cd ../src/;  make mode=deb all; cd ../tests')

  launcher = launch.pthLaunch();

  # Filter tests <could be promoted into batch execute()/verify() methods
  args = launcher.args
  subset = []
  if args.test:
    #print(registeredTests)
    tnames = args.test.split(',')
    for name in tnames:
      for t in registeredTests:
        if name == t.name:
          subset.append(t)
    if subset == []:
      raise RuntimeError('You requested to test a subset of registered tests, \n',
                         'but no registed test matched the name list provided')
    #else:
    #  print(subset)
  else:
    subset = registeredTests

  launcher = harness.pthHarness(subset)
  launcher.execute()
  launcher.verify()
  launcher.clean()

if __name__ == "__main__":
  run_unittests_example1()

  # Clean current directory after executing the tests
  os.system('make clean')


