#!/usr/bin/env python
import os,sys
sys.path.append(os.path.join(os.environ['PWD'], 'pyTestHarness'))

import argparse

import pyTestHarness.unittest as pth
import pyTestHarness.launch as launch


# Import separate tests 
sys.path.append(os.path.join(os.environ['PWD'], 't1_FB1_Direct'))
sys.path.append(os.path.join(os.environ['PWD'], 't2_FB2_MG'))

import test_1_FB1 as FB1
import test_2_FB2 as FB2

def run_unittests_example1():
  os.environ['PYTHONUNBUFFERED'] = str('1')

  if os.path.isdir('output') == False:
    os.mkdir('output')

  # Register all tests
  registeredTests = [ FB1.test_a(), FB1.test_b(), FB1.test_c(), FB1.test_d(),
                      FB2.test_a() ]

  # Force output to written somewhere else, can be invoked using -o <path>
  for test in registeredTests:
    test.setOutputPath('output')

  # Clean current directory
  os.system('make clean')
  
  # Build optimized and debug versions of LaMEM
  os.system('cd ../src/;  make mode=opt all; cd ../tests')
  os.system('cd ../src/;  make mode=deb all; cd ../tests')
 
  launcher = launch.pthLaunch()

  # Filter tests <could be promoted into batch execute()/verify() methods
  args = launcher.args
  subset = []
  if args.test:
    print(registeredTests)
    tnames = args.test.split(',')
    for name in tnames:
      for t in registeredTests:
        if name == t.name:
          subset.append(t)
    if subset == []:
      raise RuntimeError('You requested to test a subset of registered tests, \n',
                         'but no registed test matched the name list provided')
    else:
      print(subset)
  else:
    subset = registeredTests

  launcher.executeTestSuite(subset)
  launcher.verifyTestSuite(subset)


  # Clean current directory after executing the tests
  os.system('make clean')


if __name__ == "__main__":
  run_unittests_example1()
